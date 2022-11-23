import torch
import torch.nn.functional as F
import numpy as np
from Bio.PDB import Polypeptide

def res_binding(prot_res,pep_res):
    x1 = np.asarray([atom.get_coord() for atom in prot_res])
    x2 = np.asarray([atom.get_coord() for atom in pep_res])

    x1 = np.expand_dims(x1,axis=0)
    x2 = np.expand_dims(x2,axis=1)
    D = np.sqrt(np.sum((x1-x2)**2,2))
    return 1. if np.min(D) < 4 else 0.

# The following gather functions
def gather_edges(edges, neighbor_idx):
    # Features [B,N,N,C] at Neighbor indices [B,N,K] => Neighbor features [B,N,K,C]
    neighbors = neighbor_idx.unsqueeze(-1).expand(-1, -1, -1, edges.size(-1))
    edge_features = torch.gather(edges, 2, neighbors)
    return edge_features

def gather_nodes(nodes, neighbor_idx):
    # Features [B,N,C] at Neighbor indices [B,N,K] => [B,N,K,C]
    # Flatten and expand indices [N,K] => [NK] => [NK,C]
    neighbors_flat = neighbor_idx.view((-1))
    neighbors_flat = neighbors_flat.unsqueeze(-1).expand(-1, nodes.size(1))
    # Gather and re-pack
    neighbor_features = torch.gather(nodes, 0, neighbors_flat)
    neighbor_features = neighbor_features.view(list(neighbor_idx.shape)[:2] + [-1])
    return neighbor_features

def quaternions(R):
    """ Convert a batch of 3D rotations [R] to quaternions [Q]
        R [...,3,3]
        Q [...,4]
    """
    diag = torch.diagonal(R, dim1=-2, dim2=-1)
    Rxx, Ryy, Rzz = diag.unbind(-1)
    magnitudes = 0.5 * torch.sqrt(torch.abs(1 + torch.stack([
        Rxx - Ryy - Rzz,
        - Rxx + Ryy - Rzz,
        - Rxx - Ryy + Rzz
    ], -1)))
    _R = lambda i, j: R[:, :, i, j]
    signs = torch.sign(torch.stack([
        _R(2, 1) - _R(1, 2),
        _R(0, 2) - _R(2, 0),
        _R(1, 0) - _R(0, 1)
    ], -1))
    xyz = signs * magnitudes
    # The relu enforces a non-negative trace
    w = torch.sqrt(F.relu(1 + diag.sum(-1, keepdim=True))) / 2.
    Q = torch.cat((xyz, w), -1)
    Q = F.normalize(Q, dim=-1)

    return Q

def dist(X, mask,eps=1E-6, top_k=30,):
    mask_2D = torch.unsqueeze(mask, 0) * torch.unsqueeze(mask, 1)
    dX = torch.unsqueeze(X,0) - torch.unsqueeze(X,1)
    D = mask_2D * torch.sqrt(torch.sum(dX**2, 2) + eps)

    D_max, _ = torch.max(D, -1, keepdim=True)
    D_adjust = D + (1. - mask_2D) * D_max
    D_neighbors, E_idx = torch.topk(D_adjust, top_k, dim=-1, largest=False)

    # mask_neighbors = gather_edges(mask_2D.unsqueeze(-1), E_idx)

    return D_neighbors, E_idx

def orientations_coarse(X, E_idx, S, eps=1E-6):

    # Shifted slices of unit vectors
    dX = X[1:, :] - X[:-1, :] # x_i - x_(i-1)
    U = F.normalize(dX, dim=-1) #Unit vectors
    u_2 = U[:-2, :] #[x_1 - x_0,
    u_1 = U[1:-1, :] #[x_2 - x_1,
    u_0 = U[2:, :] # #[x_3 - x_2，


    # Backbone normals
    n_2 = F.normalize(torch.cross(u_2, u_1), dim=-1) # normal of x_1
    n_1 = F.normalize(torch.cross(u_1, u_0), dim=-1) # normal of x_2

    # Bond angle calculation
    cosA = -(u_1 * u_0).sum(-1) # x_2
    cosA = torch.clamp(cosA, -1 + eps, 1 - eps)
    A = torch.acos(cosA)

    # Angle between normals
    cosD = (n_2 * n_1).sum(-1)
    cosD = torch.clamp(cosD, -1 + eps, 1 - eps)
    D = torch.sign((u_2 * n_1).sum(-1)) * torch.acos(cosD)

    # Backbone features
    #AD_features = torch.stack((torch.cos(A), torch.sin(A) * torch.cos(D), torch.sin(A) * torch.sin(D)), 2)
    #AD_features = F.pad(AD_features, (0, 0, 1, 2), 'constant', 0)

    # Build relative orientations
    o_1 = F.normalize(u_2 - u_1, dim=-1)
    ori = torch.stack((o_1, n_2, torch.cross(o_1, n_2)), 2)
    ori = ori.view([ori.shape[0]] + [9])
    ori = F.pad(ori, (0, 0, 1, 2), 'constant', 0) # why pad? Pad ori to shape [L,9]

    O_neighbors = gather_nodes(ori, E_idx) # ori [L,9] E_idx [L,30]
    X_neighbors = gather_nodes(X, E_idx) # X [L,3]

    # Re-view as rotation matrices
    ori = ori.view([ori.shape[0]] + [3, 3])
    O_neighbors = O_neighbors.view(list(O_neighbors.shape[:2]) + [3, 3])

    # Rotate into local reference frames
    dX = X_neighbors - X.unsqueeze(-2)
    dU = torch.matmul(ori.unsqueeze(1), dX.unsqueeze(-1)).squeeze(-1) # rotate the dX vector and return the rotated Vector
    dU = F.normalize(dU, dim=-1)
    R = torch.matmul(ori.unsqueeze(1).transpose(-1, -2), O_neighbors)
    Q = quaternions(R)

    # Orientation features
    O_features = torch.cat((dU, Q), dim=-1)

    # Side Orientation
    dS = S - X
    dU = torch.matmul(ori, dS.unsqueeze(-1)).squeeze(-1)
    O_side = F.normalize(dU, dim=-1)

    return O_features, O_side# AD_features

def rbf(D, side=False):
    # Distance radial basis function
    if not side:
        D_min, D_max, D_count = 0., 20., 16
        D_mu = torch.linspace(D_min, D_max, D_count) # Here should +1 !!!!!!!
        D_mu = D_mu.view([1, 1, -1])
        D_sigma = (D_max - D_min) / D_count
        D_expand = torch.unsqueeze(D, -1)
        RBF = torch.exp(-((D_expand - D_mu) / D_sigma) ** 2)

    else:
        D_min, D_max, D_count = 0., 6., 3
        D_mu = torch.linspace(D_min, D_max, D_count)# Here should +1 !!!!!!!
        D_mu = D_mu.view([1, -1])
        D_sigma = (D_max - D_min) / D_count
        D_expand = torch.unsqueeze(D, -1)
        RBF = torch.exp(-((D_expand - D_mu) / D_sigma) ** 2)

    return RBF

def P_embeddings(E_idx, num_positional_embeddings = 16):

    N_nodes = E_idx.size(0)
    N_neighbors = E_idx.size(1)
    ii = torch.arange(N_nodes, dtype=torch.float32).view((-1, 1))
    d = (E_idx.float() - ii).unsqueeze(-1) # sequence level relative distance
    # Original Transformer frequencies
    frequency = torch.exp(
        torch.arange(0, num_positional_embeddings, 2, dtype=torch.float32)
        * -(np.log(10000.0) / num_positional_embeddings)
    )
    angles = d * frequency.view((1, 1, -1))
    E = torch.cat((torch.cos(angles), torch.sin(angles)), -1)
    return E

def dihedrals(X, eps = 1e-7):

    # First 3 coordinates are N, CA, C
    X = X[:, :3, :].reshape(3 * X.shape[0], 3)

    # Shifted slices of unit vectors
    dX = X[1:, :] - X[:-1, :]
    U = F.normalize(dX, dim=-1)
    u_2 = U[:-2, :]  # [x_1 - x_0,
    u_1 = U[1:-1, :]  # [x_2 - x_1,
    u_0 = U[2:, :]  # #[x_3 - x_2，

    # Backbone normals
    n_2 = F.normalize(torch.cross(u_2, u_1), dim=-1)# normal of x_1 (x0,x2) plane
    n_1 = F.normalize(torch.cross(u_1, u_0), dim=-1)# normal of x_2 (x1,x3) plane

    # Angle between normals
    cosD = (n_2 * n_1).sum(-1)
    cosD = torch.clamp(cosD, -1 + eps, 1 - eps)
    D = torch.sign((u_2 * n_1).sum(-1)) * torch.acos(cosD)

    # This scheme will remove phi[0], psi[-1], omega[-1]
    D = F.pad(D, (1, 2), 'constant', 0)
    D = D.view((int(D.size(0) / 3), 3)) # phi, psi, omega

    # Lift angle representations to the circle
    D_features = torch.cat((torch.cos(D), torch.sin(D)), 1)

    return D_features

def side_central(A):

    if len(A) >4:
        return np.mean([atom.get_coord() for atom in A[4:]], axis=0)
    else:
        return A[1].get_coord()

def side_node(B_ca,S):
    relative_S = torch.linalg.norm(B_ca - S, dim =1)
    return rbf(relative_S,side=True)

def get_graph(protein_chain):
    prot_info = {
        'chain_bb_x': [],
        'child_x': [],
        'res_num': [],
        'res_type': [],
        'atom_type': [],
    }
    res_name = []
    chain_bb_x = []
    side_C = []

    for prot_res in protein_chain:
        if not Polypeptide.is_aa(prot_res.get_resname(), standard=True):
            continue

        res_name.append(prot_res.get_resname())

        all_atoms = [atom.get_name().strip() for atom in prot_res]
        del_res = []

        # skip residues lacking a backbone atom
        if "C" not in all_atoms or "CA" not in all_atoms or "N" not in all_atoms:
            del_res.append(prot_res)
            # protein_chain.detach_child(prot_res.id)
            continue
        sorted(prot_res)# In order [<Atom N>, <Atom CA>, <Atom C>, <Atom O>]
        atom_list = prot_res.child_list
        chain_bb_x.append([atom.get_coord() for atom in atom_list[:4]])
        side_C.append(side_central(atom_list))

    onehot_prot = np.asarray([Polypeptide.three_to_index(res) for res in res_name])
    onehot_prot = torch.from_numpy(onehot_prot).to(dtype=torch.int64)
    onehot_prot = F.one_hot(onehot_prot, 20)

    X = torch.from_numpy(np.asarray(chain_bb_x)).to(dtype=torch.float32)
    X_ca = X[:, 1, :]
    side_C = torch.from_numpy(np.asarray(side_C)).to(dtype=torch.float32)

    mask = torch.isfinite(torch.sum(X_ca, 1)).to(dtype=torch.float32)

    dihedral_node = dihedrals(X)
    RBF_side_node = side_node(X_ca,side_C)

    D_neighbors, E_idx = dist(X_ca, mask)

    # Pairwise features
    O_features, O_side = orientations_coarse(X_ca, E_idx, side_C)  # AD_features
    RBF = rbf(D_neighbors)

    # Pairwise embeddings
    E_positional = P_embeddings(E_idx)

    # Use features_type = 'full'
    V = torch.cat((onehot_prot,dihedral_node,RBF_side_node,O_side), -1)
    E = torch.cat((E_positional, RBF, O_features), -1)

    return V, E, E_idx


