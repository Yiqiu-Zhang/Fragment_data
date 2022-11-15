from __future__ import print_function
import numpy as np
import utils
import torch
import torch.nn.functional as F
from Bio.PDB import PDBIO, PDBParser, Polypeptide

# Test example, will be array latter
"""
name = '1cso'
L_pep = 11
start_position = 9  - 1
pdbpath = '/home/PJLAB/zhangyiqiu/PycharmProjects/Fragment_data/6fay_PepID_26_42_0001.pdb'
pdbpath = '1cso.pdb'
receptor= 'E'
"""

'''parse a PDB file contains BOTH protein chain and peptide chain, 
Calculate binding_matrix & binding_sites
'''
def PDBtoFeature(complex_info, pdbpath):
    for info in complex_info:
        receptor, L_pep, start_position, _ = info
        start_position -= 1
        parser = PDBParser()
        base_structure = parser.get_structure('1cso', pdbpath)
        model = base_structure[0]
        prot = model[receptor]
        peptide = model['I']
        peptide = peptide.child_list[start_position: start_position+L_pep]
        pep_sequence = np.asarray([Polypeptide.three_to_index(seq.get_resname()) for seq in peptide])
        pep_sequence = torch.from_numpy(pep_sequence).to(dtype=torch.int64)

        # calculate Fragment data
        binding_matrix = np.asarray([[utils.res_binding(a,z) for a in prot] for z in peptide])
        binding_sites = np.asarray([i for i in range(len(prot)) if any(binding_matrix[:,i])>0 ])
        target_sequence = F.one_hot(pep_sequence,20)

        node, edge, neighbor_indices = utils.get_graph(prot)

        np.save('binding_matrix_4.npy', binding_matrix)
        np.save('binding_sites.npy', binding_sites)
        np.save('target_sequence.npy', target_sequence)

        np.save('nodes.npy', node)
        np.save('edges.npy', edge)
        np.save('neighbor_indices.npy', neighbor_indices)