from __future__ import print_function

import subprocess

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
def PDBtoFeature(complex_info, pdbpath, block_num, fragroot, frag_position):

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

        local_feature_path = f'/mnt/lustre/zhangyiqiu/Fragment_data/feature/feat_{block_num}'
        f_peptide_path = f'peptide/{fragroot}/{frag_position}/{receptor}_{start_position}_{L_pep}'
        f_receptor_path = f'receptor/{fragroot}/{frag_position}/{receptor}'
        # eg. fragroot = 101M_A_1_renum 
        np.save(f'{local_feature_path}/{f_peptide_path}/binding_matrix_4.npy', binding_matrix)
        np.save(f'{local_feature_path}/{f_peptide_path}/binding_sites_4.npy', binding_sites)
        np.save(f'{local_feature_path}/{f_peptide_path}/target_sequence.npy', target_sequence)

        np.save(f'{local_feature_path}/{f_receptor_path}/nodes.npy', node)
        np.save(f'{local_feature_path}/{f_receptor_path}/edges.npy', edge)
        np.save(f'{local_feature_path}/{f_receptor_path}/neighbor_indices.npy', neighbor_indices)
        bucket_feature_path = f's3://Fragment_data/feature/feat_{block_num}'
        subprocess.call(['aws', 's3', 'cp',f'{local_feature_path}/{f_peptide_path}/',
                         f'{bucket_feature_path}/{f_peptide_path}/', '--recursive'])

        subprocess.call(['aws', 's3', 'cp', f'{local_feature_path}/{f_receptor_path}/',
                         f'{bucket_feature_path}/{f_receptor_path}/', '--recursive'])
        subprocess.call(['rm', '-r', f'{local_feature_path}/{f_peptide_path}/'])
        subprocess.call(['rm', '-r', f'{local_feature_path}/{f_receptor_path}/'])