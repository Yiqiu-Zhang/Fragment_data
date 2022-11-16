'''
/mnt/lustre/zhangyiqiu/rosetta_src_2021.16.61629_bundle/main/source/bin/rosetta_scripts.default.linuxgccrelease
-in:file:s 2hle_remixed.pdb -jd2:delete_old_poses -nstruct 1 -out:chtimestamp -ex1 -ex2 -randomize_missing_coords -ignore_unrecognized_res -overwrite -parser:protocol test_peptiderive.xml -scorefile testscore

'''

import os
import subprocess
import feature_pipeline as featurepipe

rosetta_base = '/mnt/lustre/zhangyiqiu/rosetta_src_2021.16.61629_bundle/main/source/bin'
local_base = '/mnt/lustre/zhangyiqiu/Fragment_data'
fragment_base = f's3://Fragment_data/frag'
# rosetta_base = '/home/PJLAB/zhangyiqiu/Documents/rosetta_src_2021.16.61629_bundle/main/source/bin'
#base = '/home/PJLAB/zhangyiqiu/PycharmProjects/Fragment_data/peptiderive_test/pdb'
#filename = '1QHP_A_1_renum_179_184.pdb'
# base_filename = filename.split('.')[0]

mu = -49.48
std =25.52
isc_range = [mu-std,mu+std]
# isc_range = [-1000,0] # test

def parse_peptiderive(peptiderive_str):

    descrip =[]
    for line in peptiderive_str.splitlines():

        if not line or line.startswith('| 1'): # Skip comment and blank lines
            continue
        elif line.startswith('> '): # Chain_pair line
            info_list = line.split('=')
            receptor = info_list[1][1]
        elif line.startswith('>>'): # peptide length
            peptide_length = float(line.split(':')[1])
        elif line.startswith('| 0'): # Data line
            _,_,start_position,isc,_  = line.split(' ')
            if isc_range[0]< float(isc) < isc_range[1]:
                descrip.append([receptor, peptide_length, float(start_position),
                                float(isc)])
        elif line.startswith('# end chain pair'):
            return descrip

for i in range(1):
    # subprocess.call(['aws', 's3', 'cp', f'{fragment_base}/frag_{i}', f'{local_base}/frag/frag_{i}', '--recursive'])
    block = os.listdir(f'{local_base}/frag/frag_{i}')
    for fragroot in block:
        frag_files = os.listdir(f'{local_base}/frag/frag_{i}/{fragroot}')
        for file in frag_files:
            frag_position = file.split('.')[0]
            subprocess.call(['mkdir','-p', f'{local_base}/pepderive/derive_{i}/{fragroot}/{frag_position}'])
            subprocess.call([f'{rosetta_base}/rosetta_scripts.default.linuxgccrelease',
                            '-in:file:s', f'{local_base}/frag/frag_{i}/{fragroot}/{file}',
                            '-out:path:all', f'{local_base}/pepderive/derive_{i}/{fragroot}/{frag_position}',
                            '-jd2:delete_old_poses',
                            '-nstruct 1',
                            '-out:chtimestamp',
                            '-ex1',
                            '-ex2',
                            '-randomize_missing_coords',
                            '-ignore_unrecognized_res',
                            '-overwrite',
                            '-parser:protocol test_peptiderive.xml'])

            pdbpath = f'{local_base}/pepderive/derive_{i}/{fragroot}/{frag_position}/{frag_position}_0001.pdb'

            with open(f'{local_base}/pepderive/derive_{i}/{fragroot}/{frag_position}/{frag_position}_0001.peptiderive.txt') as f:
                input_peptiderive_str = f.read()
                descriptions = parse_peptiderive(input_peptiderive_str)

            featurepipe.PDBtoFeature(descriptions, pdbpath, i, fragroot, frag_position)
            subprocess.call(['rm', '-r', f'{local_base}/pepderive/derive_{i}/{fragroot}/{frag_position}/'])
        subprocess.call(['rm', '-r', f'{local_base}/frag/frag_{i}/{fragroot}/'])

            # keep receprot in the pdb file,
            # remove the third chain
            # remove the unused peptide area
            # put the receptor and the peptide in the same pdb
