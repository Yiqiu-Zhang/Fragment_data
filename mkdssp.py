import os
import subprocess


sed_cmd = '1i CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1'
bucket_base = 's3://AF_data/true_structure_dataset/pdb'
local_base = '/mnt/lustre/zhangyiqiu/Fragment_data'
fragment_base = 's3://Fragment_data'



for i in range(256):

    subprocess.call(['aws', 's3', 'cp', f'{bucket_base}/pdb_{i}/',
                     f'{local_base}/pdb/pdb_{i}', '--recursive'])
    subprocess.call(['mkdir', f'{local_base}/dssp/ds p_{i}'])
    files = os.listdir(f'{local_base}/pdb/pdb_{i}')

    for pdb_file_name in files:

        pdb_path = f'{local_base}/pdb/pdb_{i}/{pdb_file_name}'
        subprocess.call(['sed', '-i', sed_cmd, pdb_path])


        dssp_name = pdb_file_name.split('.')[0]+'.dssp'
        dssp_path = f'{local_base}/dssp/dssp_{i}/{dssp_name}'
        subprocess.call(['mkdssp', '-i', pdb_path, '-o', dssp_path])


    subprocess.call(['aws', 's3', 'cp', f'{local_base}/dssp/dssp_{i}',
                     f'{fragment_base}/dssp/dssp_{i}/', '--recursive'])
    subprocess.call(['rm', '-r', f'./pdb/pdb_{i}/*'])
    subprocess.call(['rm', '-r', f'./dssp/dssp_{i}/*'])