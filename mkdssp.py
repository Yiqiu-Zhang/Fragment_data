import os
import subprocess

#pdb_path = '/home/PJLAB/zhangyiqiu/PycharmProjects/Fragment_data/testfile/6fay_PepID_124_130.pdb'
#dssp_path = '/home/PJLAB/zhangyiqiu/PycharmProjects/Fragment_data/test.dssp'
#subprocess.call(['sed', '-i', sed_cmd, pdb_path])

sed_cmd = '1i CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1'
bucket_base = 's3://AF_data/true_structure_dataset/pdb'
local_base = '/mnt/lustre/zhangyiqiu/Fragment_data/pdb'
fragment_base = 's3://Fragment_data'



for i in range(1):
    pdb_path= f'{bucket_base}/pdb_{i}/'
    cp_cmd = f'aws s3 cp  . --recursive'
    '''
    subprocess.call(['aws', 's3', 'cp', f'{bucket_base}/pdb_{i}/',
                     f'{local_base}/pdb_{i}', '--recursive'])
    '''
    files = os.listdir(f'{local_base}/pdb_{i}')
    for pdb_file_name in files:
        print(pdb_file_name)
        pdb_path = f'{local_base}/pdb_{i}/{pdb_file_name}'
        print(pdb_path)
        subprocess.call(['sed', '-i', sed_cmd, pdb_path])
        print('Successful sed')

        dssp_name = pdb_file_name.split('.')[0]+'.dssp'
        dssp_path = f'{local_base}/dssp_{i}/{dssp_name}'
        subprocess.call(['mkdssp', '-i', pdb_path, '-o', dssp_path])
        print('Successful mkdssp')
    subprocess.call(['aws', 's3', 'cp', f'{local_base}/dssp_{i}',
                     f'{fragment_base}/dssp/dssp_{i}/', '--recursive'])