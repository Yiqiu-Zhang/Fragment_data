import os
import subprocess
import multiprocessing

PROCESS_NUM = 2

sed_cmd = '1i CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1'
bucket_base = 's3://AF_data/true_structure_dataset/pdb'
local_base = '/mnt/lustre/zhangyiqiu/Fragment_data'
fragment_base = 's3://Fragment_data'


def DSSP(pdb_names,i):
    for pdb_name in pdb_names:
        pdb_path = f'{local_base}/pdb/pdb_{i}/{pdb_name}'
        subprocess.call(['sed', '-i', sed_cmd, pdb_path])  # Adding line to .pdb make it a valid pdb file

        dssp_name = pdb_name.split('.')[0] + '.dssp'
        dssp_path = f'{local_base}/dssp/dssp_{i}/{dssp_name}'
        subprocess.call(['mkdssp', '-i', pdb_path, '-o', dssp_path])


for index in range(100,256): #

    subprocess.call(['aws', 's3', 'cp', f'{bucket_base}/pdb_{index}/',
                     f'{local_base}/pdb/pdb_{index}', '--recursive'])
    subprocess.call(['mkdir', f'{local_base}/dssp/dssp_{index}'])
    files = os.listdir(f'{local_base}/pdb/pdb_{index}')

    threads = []
    sub_l = len(files) // PROCESS_NUM
    for n in range(PROCESS_NUM + 1):
        sub_process = multiprocessing.Process(target=DSSP, args=(files[n * sub_l:(n + 1) * sub_l],index,))
        threads.append(sub_process)

    for x in threads:
        x.start()
    for x in threads:
        x.join()

    subprocess.call(['aws', 's3', 'cp', f'{local_base}/dssp/dssp_{index}',
                     f'{fragment_base}/dssp/dssp_{index}/', '--recursive'])
    subprocess.call(['rm', '-r', f'./pdb/pdb_{index}/*'])
    subprocess.call(['rm', '-r', f'./dssp/dssp_{index}/*'])


