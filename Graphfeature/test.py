

# export https_proxy="http://10.1.8.5:32680"

import os
import subprocess
import feature_pipeline as featurepipe
from multiprocessing import Pool

rosetta_base = '/mnt/lustre/zhangyiqiu/rosetta_src_2021.16.61629_bundle/main/source/bin'
local_base = '/mnt/lustre/zhangyiqiu/Fragment_data'
fragment_base = f's3://Fragment_data/frag'
mu = -49.48
std = 25.52
isc_range = [mu - std, mu + std]


def parse_peptiderive(peptiderive_str):
    descrip = []
    for line in peptiderive_str.splitlines():

        if not line or line.startswith('| 1'):  # Skip comment and blank lines
            continue
        elif line.startswith('> '):  # Chain_pair line
            info_list = line.split('=')
            receptor = info_list[1][1]
        elif line.startswith('>>'):  # peptide length
            peptide_length = float(line.split(':')[1])
        elif line.startswith('| 0'):  # Data line
            _, _, start_position, isc, _ = line.split(' ')
            if isc_range[0] < float(isc) < isc_range[1]:
                descrip.append([receptor, int(peptide_length), int(start_position),
                                float(isc)])
        elif line.startswith('# end chain pair'):
            return descrip


def Rosetta(file, i, fragroot):
    frag_position = file.split('.')[0]
    subprocess.call(['mkdir', '-p', f'{local_base}/pepderive/derive_{i}/{fragroot}/{frag_position}'])
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
                     '-parser:protocol test_peptiderive.xml'],
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL)

    pdbpath = f'{local_base}/pepderive/derive_{i}/{fragroot}/{frag_position}/{frag_position}_0001.pdb'

    with open(
            f'{local_base}/pepderive/derive_{i}/{fragroot}/{frag_position}/{frag_position}_0001.peptiderive.txt') as f:
        input_peptiderive_str = f.read()
        descriptions = parse_peptiderive(input_peptiderive_str)

    featurepipe.PDBtoFeature(descriptions, pdbpath, i, fragroot, frag_position)


for num in range(1, 2):# 0,256
    num = 1
    '''
    subprocess.call(['mkdir', f'{local_base}/frag/frag_{num}'])
    subprocess.call(['aws', 's3', 'cp', f'{fragment_base}/frag_{num}',
                     f'{local_base}/frag/frag_{num}', '--recursive'],
                    stdout=subprocess.DEVNULL)
    '''
    block = os.listdir(f'{local_base}/frag/frag_{num}/')

    for fragroot in block:
        frag_files = os.listdir(f'{local_base}/frag/frag_{num}/{fragroot}/')
        pool = Pool(2)  # 128

        for pdb_file in frag_files:
            pool.apply_async(func=Rosetta, args=(pdb_file, num, fragroot,))

        pool.close()
        pool.join()

    subprocess.call(['aws', 's3', 'cp', f'{local_base}/feature/feat_{num}/',
                     f's3://Fragment_data/feature/feat_{num}', '--recursive'],
                    stdout=subprocess.DEVNULL)
    subprocess.call(['rm', '-r', f'{local_base}/feature/feat_{num}/'])
    subprocess.call(['rm', '-r', f'{local_base}/pepderive/derive_{num}/'])
    subprocess.call(['rm', '-r', f'{local_base}/frag/frag_{num}/'])
