from Bio.PDB.Chain import Chain
from Bio.PDB import PDBIO, PDBParser
from Bio.PDB.DSSP import DSSP
import os
import subprocess
import multiprocessing

MIN_PROT_L = 50
MAX_PEP_L = 25
MIN_PEP_L = 5
MIN_TOTAL_L = 150
MAX_RSA = 0.25 # The relative solvent accessibility, below which, the residue are determined as core residue
PROCESS_NUM = 96

turePDB_base = 's3://AF_data/true_structure_dataset/pdb'
local_base = '/mnt/lustre/zhangyiqiu/Fragment_data'
fragment_base = 's3://Fragment_data'

io = PDBIO()
parser = PDBParser()

def find_pep(core):
    """Return core fragment index list. e.g. [(start, end, frag_len)] """

    pep_list = []
    s = e = i = 0
    while i < len(core):

        if core[i] is None:  # if i is None, move to the next
            l = e - s + 1
            if MIN_PEP_L <= l:
                pep_list.append((s, e, l))

            s = e = 0
            i += 1

        else:  # if i is number
            if not s:
                s = e = core[i]
            else:
                e = core[i]
            i += 1

    return pep_list
def Process_dssp(dssp_files):

    while len(dssp_files) > 0:
        file = dssp_files.pop(0)
        name = file.split('.')[0]
        full_dssp_path = f'{local_base}/dssp/dssp_{i}/' + name + '.dssp'
        full_pdb_path = f'{local_base}/pdb/pdb_{i}/' + name + '.pdb'

        # From the .pdb file, get the structures and Pick the first conformer.
        base_structure = parser.get_structure(name, full_pdb_path)
        model = base_structure[0]
        # noinspection PyBroadException
        try:
            dssp = DSSP(model, full_dssp_path, dssp='mkdssp', file_type='DSSP')

            # Return residue idx if RSA < MAX_RSA (core residue)
            res_index = list(map(lambda x: x[0] - 1 if x[3] < MAX_RSA else None, dssp.property_list))
            pep_lib = find_pep(res_index)
            subprocess.call(['mkdir', f'{local_base}/frag/frag_{i}/{name}'])

            # Add two 0 length chians, one for the peptide fragment(ID Z), two for protein fragment(ID A B)
            prot_chain = model.child_list[0]
            prot_len = len(prot_chain)
            if prot_len < MIN_TOTAL_L:
                continue
            prot_chain.id = "A"
            pep_chain = Chain("Z")
            third_chain = Chain("B")
            model.add(pep_chain)
            model.add(third_chain)

            for start, end, frag_len in pep_lib:

                cp_struc = base_structure.copy()
                cp_model = cp_struc[0]
                prot_chain, pep_chain, third_chain = cp_model.child_list

                # if frag_len <= MAX_PEP_L: # if fragment length is less than 25, do not apply slide window
                pep = prot_chain.child_list[start:end + 1]

                # ignore prot fragment if the length is less than MIN_PROT_L
                # ignore residues that next the peptide fragment
                del_res = []
                if start < MIN_PROT_L < prot_len - end:
                    del_res = prot_chain.child_list[:end + 6]

                elif start > MIN_PROT_L > prot_len - end:
                    del_res = prot_chain.child_list[start - 5:]

                elif start > MIN_PROT_L < prot_len - end:
                    del_res = prot_chain.child_list[start - 5:]
                    third_chain_res = prot_chain.child_list[end + 6:]
                    for res in third_chain_res:
                        third_chain.add(res)

                # add residues to its chain
                for res in pep:
                    pep_chain.add(res)
                for res in del_res:
                    prot_chain.detach_child(res.id)

                frag_path = f'{local_base}/frag/frag_{i}/{name}/{start}_{end}.pdb'
                bucket_frag_path = f'{fragment_base}/frag/frag_{i}/{name}/{start}_{end}.pdb'
                io.set_structure(cp_struc)
                io.save(frag_path)
                subprocess.call(['aws', 's3', 'cp', frag_path, bucket_frag_path])
                subprocess.call(['rm', '-r', frag_path])
        except:
            print('Drop the pdb file As it is not valid')
            continue

for i in range(256):
    subprocess.call(['mkdir', f'{local_base}/frag/frag_{i}'])
    files = os.listdir(f'{local_base}/dssp/dssp_{i}')

    threads = []
    sub_l = len(files) // PROCESS_NUM
    for n in range(PROCESS_NUM + 1):
        sub_process = multiprocessing.Process(target=Process_dssp, args=(files[n * sub_l:(n + 1) * sub_l],))
        threads.append(sub_process)

    for x in threads:
        x.start()
    for x in threads:
        x.join()



