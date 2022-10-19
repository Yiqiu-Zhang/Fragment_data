
'''
/mnt/lustre/zhangyiqiu/rosetta_src_2021.16.61629_bundle/main/source/bin/rosetta_scripts.default.linuxgccrelease
-in:file:s 2hle_remixed.pdb -jd2:delete_old_poses -nstruct 1 -out:chtimestamp -ex1 -ex2 -randomize_missing_coords -ignore_unrecognized_res -overwrite -parser:protocol test_peptiderive.xml -scorefile testscore

'''

from Bio.PDB.Chain import Chain
from Bio.PDB import PDBIO, PDBParser
from Bio.PDB.DSSP import DSSP

MIN_PROT_L = 50
MAX_PEP_L = 25
MIN_PEP_L = 5


io = PDBIO()
parser = PDBParser()

structure = parser.get_structure('6FAY', '6fay.pdb') # id - string, the id that will be used for the structure
res_list = []

model = structure[0] # Pick the first conformer.

# calculate DSSP, get the secondary structure information.
# mkdssp -i 1mot.pdb  -o 1mot.dssp
dssp = DSSP(model, '6fay.dssp', dssp='mkdssp', file_type='DSSP')#

res_index = list(map(lambda x: x[0]-1 if x[3]<0.25 else None,dssp.property_list))

sequence = list(map(lambda x: x[1], dssp.property_list))
sec_structure = list(map(lambda x: x[2], dssp.property_list))
rel_acc = list(map(lambda x: x[3], dssp.property_list))

def find_pep(res):

    pep_list = []
    start = end = i = 0
    while i < len(res):

        if res[i] is None:  # if i is None, move to the next
            frag_len = end - start + 1
            if MIN_PEP_L <= frag_len:
                pep_list.append((start, end, frag_len))

            start = end = 0
            i += 1

        else:  # if i is number
            if not start:
                start = end = res[i]
            else:
                end = res[i]
            i += 1

    return pep_list

pep_lib = find_pep(res_index)

prot_chain = model.child_list[0]
prot_len = len(prot_chain)

prot_chain.id = "A"
pep_chain = Chain("Z")
third_chain = Chain("B")
model.add(pep_chain)
model.add(third_chain)


for start, end, frag_len in pep_lib:

    # if frag_len <= MAX_PEP_L: # if fragment length is less than 25, do not apply slide window
    pep = prot_chain.child_list[start:end+1]  # -1 because peptide label start from 1

    if start < MIN_PROT_L < prot_len-end: # ignore residues before the peptide
        del_res = prot_chain.child_list[:end+6]  # The residues that will be deleted from the main chain

    elif start > MIN_PROT_L > prot_len - end:
        del_res = prot_chain.child_list[start-5:]

    elif start > MIN_PROT_L < prot_len - end:
        del_res = prot_chain.child_list[start-5:]
        third_chain_res = prot_chain.child_list[end+6:]
        for res in third_chain_res:
            third_chain.add(third_chain_res)

    for res in pep:
        pep_chain.add(res)
    # noinspection PyUnboundLocalVariable
    for res in del_res:
        prot_chain.detach_child(res.id)

    io.set_structure(structure)
    io.save(f'testfile/6fay_PepID_{start}_{end}.pdb')

