
'''
/mnt/lustre/zhangyiqiu/rosetta_src_2021.16.61629_bundle/main/source/bin/rosetta_scripts.default.linuxgccrelease
-in:file:s 2hle_remixed.pdb -jd2:delete_old_poses -nstruct 1 -out:chtimestamp -ex1 -ex2 -randomize_missing_coords
-ignore_unrecognized_res -overwrite
-parser:protocol test_peptiderive.xml -scorefile testscore

'''

from Bio.PDB.Chain import Chain
from Bio.PDB import PDBIO, PDBParser

my_chain = Chain("Z")

io = PDBIO()
parser = PDBParser()

structure = parser.get_structure('6fay', '6fay.pdb') # id - string, the id that will be used for the structure
reslist = []
for model in structure:
    model.add(my_chain)
    for chain in model:
        for res in chain:
            res_id = res.id[1]
            if  900< res_id < 1000:
                reslist.append(res)
        for res in reslist:
            chain.detach_child(res.id)
            my_chain.add(res)
io.set_structure(structure)
io.save('6fay_3.pdb')
