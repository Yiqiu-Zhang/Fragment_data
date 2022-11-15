import Graphfeature.feature_pipeline as featurepipe
rosetta_base = 'mnt/lustre/zhangyiqiu/rosetta_src_2021.16.61629_bundle/main/source/bin'
rosetta_base = '/home/PJLAB/zhangyiqiu/Documents/rosetta_src_2021.16.61629_bundle/main/source/bin'
base = '/home/PJLAB/zhangyiqiu/PycharmProjects/Fragment_data/peptiderive_test/pdb'
filename = '1QHP_A_1_renum_179_184.pdb'
base_filename = filename.split('.')[0]

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

with open(f'{base}/{base_filename}_0001.peptiderive.txt') as f:
    pdbpath = f'{base}/{base_filename}_0001.pdb'
    input_peptiderive_str = f.read()
    descriptions = parse_peptiderive(input_peptiderive_str)
    featurepipe.PDBtoFeature(descriptions, pdbpath)

    # keep receprot in the pdb file,
    # remove the third chain
    # remove the unused peptide area
    # put the receptor and the peptide in the same pdb








