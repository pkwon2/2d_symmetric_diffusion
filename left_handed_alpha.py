from Bio import PDB
from Bio.PDB import DSSP
import os,sys,glob
import warnings
warnings.filterwarnings("ignore")

parser = PDB.PDBParser()

def is_left_handed(pdb):
    structure = parser.get_structure("structure", pdb)
    model = structure[0]
    dssp = DSSP(model, pdb)
    for resi, data in dssp.property_dict.items():
        #phi, psi = data[2], data[3]
        phi = data['PHI']
        if phi < 0:# and psi < 0:
            return True
    return False

pdbs = glob.glob("/net/scratch/lhtran/repeat_dev/chain_break_debug/*.pdb")
# #is_left_handed("/home/lhtran/projects/sm_binders/ef2_binder/EF2_0696_MET.pdb")
# for pdb in pdbs:
#     if is_left_handed(pdb):
#         print('left', pdb.split('/')[-1])

#pdb='/net/scratch/lhtran/repeat_dev/chain_break_debug/EF2_160_C2_dimmer_user_T_25_break_sym_16.pdb'
#is_left_handed(pdb)
####

pdb = "/net/scratch/lhtran/repeat_dev/chain_break_debug/EF2_160_C4_dimmer_14.pdb"

def calculate_phi_psi(residue):
    #try:
        phi = residue.get_phi()
        psi = residue.get_psi()
        return phi, psi
    #     if phi and psi:
    #         return phi, psi
    #     else:
    #         return None, None
    # except:
    #     return None, None

def is_left_handed_alpha_helix(residue):
    phi, psi = calculate_phi_psi(residue)
    if phi is not None and psi is not None:
        if phi < 0:# and psi < 0:
            return True
    return False

def find_left_handed_helices(structure):
    left_handed_helices = []
    for model in structure:
        for chain in model:
            phi_psi_values = []
            for residue in chain:
                if is_left_handed_alpha_helix(residue):
                    phi_psi_values.append((residue.get_id(), calculate_phi_psi(residue)))
                else:
                    if len(phi_psi_values) >= 4:
                        left_handed_helices.append(phi_psi_values)
                    phi_psi_values = []
            if len(phi_psi_values) >= 3:
                left_handed_helices.append(phi_psi_values)
    return left_handed_helices

# Load the protein structure
parser = PDB.PDBParser()
structure = parser.get_structure("protein", pdb)
left_handed_helices = find_left_handed_helices(structure)
print(left_handed_helices)

# for pdb in pdbs:
#     structure = parser.get_structure("protein", pdb)

#     left_handed_helices = find_left_handed_helices(structure)
#     if len(left_handed_helices) > 0:
#         print("Left-handed alpha-helical bundles found:")
#         for bundle in left_handed_helices:
#             print(bundle)
#     else:
#         print("No left-handed alpha-helical bundles found.")