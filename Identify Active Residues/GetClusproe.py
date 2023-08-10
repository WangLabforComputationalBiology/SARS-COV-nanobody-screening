#------------------------------
# This is a code for searching the active residues in antigen chain.
# The mark with”-” for non-active residues.
# Output the ID of the active residues in [...], followed by their total number.
#------------------------------


import Bio.PDB
import numpy
input_dir = ""
# output_dir = ""

def analyse_PDB(file):  
    print()
    parser = Bio.PDB.PDBParser()
    pdb_name = file[0:4]
    structure = parser.get_structure(pdb_name, file)

    chains = list(structure.get_chains())
    antigen_output = []
    nanobody_list = []

 # Judging by the chain length: the longer one is the antigen chain,and another is the nanobody.
    if len(chains[0]) < len(chains[1]):
        nanobody_list.append(chains[0].get_id())
        antigen_output.append(chains[1].get_id())
    else:
        nanobody_list.append(chains[1].get_id())
        antigen_output.append(chains[0].get_id())


    sequences = []
    label_lists = []
    temp_combined=[]
    for antigen_char in nanobody_list:
        for nanobody_char in antigen_output:
            # Acquire nanobodies and antigen chains
            try:
                nanobody = structure[0][nanobody_char.strip()]
                antigen = structure[0][antigen_char.strip()]
            except Exception as e:
                print(f"An error occurred: {e}")
                print(antigen_output,nanobody_list)
                print(f"The current parameters are: {pdb_name}-{antigen_char}-{nanobody_char}")
                exit()
            # Find the binding site between nanobody and antigen
            ns = Bio.PDB.NeighborSearch(list(antigen.get_atoms()))
            contact_residues = []
            paratope_seq = ''
            sequence = ''
            for residue in nanobody:
                if Bio.PDB.is_aa(residue.get_resname(), standard=True):  
                    sequence += Bio.PDB.Polypeptide.three_to_one(residue.get_resname())
                    for atom in residue.get_atoms():
                        close_atoms = ns.search(atom.coord, 4.5, level='A')   # Set the threshold to 4.5 Å
                        if len(close_atoms) > 0:
                            paratope_seq += Bio.PDB.Polypeptide.three_to_one(residue.get_resname())
                            contact_residues.append(residue)
                            break
                    else:
                        paratope_seq += '-'
            contact_positions = [residue.get_id()[1] for residue in contact_residues]
            # print("#" * 20, nanobody_char, 'against', antigen_char, "#" * 20)
            if len(contact_positions):
                # print('Y')
                print(paratope_seq)
                print(contact_positions, len(contact_positions))
                label_list_temp = ["N" if c == "-" else "P" for c in paratope_seq]
                sequences.append(sequence)
                label_lists.append(label_list_temp)
                temp_combined.append(pdb_name + ' ' + nanobody_char + ' anti ' + antigen_char)
                pass
            else:
                # print("NONE")
                pass

    return sequences, label_lists, temp_combined

if __name__ == '__main__':

    print(analyse_PDB('modl.pdb'))
    pass
