import Bio.PDB
import numpy
input_dir = ""
# output_dir = ""



def analyse_PDB(file):  # Calculation of binding sites for PDB
    # check filename
    print()
    parser = Bio.PDB.PDBParser()
    pdb_name = file[0:4]
    structure = parser.get_structure(pdb_name, file)
    compnds = []
    for record in structure.header['compound']:
        compnd = {}
        compnd['mol_id'] = record
        compnd['molecule'] = structure.header['compound'][str(record)]['molecule']
        compnd['chain'] = structure.header['compound'][str(record)]['chain']
        compnds.append(compnd)

    nanobody_list = []
    antigen_list = []
    # Simple filtering of COMPND records
    for item in compnds:
        # print(item)
        if 'spike' in item['molecule']:
            nanobody_list.extend([s[0] if len(s) > 1 else s for s in item['chain'].upper().split(', ')])
            pass
        else:
            antigen_list.extend([s[0] if len(s) > 1 else s for s in item['chain'].upper().split(', ')])
            pass
    print(nanobody_list)
    print(antigen_list)
    if len(nanobody_list) == 0 or len(antigen_list) == 0:
        return
        pass
    sequences = []
    label_lists = []
    temp_combined=[]
    for antigen_char in antigen_list:
        for nanobody_char in nanobody_list:
            # 获取纳米抗体和抗原的链
            try:
                nanobody = structure[0][nanobody_char.strip()]
                antigen = structure[0][antigen_char.strip()]
            except Exception as e:
                print(f"An error occurred: {e}")
                print(nanobody_list,antigen_list)
                print(f"The current parameters are: {pdb_name}-{antigen_char}-{nanobody_char}")
                exit()
            # Find the contact area between antibody and antigen
            ns = Bio.PDB.NeighborSearch(list(antigen.get_atoms()))
            contact_residues = []
            paratope_seq = ''
            sequence = ''
            for residue in nanobody:
                if Bio.PDB.is_aa(residue.get_resname(), standard=True):  # Determine whether it is an amino acid
                    sequence += Bio.PDB.Polypeptide.three_to_one(residue.get_resname())
                    for atom in residue.get_atoms():
                        close_atoms = ns.search(atom.coord, 4.5, level='A')  # The threshold defining the contact area is 4.5 Å
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
    # random_cut((0.8, 0.1, 0.1))
    print(analyse_PDB('8gz5.pdb'))  # sample file
    pass
