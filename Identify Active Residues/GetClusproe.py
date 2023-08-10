#!wget https://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/archive/all_nano/
import Bio.PDB
import numpy
input_dir = ""
# output_dir = ""

#？改：原子量更大的--设为抗原链

def analyse_PDB(file):  # 计算PDB的结合位点
    # 检查文件名
    print()
    parser = Bio.PDB.PDBParser()
    pdb_name = file[0:4]
    structure = parser.get_structure(pdb_name, file)

#########################################
####获取两条链-原子量更大的--设为抗原链 #######
#########################################
    chains = list(structure.get_chains())
    antigen_output = []
    nanobody_list = []

    # antigen_output 输出链    另一条结合链
    if len(chains[0]) < len(chains[1]):
        nanobody_list.append(chains[0].get_id())# 原子量更大的链
        antigen_output.append(chains[1].get_id())# 原子量更小的链
    else:
        nanobody_list.append(chains[1].get_id())# 原子量更大的链
        antigen_output.append(chains[0].get_id())# 原子量更小的


    sequences = []
    label_lists = []
    temp_combined=[]
    for antigen_char in nanobody_list:
        for nanobody_char in antigen_output:
            # 获取纳米抗体和抗原的链
            try:
                nanobody = structure[0][nanobody_char.strip()]
                antigen = structure[0][antigen_char.strip()]
            except Exception as e:
                print(f"An error occurred: {e}")
                print(antigen_output,nanobody_list)
                print(f"The current parameters are: {pdb_name}-{antigen_char}-{nanobody_char}")
                exit()
            # 查找抗体和抗原之间的接触面积
            ns = Bio.PDB.NeighborSearch(list(antigen.get_atoms()))
            contact_residues = []
            paratope_seq = ''
            sequence = ''
            for residue in nanobody:
                if Bio.PDB.is_aa(residue.get_resname(), standard=True):  # 判断是否是氨基酸
                    sequence += Bio.PDB.Polypeptide.three_to_one(residue.get_resname())
                    for atom in residue.get_atoms():
                        close_atoms = ns.search(atom.coord, 4.5, level='A')  # 定义接触面积的阈值为4.5埃
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
    print(analyse_PDB('model.000.00.pdb'))
    pass
