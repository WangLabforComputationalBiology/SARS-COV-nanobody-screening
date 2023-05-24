import os
import shutil

def copyfile_base_keyword_in_file(file_pathname):
    # 遍历该目录下的所有文件
    for filename in os.listdir(file_pathname):
        path = os.path.join(file_pathname, filename)
        with open(path, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if 'NANOBODY' in line:
                    src_file_path = os.path.join(file_pathname, filename)
                    dst_file_path = os.path.join(to_dir_path, filename)
                    shutil.copy(src_file_path, dst_file_path)

# 指定要查找的文件夹路径
file_pathname = r'E:\Download\all_structures\SARS-COV'

# 指定要复制到的目标文件夹路径
to_dir_path = r'E:\Download\all_structures\SARS-COV nanobody'

# 调用函数，查找含有关键字的文件并复制到目标文件夹中
copyfile_base_keyword_in_file(file_pathname)