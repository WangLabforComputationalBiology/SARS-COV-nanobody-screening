import os
import shutil

def copyfile_base_keyword_in_file(file_pathname):
    # Traverse all files in this directory
    for filename in os.listdir(file_pathname):
        path = os.path.join(file_pathname, filename)
        with open(path, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if 'NANOBODY' in line:
                    src_file_path = os.path.join(file_pathname, filename)
                    dst_file_path = os.path.join(to_dir_path, filename)
                    shutil.copy(src_file_path, dst_file_path)

# Specify the folder path to look for
file_pathname = r'E:\Download\all_structures\SARS-COV'

# Specify the destination folder path to copy to
to_dir_path = r'E:\Download\all_structures\SARS-COV nanobody'

# Call the function to find the file containing the keyword and copy it to the target folder
copyfile_base_keyword_in_file(file_pathname)
