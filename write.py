# coding=utf-8
import os
import xlwt  # Operate the excel module
import sys

file_path = sys.path[0] + '\\filenamelist.xls'  # sys.path[0] is to get the current path, filenamelist is the file to be written
f = xlwt.Workbook(encoding='utf-8', style_compression=0)  # Create a new excel
sheet = f.add_sheet('sheet1')  # Create a new sheet
pathDir = os.listdir(sys.path[0])  # The file is placed in the current folder, used to get all the file directories in the current folder

i = 0  # Write the list of files to test.xls
for s in pathDir:
    sheet.write(i, 0, s)  # Parameters i, 0, s represent row, column, write value respectively
    i = i + 1

print(file_path)
print(i)  # Show the number of filenames
f.save(file_path)
