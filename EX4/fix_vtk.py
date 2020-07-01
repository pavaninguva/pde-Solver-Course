import glob

from tempfile import mkstemp
from shutil import move
from os import remove
import os
# Obtain list of VTK files in directory
vtk_list = glob.glob("*.vtk")

# Rename the extensions to .txt files to process
for vtk in vtk_list: 
    # Access the last four characters in the path and replace them with .txt
    os.rename(vtk, vtk[:-4] + ".txt")

txt_list = glob.glob("*output.*.txt")

def replace(source_file_path, pattern, substring):
    fh, target_file_path = mkstemp()
    with open(target_file_path, 'w') as target_file:
        with open(source_file_path, 'r') as source_file:
            for line in source_file:
                if line.strip() == pattern:
                    lineout = f"{substring}\n"
                    target_file.write(lineout)
                else:
                    target_file.write(line)
    remove(source_file_path)
    move(target_file_path, source_file_path)

for text in txt_list:
    replace(text, "41", "9")
    os.rename(text, text[:-4] + ".vtk")
            
