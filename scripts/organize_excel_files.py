# Lingjun Xie
# 08/29/2023
# Rutgers, the State University of New Jersey
# Khare Lab

# import library

import os
import shutil

# Define the source and destination directories
source_dir = './'  # Replace with your source folder path
dest_dir = './excel_files'  # Replace with your destination folder path

# Create destination directory if it doesn't exist
if not os.path.exists(dest_dir):
    os.makedirs(dest_dir)

# Walk through the source directory, including subdirectories
for dirpath, dirnames, filenames in os.walk(source_dir):
    for filename in filenames:
        if filename.endswith('.xlsx') or filename.endswith('.xls'):
            # Construct the full path of the file
            full_file_path = os.path.join(dirpath, filename)
            
            # Construct the destination path
            dest_file_path = os.path.join(dest_dir, filename)
            
            # Check if a file with the same name already exists in the destination folder
            counter = 1
            while os.path.exists(dest_file_path):
                dest_file_path = os.path.join(dest_dir, f"{filename.split('.')[0]}_{counter}.xlsx")
                counter += 1
            
            # Move the file
            shutil.move(full_file_path, dest_file_path)
            print(f"Moved {full_file_path} to {dest_file_path}")
