# Lingjun Xie
# 08/29/2023
# Rutgers, the State University of New Jersey
# Khare Lab

# import library

import pandas as pd
import os

# List all Excel files in the current directory
excel_files = [f for f in os.listdir() if f.endswith('.xlsx') or f.endswith('.xls')]

# Initialize an empty list to store the rows with the smallest numbers
min_rows = []

for file in excel_files:
    # Read the Excel file into a DataFrame
    df = pd.read_excel(file)
    
    # Find the row with the smallest number
    # Assuming the number can be in any column, we'll first stack the DataFrame to get all numbers
    # Then, we'll find the index of the smallest number
    min_idx = df.stack().idxmin()
    
    # Get the row with the smallest number
    min_row = df.loc[min_idx[0]].to_frame().transpose()
    
    # Add a new column to store the file name
    min_row['Source_File'] = file
    
    # Reorder the columns to make 'Source_File' the first column
    cols = ['Source_File'] + [col for col in min_row.columns if col != 'Source_File']
    min_row = min_row[cols]
    
    # Append the row to our list
    min_rows.append(min_row)

# Combine the rows into a new DataFrame
result_df = pd.concat(min_rows, ignore_index=True)

# Write the combined rows to a new Excel file
result_df.to_excel('minimum_distance_for_each_rotamer.xlsx', index=False)
