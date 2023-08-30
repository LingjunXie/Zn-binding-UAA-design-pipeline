# Lingjun Xie
# 08/28/2023
# Rutgers, the State University of New Jersey
# Khare Lab

# import library

from pymol import cmd
import pandas as pd
import argparse, sys, subprocess, os

# flags

def print_notice(outstring):
    print("lx110::"+str(outstring))

def parse_args():
    info = """
        This script should calculate the distance between AAPF's end moity N atom and lignad's Zn binding atom.
        """
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-l', '--ligand', type=str, required=True,
                        help='Three letter name of substrate.')
    parser.add_argument('-u', '--uaa', type=str, required=True,
                        help='Three letter name of UAA.')
    parser.add_argument('-t', '--target_atoms', type=str, required=True,
                        help='Zn binding atom names in original substrate, separate with comma.')
    parser.add_argument('-o', '--out_name', type=str, required=True,
                        help='Output file name.')

    return parser.parse_args()
    
def generate_distance_list(pdb_file, uaa, ligand, target_atoms):

    cmd.load(pdb_file)
    print_notice(pdb_file + ' loaded')
    
    dis_list = []
    for atom in target_atoms:
        distance = cmd.distance('resn ' + uaa + ' and name NM1', 'resn ' + ligand + ' and name ' + atom)
        dis_list.append(distance)
        print_notice('NM1-' + atom + ' = ' + str(distance))
    for atom in target_atoms:
        distance = cmd.distance('resn ' + uaa + ' and name NM2', 'resn ' + ligand + ' and name ' + atom)
        dis_list.append(distance)
        print_notice('NM1-' + atom + ' = ' + str(distance))
        
    return dis_list
    
    file_list = []
    for file in os.listdir():
        if file.endswith('.pdb'):
            file_list.append(file)

# Main function starts here.

def main(args):

    print(args.target_atoms)
    target_atoms_comma = args.target_atoms
    target_atoms = target_atoms_comma.split(',')
    print(target_atoms)

    # make a list that contains all file names end with .pdb
    
    pdb_files = []
    for file in os.listdir():
        if file.endswith('.pdb'):
            pdb_files.append(file)
            
    # generate column list for dataframe, e.g. NM1-NAN

    column_list = []
    for atom in target_atoms:
        column_name = "NM1-" + atom
        column_list.append(column_name)
    for atom in target_atoms:
        column_name = "NM2-" + atom
        column_list.append(column_name)
        
    print(column_list)
        
    # generate all distance and put them into one dataframe

    df = pd.DataFrame(columns = column_list)
    for pdb_file in pdb_files:
        dis_list = generate_distance_list(pdb_file, args.uaa, args.ligand, target_atoms)
        df.loc[len(df)] = dis_list
        print_notice('All distances in file ' + pdb_file + ' were obtained')
        
    # dump all data into excel file
    
    df.to_excel(args.out_name + '.xls', index = False, header=True)

if __name__ == '__main__':
    args = parse_args()
    main(args)

