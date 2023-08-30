# Lingjun Xie
# 08/22/2023
# Rutgers, the State University of New Jersey
# Khare Lab

# import library

from pyrosetta import *
init()
import argparse, sys, subprocess
from pyrosetta.rosetta.core.select import get_residues_from_subset

# flags

def print_notice(outstring):
    print("lx110::"+str(outstring))

def parse_args():
    info = """
        This script should take a protein file and a few other parameters
        and output a single file where residues around Zn mutated to Ala.
        """
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-p', '--input_pdb', type=str, required=True,
                        help='Input starting PDB file of Zn metalloenzyme.')
    parser.add_argument('-r', '--radius', type=float, required=False, default=8.0,
                        help='Radius around Zn to identify residues to mutate.')
    parser.add_argument('-o', '--out_file', type=str, required=False,
                        help='Name of the output file.')
    parser.add_argument('-c', '--catalytic_residues', nargs='+', required=True,
                        help='PDB ID of catalytic residues, space between PDB IDs')

    return parser.parse_args()

# Main function starts here.

def main(args):

    # load pose

    start_pose = pose_from_pdb(args.input_pdb)
    pose = start_pose.clone()
    scorefxn = pyrosetta.create_score_function("ref2015_cart.wts")
    print_notice("Scaffold protein Loaded Successfully!")
    print_notice("Scaffold protein has"+str(pose.total_residue())+"residues.")

    # select Zn neighboring residues

    Zn_res = pyrosetta.rosetta.core.select.residue_selector.ResidueNameSelector("ZN")

    # select residues within defined range from Zn

    Zn_neighbor = pyrosetta.rosetta.core.select.residue_selector.CloseContactResidueSelector()
    Zn_neighbor.central_residue_group_selector(Zn_res)
    Zn_neighbor.threshold(args.radius)

    # selet Not_Zn_neighbor

    Not_Zn_neighbor = pyrosetta.rosetta.core.select.residue_selector.NotResidueSelector(Zn_neighbor)
    
    # select catalytic residue
    
    catalytic_pdb = args.catalytic_residues
    catalytic_res = pyrosetta.rosetta.core.select.residue_selector.OrResidueSelector()
    for i in catalytic_pdb:
        resi = pose.pdb_info().pdb2pose("A", int(i))
        select_res = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(str(resi))
        catalytic_res.add_residue_selector(select_res)
        
    # selet residues not to design
    
    no_design_res = pyrosetta.rosetta.core.select.residue_selector.OrResidueSelector()
    no_design_res.add_residue_selector(Not_Zn_neighbor)
    no_design_res.add_residue_selector(catalytic_res)
    
    # select residues to design
    
    design_res = pyrosetta.rosetta.core.select.residue_selector.NotResidueSelector()
    design_res.set_residue_selector(no_design_res)
    
    # put residues to design in a list
    
    res_design = str(pyrosetta.rosetta.core.select.get_residues_from_subset(design_res.apply(pose)))
    res_design = res_design.replace('vector1_unsigned_long[' , '')
    res_design = res_design.replace(']' , '')
    res_list = res_design.split(', ')
    
    # output pos file contains positions to design
    
    out_name = args.input_pdb.strip('.pdb') + '.pos'
    f = open(out_name, "w")
    for x in res_list:
        f.write(x+'\n')
    f.close

    # set up task factory

    tf = pyrosetta.rosetta.core.pack.task.TaskFactory()

    # these are pretty standard
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())

    # disable design on not_Zn_neighbor
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
        pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(), no_design_res))

    # enable design on Zn_neighbor to Ala
    aa_to_design = pyrosetta.rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
    aa_to_design.aas_to_keep("A")
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
        aa_to_design, Zn_neighbor))

    # convert the task factory into a PackerTask
    packer_task = tf.create_task_and_apply_taskoperations(pose)

    # set up MoveMap

    mm = pyrosetta.rosetta.core.kinematics.MoveMap()
    mm.set_bb(False)
    mm.set_chi(False)
    mm.set_jump(False)

    # set up FastDesign

    rel_design = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn_in=scorefxn, standard_repeats=1, script_file="MonomerDesign2019")
    rel_design.cartesian(True)
    rel_design.set_task_factory(tf)
    rel_design.set_movemap(mm)
    rel_design.minimize_bond_angles(True)
    rel_design.minimize_bond_lengths(True)

    # run FastDesign

    rel_design.apply(pose)

    # dump designed pose to pdb

    if args.out_file:
        out_name = args.out_file
    else:
        out_name = args.input_pdb.strip('_relaxed.pdb') + '_Zn_neighbor_to_Ala.pdb'

    pose.dump_pdb(out_name)

    print_notice("PDB file was dumped to the folder")

if __name__ == '__main__':
    args = parse_args()
    main(args)

