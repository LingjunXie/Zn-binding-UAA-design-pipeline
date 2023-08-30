# codes originally from Elliott
# https://github.com/emd182/uaa_azo_proteins/blob/main/uaa_chi_sample.py
# Modified by Lingjun Xie
# 08/25/2023
# Rutgers, the State University of New Jersey
# Khare Lab

#import library

import pyrosetta as pr
from pyrosetta.rosetta.protocols.minimization_packing import \
    PackRotamersMover
from pyrosetta.rosetta.core.pack.rotamer_set import \
    RotamerSetFactory
import argparse, sys, subprocess

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
                        help='PDB code for protein to design, lower case')
    parser.add_argument('-u', '--uaa', type=str, required=True,
                        help='Three letter name of UAA used for design.')
    parser.add_argument('-r', '--res_id', type=int, required=True,
                        help='res id of position to mutate')

    return parser.parse_args()

# Main function starts here.

def main(args):

    # load pose and UAA

    pr.init('-extra_res_fa ../../../UAAs/' + args.uaa + '/' + args.uaa + '.params -ex1 -ex2')
    start_pose = pr.pose_from_pdb('../' + args.input_pdb + '_Zn_neighbor_to_Ala.pdb')
    pose = start_pose.clone()

    # set up and make a mutation

    sf = pr.rosetta.core.scoring.ScoreFunction()
    sf.add_weights_from_file('ref2015')

    res_mut = args.res_id
    mutater = pr.rosetta.protocols.simple_moves.MutateResidue()
    mutater.set_target(res_mut)
    mutater.set_res_name(args.uaa)
    mutater.apply(pose)

    # make Packertask and restrict to repacking

    packer_task = pr.standard_packer_task(pose)
    packer_task.restrict_to_repacking()
    packer_task.set_bump_check(False)
    pack_mover = PackRotamersMover(sf, packer_task)

    rsf = RotamerSetFactory()
    rs = rsf.create_rotamer_set(pose)
    rs.set_resid(res_mut)
    sf(pose)
    packer_graph = pr.rosetta.core.pack.create_packer_graph(pose, sf, packer_task)
    rs.build_rotamers(pose, sf, packer_task, packer_graph)

    print(rs.rotamer(1).name(), rs.num_rotamers())
    short_pose = pr.rosetta.core.pose.Pose()
    short_pose.detached_copy(pose)
    for i in range(1, rs.num_rotamers()+1):
        short_pose.residue(res_mut).set_all_chi(rs.rotamer(i).chi())
        short_pose.dump_pdb('mutated_pos_'+str(i)+'.pdb')#, res_mut_loc)
        
if __name__ == '__main__':
    args = parse_args()
    main(args)
