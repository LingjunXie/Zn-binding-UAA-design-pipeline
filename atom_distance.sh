#!/bin/bash
while (( $# > 1 ))
do
    case $1 in
        -pdb) pdb_id="$2";;
        -uaa) uaa_name="$2";;
        -ligand) ligand_name="$2";;
        -target) target_atoms="$2";;
        *) break;
    esac
    shift 2
done

cat ../${pdb_id}.pos | while read line
do
    cd ${line}_${uaa_name}
    echo lx110::Now entering ${line}_${uaa_name}
    python /Users/lingjunxie/Desktop/Lab/UAA_project/Computational/Zn_binding/scripts/atom_distance.py -l ${ligand_name} -t ${target_atoms}  -o ${line}_distance -u ${uaa_name}
    cd ..
    sleep 0.05
done
