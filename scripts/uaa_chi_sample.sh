#!/bin/bash
while (( $# > 1 ))
do
    case $1 in
        -pdb) pdb_id="$2";;
        -uaa) uaa_name="$2";;
        *) break;
    esac
    shift 2
done

cat ${pdb_id}.pos | while read line
do
    mkdir ${line}_${uaa_name}
    cd ${line}_${uaa_name}
    python /Users/lingjunxie/Desktop/Lab/UAA_project/Computational/Zn_binding/scripts/uaa_chi_sample.py -p $pdb_id -u $uaa_name -r $line
    cd ..
    sleep 0.05
done
