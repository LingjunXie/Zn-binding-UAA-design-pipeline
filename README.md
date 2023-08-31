# Zn binding UAA design pipeline
## folder structure ##
    /scripts
        atom_distance.py
        find_closest_distance.py
        organize_excel_files.py
        uaa_chi_sample.py
        Zn_neighbor_to_Ala.py
        atom_distance.sh
        organize_folders.sh
        uaa_chi_sample.sh
        
    /protins
        /$protein
            ${protein}_relaxed.pdb
            ${protein}_Zn_neighbor_to_Ala.pdb // genarated by Zn_neighbor_to_Ala.py
            ${protein}_no_lig.pdb // remove lig from protein_Zn_neighbor_to_Ala.pdb
            ${protein}.pos // genarated by Zn_neighbor_to_Ala.py
            /$pos
            
    /UAAs
        /UAA
            UAA.params
            UAA.rotlib
            
### Make rotamer library for UAA ###

**make mol2 file**

Draw unnatural amino acid with Acetylated, N-methylated capping groups on backbone
- Optimize geometry with Avogadro
- Name every atom a different name 
- Save as mol2 file

**Add instruction lines at bottom of mol2 file**

```rust
@
M  ROOT \[the number of the N-terminal atom\]  
M  POLY_N_BB \[ditto, unless for some reason the most N terminal atom in your residue type is not N\]  
M  POLY_CA_BB \[the CA atom number\]  
M  POLY_C_BB \[the C atom number\]  
M  POLY_O_BB \[the O atom number\]  
M  POLY_IGNORE \[all the atoms of the capping groups except UPPER and LOWER\]   (10 atoms)
M  POLY_UPPER \[the nitrogen atom number of the C-terminal methyl amide\]  
M  POLY_LOWER \[the carbonyl atom number of the N-terminal acetyl\]  
M  POLY_CHG \[the charge\]  
M  POLY_PROPERTIES \[any properties, like PROTEIN and CHARGED, etc.\]  
M  END
```

**Convert mol2 file to params file**

put mol2 file under /home/lx110/mol_file
```rust
cd /home/lx110/mol_file
/home/lx110/Rosetta/main/source/scripts/python/public/molfile_to_params_polymer.py -i APC.mol2 --name APC --clobber --polymer
```

**make input file**

make iput file in the format showed below

```rust
AA_NAME APT
OMG_RANGE 180 180 1
PHI_RANGE 0 350 10
PSI_RANGE 0 30 10
EPS_RANGE 180 180 1
NUM_CHI 4
NUM_BB 2
CHI_RANGE 1 0  330  30
CHI_RANGE 2 0  300  60
CHI_RANGE 3 0  330  30
CHI_RANGE 4 0  330  30
ROTWELLS 1 3  60 180 300
ROTWELLS 2 2  90 270
ROTWELLS 3 2  0 180
ROTWELLS 4 2  0 180
TEMPERATURE 1
```
OMG, PHI,PSI, EPS are pretty standard

**Notice** PSI_RANGE should also be 0 350 10, however, to prevent the job being killed for using too much memory, divide it into couple pieces, in this case, we'll have 0 30 10, 40 70 10, 80 110 10... all the way until 320 350 10, so we would have 9 input files in total

NUM_CHI should align with number of chi in params file

CHI_RANGE and ROTWELLS for chi 1 and chi 2 are pretty standard for Phe based UAA

if any chi are not independent to others, for example, chi 4 can only be certain value if chi 3 is set to a certain value, use CENTROID instead of ROTWELLS, for example:

```rust
AA_NAME APC
OMG_RANGE 180 180 1
PHI_RANGE 0 350 10
PSI_RANGE 0 30 10
EPS_RANGE 180 180 1
NUM_CHI 4
NUM_BB 2
CHI_RANGE 1 0  330  30
CHI_RANGE 2 0  300  60
CHI_RANGE 3 0  330  30
CHI_RANGE 4 0  330  30
CENTROID  60 1  90 1    46.8 1  -133.2 1
CENTROID  60 1  90 1  -133.2 2    46.8 2
CENTROID  60 1 270 2    46.8 1  -133.2 1
CENTROID  60 1 270 2  -133.2 2    46.8 2
CENTROID 180 2  90 1    46.8 1  -133.2 1
CENTROID 180 2  90 1  -133.2 2    46.8 2
CENTROID 180 2 270 2    46.8 1  -133.2 1
CENTROID 180 2 270 2  -133.2 2    46.8 2
CENTROID 300 3  90 1    46.8 1  -133.2 1
CENTROID 300 3  90 1  -133.2 2    46.8 2
CENTROID 300 3 270 2    46.8 1  -133.2 1
CENTROID 300 3 270 2  -133.2 2    46.8 2
TEMPERATURE 1
```

**Run MakeRotLib protocol**

```rust
cd /home/lx110/Rosetta/main/demos/public/make_rot_lib/inputs
mkdir APC
cd APC
```
put all input file in this folder

```rust
cd ../../outputs
mkdir APC
cd APC
for i in {1..9}; do slurmit.py  --job APC_${i} --command "/home/lx110/Rosetta/main/source/bin/MakeRotLib.default.linuxgccrelease -extra_res_fa ../../inputs/APC/APC.params -options_file ../../inputs/APC/APC_${i}.in -make_rot_lib:output_logging false"; done
```

**Assemble libraries into one file**
```rust
../../scripts/make_final_from_phi_psi.pl APC
```
looking for APC.rotlib in output/APC folder, copy it into PATH/TO/Zn_binding_UAA_pipeline/UAAs/APC

**Modifying params file**

add those lines at the very bottom of original params file:

```rust
NCAA_ROTLIB_PATH ABSOLUTE_PATH_TO ROTLIB/APC.rotlib
NCAA_ROTLIB_NUM_ROTAMER_BINS (# chi angles) (# rotamer number for each angle)
```

for example:

```rust
NCAA_ROTLIB_PATH /Users/lingjunxie/Desktop/Lab/UAA_project/Computational/Zn_binding/UAAs/APC/rotlib/APC.rotlib
NCAA_ROTLIB_NUM_ROTAMER_BINS 4 3 2 2 2
```


**Mutate Zn neighbors to Ala without touching catalytic important residues and genarate position file**
conda activate pr
> cd proteins/4eyl
> python ../../scripts/Zn_neighbor_to_Ala.py -p 4eyl_relaxed.pdb -c 120 122 124 189 208 250 211 220 93 67

**sample UAA chi**
bash ../../scripts/uaa_chi_sample.sh -pdb 4eyl -uaa APT
bash ../../scripts/organize_folders.sh -uaa APT

**calculate distance between AAPF's end N and Zn binding atoms on original ligand**
cd APT
conda activate pymol
bash ../../../scripts/atom_distance.sh -pdb 4eyl -uaa APT -ligand 0RV -target NAN,OAE,OAI,OAF

**extract excel files from each folder**
mkdir excel_files
python ../../../scripts/orgainze_excel_files.py

**extract rows with closest distance**
cd excel_files
python ../../../scripts/find_closest_distance.py
