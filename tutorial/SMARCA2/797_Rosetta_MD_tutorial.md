# Rosetta

###  1. Preparation

Open AutoPROTAC ternaries.

```python
delete SMARCA2-FWZ
delete VHL-4YYA
alter resn 4YYA, resn="4YY"  # Rosetta only accept 3-letter ligand
alter resn HSD, resn="HIS"

set pdb_retain_ids,0
cd C:\Data\OneDrive - Baylor College of Medicine\Hanfeng Modeling PC\SMARCA2 VHL\20250903_797test

python
import os
os.makedirs("autoPROTAC_ter", exist_ok=True)
for obj in cmd.get_object_list():
    cmd.save(f"autoPROTAC_ter/{obj}.pdb", obj)
python end
```

run `get_atom_distance_for_objects.py` in PyMOL with bridging atom id, to understand the **linker distance distribution**:

```
get_atom_distance_for_objects C16,O19

# output between 6.4-7.91
```



Upload all ternaries pdbs to server

Prepare Rosetta params files.

```bash
# Server bash

pymol "$(ls *.pdb | head -n 1)" 

# PyMOL command
h_add chain X
save FWZ.mol2, chain X
h_add chain Y
save 4YY.mol2, chain Y

# bash
molfile_to_params.py FWZ.mol2 -n FWZ --keep-names --centroid
molfile_to_params.py 4YY.mol2 -n 4YY --keep-names --centroid
```

Prepare cst file (according to **linker distance distribution**)

avg(6.4, 7.91)=7.155

7.155-6.4=0.755

sd=0.2 means penalty = ((x-tolerance)/sd)^2

Bridging atom ids are on the FWZ/4YYA ligand (from `PROTAC.xml`)

```bash
echo "AtomPair CAJ 1Y N41 1X FLAT_HARMONIC 7.155 0.2 0.755">cst
```



### 2. Rosetta Prepacking

```bash
# Manual MPI, since SBGrid Rosetta 3.14 doesn't support mpi

readlink -f *ter.pdb > pdblist

CPU_nodes=$(nproc)
N=$(wc -l < pdblist)
CHUNK=$(( (N + CPU_nodes - 1) / CPU_nodes ))
split -l $CHUNK pdblist pdblist.

for pdblist in pdblist.*;
do nohup docking_prepack_protocol.static.linuxgccrelease -in:file:l $pdblist -out:suffix _prepacked -run:constant_seed -docking:partners PX_ABY -docking:sc_min -extra_res_fa FWZ.fa.params 4YY.fa.params -extra_res_cen FWZ.cen.params 4YY.cen.params -ignore_unrecognized_res -load_PDB_components False & 
rm $pdblist 
done


# Alternatively, configure SBGrid to use Rosetta 3.13 mpi
readlink -f *ter.pdb > pdblist
mpirun -np 64 docking_prepack_protocol.mpi.linuxgccrelease -in:file:l pdblist -out:suffix _prepacked -run:constant_seed -docking:partners PX_ABY -docking:sc_min -extra_res_fa FWZ.fa.params 4YY.fa.params -extra_res_cen FWZ.cen.params 4YY.cen.params -ignore_unrecognized_res -load_PDB_components False
```

```bash
#Appendix: configuring SBGrid to use Rosetta 3.13
cat > ~/.sbgrid.conf <<EOF
ROSETTA_X=3.13
EOF
source /programs/sbgrid.shrc
```



### 3. Rosetta Local refine constraint docking

create docking opinion file `flag_docking_local_refine_0.3_0.8_constraint`

```bash
cat > flag_docking_local_refine_0.3_0.8_constraint <<EOF
-nstruct 50

-partners PX_ABY
-dock_pert 0 0
-docking:docking_local_refine
-dock_mcm_trans_magnitude 0.3
-dock_mcm_rot_magnitude 0.8
-docking:sc_min

-extra_res_fa FWZ.fa.params 4YY.fa.params
-extra_res_cen FWZ.cen.params 4YY.cen.params

-ex1
-ex2aro
-ex3
-ex4
-ex3::level 4
-ex4::level 4

-constraints:cst_file cst

-out:suffix _local_dock_cst
-out:path:all ./docking_0.3_0.8_cst
EOF

```

```bash
# Manual MPI, since SBGrid Rosetta 3.14 doesn't support mpi
readlink -f *prepack*.pdb > pdblist
mkdir docking_0.3_0.8_cst
CPU_nodes=$(nproc)
N=$(wc -l < pdblist)
CHUNK=$(( (N + CPU_nodes - 1) / CPU_nodes ))
split -l $CHUNK pdblist pdblist.

for pdblist in pdblist.*;
do nohup docking_protocol.static.linuxgccrelease -in:file:l $pdblist @flag_docking_local_refine_0.3_0.8_constraint & done

# Alternatively, configure SBGrid to use Rosetta 3.13 mpi
readlink -f *prepack*.pdb > pdblist
mkdir docking_0.3_0.8_cst
mpirun -np 64 docking_protocol.mpi.linuxgccrelease -in:file:l pdblist @flag_docking_local_refine_0.3_0.8_constraint

```

#### TODO: Is top100 better or top1% better for next step?

top100:

```bash
# View Top 100 Result
sort -nk 6 score_local_dock_cst.sc | head -n 100 
# Save Top 100 Result
sort -nk 6 score_local_dock_cst.sc | head -n 100 | awk '{print $NF".pdb"}' > top100list.txt

```

```bash
# Cluster top 100 decoys
mkdir top100
while IFS= read -r filename; do cp "$filename" top100; done < top100list.txt

cd top100
cp ../../*.params ./
readlink -f *.pdb > top100list.txt

time energy_based_clustering.static.linuxgccrelease -in:file:l top100list.txt -cluster:energy_based_clustering:cluster_radius 3.5 -extra_res_fa FWZ.fa.params 4YY.fa.params -extra_res_cen FWZ.cen.params 4YY.cen.params   > cluster.log

```

Top1%

```bash
# View Top 1% Result
sort -nk 6 score_local_dock_cst.sc | head -n $(( ( $(wc -l < score_local_dock_cst.sc) - 2 ) / 100 )) 
# Save Top 1% Result
sort -nk 6 score_local_dock_cst.sc | head -n $(( ( $(wc -l < score_local_dock_cst.sc) - 2 ) / 100 )) | awk '{print $NF".pdb"}' > top1percent.txt

# Cluster top 1% decoys
mkdir top1percent
while IFS= read -r filename; do cp "$filename" top1percent; done < top1percent.txt

cd top100
cp ../../*.params ./
readlink -f *.pdb > top1percent.txt

time energy_based_clustering.static.linuxgccrelease -in:file:l top1percent.txt -cluster:energy_based_clustering:cluster_radius 3.5 -extra_res_fa FWZ.fa.params 4YY.fa.params -extra_res_cen FWZ.cen.params 4YY.cen.params   > cluster.log
```

Output: a total of 14 clusters generated (The cluster 9 looks similar to crystal structure).

### 4. run linker_builder.py in the autoPROTAC

```bash
mkdir cluster_rep
cp c.*.1.pdb cluster_rep

# If 4YY was used for rosetta, always rename `4YY` back to `4YYA`, otherwise `linker_builder.py` will use `4YY` xml.
# rename `4YY` back to `4YYA`
for f in cluster_rep/*.pdb; do sed -i 's/4YY /4YYA/g' "$f"; done

```

Load all the c.*.1.pdb and omega linker conformers into pymol

```bash
# PyMOL command
run linker_builder.py
linker_builder

run merge_PROTAC.py
merge_PROTAC

# save results as pdb, one object per file.
```

Only 5 clusters survived.

Use Maestro to assign bond order, add hydrogen, minimize small molecule

save objects as `minimized__c_*_1_with_linker_rmsd_*.pdb`



# GROMACS

## Vanilla MD for each conformation

### 1. File preparation



Linux Bash:

we have 5 clusters from last step, so

#### rename pdb file with suffix number 1-5:



```bash
chmod 754 *

suffix=1
for file in *.pdb
do
    # Get the filename without the extension
    base_name=$(basename "$file" .pdb)

    # Form the new filename, appending the suffix to the original name
    new_file="${base_name}_${suffix}.pdb"

    # Rename the file
    mv "$file" "$new_file"

    # Increment the suffix for the next file
    ((suffix++))
done
```



#### Separate molecule and protein



```bash
( for i in {1..5};
do grep "797 Z" *c_*_$i.pdb > 797_$i.pdb;
grep CONECT *c_*_$i.pdb >> 797_$i.pdb;
sed /"797 Z"/d *c_*_$i.pdb > binary_$i.pdb;
sed -i /CONECT/d binary_$i.pdb
done)

# Make sure always use lowercase in LIG name '797', otherwise cgenff_charmm2gmx_py3_nx1.py will give you trouble
```



#### CGenFF for ligand parameterization

`pymol 797*`

```python
# Save as mol2
# pymol command input
python
for obj in cmd.get_object_list("all"):
    cmd.save(obj + ".mol2", obj)
python end
```

Edit the Resi `797_001` into `797`,

```bash
chmod 754 *
sed -i s/797_[0-9]*/797/g *.mol2
```

upload one of them to CGenFF: https://cgenff.com, get `797.str` based on `CGenFF version 4.6`

Run conversion using `cgenff_charmm2gmx_py3_nx1.py` and `charmm36-jul2021.ff` 

```bash
#cp ~/cgenff_charmm2gmx_py3_nx1.py .
#ln -s ~/charmm36-jul2021.ff ./

#python3 -m pip3 install networkx==1.11
#python3 -m pip3 install numpy

chmod 754 *

for i in {1..5};
do python3 ~/cgenff_charmm2gmx_py3_nx1.py 797 797_$i.mol2 797.str ~/charmm36-jul2021.ff;
	for f in 797.[itp]*; 
	do mv $f "$(echo "$f"|sed s/797/797_$i/)";
	# rename ligand topology files with numbers;
    done;
    mv 797_ini.pdb 797_$i'_ini'.pdb;    
done
```

#### Upload everything to `$i` folder

```bash
for i in {1..5};
do mkdir $i;
mv *_$i.* ./$i;
mv *_$i'_ini.pdb' ./$i
ln -s ~/797ter/charmm36-jul2021.ff ./$i
done
```



### 2. Small molecule modeling script (FYI, can directly go to step 3)

Usage example:`~/protein_ligand_MD_script_100ns.sh 797 1`. This script will be called in step 3. No need to run now.

```bash
#!/bin/bash
# Check if this directory contains: 797_1.mol2/.pdb/.itp/.prm/.top  binary_1.pdb 797_1_ini.pdb
# Check if both small molecule name and $i value are provided as arguments
set -e

if [ -z "$1" ] || [ -z "$2" ]; then
    echo 'Please provide both the small molecule name and the value for $i as arguments.'
    exit 1
fi

# Assign the small molecule name and $i value to variables
molecule_name="$1"
i="$2"

if [ ! -d "charmm36-jul2021.ff" ]; then
        ln -s ../charmm36-jul2021.ff ./
fi
gmx editconf -f "${molecule_name}_${i}_ini.pdb" -o "${molecule_name}_${i}.gro"
sed -i "s/ZN    ZN/ZN   ZN2/g" binary_$i.pdb
cp /data/apps/sbgrid/x86_64-linux/gromacs/2022.1_cu11.5.2/share/gromacs/top/residuetypes.dat ./
sed -i "s/ZN2/ZN/g" charmm36-jul2021.ff/ions.itp
sed -i '$aZN2     Ion\
SOD     Ion\
CLA     Ion' residuetypes.dat
gmx pdb2gmx -f binary_$i.pdb -o binary_$i.gro -ignh -ff charmm36-jul2021 -water tip3p
sleep 5
cp binary_$i.gro ternary_$i.gro
tail -n +3 "${molecule_name}_$i.gro" | head -n -1 > lig1coord
sed -i "/^$(tail -n2 ternary_$i.gro | head -n1)/ r lig1coord" ternary_$i.gro
totalatom=$(expr $(sed -n 2p ternary_$i.gro) + $(sed -n 2p "${molecule_name}_$i.gro"))
sed -i "2c $totalatom" ternary_$i.gro
sed -i 's/; Include chain topologies/; Include ligand parameters\
#include "'"$molecule_name"'_'"$i"'.prm"\
\
; Include chain topologies/g' topol.top
sed -i 's/; Include water topology/; Include ligand topology\
#include "'"$molecule_name"'_'"$i"'.itp"\
\
; Include water topology/g' topol.top
sed -i '$a'"$molecule_name"'             1' topol.top
cp ~/mdp_protein_ligand/*.mdp ./
sed -i "s/mz1/$molecule_name/g" *.mdp
gmx editconf -f ternary_$i.gro -o newbox.gro -bt cubic -d 1.0
gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro
gmx editconf -f solv.gro -o solv.pdb
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr
echo "SOL" | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -conc 0.15
gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
sleep 5
gmx_mpi mdrun -v -s em.tpr -deffnm em -pin on
echo non-Water | gmx trjconv -f em.gro -s em.tpr -o em.pdb
( echo '0&!aH*'; echo q ) | gmx make_ndx -f "${molecule_name}_$i.gro" -o index_"$molecule_name"_$i.ndx
echo 3 | gmx genrestr -f "${molecule_name}_$i.gro" -n index_"$molecule_name"_$i.ndx -o posre_"$molecule_name".itp -fc 1000 1000 1000
sed -i 's/; Include water topology/; Ligand position restraints\
#ifdef POSRES\
#include "posre_'"$molecule_name"'.itp"\
#endif\
\
; Include water topology/g' topol.top
echo q | gmx make_ndx -f em.gro -o index.ndx
gmx select -f em.gro -s em.tpr -on index_protein_"$molecule_name".ndx -select 'resname '"$molecule_name"' or group 1'
sed -i 's/resname_'"$molecule_name"'_or_group_1/Protein_'"$molecule_name"'/g' index_protein_"$molecule_name".ndx
cat index_protein_"$molecule_name".ndx >> index.ndx
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr
sleep 5
export OMP_NUM_THREADS=128
gmx_mpi mdrun -v -s nvt.tpr -deffnm nvt -pin on
sed -i 's/Berendsen/C-rescale/g' npt.mdp
gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr
sleep 5
gmx_mpi mdrun -v -s npt.tpr -deffnm npt
sed -i -e 's/500000 /50000000 /g' -e 's/1000 ps/100,000 ps/g' -e 's/1 ns/100 ns/g' md.mdp
sed -i -e 's/5000 /50000 /g' -e 's/10.0 ps/100 ps/g' md.mdp
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_"$molecule_name".tpr
```



### 3. Prepare MD files (NVT, NPT)

```
cd 1   # This folder should contain 797_1.mol2/.pdb/.itp/.prm/.top  binary_1.pdb 797_1_ini.pdb
~/protein_ligand_MD_script_100ns.sh 797 1
```



### 4. Run production MD

#### Server

```bash
export OMP_NUM_THREADS=8
# 1 job on multiple GPUs
mpirun -n 8 gmx_mpi mdrun -v -s md_797.tpr -deffnm md_797 -ntomp 8 -gpu_id 0123 -pin on -pme gpu -npme 1 -nb gpu -bonded gpu

# multiple jobs on multiple GPUs
export OMP_NUM_THREADS=10
gmx mdrun -v -s md_797.tpr -deffnm md_797 -ntomp 10 -ntmpi 1 -gpu_id 0

#continue simulation
export OMP_NUM_THREADS=10
gmx mdrun -v -s md_797.tpr -deffnm md_797 -ntomp 10 -gpu_id 0 -ntmpi 1 -cpi md_797.cpt
```



#### Post md visualization (optional)



```bash
# Find ligand atom number 4281
less em.gro | grep "797    "


( echo 'a 4281'; echo 'name 25 center'; echo q)|gmx make_ndx -n index.ndx -o center.ndx -f em.gro
~/visualization_center_dt100.sh md_797
```



```bash
# create rms matrix
( echo 'ri 1-381'; echo 'name 24 CRBN'; echo 'ri 382-639'; echo 'name 25 BTK'; echo ''; echo q )| gmx make_ndx -f em.gro -o CRBN_BTK.ndx;
( echo 'CRBN'; echo 'BTK')|gmx rms -f md_797_fit.xtc -s md_797.tpr -n CRBN_BTK.ndx -m rms.xpm -o rmsd.xvg;
# cluster trajectory
gmx cluster -f md_797_fit.xtc -s md_797.tpr -dm rms.xpm -o rmsd.xpm -n index.ndx -cutoff 0.3 -nofit -cl -clid cluster.xvg
# lsq choose protein_797, output choose protein_797
```



```
# Use this script as "visualization_100ns.sh filename"
filename="${1%.*}"
echo "This script use group 'center' in center.ndx to center" 
echo $filename
sleep 2

if [ ! -f center.ndx ];
then echo "ERROR: Please choose one atom as the center using gmx make_ndx, and output as center.ndx";
echo "Exiting..."
exit 1
fi

current_folder=$(basename "$PWD")

( echo system)| gmx trjconv -s "$filename".tpr -f "$filename".xtc -o "$filename"_whole.xtc -pbc whole;
( echo center; sleep 1; echo system)| gmx trjconv -s "$filename".tpr -f "$filename"_whole.xtc -o "$filename"_center.xtc -center -pbc nojump -ur rect -n center.ndx;
( echo backbone; sleep 1; echo system)| gmx trjconv -s "$filename".tpr -f "$filename"_center.xtc -o "$filename"_fit.xtc -fit rot+trans # rot doesn't compatible with -pbc;
( echo non-water)| gmx trjconv -s "$filename".tpr -f "$filename"_fit.xtc -o "$filename"_visual_"$current_folder".pdb -pbc nojump -dt 100
```

### 5. Optional: MMPBSA

##### Concatenate trajectories and RMSD matrix output



for 100ns md for mmpbsa:

```
# remove water
for folder in [1-5]; do
    cd "$folder";
    ( echo '24|13'; echo q )|gmx make_ndx -n index.ndx -o index_nowater.ndx -f em.gro; #choose Protein_797|Zn2;
    echo Protein_797_ZN2|gmx trjconv -s md_797.tpr -f md_797_fit.xtc -o md_797_nowater_fit_`$folder`.xtc -n index_nowater.ndx;
    cd ..
done

gmx trjcat -f [1-5]/md_797_nowater_fit*.xtc -o md_797_nowater_100ns_cat.xtc -cat
```



##### gmx_MMPBSA



```
# At lab server
conda activate gmxMMPBSA
```



Conside 'CRBN' as A, and target protein as B

```
vim mmpbsa.in
```



```
Sample input file for PB calculation
This input file is meant to show only that gmx_MMPBSA works. Althought,
we tried to used the input files as recommended in the Amber manual,
some parameters have been changed to perform more expensive calculations
in a reasonable amount of time. Feel free to change the parameters 
according to what is better for your system.

&general
sys_name="Prot-Lig-CHARMM",
startframe=1,
# 0ns
endframe=5004,
# 100ns
#forcefields="oldff/leaprc.ff99SB,leaprc.gaff"
#forcefields="leaprc.protein.ff14SB,leaprc.gaff2" 
# NOt setting forcefield will make gmx_MMPBSA use the charmm forcefield in the topology
# In gmx_MMPBSA v1.5.0 we have added a new PB radii set named charmm_radii. 
# This radii set should be used only with systems prepared with CHARMM force fields. 
# Uncomment the line below to use charmm_radii set
PBRadii=7,
/
&pb
# radiopt=0 is recommended which means using radii from the prmtop file for both the PB calculation and for the NP
# calculation
istrng=0.15, fillratio=4.0, radiopt=0, inp=1
/
&decomp
idecomp=2, dec_verbose=3
print_res="within 4"
```



```
# Trying to use concatenated traj...
cp 1/em.gro .
cp 1/annealing_1/CRBN_BTK.ndx ./
cp 1/md_797.tpr ./
 gmx make_ndx -f em.gro -n index.ndx -o index.ndx;
# 24 CRBN  as receptor
# 25 BTK  as ligand


conda activate gmxMMPBSA
ln -s ~/charmm36-jul2021.ff/ ./
cp */797_*.prm ./
cp 1/*.itp ./
cp 1/topol.top ./
# make sure MPI is using /home/hlin/.conda/envs/gmxMMPBSA/bin/mpirun, otherwise it has parellel problems
# conda environment gmxMMPBSA has set the PATH 
conda env config vars set PATH=/data/apps/sbgrid/x86_64-linux/xv/3.10a/bin:/data/apps/sbgrid/x86_64-linux/gromacs/2022.1_cu11.5.2/bin:/home/hlin/.conda/envs/gmxMMPBSA/bin:/data/apps/sbgrid/x86_64-linux/pymol/2.5.3_386:/usr/local/cuda-11.8/bin:/usr/bin:/home/hlin/rosetta_bin_linux_2021.16.61629_bundle/main/source/scripts/python/public:/home/hlin/rosetta_bin_linux_2021.16.61629_bundle/main/source:/home/hlin/rosetta_bin_linux_2021.16.61629_bundle/main/source/bin:/programs/x86_64-linux/system/sbgrid_bin:/usr/bin:/home/hlin/.local/bin:/home/hlin/bin:/usr/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/var/lib/snapd/snap/bin:/programs/share/bin:/programs/share/sbgrid/bin:/programs/x86_64-linux/sbgrid_installer/latest



mpirun -np 80 gmx_MMPBSA MPI -O -i mmpbsa.in -cs md_797.tpr -ci CRBN_BTK.ndx -cg 24 25 -ct md_797_nowater_100ns_cat.xtc -cp topol.top -o FINAL_RESULTS_MMPBSA.dat -eo FINAL_RESULTS_MMPBSA.csv -nogui;
# CRBN as receptor, BTK as ligand

gmx_MMPBSA_ana
# if gmx_MMPBSA_ana has BDC keyvalue problem, try to lower the threshold of energy kcal/mol
```





## HAPOD-MD

Continue from 2.

### 3. AnnealingMD with shorter gradient:



**starting from rosetta output pdb** minimized with maestro

 0-10ns 300K

 10-30ns 300-400K

 30-50ns 400K

vim `annealing_md.mdp ` in the `~/` folder 

```
title                   = Protein-ligand complex MD simulation
; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = 25000000   ; 2 fs * 25000000 = 50000 ps (50 ns)
dt                      = 0.002     ; 2 fs
; Output control
nstenergy               = 50000      ; save energies every 100.0 ps
nstlog                  = 50000      ; update log file every 100.0 ps
nstxout-compressed      = 50000      ; save coordinates every 100.0 ps
; Bond parameters
continuation            = yes       ; continuing from NPT
constraint_algorithm    = lincs     ; holonomic constraints
constraints             = h-bonds   ; bonds to H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Neighbor searching and vdW
cutoff-scheme           = Verlet
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 20        ; largely irrelevant with Verlet
rlist                   = 1.2
vdwtype                 = cutoff
vdw-modifier            = force-switch
rvdw-switch             = 1.0
rvdw                    = 1.2       ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
rcoulomb                = 1.2
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling
tcoupl                  = V-rescale                     ; modified Berendsen thermostat
tc-grps                 = Protein_797 Water_and_ions    ; two coupling groups - more accurate
tau_t                   = 0.1   0.1                     ; time constant, in ps
ref_t                   = 300   300                     ; reference temperature, one for each group, in K
; Pressure coupling
pcoupl                  = Parrinello-Rahman             ; pressure coupling is on for NPT
pcoupltype              = isotropic                     ; uniform scaling of box vectors
tau_p                   = 2.0                           ; time constant, in ps
ref_p                   = 1.0                           ; reference pressure, in bar
compressibility         = 4.5e-5                        ; isothermal compressibility of water, bar^-1
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Dispersion correction is not used for proteins with the C36 additive FF
DispCorr                = no
; Velocity generation
gen_vel                 = no        ; continuing from NPT equilibration

; SIMULATED ANNEALING CONTROL
; If you are running protein- RNA/DNA/metal ions/ligands than you have to use what you have given and do little change in annealing protocol i.e.,
annealing = single single
annealing-npoints = 4 4
annealing-time = 0 10000 30000 50000 0 10000 30000 50000 ; ps
annealing-temp = 300 300 400 400 300 300 400 400
```



#### Annealing MD preparation script



make a series of directories like this using the code snippet:

```
1
|-797_1.mol2  797_1.pdb  binary_1.pdb  797_1.itp  797_1.prm  797_1.top  797_1_ini.pdb
2
|-797_2.mol2  797_2.pdb  binary_2.pdb  797_2.itp  797_2.prm  797_2.top  797_2_ini.pdb
......
5
|-797_5.mol2  797_5.pdb  binary_5.pdb  797_5.itp  797_5.prm  797_5.top  797_5_ini.pdb
```

```bash
mkdir annealing_md
for i in {1..5};
do mkdir annealing_md/$i;
mv $i/797_$i.mol2 ./annealing_md/$i;
mv $i/797_$i.pdb ./annealing_md/$i;
mv $i/binary_$i.pdb ./annealing_md/$i;
mv $i/797_$i.itp ./annealing_md/$i;
mv $i/797_$i.prm ./annealing_md/$i;
mv $i/797_$i.top ./annealing_md/$i;
mv $i/797_$i'_ini.pdb' ./annealing_md/$i;
ln -s ~/charmm36-jul2021.ff ./annealing_md/$i
done

cd annealing_md
ls -R

# In this annealing_md folder, you should have folder 1-5, prepare_annealing.sh, annealing_md.mdp
./prepare_annealing.sh 797 1 2 3 4 5 #each annealing conformation run k=4 replicates
```



./prepare_annealing.sh :

```bash
#!/bin/bash
set -e

if [ -z "$1" ]; then
    echo 'Please provide the small molecule name for $name as arguments.'
    exit 1
fi
name="$1"

if [ -z "$2" ]; then
    echo 'Please provide at least one value for $i as arguments.'
    exit 1
fi

# Shift arguments to get rid of $name
shift

# Iterate over each $i value provided
for i in "$@"; do
    # Ensure the provided value for $i is a valid directory
    if [ ! -d "$i" ]; then
        echo "Directory '$i' does not exist."
        continue
    fi

    # Create 4 replicate directories dynamically
    for ((k=1; k<=4; k++)); do
        mkdir -p "$i/annealing_$k"
    done

    # Iterate over annealing folders
    for annealing_folder in "$i"/annealing_*; do
        (
        cd "$annealing_folder"
        
        # Copy required files
        cp "../binary_$i".* ./
        cp "../${name}_${i}"* ./
        cp "../../annealing_md.mdp" md.mdp
        
        # Execute script
        ~/protein_ligand_MD_script_100ns.sh "$name" "$i"
        )
    done
done
```



### 4. Production MD Run

#### TACC

```bash
# Lonestar 6
module load gromacs/2022.1
idev -p gpu-a100 -N 1 -n 6
export
ibrun gmx_mpi_gpu mdrun -v -s md_797.tpr -deffnm md_797 -ntomp 6 -gpu_id 012


#!/bin/bash
#SBATCH -J md           # Job name
#SBATCH -o md.o%j       # Name of stdout output file
#SBATCH -e md.e%j       # Name of stderr error file
#SBATCH -p gpu-a100          # Queue (partition) name
#SBATCH -N 1              # Total # of nodes (must be 1 for serial)
#SBATCH -n 6               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 12:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=hanfeng.lin@bcm.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A CHE23001        # Allocation name (req'd if you have more than 1)

# Other commands must follow all #SBATCH directives...

module list
pwd
date

# Launch serial code...
source ~/.bash_profile
module load intel/19.1.1 
module load impi/19.0.9
module load cuda/11.4
module load gromacs/2022.1
export OMP_NUM_THREADS=6

ibrun gmx_mpi_gpu mdrun -v -s md_797.tpr -deffnm md_797 -ntomp 6 -gpu_id 012
```



#### Server

```bash
export OMP_NUM_THREADS=8
# 1 job on multiple GPUs
mpirun -n 8 gmx_mpi mdrun -v -s md_797.tpr -deffnm md_797 -ntomp 8 -gpu_id 0123 -pin on -pme gpu -npme 1 -nb gpu -bonded gpu

# multiple jobs on multiple GPUs (Recommended)
export OMP_NUM_THREADS=10
gmx mdrun -v -s md_797.tpr -deffnm md_797 -ntomp 10 -ntmpi 1 -gpu_id 0

#continue simulation
export OMP_NUM_THREADS=10
gmx mdrun -v -s md_797.tpr -deffnm md_797 -ntomp 10 -gpu_id 0 -ntmpi 1 -cpi md_797.cpt
```

```bash
# loop over all the replicates in one cluster

export OMP_NUM_THREADS=10
# k=4 replicates
for i in {1..4};
do cd annealing_$i;
gmx mdrun -v -s md_797.tpr -deffnm md_797 -ntomp 10 -ntmpi 1 -gpu_id 3; # change -gpu_id
cd ..;
done


export OMP_NUM_THREADS=10
# k=4 replicates
for i in {1..4};
do cd annealing_$i;
gmx mdrun -v -s md_p10.tpr -deffnm md_p10 -ntomp 10 -ntmpi 1 -gpu_id 0; # change -gpu_id
cd ..;
done

```



#### post md visualization & rmsd



```bash


#visualization
( echo 'a 4281'; echo 'name 25 center'; echo q)|gmx make_ndx -n index.ndx -o center.ndx -f em.gro
 ~/visualization_center_dt100.sh md_797

# RMSD calculation for target protein when align for E3 
# ref = 0ns
for annealing_folder in */annealing_*; do
cd "$annealing_folder";
( echo 'a 4281'; echo 'name 25 center'; echo q)|gmx make_ndx -n index.ndx -o center.ndx -f em.gro
 ~/visualization_center_dt100.sh md_797;
( echo 'chain A'; echo 'chain B'; echo q )| gmx make_ndx -f em.pdb -o E3_target.ndx;
( echo 'chA'; echo 'chB')|gmx rms -f md_797_fit.xtc -s md_797.tpr -n E3_target.ndx -m rms.xpm -o rmsd.xvg;
cd ../..
done

paste */annealing_*/rmsd.xvg > merged_rmsd.xvg


# RMSD with ref=10ns frame
for annealing_folder in */annealing_*; do
cd "$annealing_folder";
echo > blank.ndx
( echo 'chain A'; echo 'chain B'; echo 'chain C'; echo '0|1|2'; echo q )| gmx make_ndx -f em.pdb -n blank.ndx -o E3_target.ndx;
( echo 'chA_chB_chC')|gmx trjconv -s md_797.tpr -f md_797_fit.xtc -n E3_target.ndx -o ${annealing_folder//\//_}_ref_10ns.pdb -dump 10000
( echo 'chA'; echo 'chB')|gmx rms -f md_797_fit.xtc -s ${annealing_folder//\//_}_ref_10ns.pdb -n E3_target.ndx -m rms.xpm -o rmsd_10ns.xvg;
cd ../..
done
#xv rms.xpm
#grace rmsd.xvg

paste */annealing_*/rmsd_5ns.xvg > merged_rmsd_5ns.xvg

eval paste $(for f in */annealing_*/rmsd_5ns.xvg; do
    echo "<(awk '{\$1=\"\"; sub(/^ /,\"\"); print}' $f)"
done) > merged_rmsd_5ns.xvg


```





### Concatenate trajectories and RMSD matrix output



for annealing MD:

```bash
# Remove water from fit trajectory to make atom number consistent
for annealing_folder in */annealing_*; do
    cd "$annealing_folder";
    echo Protein_797|gmx trjconv -s md_797.tpr -f md_797_fit.xtc -o md_797_nowater_fit_`$annealing_folder`.xtc -n index.ndx;
    cd ../..
done

gmx trjcat -f [1-5]/anneal*/md_797_nowater_fit*.xtc -o md_797_nowater_anneal_cat.xtc -cat

#First do fitting for CRBN, then cluster for BTK
cp 1/annealing_1/E3_target.ndx ./
cp 1/annealing_1/md_797.tpr ./
( echo 'chA'; echo 'chA_chB_chC') | gmx trjconv -s md_797.tpr -f md_797_nowater_anneal_cat.xtc -o md_797_nowater_anneal_cat_fit.xtc -n E3_target.ndx -fit rot+trans
( echo 'chA_chB_chC') |gmx trjconv -s md_797.tpr -f md_797_nowater_anneal_cat_fit.xtc -n E3_target.ndx -o md_797_nowater_anneal_cat_fit.pdb -dt 1000
# choose chA for lsq, then chA_chB_chC for output
# gmx cluster -f md_797_cat_fit.xtc -s md_797.tpr -o rmsd.xpm -n CRBN_BTK.ndx -cutoff 0.3 -fit no -cl
# choose BTK for RMSD (-fit no), Protein_797 as output


# should be the same matrix as the following rms calculation
( echo 'chA'; echo 'chB') | gmx rms -f md_797_nowater_anneal_cat_fit.xtc -s md_797.tpr -n E3_target.ndx -m rms.xpm -dt 1000 #ps
least square fit: VHL
RMSD: SMARCA2

gmx xpm2ps -f rms.xpm -o rmsd.eps -rainbow blue -size 800
xv rmsd.eps
```





