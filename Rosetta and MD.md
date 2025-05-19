# Rosetta

p10 rosetta packing and local docking

C:\Users\Hanfeng\OneDrive - Baylor College of Medicine\Hanfeng Modeling PC\BTK\PS10 modeling\BTK_fixed\PS10_e10_7conf_merged.pse

**Did not found PS10 cryoEM like conformation for e10, e11**

**e12 gave 3 cryoEM like conformations**

Start from Before merging, Rename ternary chain X from PS2 into P2W  (PS10 has the same warhead as PS2)

~~~
alter chain X, resn="P2W"
remove resn P2W and (name N1 or name N2 or name C7 or name C8 or name C9 or name C10 or name C11)
~~~

Remove P2W piperazine ring...

Generate Binary backbone complex, `alter` HSD to HIS. dump backbone and ligands from autoProtac (don’t retain atom ids) one object per file

Upload to server, create pdb list

```
readlink -f *.pdb > pdblist
```



#### Generate param file

**When generating mol2 files in PyMOL, always add new hydrogens to no-hydrogen structure. CHOOSE the "retain atom id" option** in pymol export (Can only work for **3-letter** name)

Export .mol2 files from autoPROTAC ternary output -> `Y70.mol2`, `P2W.mol2`, upload to server, run script:

```
molfile_to_params.py Y70.mol2 -n Y70 --keep-names --centroid
molfile_to_params.py P2W.mol2 -n P2W --keep-names --centroid
```



#### Docking Prepack Protocol

https://www.rosettacommons.org/docs/latest/application_documentation/docking/docking-prepack-protocol

~~~bash
docking_prepack_protocol.static.linuxgccrelease -in:file:l pdblist -out:suffix _prepacked -run:constant_seed -docking:partners PX_ABY -docking:sc_min -extra_res_fa P2W.fa.params Y70.fa.params -extra_res_cen P2W.cen.params Y70.cen.params -ignore_unrecognized_res -load_PDB_components False

~~~

This will also generate hydrogen for ligands.

#### Docking_local_refine

this part can add distance constaint between two bridging atoms. Allows for deeper docking sampling



##### High Resolution Docking Only Flags (no constraint)

```flag_docking_local_refine
-nstruct 50

-partners PX_ABY
-dock_pert 0 0
-docking:docking_local_refine
-dock_mcm_trans_magnitude 0
-dock_mcm_rot_magnitude 0	
-docking:sc_min


-extra_res_fa P2W.fa.params Y70.fa.params
-extra_res_cen P2W.cen.params Y70.cen.params 

-ex1
-ex2aro
-ex3
-ex4
-ex3::level 4
-ex4::level 4

-out:suffix _local_dock
-out:path:all ./docking_0_0
```

```bash
readlink -f *prepack*.pdb > pdblist
mkdir docking_0_0

#nohup docking_protocol.static.linuxgccrelease -in:file:l pdblist @flag_docking_local_refine &
```

* Without `-docking:sc_min` it will give some very large moving. 

```bash
wc -l pdblist  
split -l 1 pdblist pdblist.

for pdblist in pdblist.*;
do nohup docking_protocol.static.linuxgccrelease -in:file:l $pdblist @flag_docking_local_refine & done
```

```bash

sort -nk 6 score_local_dock.sc | head -n 100 | awk '{print $28".pdb"}' > top100list.txt

mkdir top100
while IFS= read -r filename; do cp "$filename" top100; done < top100list.txt



```

found 2 cryoEM like structure 

2.6_P10_e12_7126ter_prepacked_0001_local_dock_0002.pdb 2.6_P10_e12_7126ter_prepacked_0001_local_dock_0036.pdb



##### Try constraint file for bridging atoms (Indeed perform better than no constraint)

cst file:

```cst
AtomPair C6 1Y C5 1X FLAT_HARMONIC 7.46 0.2 1.46
```

No penalty for **6.0-8.92A (according to oeomega bridging atom distance distribution)**

7.46 (center x0) +-1.46 (tolerance)

sd=0.5 means penalty = ((x-tolerance)/sd)^2

**flag_docking_local_refine_0.3_0.8_constraint :**

```bash
-nstruct 50

-partners PX_ABY
-dock_pert 0 0
-docking:docking_local_refine
-dock_mcm_trans_magnitude 0.3
-dock_mcm_rot_magnitude 0.8
-docking:sc_min

-extra_res_fa P2W.fa.params Y70.fa.params
-extra_res_cen P2W.cen.params Y70.cen.params

-ex1
-ex2aro
-ex3
-ex4
-ex3::level 4
-ex4::level 4

-constraints:cst_file cst

-out:suffix _local_dock_cst
-out:path:all ./docking_0.3_0.8_cst                    
```

```bash

for pdblist in pdblist.*;
do nohup docking_protocol.static.linuxgccrelease -in:file:l $pdblist @flag_docking_local_refine_0.3_0.8_constraint & done

sort -nk 6 score_local_dock_cst.sc | head -n 100 | awk '{print $28".pdb"}' > top100list.txt
```



#### Cluster top 100 decoys

~~~bash
sort -nk 6 score_local_dock.sc | head -n 100 | awk '{print $29".pdb"}' > top100list.txt

mkdir top100
while IFS= read -r filename; do cp "$filename" top100; done < top100list.txt

cd top100

cp ../../*.params ./
readlink -f *.pdb > top100list.txt

cluster.static.linuxgccrelease  -in:file:l top100list.txt -in:file:fullatom -cluster:radius 3 -extra_res_fa P2W.fa.params Y70.fa.params -extra_res_cen P2W.cen.params Y70.cen.params -cluster:sort_groups_by_energy  > cluster.log 

# recommend, could not use MPI, ~5 min
time energy_based_clustering.static.linuxgccrelease -in:file:l top100list.txt -cluster:energy_based_clustering:cluster_radius 3.5 -extra_res_fa P2W.fa.params Y70.fa.params -extra_res_cen P2W.cen.params Y70.cen.params  > cluster.log
~~~

got 6 clusters

#### Linker_builder.py

keep c.*.1 (because they are the lowest energy cluster center), import in PyMOL, also import PS10_e12 oeomega conformers

```
alter resn P2W, resn="PS2"

linker_builder c.4.1

#copy linker to c.1.23

merge_protac

```

* c.6.1 has no linker with RMS<0.4. Try 0.8. Get 0.501



Use maestro to assign bond order, add hydrogen, and minimize small molecule



# Annealing MD for each conformation

### 1. File preparation

Linux Bash:

we have 6 clusters from last step, so

#### rename pdb file with suffix number 1-6:

~~~bash
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
~~~



#### Separate molecule and protein

```bash
( for i in {1..6};
do grep P10 *c_*_$i.pdb > p10_$i.pdb;
grep CONECT *c_*_$i.pdb >> p10_$i.pdb;
sed /P10/d *c_*_$i.pdb > binary_$i.pdb;
sed -i /CONECT/d binary_$i.pdb
done)

# Make sure always use lowercase 'p10', otherwise cgenff_charmm2gmx_py3_nx1.py will give you trouble
```



#### CGenFF for ligand parameterization

In pymol, Add hydrogen for `P10_$i.pdb`, do `rename *`, export by **`p10_$i.mol2`** format.

Edit the Resi `P10_001` into `P10`,

~~~
chmod 754 *
sed -i s/p10_[0-9]*/p10/g *.mol2
~~~

 upload to CGenFF, get `p10.str` based on `CGenFF version 4.6`

No need to copy `cgenff_charmm2gmx_py3_nx1.py` and `charmm36-jul2021.ff` from other folder to here

~~~bash
#cp ~/cgenff_charmm2gmx_py3_nx1.py .
#ln -s ~/charmm36-jul2021.ff ./

#python3 -m pip3 install networkx==1.11
#python3 -m pip3 install numpy

chmod 754 *

for i in {1..6};
do python3 ~/cgenff_charmm2gmx_py3_nx1.py p10 p10_$i.mol2 p10.str ~/charmm36-jul2021.ff;
	for f in p10.[itp]*; 
	do mv $f "$(echo "$f"|sed s/p10/p10_$i/)";
	# rename ligand topology files with numbers;
    done;
    mv p10_ini.pdb p10_$i'_ini'.pdb;    
done
~~~

## 

#### Upload everything to `~/p10ter/$i` folder

```bash
mkdir ~/p10ter
```

```bash
cd ~/p10ter
for i in {1..6};
do mkdir $i;
mv *_$i.* ./$i;
mv *_$i'_ini.pdb' ./$i
done
```

```bash
 # forcefield
 cd p10ter # Forcefield in this root folder
 for i in {1..6};
 do ln -s ~/p10ter/charmm36-jul2021.ff ./$i
 done
 
 
```



#### 

### 2. Small molecule modeling script!! 

Usage example:`~/protein_ligand_MD_script_100ns.sh p10 1`. This script will be called in step 3. No need to run now.

~~~bash
#!/bin/bash
# Check if this directory contains: p10_1.mol2/.pdb/.itp/.prm/.top  binary_1.pdb p10_1_ini.pdb
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

~~~



### 3. AnnealingMD with shorter gradient: 

**starting from rosetta output pdb** minimized with maestro

​	0-10ns 300K

​	10-30ns 300-400K

​	30-50ns 400K

vim `annealing_md.mdp ` in the root folder `~/ps10ter/minimized`

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
tc-grps                 = Protein_p10 Water_and_ions    ; two coupling groups - more accurate
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

./prepare_annealing.sh p10 1 2 3 4 5 6 (each annealing conformation run k=4 replicates)

make directory contains:

`p10_1.mol2  p10_1.pdb  binary_1.pdb  p10_1.itp  p10_1.prm  p10_1.top  p10_1_ini.pdb`

...

`p10_6.mol2  p10_6.pdb  binary_6.pdb  p10_6.itp  p10_6.prm  p10_6.top  p10_6_ini.pdb`



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

#### 

### 4. Production MD Run



#### TACC

~~~bash
# Lonestar 6
module load gromacs/2022.1
idev -p gpu-a100 -N 1 -n 6
export
ibrun gmx_mpi_gpu mdrun -v -s md_p10.tpr -deffnm md_p10 -ntomp 6 -gpu_id 012


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

ibrun gmx_mpi_gpu mdrun -v -s md_p10.tpr -deffnm md_p10 -ntomp 6 -gpu_id 012

~~~

#### Server

```bash
export OMP_NUM_THREADS=8
# 1 job on multiple GPUs
mpirun -n 8 gmx_mpi mdrun -v -s md_p10.tpr -deffnm md_p10 -ntomp 8 -gpu_id 0123 -pin on -pme gpu -npme 1 -nb gpu -bonded gpu

# multiple jobs on multiple GPUs
export OMP_NUM_THREADS=10
gmx mdrun -v -s md_p10.tpr -deffnm md_p10 -ntomp 10 -ntmpi 1 -gpu_id 0

#continue simulation
export OMP_NUM_THREADS=10
gmx mdrun -v -s md_p10.tpr -deffnm md_p10 -ntomp 10 -gpu_id 0 -ntmpi 1 -cpi md_p10.cpt
```



### Post md visualization

```bash
# Find ligand atom number
less em.gro | grep "P10    "


( echo 'a 10350'; echo 'name 25 center'; echo q)|gmx make_ndx -n index.ndx -o center.ndx -f em.gro
~/visualization_center_dt100.sh md_p10

```



~~~bash
# create rms matrix
( echo 'ri 1-381'; echo 'name 24 CRBN'; echo 'ri 382-639'; echo 'name 25 BTK'; echo ''; echo q )| gmx make_ndx -f em.gro -o CRBN_BTK.ndx;
( echo 'CRBN'; echo 'BTK')|gmx rms -f md_p10_fit.xtc -s md_p10.tpr -n CRBN_BTK.ndx -m rms.xpm -o rmsd.xvg;
# cluster trajectory
gmx cluster -f md_p10_fit.xtc -s md_p10.tpr -dm rms.xpm -o rmsd.xpm -n index.ndx -cutoff 0.3 -nofit -cl -clid cluster.xvg
# lsq choose protein_p10, output choose protein_p10
~~~



~~~bash
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


~~~

### Hydrogen bond analysis

~~~bash
gmx hbond -f md_p10_fit.xtc -s md_p10.tpr -n index.ndx -num -hbn -hbm
# Choose CRBN and BTK

bash ~/gmx_hbdat.bsh -s md_p10.tpr -n hbond.ndx -m hbmap.xpm
~~~

`gmx_hbdat.bsh`

~~~bash
echo -e "\
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   gmx_hbdat    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
>>>>>>>>>>>>>>>>>>>>>>>>>>     2023-06-05 00:54:15     <<<<<<<<<<<<<<<<<<<<<<<<<\n
>>   Usage: gmx_hbdat -s     *.tpr  -n     *.ndx  -m     *.xpm
>> Default: gmx_hbdat -s topol.tpr  -n hbond.ndx  -m hbmap.xpm

--------------------------------------------------------------------------------
>> Log:
   2022-06-25: add #res for atoms
   2023-06-05: revise output
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"

tpr=topol.tpr
ndx=hbond.ndx
xpm=hbmap.xpm

gmx='gmx'        # /path/to/GMX/bin/gmx_mpi
dump="$gmx dump" # gmx dump

opt=($*); N=${#opt[@]}
for((i=0; i<N; i++)); do
	arg=${opt[$i]}; j=$((i+1)); val=${opt[$j]}
	[[ $arg =~ -s   ]] && { tpr=$val; }
	[[ $arg =~ -n   ]] && { ndx=$val; }
	[[ $arg =~ -m   ]] && { xpm=$val; }
done

$dump -s $tpr -quiet 2>>/dev/null | awk -v ndx=$ndx -v xpm=$xpm '
	BEGIN {
		isHB=0; nhb=0
		while(getline < ndx ) {
			if(index($0, "hbonds_")) isHB=1
			if(isHB && !index($1,"[")) {
				nhb++; don[nhb]=$1; hyd[nhb]=$2; acc[nhb]=$3
			}
		}
		close(ndx)

		isHB=0; nhb=0
		while(getline < xpm) {
			if(index($0, "y-axis")) isHB=1
			if(isHB && index($0, "\"")) {
				nhb++
				n=0; gsub(/[,\"]/, "")
				for(i=1; i<=length($0); i++) if(substr($0, i, 1)=="o") n++
				occ[nhb]=n*100/length($0)
			}
		}
		close(xpm)
	}

	/#molblock/  { Ntyp=$3 }
	/moltype.+=/ { Imol=$3; getline; Nmol[Imol]=$3 }
	/moltype.+\(/ { Imol=$0; gsub(/[^0-9]/,"",Imol)
		getline txt; sub(/.*=/,"",txt); gsub(" ","_",txt)
		Name[Imol]=txt
		getline; getline txt;       gsub(/[^0-9]/,"",txt); Natm[Imol]=txt+0
		for(i=0; i<Natm[Imol]; i++) {
			getline; txt=$0; idx=$3; resID[Imol, i]=$(NF-2)+1+nres
		}
		getline
		for(i=0; i<Natm[Imol]; i++) {
			getline txt
			sub(/.+=./, "", txt); sub(/..$/, "", txt)
			Tatm[Imol, i]=txt
		}
	}

	/residue\[/ { nres++
		sub(/.*="/,"",$0); sub(/".*/,"",$0);
		resName[nres]=sprintf("%4d%s", nres, $0)
	}

	END {

	n=split("ALA A CYS C ASP D GLU E PHE F GLY G HIS H ILE I LYS K LEU L MET M ASN N PRO P GLN Q ARG R SER S THR T VAL V TRP W TYR Y UNK X", arr)
	for(i=1; i<=n/2; i++) sname[arr[2*i-1]]=arr[2*i]

		Ntot=0; maxlen=0
		for(i=0; i<Ntyp; i++) {
			if(length(Name[i])>maxlen) maxlen=length(Name[i])
			for(n=0; n<Nmol[i]; n++) {
				for(j=0; j<Natm[i]; j++) {
					Ntot++
					Label[Ntot]=Ntot" "Name[i]" "resName[resID[i,j]]" "Tatm[i, j]
				}
			}
		}

		print "#HBond #atom                 mol.         #res Donor-H     mol.         #res Acceptor  Occupancy%   inter?  Label"
		fmt="%4d   %-20s %-"maxlen"s %8s %-9s  %-"maxlen"s %8s %-9s %9.3f    %4d    %-s\n"
		for(i=1; i<=nhb; i++) {
			split(Label[don[i]]" "Label[hyd[i]]" "Label[acc[i]], arr)

			n=arr[3]; sub(/[^0-9]+/, "", n)
			s=arr[3]; sub(/[0-9]+/,  "", s)
			split(arr[2], donor_chain_arr, "_")
			donor_chain_id = donor_chain_arr[3]
			gsub("\"", "", donor_chain_id)
			split(arr[10], acceptor_chain_arr, "_")
			acceptor_chain_id = acceptor_chain_arr[3]
			gsub("\"", "", acceptor_chain_id)

			tag=donor_chain_id":"sname[s]""n"@"arr[4]"-"arr[8]
			n=arr[11]; sub(/[^0-9]+/, "", n)
			s=arr[11]; sub(/[0-9]+/,  "", s)
			tag=tag"..."acceptor_chain_id":"sname[s]""n"@"arr[12]

			printf fmt, i,
				arr[1]"-"arr[5]"..."arr[9],
				arr[2],  arr[3], arr[4]"-"arr[8],
				arr[10], arr[11], arr[12], occ[nhb-i+1], 
				arr[3]==arr[11]? 0:1, tag
		}
	}
' > hbdat.dat
~~~

### Concatenate trajectories and RMSD matrix output

for annealing MD:

~~~bash
# Remove water from fit trajectory to make atom number consistent
for annealing_folder in */annealing_*; do
    cd "$annealing_folder";
    ( echo '24|13'; echo q )|gmx make_ndx -n index.ndx -o index_nowater.ndx -f em.gro; #choose Protein_p10|Zn2;
    echo Protein_p10_ZN2|gmx trjconv -s md_p10.tpr -f md_p10_fit.xtc -o md_p10_nowater_fit_`$annealing_folder`.xtc -n index_nowater.ndx;
    cd ../..
done

gmx trjcat -f [1-5]/anneal*/md_p10_nowater_fit*.xtc -o md_p10_nowater_anneal_cat.xtc -cat

#First do fitting for CRBN, then cluster for BTK
cp 1/annealing_1/CRBN_BTK.ndx ./
cp 1/annealing_1/md_p10.tpr ./
gmx trjconv -s md_p10.tpr -f md_p10_nowater_anneal_cat.xtc -o md_p10_nowater_anneal_cat_fit.xtc -n CRBN_BTK.ndx -fit rot+trans
gmx trjconv -s md_p10.tpr -f md_p10_nowater_anneal_cat_fit.xtc -n CRBN_BTK.ndx -o md_p10_nowater_anneal_cat_fit.pdb -dt 1000
# choose CRBN for lsq, then Protein|p10 for output
# gmx cluster -f md_p10_cat_fit.xtc -s md_p10.tpr -o rmsd.xpm -n CRBN_BTK.ndx -cutoff 0.3 -fit no -cl
# choose BTK for RMSD (-fit no), Protein_p10_Zn2 as output


# should be the same matrix as the following rms calculation
gmx rms -f md_p10_nowater_anneal_cat_fit.xtc -s md_p10.tpr -n CRBN_BTK.ndx -m rms.xpm -dt 1000 #ps
least square fit: CRBN
RMSD: BTK

gmx xpm2ps -f rms.xpm -o rmsd.eps -rainbow blue -size 800
xv rmsd.eps
~~~







# Vanilla MD for each conformation

### 3. Prepare MD files (NVT, NPT)

`cd ~/p10ter/1   # This folder should contain p10_1.mol2/.pdb/.itp/.prm/.top  binary_1.pdb p10_1_ini.pdb` 

`~/protein_ligand_MD_script_100ns.sh p10 1`



### 4. Run production MD

### Concatenate trajectories and RMSD matrix output

for 100ns md for mmpbsa:

```bash
# remove water
for folder in [1-5]; do
    cd "$folder";
    ( echo '24|13'; echo q )|gmx make_ndx -n index.ndx -o index_nowater.ndx -f em.gro; #choose Protein_p10|Zn2;
    echo Protein_p10_ZN2|gmx trjconv -s md_p10.tpr -f md_p10_fit.xtc -o md_p10_nowater_fit_`$folder`.xtc -n index_nowater.ndx;
    cd ..
done

gmx trjcat -f [1-5]/md_p10_nowater_fit*.xtc -o md_p10_nowater_100ns_cat.xtc -cat

```



### gmx_MMPBSA

~~~bash
# At lab server
conda activate gmxMMPBSA
~~~

Conside 'CRBN' as A, and target protein as B

```
vim mmpbsa.in
```

~~~python
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
~~~



~~~bash
# Trying to use concatenated traj...
cp 1/em.gro .
cp 1/annealing_1/CRBN_BTK.ndx ./
cp 1/md_p10.tpr ./
 gmx make_ndx -f em.gro -n index.ndx -o index.ndx;
# 24 CRBN  as receptor
# 25 BTK  as ligand


conda activate gmxMMPBSA
ln -s ~/charmm36-jul2021.ff/ ./
cp */p10_*.prm ./
cp 1/*.itp ./
cp 1/topol.top ./
# make sure MPI is using /home/hlin/.conda/envs/gmxMMPBSA/bin/mpirun, otherwise it has parellel problems
# conda environment gmxMMPBSA has set the PATH 
conda env config vars set PATH=/data/apps/sbgrid/x86_64-linux/xv/3.10a/bin:/data/apps/sbgrid/x86_64-linux/gromacs/2022.1_cu11.5.2/bin:/home/hlin/.conda/envs/gmxMMPBSA/bin:/data/apps/sbgrid/x86_64-linux/pymol/2.5.3_386:/usr/local/cuda-11.8/bin:/usr/bin:/home/hlin/rosetta_bin_linux_2021.16.61629_bundle/main/source/scripts/python/public:/home/hlin/rosetta_bin_linux_2021.16.61629_bundle/main/source:/home/hlin/rosetta_bin_linux_2021.16.61629_bundle/main/source/bin:/programs/x86_64-linux/system/sbgrid_bin:/usr/bin:/home/hlin/.local/bin:/home/hlin/bin:/usr/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/var/lib/snapd/snap/bin:/programs/share/bin:/programs/share/sbgrid/bin:/programs/x86_64-linux/sbgrid_installer/latest



mpirun -np 80 gmx_MMPBSA MPI -O -i mmpbsa.in -cs md_p10.tpr -ci CRBN_BTK.ndx -cg 24 25 -ct md_p10_nowater_100ns_cat.xtc -cp topol.top -o FINAL_RESULTS_MMPBSA.dat -eo FINAL_RESULTS_MMPBSA.csv -nogui;
# CRBN as receptor, BTK as ligand

gmx_MMPBSA_ana
# if gmx_MMPBSA_ana has BDC keyvalue problem, try to lower the threshold of energy kcal/mol
~~~

### #RMSD

Compare the RMSD for ternary and binary trajectories (reference: 10ns frame)

~~~bash
# Binary

( echo 'ri 1-381'; echo 'name 20 CRBN'; echo 'ri 382-625'; echo 'name 21 BTK'; echo ''; echo q )| gmx make_ndx -f em.gro -o index.ndx
# ri 1-381 
# name 20 CRBN
# ri 382-625
# name 21 BTK
# Choose BTK and CRBN

(echo '0')|gmx trjconv -f md_p10_binary_fit.xtc -s md_p10_binary.tpr -b 10000 -e 10000 -o md_p10_binary_fit_10ns.gro
# choose system
(echo 'CRBN'; echo 'BTK')|gmx rms -f md_p10_binary_fit.xtc -s md_p10_binary_fit_10ns.gro -n index.ndx -m rms_binary.xpm -o rms_binary.xvg -b 10000 -e 100000 -nomw
# LSQ: CRBN  RMSD: BTK
gmx xpm2ps -f rms_binary.xpm -o rms_binary.eps -rainbow blue

~~~

~~~bash
# Ternary

( echo 'ri 1-381'; echo 'name 22 CRBN'; echo 'ri 382-625'; echo 'name 23 BTK'; echo ''; echo q )| gmx make_ndx -f em.gro -o index.ndx
# Choose BTK and CRBN

(echo '0')|gmx trjconv -f md_p10_fit.xtc -s md_p10.tpr -b 10000 -e 10000 -o md_p10_fit_10ns.gro
# choose system
(echo 'CRBN'; echo 'BTK')|gmx rms -f md_p10_fit.xtc -s md_p10_fit_10ns.gro -n index.ndx -m rms_ternary.xpm -o rms_ternary.xvg -b 10000 -e 100000 -nomw
# LSQ: CRBN  RMSD: BTK
gmx xpm2ps -f rms_ternary.xpm -o rms_ternary.eps -rainbow blue

~~~

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

#### #post md visualization & rmsd

~~~bash
export OMP_NUM_THREADS=12
gmx mdrun -v -s md_p10.tpr -deffnm md_p10 -ntomp 12 -ntmpi 1 -gpu_id 3


#visualization
( echo 'a 10350'; echo 'name 25 center'; echo q)|gmx make_ndx -n index.ndx -o center.ndx -f em.gro
 ~/visualization_center_dt100.sh md_p10

# rmsd calculation for BTK when align for CRBN 
for annealing_folder in */annealing_*; do
cd "$annealing_folder";
( echo 'ri 1-381'; echo 'name 24 CRBN'; echo 'ri 382-639'; echo 'name 25 BTK'; echo ''; echo q )| gmx make_ndx -f em.gro -o CRBN_BTK.ndx;
( echo 'CRBN'; echo 'BTK')|gmx rms -f md_p10_fit.xtc -s md_p10.tpr -n CRBN_BTK.ndx -m rms.xpm -o rmsd.xvg;
cd ../..
done

paste */anneal_*/rmsd.xvg > merged_rmsd.xvg

~~~



# Update: Random try: flexible PPI docking using Rosetta 4.0 Protocol (Bullshit)

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6184633/

**Ensemble Generation** 

Start from BTK-PS2.pdb(PX) and CRBN-Y70.pdb(ABY) from prepacked AutoProtac conformation. (ligands have hydrogen) 



- Relax 

Rosetta FastRelax protocol (42) is a refinement algorithm relying on iterations of side-chain packing and energy minimization in torsion space (ϕ, ψ and χi). Five cycles of refinement are carried out while ramping the repulsive part of the van der Waals score term. The command line was: 

```bash
mpirun -n 30 relax.mpi.linuxgccrelease -in:file:s CRBN_Y70.pdb -out:suffix _relax -nstruct 30 -relax:thorough  -extra_res_fa P2W.fa.params Y70.fa.params -extra_res_cen P2W.cen.params Y70.cen.params -ignore_unrecognized_res

```



- Normal Mode Analysis 

Normal Mode Analysis is implemented as the NormalModeRelaxMover in Rosetta. We accessed this through an XML interface called RosettaScripts (56). This protocol mixes motion along the first 5 normal modes, with perturbation of 1 Å. This is iterated with the relax protocol described previously. To prevent non-physical bond angles and bond lengths, we added a term to the score function to penalize deviations from ideal bond angles and lengths. The command line was: 

```bash
mpirun -n 40 rosetta_scripts.mpi.linuxgccrelease -in:file:s CRBN_Y70.pdb -out:suffix _nma -nstruct 40 -parser:protocol nma.xml -extra_res_fa P2W.fa.params Y70.fa.params -extra_res_cen P2W.cen.params Y70.cen.params -ignore_unrecognized_res
```

where nma.xml is:

```xml
<ROSETTASCRIPTS>
 <SCOREFXNS>
 <ScoreFunction name="bn15_cart" weights="beta_nov15_cart" />
 </SCOREFXNS>
 <RESIDUE_SELECTORS>
 </RESIDUE_SELECTORS>
 <TASKOPERATIONS>
 </TASKOPERATIONS>
 <FILTERS>
 </FILTERS>
 <MOVERS>
 <NormalModeRelax name="nma" cartesian="true" centroid="false"
scorefxn="bn15_cart" nmodes="5" mix_modes="true" pertscale="1.0"
randomselect="false" relaxmode="relax" nsample="20"
cartesian_minimize="false" />
 </MOVERS>
 <APPLY_TO_POSE>
 </APPLY_TO_POSE>
 <PROTOCOLS>
 <Add mover="nma" />
 </PROTOCOLS>
 <OUTPUT scorefxn="bn15_cart" />
</ROSETTASCRIPTS>
```

- Backrub 

Rosetta Backrub protocol (43) rotates segments of the protein backbone about an axis defined by the starting and ending atoms of the segment. This is followed by side chain packing. The command line was: 

```bash
mpirun -n 30 backrub.mpi.linuxgccrelease -in:file:s CRBN_Y70.pdb -out:suffix _backrub -nstruct 30 -backrub:ntrials 20000 -backrub:mc_kt 0.6 -extra_res_fa P2W.fa.params Y70.fa.params -extra_res_cen P2W.cen.params Y70.cen.params -ignore_unrecognized_res

rm *last.pdb
rm *low.pdb
```

**After generating ensemble for CRBN, need to remove all the Zn atom, because the following ensemble prepack doesn't support Zn.**

```bash
for file in CRBN_Y70_*.pdb; do
  if [ -f "$file" ]; then
    sed -i '/ZN/d' "$file"
  fi
done
```

Generate ensemble files:

```bash
ls | grep -e "BTK_PS2_.*.pdb$" > ensemble_BTK 

ls | grep -e "CRBN_Y70_.*.pdb$" > ensemble_CRBN
```



Prior to docking simulations, the side chains of all backbone conformers, including the unbound state were optimized in isolation using the following command line: 

**modify 2.6_P10_e12_7126ter.pdb chain sequence as A,Y,P,X**

```bash
# in pymol, load all pdb files
alter chain Y, chain="C"
remove chain B
export molecule (don't tick Original atom order)
# then in bash
sed -i s/"Y70 C"/"Y70 Y"/g test.pdb
```



```bash
# This command doesn't generate ligands
docking_prepack_protocol.static.linuxgccrelease -in:file:s 2.6_P10_e12_7126ter_prepacked_0001.pdb -out:suffix _prepacked -nstruct 1 -ensemble1 ensemble2 -ensemble2 ensemble1 -run:constant_seed -docking:partners AB_P -docking:sc_min -ignore_unrecognized_res -load_PDB_components False

# modify 2.6_P10_e12_7126ter.pdb chain sequence as A,B,Y,P,X using vscode

docking_prepack_protocol.static.linuxgccrelease -in:file:s 2.6_P10_e12_7126ter.pdb -out:suffix _prepacked_test -nstruct 1 -ensemble1 ensemble_CRBN -ensemble2 ensemble_BTK -run:constant_seed -docking:partners ABY_PX -docking:sc_min -extra_res_fa P2W.fa.params Y70.fa.params -extra_res_cen P2W.cen.params Y70.cen.params -ignore_unrecognized_res -load_PDB_components False
# Failed. Prepack always exit in error
# Tried deleting small molecules from pdb, chain sequence as B,A,P., not working

#Error:
protocols.moves.RigidBodyMover: (0) Translate: Jump (before): RT 0.0919571 -0.995137 -0.0353039 -0.158495 -0.0496297 0.986112 -0.983068 -0.0850844 -0.162289 964.489 279.885 281.969  
protocols.moves.RigidBodyMover: (0) Translate: Jump (after):  RT 0.0919571 -0.995137 -0.0353039 -0.158495 -0.0496297 0.986112 -0.983068 -0.0850844 -0.162289 40.2189 8.73937 13.2604  
core.conformation.Residue: (0) Building atom  H   based on standard internal coordinates.
core.util.switchresiduetypeset: (0) trying to preserve existing coords for non-protein residue: 382 ZN
core.conformation.Residue: (0) Building atom  H   based on standard internal coordinates.
core.conformation.Interface: (0) Calculating protein-protein interface
core.conformation.Interface: (0) Calculating protein-protein interface


AN INTERNAL ERROR HAS OCCURED. PLEASE SEE THE CONTENTS OF ROSETTA_CRASH.log FOR DETAILS.


terminate called after throwing an instance of 'std::out_of_range'
  what():  map::at
Got some signal... It is:6
Signal 6 (SIGABRT) means that the process was aborted.  This usually means an internal Rosetta error caused by (often) bad inputs, (sometimes) developer error, or (rarely) hardware problems.

# Try ensemble without Zn, and ternary without zinc, ternary chain order: AYPX
# success!!! no, in some case it might fail, when ligand is at interface, the script cannot find CA atom in the ligand

docking_prepack_protocol.mpi.linuxgccrelease -in:file:s 2.6_P10_e12_7126ter.pdb -out:suffix _prepacked_test -nstruct 1 -ensemble1 ensemble_CRBN -ensemble2 ensemble_BTK  -docking:partners AY_PX -ignore_unrecognized_res -load_PDB_components False -extra_res_fa P2W.fa.params Y70.fa.params -extra_res_cen P2W.cen.params Y70.cen.params

# Have an atom named CA in the ligands before autoPROTAC!
sed s/"C1 "/"CA "/g 2.6_P10_e12_7126ter.pdb > 2.6_P10_e12_7126ter_CA.pdb
sed s/"C1 "/"CA "/g P2W.cen.params > P2W_CA.cen.params
sed s/"C1 "/"CA "/g P2W.fa.params > P2W_CA.fa.params
sed s/"C1 "/"CA "/g Y70.cen.params > Y70_CA.cen.params
sed s/"C1 "/"CA "/g Y70.fa.params > Y70_CA.fa.params
for file in CRBN_Y70_*.pdb; do
  if [ -f "$file" ]; then
    sed -i s/"C1 "/"CA "/g "$file"
  fi
done
for file in BTK_PS2_*.pdb; do
  if [ -f "$file" ]; then
    sed -i s/"C1 "/"CA "/g "$file"
  fi
done
docking_prepack_protocol.mpi.linuxgccrelease -in:file:s 2.6_P10_e12_7126ter_CA.pdb -out:suffix _prepacked_test -nstruct 1 -ensemble1 ensemble_CRBN -ensemble2 ensemble_BTK  -docking:partners AY_PX -ignore_unrecognized_res -load_PDB_components False -extra_res_fa P2W_CA.fa.params Y70_CA.fa.params -extra_res_cen P2W_CA.cen.params Y70_CA.cen.params -out:level 300 -in:detect_disulf false -ex1 -ex2aro
# Worked!
```

Using this prepacked structure, we then performed the docking simulations using the command line: 

##### Ensemble Docking Test

```bash
mkdir flexdocking_0.3_0.8_cst

docking_protocol.mpi.linuxgccrelease -in:file:s 2.6_P10_e12_7126ter_CA_prepacked_test_0001.pdb -docking_centroid_outer_cycles 2 -docking_centroid_inner_cycles 10 -ensemble1 ensemble_CRBN -ensemble2 ensemble_BTK @flag_docking_local_refine_0.3_0.8_cst

#cycles must be called before ensemble since the cycle number is overwritten if in ensemble mode

echo 'AtomPair C6 1Y C5 1X FLAT_HARMONIC 7.46 0.2 1.46' > cst
```



```flag_docking_local_refine_0.3_0.8_cst
-nstruct 1000

-partners AY_PX
-dock_pert 0 0
-docklowres_trans_magnitude 0
-docklowres_rot_magnitude 0
#-docking:docking_local_refine #ensemble_docking must undergo lowres docking
-dock_mcm_trans_magnitude 0.3
-dock_mcm_rot_magnitude 0.8
-docking:sc_min

-extra_res_fa P2W_CA.fa.params Y70_CA.fa.params
-extra_res_cen P2W_CA.cen.params Y70_CA.cen.params

-ex1
-ex2aro
-ex3
-ex4
#-ex3::level 4
#-ex4::level 4

-constraints:cst_file cst

-out:suffix _local_flexdock_cst
-out:path:all ./flexdocking_0.3_0.8_cst
```

### Complete workflow

```bash
mkdir EnsembleDocking
cp *ter.pdb EnsembleDocking
cp BTK_PS2.pdb EnsembleDocking/ #from prepacked AutoProtac conformation. (ligands have hydrogen) 
cp CRBN_Y70.pdb EnsembleDocking/ #from prepacked AutoProtac conformation. (ligands have hydrogen and Zn) 
cp *.params EnsembleDocking

# generate ensembles as described above
...

# remove Zn from all the ensemble
for file in CRBN_Y70_*.pdb; do
  if [ -f "$file" ]; then
    sed -i '/ZN/d' "$file"
  fi
done

# modify autoPROTAC pdb files, including remove Zn, reorder chain order, convert C1 to CA atom
pymol *ter.pdb
# pymol command line:
remove name ZN
alter chain Y, chain="C"
# export molecule (don't tick Original atom order, choose one file per object)
# then in bash
sed -i s/"Y70 C"/"Y70 Y"/g *ter.pdb 
# Order should look like AYPXZ (chain Z won't get read by rosetta)
# convert C1 to CA
sed s/"C1 "/"CA "/g P2W.cen.params > P2W_CA.cen.params
sed s/"C1 "/"CA "/g P2W.fa.params > P2W_CA.fa.params
sed s/"C1 "/"CA "/g Y70.cen.params > Y70_CA.cen.params
sed s/"C1 "/"CA "/g Y70.fa.params > Y70_CA.fa.params
for file in CRBN_Y70_*.pdb; do
  if [ -f "$file" ]; then
    sed -i s/"C1 "/"CA "/g "$file"
  fi
done
for file in BTK_PS2_*.pdb; do
  if [ -f "$file" ]; then
    sed -i s/"C1 "/"CA "/g "$file"
  fi
done
# Create ensemble list
ls | grep -e "BTK_PS2_.*.pdb$" > ensemble_BTK  
ls | grep -e "CRBN_Y70_.*.pdb$" > ensemble_CRBN

# Prepare all the files needed for prepack, then prepack
for model in `ls|grep ter.pdb`; do
  mkdir ${model%.pdb}
  cp $model ${model%.pdb}/
  sed s/"C1 "/"CA "/g ${model%.pdb}/$model > ${model%.pdb}/${model%.pdb}_CA.pdb 
  cp CRBN_Y70_*.pdb ${model%.pdb}
  cp BTK_PS2_*.pdb ${model%.pdb}
  cp *.params ${model%.pdb}
  cp ensemble_* ${model%.pdb}
  cd ${model%.pdb}
  nohup docking_prepack_protocol.mpi.linuxgccrelease -in:file:s ${model%.pdb}_CA.pdb -out:suffix _prepacked -nstruct 1 -ensemble1 ensemble_CRBN -ensemble2 ensemble_BTK -docking:partners AY_PX -ignore_unrecognized_res -load_PDB_components False -extra_res_fa P2W_CA.fa.params Y70_CA.fa.params -extra_res_cen P2W_CA.cen.params Y70_CA.cen.params -out:level 300 -in:detect_disulf false -ex1 -ex2aro &
  cd ..
done

# prepare docking files
mkdir flexdocking_0.3_0.8_cst # for output

echo 'AtomPair C6 1Y C5 1X FLAT_HARMONIC 7.46 0.2 1.46' > cst
echo "-nstruct 1000

-partners AY_PX
-dock_pert 0 0
-docklowres_trans_magnitude 0
-docklowres_rot_magnitude 0
#-docking:docking_local_refine #ensemble_docking must undergo lowres docking
-dock_mcm_trans_magnitude 0.3
-dock_mcm_rot_magnitude 0.8
-docking:sc_min

-extra_res_fa P2W_CA.fa.params Y70_CA.fa.params
-extra_res_cen P2W_CA.cen.params Y70_CA.cen.params

-ex1
-ex2aro
-ex3
-ex4
#-ex3::level 4
#-ex4::level 4

-constraints:cst_file cst

-out:suffix _local_flexdock_cst
-out:path:all ../flexdocking_0.3_0.8_cst" >flag_docking_local_refine_0.3_0.8_cst


# Docking Time!!
for model in `ls|grep ter.pdb`; do
  cd ${model%.pdb}
  cp ../cst .
  cp ../flag_docking_local_refine_0.3_0.8_cst .
  nohup mpirun -n 4 docking_protocol.mpi.linuxgccrelease -in:file:s ${model%.pdb}_CA_prepacked_0001.pdb -docking_centroid_outer_cycles 2 -docking_centroid_inner_cycles 10 -ensemble1 ensemble_CRBN -ensemble2 ensemble_BTK @flag_docking_local_refine_0.3_0.8_cst > docking_nohup & # cycle problem not solved yet
  cd ..
done
```

### Flexible Docking result

```bash
sort -nk 6 score_local_flexdock_cst.sc | head -n 100 | awk '{print $41".pdb"}' > top100list.txt

mkdir top100
while IFS= read -r filename; do cp "$filename" top100; done < top100list.txt

```

Not as great as expected. CRBN loop conformation still wrong...



# CSK modeling

3d7t chain A: CST, align to BTK_PS2 docking model, minimize by maestro. Protein preparation. Build activation loop by prime. Mutate 361-362 back to KK.

### AutoPROTAC (ternary reduction ON):

use P10_e12 conformer

```
autoProtac PS2_CSK_alignfromBTK, CRBN-Poma-8D81
```

**e12 gave 24confs, 2 BTK cryoEM like conformations**

Start from Before merging, Rename ternary chain X from PS2 into P2W  (PS10 has the same warhead as PS2)

~~~
alter chain X, resn="P2W"
remove resn P2W and (name N1 or name N2 or name C7 or name C8 or name C9 or name C10 or name C11)
~~~

Remove P2W piperazine ring...

Generate Binary backbone complex, `alter` HSD to HIS. dump backbone and ligands from autoProtac (don’t retain atom ids) one object per file

### Rosetta prepack

Upload to server, create pdb list

```
readlink -f *.pdb > pdblist
```

get P2W and Y70 params

```
cp ~/ps10rosetta/*.params .
```

https://www.rosettacommons.org/docs/latest/application_documentation/docking/docking-prepack-protocol

~~~bash
docking_prepack_protocol.static.linuxgccrelease -in:file:l pdblist -out:suffix _prepacked -run:constant_seed -docking:partners PX_ABY -docking:sc_min -extra_res_fa P2W.fa.params Y70.fa.params -extra_res_cen P2W.cen.params Y70.cen.params -ignore_unrecognized_res -load_PDB_components False

~~~

This will also generate hydrogen for ligands.

### constraint local refine

```cst
AtomPair C6 1Y C5 1X FLAT_HARMONIC 7.46 0.2 1.46
```

No penalty for 6.0-8.92A (according to oeomega bridging atom distance distribution)

7.46 (center x0) +-1.46 (tolerance)

sd=0.5 means penalty = ((x-tolerance)/sd)^2

```bash
-constrains cst
-dock_mcm_trans_magnitude 0.3
-dock_mcm_rot_magnitude 0.8

readlink -f *prepacked*.pdb > pdblist

mkdir docking_0.3_0.8_cst

mpirun -n 60 docking_protocol.mpi.linuxgccrelease -in:file:l pdblist @flag_docking_local_refine_0.3_0.8_cst

sort -nk 6 score_local_dock_cst.sc | head -n 100 | awk '{print $29".pdb"}' > top100list.txt

mkdir top100
while IFS= read -r filename; do cp "$filename" top100; done < top100list.txt

cd top100

cp ../../*.params ./
readlink -f *.pdb > top100list.txt


# recommend, could not use MPI, ~5 min
time energy_based_clustering.static.linuxgccrelease -in:file:l top100list.txt -cluster:energy_based_clustering:cluster_radius 3.5 -extra_res_fa P2W.fa.params Y70.fa.params -extra_res_cen P2W.cen.params Y70.cen.params  > cluster.log
```

```
-nstruct 50

-partners PX_ABY
-dock_pert 0 0
-docking:docking_local_refine
-dock_mcm_trans_magnitude 0.3
-dock_mcm_rot_magnitude 0.8
-docking:sc_min

-extra_res_fa P2W.fa.params Y70.fa.params
-extra_res_cen P2W.cen.params Y70.cen.params

-ex1
-ex2aro
-ex3
-ex4
-ex3::level 4
-ex4::level 4

-constraints:cst_file cst

-out:suffix _local_dock_cst
-out:path:all ./docking_0.3_0.8_cst
```

#### Linker_builder.py

keep c.*.1 (because they are the lowest energy cluster center), import in PyMOL, also import PS10_e12 oeomega conformers

```
alter resn P2W, resn="PS2"

linker_builder c.1.1

#copy linker to c.1.1

merge_protac

```

CSK cluster 3 is the closest to the BTK cryoEM, but still significantly different

### MD 100ns

#### rename pdb file with suffix number 1-4:

~~~bash
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
~~~



#### Separate molecule and protein

```bash
( for i in {1..4};
do grep P10 *c_*_$i.pdb > p10_$i.pdb;
grep CONECT *c_*_$i.pdb >> p10_$i.pdb;
sed /P10/d *c_*_$i.pdb > binary_$i.pdb;
sed -i /CONECT/d binary_$i.pdb
done)

# Make sure always use lowercase 'p10', otherwise cgenff_charmm2gmx_py3_nx1.py will give you trouble
```

#### CGenFF

In pymol, Add hydrogen for `P10_$i.pdb`, do `rename *`, export by `p10_$i.mol2` format.

Edit the Resi `P10_001` into `P10`,

~~~
chmod 754 *
sed -i s/p10_[0-9]*/p10/g *.mol2
~~~

 upload to CGenFF, get `p10.str` based on `CGenFF version 4.6`

No need to copy `cgenff_charmm2gmx_py3_nx1.py` and `charmm36-jul2021.ff` from other folder to here

~~~bash
#cp ~/cgenff_charmm2gmx_py3_nx1.py .
#ln -s ~/charmm36-jul2021.ff ./

#python3 -m pip3 install networkx==1.11
#python3 -m pip3 install numpy

chmod 754 *

for i in {1..4};
do python3 ~/cgenff_charmm2gmx_py3_nx1.py p10 p10_$i.mol2 p10.str ~/charmm36-jul2021.ff;
	for f in p10.[itp]*; 
	do mv $f "$(echo "$f"|sed s/p10/p10_$i/)";
	# rename ligand topology files with numbers;
    done;
    mv p10_ini.pdb p10_$i'_ini'.pdb;    
done
~~~



#### Upload everything to `~/p10ter/$i` folder

```bash
# At TACC
mkdir $SCRATCH/p10ter
```

```bash
# At Local
scp -r * hanfeng@ls6.tacc.utexas.edu:\$SCRATCH/p10ter;
```

```bash
# At TACC
cd $SCRATCH/p10ter
for i in {3..3};
do mkdir $i;
mv *_$i.* ./$i;
mv *_$i'_ini.pdb' ./$i
done
```

```bash
 # forcefield
 cd p10ter # Forcefield in this root folder
 for i in {1..4};
 do ln -s ~/charmm36-jul2021.ff ./$i
 done
 
 
 for i in {1..6}; do i=$(expr $i + $i); echo $i; done
```

#### Write a small molecule modeling script!!

make directory contains:

`p10_1.mol2  p10_1.pdb  binary_1.pdb  p10_1.itp  p10_1.prm  p10_1.top  p10_1_ini.pdb`

```bash
for file in ./*; do
    # Check if the file is a regular file
    if [ -f "$file" ]; then
        # Extract the suffix from the file name
        filename=$(basename "$file")
        if [[ $filename == *_ini.pdb ]]; then
            i=$(echo "$filename" | cut -d '_' -f 2)
        elif [[ "$filename" == *_* ]]; then
            i=$(echo "$filename" | cut -d '.' -f 1 | rev | cut -d '_' -f 1 | rev)
        else
            continue  # Skip files without underscores
        fi
        
        # Create a directory if it doesn't exist
        mkdir -p "./$i"
        
        # Move the file to the corresponding directory
        mv "$file" "./$i/"
    fi
done
```



Forcefield is at: `$SCRATCH/charmm36-jul2021.ff`



`~/protein_ligand_MD_script_100ns.sh p10 1`

~~~bash
#!/bin/bash
# Check if this directory contains: p10_1.mol2/.pdb/.itp/.prm/.top  binary_1.pdb p10_1_ini.pdb
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
echo "SOL" | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral
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

~~~

#### Server

```bash
export OMP_NUM_THREADS=8
# 1 job on multiple GPUs
mpirun -n 8 gmx_mpi mdrun -v -s md_p10.tpr -deffnm md_p10 -ntomp 8 -gpu_id 0123 -pin on -pme gpu -npme 1 -nb gpu -bonded gpu

# multiple jobs on multiple GPUs
export OMP_NUM_THREADS=10
gmx mdrun -v -s md_p10.tpr -deffnm md_p10 -ntomp 10 -ntmpi 1 -gpu_id 0

#continue simulation
export OMP_NUM_THREADS=10
gmx mdrun -v -s md_p10.tpr -deffnm md_p10 -ntomp 10 -gpu_id 0 -ntmpi 1 -cpi md_p10.cpt
```



#### post md visualization

```
( echo 'a 10350'; echo 'name 25 center'; echo q)|gmx make_ndx -n index.ndx -o center.ndx -f em.gro
~/visualization_center_dt100.sh md_p10
```

#### Concatenate trajectories and Clustering

for 100ns md for mmpbsa:

```bash
# remove water
for folder in [1-4]; do
    cd "$folder";
    ( echo '24|13'; echo q )|gmx make_ndx -n index.ndx -o index_nowater.ndx -f em.gro; #choose Protein_p10|Zn2;
    echo Protein_p10_ZN2|gmx trjconv -s md_p10.tpr -f md_p10_fit.xtc -o md_p10_nowater_fit_`$folder`.xtc -n index_nowater.ndx;
    cd ..
done

gmx trjcat -f [1-4]/md_p10_nowater_fit*.xtc -o md_p10_nowater_100ns_cat.xtc -cat

```

#### gmx_MMPBSA

~~~bash
# At lab server
conda activate gmxMMPBSA
~~~

Conside 'CRBN' as A, and target protein as B

```
vim mmpbsa.in
```

~~~python
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
endframe=4003,
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
~~~



~~~bash
# Trying to use concatenated traj...
cp 1/em.gro .
cp 1/index.ndx ./
cp 1/md_p10.tpr ./
 gmx make_ndx -f em.gro -n index.ndx -o CRBN_CSK.ndx;
# 25 CRBN  as receptor: ri 1-381
# 26 CSK  as ligand: ri 382-644


conda activate gmxMMPBSA
ln -s ~/charmm36-jul2021.ff/ ./
cp */p10_*.prm ./
cp 1/*.itp ./
cp 1/topol.top ./
# make sure MPI is using /home/hlin/.conda/envs/gmxMMPBSA/bin/mpirun, otherwise it has parellel problems
# conda environment gmxMMPBSA has set the PATH 
conda env config vars set PATH=/data/apps/sbgrid/x86_64-linux/xv/3.10a/bin:/data/apps/sbgrid/x86_64-linux/gromacs/2022.1_cu11.5.2/bin:/home/hlin/.conda/envs/gmxMMPBSA/bin:/data/apps/sbgrid/x86_64-linux/pymol/2.5.3_386:/usr/local/cuda-11.8/bin:/usr/bin:/home/hlin/rosetta_bin_linux_2021.16.61629_bundle/main/source/scripts/python/public:/home/hlin/rosetta_bin_linux_2021.16.61629_bundle/main/source:/home/hlin/rosetta_bin_linux_2021.16.61629_bundle/main/source/bin:/programs/x86_64-linux/system/sbgrid_bin:/usr/bin:/home/hlin/.local/bin:/home/hlin/bin:/usr/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/var/lib/snapd/snap/bin:/programs/share/bin:/programs/share/sbgrid/bin:/programs/x86_64-linux/sbgrid_installer/latest



mpirun -np 80 gmx_MMPBSA MPI -O -i mmpbsa.in -cs md_p10.tpr -ci CRBN_CSK.ndx -cg 25 26 -ct md_p10_nowater_100ns_cat.xtc -cp topol.top -o FINAL_RESULTS_MMPBSA.dat -eo FINAL_RESULTS_MMPBSA.csv -nogui;
# CRBN as receptor, CSK as ligand

gmx_MMPBSA_ana
# if gmx_MMPBSA_ana has BDC keyvalue problem, try to lower the threshold of energy kcal/mol
~~~

### 