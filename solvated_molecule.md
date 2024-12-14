
## step 1
submit molecular structure at
https://atb.uq.edu.au/

## step 2
after calculation is completed 
![image](https://github.com/user-attachments/assets/b9043ed8-3676-4351-8aa3-23c0ce3a96b1)
![image](https://github.com/user-attachments/assets/80f63da2-b7d1-4afd-8fb6-9f33f22e878d)



## step 3
send the optimized.pdb file into terachem to further optimized the structure

## step 4
used optimised structure to update the resp file in the molecule.ito downloaed in step 2. The molecular structure from DFT can be updated into the optimized.pdb obtained in step 2

## step 5
torsional dihedral angle can be further quantified using DFT and the potential extracted and can be quantified using this file and updated the molecule.itp file

## step 6
gmx_mpi editconf -f molecule.pdb -o molecule.gro

## step 7
make sure the index start as 1 in molecule.gro

## step 8
determine how many molecule you want to add into for 10nm x10nm x10nm box. gmx_mpi insert-molecules -ci molecule.gro -nmol 25 -box 10 10 10 -o box_with_molecule.gro

## step 9
then solvate the box with solvent.gro (you moight need to get the itp and pdb from atb website)
gmx_mpi insert-molecules -f box_with_molecule.gro -ci solvent.gro -nmol 20000 -o combined_box.gro
check the putput file and see how many solvents are inserted.

## step 10
update your comobined topology for example combined_box.itp
```plaintext
#include "gromos54a7_atb.ff/forcefield.itp"

; Include the parameters for the molecule
#include "solvent.itp"
#include "molecule.itp"

[ system ]
; Name of the system
mol only system

[ molecules ]
; Compound    #mols
molecule         25
solvent       4987
```
Beware that under [ molecules ], the order must be the same with combined_box.gro
step 11
minimize the box
gmx_mpi grompp -f minim.mdp -c combined_box.gro -r combined_box.gro -p combined_box.top -o em.tpr -maxwarn 4


minim.mdp can be

```plaintext
; minim.mdp - used as input into grompp to generate em.tpr
; Parameters describing what to do, when to stop and what to save
periodic_molecules = yes
integrator = steep          ; Algorithm (steep = steepest descent minimization)
emtol       = 1        ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.00001          ; Minimization step size
nsteps      = 50000         ; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist         = 20         ; Frequency to update the neighbor list and long range forces
ns_type         = grid      ; Method to determine neighbor list (simple, grid)
coulombtype     = PME       ; Treatment of long range electrostatic interactions
rcoulomb        = 1.2       ; Short-range electrostatic cut-off
rvdw            = 1.2      ; Short-range Van der Waals cut-off
pbc             = xyz       ; Periodic Boundary Conditions in all 3 dimensions
```

## Step 11
then you can run minimized and low-grade npt in CPU

using the command below

```plaintext
mpirun gmx_mpi mdrun -v -deffnm em -ntomp 1
sleep 5
gmx_mpi grompp -f npt.mdp -c em.gro -r em.gro -p combined_boxtop -o npt.tpr -maxwarn 3
sleep 5
mpirun gmx_mpi mdrun -v -deffnm npt -ntomp 1
```


the npt.mdp is (where molecule and solvent are the residue name)
```plaintext
title                   = NPT equilibration 
; Run parameters
periodic_molecules = yes
integrator              = md        ; leap-frog integrator
nsteps                  = 50000000    ; 10 ns
dt                      = 0.001     ;  1 fs
; Output control
nstxout                 = 2500000      ; save coordinates every 1.0 ns
nstvout                 = 2500000       ; save velocities every 1.0 ns
nstenergy               = 2500000      ; save energies every 1.0 ns
nstlog                  = 2500000      ; update log file every 1.0 ns
; Bond parameters
continuation            = no       ; Restarting after NVT 
constraint_algorithm    = lincs     ; holonomic constraints 
constraints             = H-bonds    ; bonds involving H are constrained
lincs_iter              = 4         ; accuracy of LINCS
lincs_order             = 8         ; also related to accuracy
; Nonbonded settings 
cutoff-scheme           = Verlet    ; Buffered neighbor searching
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb                = 1.2       ; short-range electrostatic cutoff (in nm)
rvdw                    = 1.2       ; short-range van der Waals cutoff (in nm)
coulomb-modifier 	= Potential-shift-verlet
vdw_type 		= cutoff
vdw-modifier 		= Potential-shift-verlet
epsilon_r 		=4
DispCorr                = EnerPres  ; account for cut-off vdW scheme
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.12      ; grid spacing for FFT
; Temperature coupling is on
tcoupl                  = V-rescale      ; modified Berendsen thermostat
tc-grps                 = molecule solvent    ; two coupling groups - more accurate
tau_t                   = 0.5  0.5      ; time constant, in ps
ref_t                   = 300  300     ; reference temperature, one for each group, in K

pcoupl                  = Berendsen   ; Pressure coupling on in NPT
pcoupltype              = semiisotropic
tau_p                   = 2.0                ; Time constant, in ps
ref_p                   = 1.0 1.0           ; Reference pressure for x/y and z directions
compressibility         = 4.5e-5 4.5e-5     ; Compressibility for x/y and z directions
refcoord_scaling        = com
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Velocity generation
gen_vel                 = yes       ; assign velocities from Maxwell distribution
```

## step 12

