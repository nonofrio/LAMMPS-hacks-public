This package enables the electrochemical dynamics with implicit degrees of freedom (EChemDID) with any potential in LAMMPS which is coupled to QEq.
Origninally, this was developed to work with ReaxFF but any potential should work as long as it uses QEq.
The method describes the equilibration of external electrochemical potentials (voltage) within metallic structures and their effect on the self-consistent partial atomic charges used in reactive molecular dynamics. 
An additional variable assigned to each atom denotes the local potential in its vicinity and we use fictitious, but computationally convenient, dynamics to describe its equilibration within connected metallic structures on-the-fly during the molecular dynamics simulation. 
This local electrostatic potential is used to dynamically modify the atomic electronegativities used to compute partial atomic changes via charge equilibration. 

Please cite the following papers if you are using this package in a publication:
- Voltage equilibration for reactive atomistic simulations of electrochemical processes. Nicolas Onofrio and Alejandro Strachan, The Journal of Chemical Physics, 143, 054109, 2015 https://doi.org/10.1063/1.4927562
- Atomic origin of ultrafast resistance-switching in nanoscale electrometallization cells. Nicolas Onofrio, David Guzman and Alejandro Strachan, Nature Materials, vol 14, 4, 2015 https://doi.org/10.1038/NMAT4221

To install just copy this directory into LAMMPS/src/ECHEMDID and from src type:

`make yes-qeq`  
`make yes-echemdid` // In this order

Note with cmake: 
- create src/ECHEMDID directory and copy fix_echemdid.* in there
- backup src/QEQ/fix_qeq.* and src/QEQ/fix_qeq_shielded.* (i.e. cp src/QEQ/fix_qeq.* src/QEQ/fix_qeq.*_bak, etc.)
- overwrite src/QEQ/fix_qeq.* and src/QEQ/fix_qeq_shielded.* with those provided in this directory

cmake ../cmake -D PKG_QEQ=on -D PKG_ECHEMDID=on && make -j4
