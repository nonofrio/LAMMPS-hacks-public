# Description
This example describes the application of an exteral voltage to a capacitor like structure using EChemDID.
For more details refer to Ref. https://doi.org/10.1063/1.4927562
The input file is self explenatory.
According to the input file, LAMMPS will return zeros.. it is a static simulation with no potential.
However the electrochemical potential will propagate as you can see from the dump file.
One can use the tcl script to map the local potential to the user field of vmd for visualization purpose.
You can plot/compute the error between the electrochemical diffusion versus the analytical relation with plot.py/get_error.py

lmp -in capa.in
python get_error.py
