MOFF Guide
By Andrew Latham and Bin Zhang
Please cite Latham, A; Zhang, B BioRxiv, 2021. (https://www.biorxiv.org/content/10.1101/2021.01.06.425600v1.abstract)

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
You will observe that there are 2 nested folders: Ready and To_setup. Ready is a ready to run simulation, with all the necessary files ready to go. To_setup is how to start a simulation from scratch with a new PDB.

In this guide, we will start with how to run a simulation from a ready to run trajectory. We will then walk you through setting up you own simulation conditions.

-------------------------------------------------------------------------------------------   Run Simulation         ----------------------------------------------------------------------------------------
1. Go to the folder Ready
2. Place the PDB of alpha-carbons in a simulation box.
	gmx_mpi editconf -f IBB_CA.pdb -box 50 50 50 -c -o IBB_CA.pdb
3. Run energy minimization
	gmx_mpi grompp -f min.mdp -c IBB_CA.pdb -p IBB.top -o min.tpr
	gmx_mpi mdrun -deffnm min -table table_MOFF.xvg -tablep table_smog.xvg
	(Notice the usage of the two tables. Use the MOFF table (table_MOFF.xvg) for non-bonded interactions (-table), and use the SMOG table (table_smog.xvg) for paired interactions (-tablep)
4. Run replica exchange simulation in langevin dynamics at Temperatures ranging from 300 to 400 K. 
	gmx_mpi grompp -f prod_0.mdp -c min.gro -p IBB.top -o prod_0.tpr
	gmx_mpi grompp -f prod_1.mdp -c min.gro -p IBB.top -o prod_1.tpr
	gmx_mpi grompp -f prod_2.mdp -c min.gro -p IBB.top -o prod_2.tpr
	gmx_mpi grompp -f prod_3.mdp -c min.gro -p IBB.top -o prod_3.tpr
	gmx_mpi grompp -f prod_5.mdp -c min.gro -p IBB.top -o prod_5.tpr
	mpirun -np 6 gmx_mpi mdrun -deffnm prod_ -table table_MOFF.xvg -tablep table_smog.xvg -multi 6 -replex 100
5.	To run biased simpulation (Optional):
	If you want to run a biased simulation instead of an unbiased simulation use:
	mpirun -np 6 gmx_mpi mdrun -deffnm prod_ -table table_MOFF.xvg -tablep table_smog.xvg -multi 6 -replex 100 -plumed bias.dat
	(plumed.dat is a plumed file that biases on the radius of gyration. The SLOPE is optimized at each itraction of our algorithm)


-------------------------------------------------------------------------------------------   Starting from Scratch        ----------------------------------------------------------------------------------------
1. Go to folder To_Setup
2. Run STRIDE at http://webclu.bio.wzw.tum.de/cgi-bin/stride/stridecgi.py. Download the output and name it IBB_stride.dat.
3. Run SMOG calculation
	smog2 -i IBB.pdb -t /path_to_smog/smog-2.2/SBM_AA -tCG /path_to_smog/smog-2.2/SBM_calpha
4. Generate topology and C-alpha PDB
	python ../../../Scripts/write_MOFF.py IBB ../../../Scripts/template_MOFF.top 
5. Generate a table potentials
	python ../../../Scripts/write_table.py 162
6. You have successfully set up you simulation! Copy over the *.mdp files from the Ready folder and follow the steps above.
