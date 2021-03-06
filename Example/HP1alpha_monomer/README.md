MOFF Guide
By Andrew Latham and Bin Zhang
Please cite Latham, A; Zhang, B J. Chem. Theory Comput. 2021, 17, 3134. (https://pubs.acs.org/doi/abs/10.1021/acs.jctc.0c01220)

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

This portion of the tutorial will demonstrate how to handle ordred regions in large proteins. While the majority of the framework is the same, there is one additional step in setting up the simulation.
This step (step 6 below), adds the structure based folded potential, U_fold. In our study, we utilize this potential to model HP1 dimers, and here we walk through how to set up such a simulation for a monomer of HP1. We recommend tuning the strength of this potential to all atom simulations.
More information on this type of modeling is available through the SMOG websites (https://smog-server.org).

python write_Ufold.py jobid start1 end1 smog_file pair_eps

jobid - name of the pdb / topology you want to run

start1 - starting index for including the folded potential

end1 - ending index for including the folded potential

smog_file - name of smog topology file (default = smog.top)

pair_eps - strength of folded potential (default = 3 kJ/mol, this should be fit to all atom simulations)

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
You will observe that there are 2 nested folders: Ready and To_setup. Ready is a ready to run simulation, with all the necessary files ready to go. To_setup is how to start a simulation from scratch with a new PDB.
In this guide, we will start with how to run a simulation from a ready to run trajectory. We will then walk you through setting up you own simulation conditions.

-------------------------------------------------------------------------------------------   Run Simulation         ----------------------------------------------------------------------------------------
1. Go to the folder Ready
2. Place the PDB of alpha-carbons in a simulation box.
	gmx_mpi editconf -f hp1a_mono_CA.pdb -box 50 50 50 -c -o hp1a_mono_CA.pdb
3. Run energy minimization
	gmx_mpi grompp -f min.mdp -c hp1a_mono_CA.pdb -p hp1a_mono.top -o min.tpr
	gmx_mpi mdrun -deffnm min -table table_MOFF.xvg -tablep table_smog.xvg
	(Notice the usage of the two tables. Use the MOFF table (table_MOFF.xvg) for non-bonded interactions (-table), and use the SMOG table (table_smog.xvg) for paired interactions (-tablep)
4. Run replica exchange simulation in langevin dynamics at Temperatures ranging from 300 to 400 K. 
	gmx_mpi grompp -f prod_0.mdp -c min.gro -p hp1a_mono.top -o prod_0.tpr
	gmx_mpi grompp -f prod_1.mdp -c min.gro -p hp1a_mono.top -o prod_1.tpr
	gmx_mpi grompp -f prod_2.mdp -c min.gro -p hp1a_mono.top -o prod_2.tpr
	gmx_mpi grompp -f prod_3.mdp -c min.gro -p hp1a_mono.top -o prod_3.tpr
	gmx_mpi grompp -f prod_5.mdp -c min.gro -p hp1a_mono.top -o prod_5.tpr
	mpirun -np 6 gmx_mpi mdrun -deffnm prod_ -table table_MOFF.xvg -tablep table_smog.xvg -multi 6 -replex 100


-------------------------------------------------------------------------------------------   Starting from Scratch        ----------------------------------------------------------------------------------------
1. Go to folder To_Setup
2. Run STRIDE at http://webclu.bio.wzw.tum.de/cgi-bin/stride/stridecgi.py. Download the output and name it hp1a_mono_stride.dat.
3. Run SMOG calculation
	smog2 -i hp1a_mono.pdb -t /path_to_smog/smog-2.2/SBM_AA -tCG /path_to_smog/smog-2.2/SBM_calpha
4. Generate topology and C-alpha PDB
	python ../../../Scripts/write_MOFF.py hp1a_mono ../../../Scripts/template_MOFF.top 
5. Generate a table potentials
	python ../../../Scripts/write_table.py 82
6. Add folded potential
	Add U_fold for the CD, residues 20 to 78:
	python write_Ufold.py hp1a_mono 20 78 smog.top 6.0 
	Add U_fold for the CSD, residues 121 to 179:
	python write_Ufold.py hp1a_mono 121 179 smog.top 6.0
	
