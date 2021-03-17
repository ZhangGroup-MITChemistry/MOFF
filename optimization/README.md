# MOFF Guide
By Andrew Latham and Bin Zhang
Please cite Latham, A; Zhang, B BioRx, 2021. (https://www.biorxiv.org/content/10.1101/2021.01.06.425600v1.abstract)

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
This portion walks you through our force field optimization procedure. The process is done iteratively, and data from iteration0 and iteration15 are in the corresponding folders. While this procedure for optimizing the tertiary potential is the same, the secondary structure is different, as highlighted in the subsection "Varying the secondary structure potential during optimization".
Note that the numbering schemes are different in the simulation folders and the folders in make_ff. The data from simulations at iteration0 are used to make the force field that will be used for iteration1, and the data from simulations at iteration15 are used to make the force field that will be used for iteration16, which inspires the differences in numbering. In this case, iteration0 is the first iteration of our force field (or MJ as defined in our paper), and iteration15 is the final iteration of our force field (or MOFF as defined in our paper).

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Do steps 1-5 on all proteins to generate the biased and unbiased ensemble:
1. Run simulation in unbiased ensemble. (Note IBB is used here but can be any protein)
        A. Place the PDB of alpha-carbons in a simulation box.
            gmx_mpi editconf -f IBB_CA.pdb -box 50 50 50 -c -o IBB_CA.pdb
        B. Run energy minimization
            gmx_mpi grompp -f min.mdp -c IBB_CA.pdb -p IBB.top -o min.tpr
            gmx_mpi mdrun -deffnm min -table table_MOFF.xvg -tablep table_smog.xvg
            (Notice the usage of the two tables. Use the MOFF table (table_MOFF.xvg) for non-bonded interactions (-table), and use the SMOG table (table_smog.xvg) for paired interactions (-tablep). Note that only proteins with alpha helices will require a second table.
        C. Run replica exchange simulation in langevin dynamics at Temperatures ranging from 300 to 400 K. 
            gmx_mpi grompp -f prod_0.mdp -c min.gro -p IBB.top -o prod_0.tpr
            gmx_mpi grompp -f prod_1.mdp -c min.gro -p IBB.top -o prod_1.tpr
            gmx_mpi grompp -f prod_2.mdp -c min.gro -p IBB.top -o prod_2.tpr
            gmx_mpi grompp -f prod_3.mdp -c min.gro -p IBB.top -o prod_3.tpr
            gmx_mpi grompp -f prod_5.mdp -c min.gro -p IBB.top -o prod_5.tpr
            mpirun -np 6 gmx_mpi mdrun -deffnm prod_ -table table_MOFF.xvg -tablep table_smog.xvg -multi 6 -replex 100
2. Run simulation in biased ensemble
            A. Repeat steps (A-B) above.
            B. Prepare .tpr files as done for the unbiased ensemble.
            C. Run simulation. Use -plumed extension to add in bias.
            mpirun -np 6 gmx_mpi mdrun -deffnm prod_ -table table_MOFF.xvg -tablep table_smog.xvg -multi 6 -replex 100 -plumed bias.dat
3. Optimize the biased ensemble at 300K. 
            A. Use rg_ave.py. The code will output 2 numbers. One is the mean Rg, after equilibration. The other is the standard deviation.Use rg_ave.py.
            rg_ave.py rg.0.dat
            B. If the Rg is within 0.05 nm of the experimental value, proceed to step 4. Or else, do step 3C.
            C. Change the value of SLOPE in the file bias.dat to optimize the Rg. Generally, try changing the Rg by 0.5 kJ/mol / nm. Larger (more positive) values make the Rg smaller, and smaller (more negative) values make the Rg bigger. Afterwards, proceed to step 3B.
4. After all of the simulations are completed, collect the necessary data in the correct format. We need a pdb file and radius of gyration file of both the biased and unbiased trajectory with the correct timestep and equilibration time.
        echo 0 | /home/aplatham2/Programs/gromacs-2018.4/bin/gmx_mpi gyrate -f prod_0.xtc -s prod_0.tpr -b 100000 -o gyrate.xvg  -dt 200
        echo 0 | /home/aplatham2/Programs/gromacs-2018.4/bin/gmx_mpi trjconv -f prod_0.xtc -s prod_0.tpr -b 100000 -o prod_0.pdb  -dt 200
5. Move simulations to correct folder ( or bias_dat).
        A. Move the unbiased ensemble to the nb_dat folder, call the files *_unbiased.xvg (for gyrate.xvg) and *_unbiased.pdb (for prod_0.pdb)
        B. Move the biased ensemble to the bias_dat folder, call the files *_bias.xvg (for gyrate.xvg) and *_bias.pdb (for prod_0.pdb)

Do steps 6-9 to generate the force field from the trajectories:
6. Update the presents and parameters section of prepare_opt.py. The important paramters are highlighted below. From iteration to iteration, only alpha should change.
        A. type - the number of types of atoms in your force field (20 for all amino acids)
        B. timesteps - the timesteps in each configuration
        C. ordered- number of ordered proteins in your training set.
        D. rg_exp- the experimental value of Rg for each protein in the training set
        E. alpha- the slope of the biasing potential for each Rg in the training set. IMPORTANT: this needs to be updated to match the value of SLOPE after each optimization
        F. The name of all training proteins is taken in as an argument. The order is improtant. The ordered proteins should come first. rg_exp and alpha values should be in the same order as the proteins in the argument.
7. Run prepare_opt.py - This code calculates the contacts and biasing energies in preparation from fitting. IMPORTANT: make sure order of proteins given matches rg_exp and alpha (step 6F)
        python prepare_opt.py 1soy  1ubq  1wla  3mzq  5tvz  6eez  6h8m  ACTR  An16 asynuclein  ERMTADn hNHE1cdt IBB N49 N98 NLS NSP NUL NUS P53 ProTa sh4ud Sic
        Outputs:
        A. D1.txt, D.txt - values from diagonal matrix generated by SVD. Used for debugging.
        B. contact_list.txt - list of contacts in the sampled structures minus contacts from the structures with Rg within 0.05 nm (C in Eq. 2). Noise been removed using SVD (Eq. S8).
        C. bias_tot.txt - Biasing energy of all structures (alpha*Rg in Eq 2)
        D. pdb_list.txt - Contacts of the PDB structure (C_PDB in Eq. 3)
        E. pdb_list_contacts.txt - list of contacts in the sampled structures (C_sim in Eq. 3). Contacts from Rg within 0.05 have not been removed, and noise has not been removed.
8. Run calc_eps.m
    This is the code that does the necessary fitting by simulataneously solving Eq. 2 and Eq. 3. 
    matlab nodisplay -nosplash - nodesktop -r "calc_eps; exit"
    Requires old contact energy in matrix form, called ("iter*.dat")
    Inputs from prepare_opt.py: contact_list.txt, bias_tot.txt, pdb_list.txt, and pdb_list_contacts.txt
    Note definitions of lb and ub. These are the upper and lower bounds on changes in the radius of gyration.
    Outputs:
        A. E_pdb.txt, E_sim.txt - new energies of the PDB and simulated energies. For debugging.
        B. delta_eps.txt - matrix representing the change in contact energy from solving Eq. 2 and Eq. 3.
9. Run finish_opt.py
    This script converts the change in contact energy to the non_bonded parameters read by gromacs.
    Inputs:
    A.  delta_eps.txt - matrix representing the change in contact energy from solving Eq. 2 and Eq. 3., generated by calc_eps.m
    B. sigma.dat - matrix of amino acid excluded volumes
    C. IMPORTANT: iter*.dat - matrix of contact energy from the previous iteration. The wild card needs to be replaced with the correct file / energy matrix on each iteration.
    Outputs:
    A. new_eps.dat - new contact energies. This is the same file as iter*.dat, but now for the current iteration instead of the previous iteration.
    B. new_ff.dat - new force field. This is the "nonbond_params" section of GROMACS topology. To run a simulation in the new ensemble, replace that section of the topology files from the previous iteration with the contents of this file.
            
