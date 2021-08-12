# MOFF - MRG-CG-DNA Guide
By Andrew Latham and Bin Zhang
Please cite: Latham, A; Zhang, B BioRxiv, 2021. ()
Latham, A; Zhang, B J. Chem. Theory Comput. 2021, 17, 3134. (https://pubs.acs.org/doi/abs/10.1021/acs.jctc.0c01220)
Savelyev, A; Papoian, GA Proc. Natl. Acad. Sci. U.S.A. 20210, 107, 20340.

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
In this folder, we have the scripts necessary to run simulations with MRG-CG-DNA. You will observe that there are 2 nested folders: Ready and To_setup. Ready is a ready to run simulation, with all the necessary files ready to go. To_setup is how to start a simulation from scratch with a new PDB.

In this guide, we will start with how to run a simulation from a ready to run trajectory. We will then walk you through setting up a simulation from a PDB structure.

---------------------------------------------------------------------------------- Run simulation -------------------------------------------------------------------------------------------------------------------
1. Go to the folder Ready
2. Place the coarse grained structure in a simulation box.
gmx_mpi editconf -f cgdna.gro -box 200 200 200 -c -o cgdna.gro
3. Run energy minimization (Note usage of tables. table_MOFF.xvg is for the non-bonded interactions. table_b*.xvg accounts for the quartic bonded interactions in the fan potential. These bonded tables need to be ordered correctly from 0 to 11)
gmx_mpi grompp -f min.mdp -c cgdna.gro -p cgdna.top -o min.tpr
gmx_mpi mdrun -deffnm min -table table_MOFF.xvg -tableb table_b0.xvg table_b1.xvg table_b2.xvg table_b3.xvg table_b4.xvg table_b5.xvg table_b6.xvg table_b7.xvg table_b8.xvg table_b9.xvg table_b10.xvg table_b11.xvg
4. Run replica exchange simulation at temperatures ranging from 300 K to 400 K
gmx_mpi grompp -f prod_0.mdp -c min.gro -p cgdna.top -o prod_0.tpr
gmx_mpi grompp -f prod_1.mdp -c min.gro -p cgdna.top -o prod_1.tpr
gmx_mpi grompp -f prod_2.mdp -c min.gro -p cgdna.top -o prod_2.tpr
gmx_mpi grompp -f prod_3.mdp -c min.gro -p cgdna.top -o prod_3.tpr
gmx_mpi grompp -f prod_4.mdp -c min.gro -p cgdna.top -o prod_4.tpr
gmx_mpi grompp -f prod_5.mdp -c min.gro -p cgdna.top -o prod_5.tpr
mpirun -np 6 gmx_mpi mdrun -deffnm prod_ -multi 6 -replex 100  -table table_MOFF.xvg -tableb table_b0.xvg table_b1.xvg table_b2.xvg table_b3.xvg table_b4.xvg table_b5.xvg table_b6.xvg table_b7.xvg table_b8.xvg table_b9.xvg table_b10.xvg table_b11.xvg



---------------------------------------------------------------------------------- Starting from Scratch -------------------------------------------------------------------------------------------------------------------
1. Go to folder Setup
2. Generate coarse grained structure from all atom structure.
python  ../../Scripts/pdb2gro_MRGcgdna.py all_atom_200bpDNA.pdb cgdna.gro
3. Generate coarse grained topology from all atom structure.
python ../../Scripts/pdb2top_MRGcgdna.py all_atom_200bpDNA.pdb cgdna.top
4. Generate a non-bonded table potentials at the appropriate ionic strength (Note: this barrows code from the MOFF protein force field and only the table_MOFF is needed. Ionic strengths are in mM)
python  ../../../Scripts/write_table.py 100
5. Copy bonded tables. These can be found in Scripts/Tables/ and are used for quartic bonded interactions in the fan potential of MRG-CG-DNA.
cp ../../Scripts/Tables/* .

        
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Email aplatham@mit.edu with any questions
