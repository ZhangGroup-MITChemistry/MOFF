# MOFF - MRG-CG-DNA Guide
By Andrew Latham and Bin Zhang
Please cite: Latham, A; Zhang, B BioRxiv, 2021. ()
Latham, A; Zhang, B J. Chem. Theory Comput. 2021, 17, 3134. (https://pubs.acs.org/doi/abs/10.1021/acs.jctc.0c01220)
Savelyev, A; Papoian, GA Proc. Natl. Acad. Sci. U.S.A. 20210, 107, 20340.

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
In this folder, we have the scripts necessary to run simulations combining the MOFF protein model and MRG-CG DNA model. You will observe that there are 2 nested folders: Ready and To_setup. Ready is a ready to run simulation, with all the necessary files ready to go. To_setup is how to start a simulation from scratch with seperate PDBs for each the protein and DNA. The setup is done in three parts: A) Part 1: Setup Protein, B) Part 2: Setup DNA, and C) Part 3: Combine Protein and DNA

---------------------------------------------------------------------------------- Run simulation -------------------------------------------------------------------------------------------------------------------
1. Go to the folder Ready
2. Run energy minimization (Note usage of tables. table_MOFF.xvg is for the non-bonded interactions. table_smog.xvg is for the pair (structure-based) interactions. table_b*.xvg accounts for the quartic bonded interactions in the fan potential. These bonded tables need to be ordered correctly from 0 to 11)
gmx_mpi grompp -f min.mdp -c start.gro -p hp1_dna.top -o min.tpr
gmx_mpi -deffnm min -table table_MOFF.xvg -tablep table_smog.xvg -tableb table_b0.xvg table_b1.xvg table_b2.xvg table_b3.xvg table_b4.xvg table_b5.xvg table_b6.xvg table_b7.xvg table_b8.xvg table_b9.xvg table_b10.xvg table_b11.xvg
3. Run replica exchange simulation at temperatures ranging from 300 K to 400 K
gmx_mpi grompp -f prod_0.mdp -c min.gro -p hp1_dna.top -o prod_0.tpr
gmx_mpi grompp -f prod_1.mdp -c min.gro -p hp1_dna.top -o prod_1.tpr
gmx_mpi grompp -f prod_2.mdp -c min.gro -p hp1_dna.top -o prod_2.tpr
gmx_mpi grompp -f prod_3.mdp -c min.gro -p hp1_dna.top -o prod_3.tpr
gmx_mpi grompp -f prod_4.mdp -c min.gro -p hp1_dna.top -o prod_4.tpr
gmx_mpi grompp -f prod_5.mdp -c min.gro -p hp1_dna.top -o prod_5.tpr
mpirun -np 6 gmx_mpi mdrun -deffnm prod_ -multi 6 -replex 100  -table table_MOFF.xvg -tablep table_smog.xvg -tableb table_b0.xvg table_b1.xvg table_b2.xvg table_b3.xvg table_b4.xvg table_b5.xvg table_b6.xvg table_b7.xvg table_b8.xvg table_b9.xvg table_b10.xvg table_b11.xvg



---------------------------------------------------------------------------------- Starting from Scratch -------------------------------------------------------------------------------------------------------------------
A) Part 1: Setup Protein
1. Go to folder To_Setup
2. Run STRIDE at http://webclu.bio.wzw.tum.de/cgi-bin/stride/stridecgi.py.
Download the output and name it hp1a_mono_stride.dat.
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
B) Part 2: Setup DNA
7. Generate coarse grained structure from all atom DNA structure.
python  ../../Scripts/pdb2gro_MRGcgdna.py all_atom_200bpDNA.pdb cgdna.gro
8. Generate coarse grained topology from all atom DNA structure.
python ../../Scripts/pdb2top_MRGcgdna.py all_atom_200bpDNA.pdb cgdna.top
9. Copy bonded tables. These can be found in Scripts/Tables/ and are used for quartic bonded interactions in the fan potential of MRG-CG-DNA.
cp ../../Scripts/Tables/* .
C) Part 3: Combine Protein and DNA
10. Convert *.top files to *.itp files. (This way, the bonded terms are described seperately for each molecule, but the nonboned terms are described in 1 location)
python ../../Scripts/top2itp.py hp1a_mono
python ../../Scripts/top2itp.py cgdna
11. Copy tertiary potential from Scripts to the active folder
cp ../../Scripts/prot_dna.itp .
12. Write new topology. This script adds the protein molecules and DNA molecules to your topology. The names of each molecule must be consistent with the name of the itp file. The same commands can be used to add more than 2 types of molecules.
python ../../Scripts/write_top.py hp1_dna.top hp1a_mono 1 cgdna 1
13. Add molecules simulation box with GROMACS insert-molecules command. The numbers here must match the numbers used in the write_top.py
gmx_mpi insert-molecules -nmol 1 -ci hp1a_mono_CA.pdb -scale 5 -radius 2 -box 200 200 200 -o temp.gro
gmx_mpi insert-molecules -nmol 1 -ci cgdna.gro -scale 5 -radius 2 -f temp.gro -o start.gro

        
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Email aplatham@mit.edu with any questions
