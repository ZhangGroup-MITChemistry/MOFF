# MOFF Guide
By Andrew Latham and Bin Zhang
Please cite Latham, A; Zhang, B J. Chem. Theory Comput. 2021, 17, 3134. (https://pubs.acs.org/doi/abs/10.1021/acs.jctc.0c01220)

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
In this folder, we have the scripts necessary to run MOFF simulations (Scripts), 2 examples of simulating proteins in MOFF (Examples), and an example of how we optimized the force field, along with its corresponding code (Optimization).

IBB is an example of our default potential on an IDP in our training set, and HP1_alpha_monomer is an example of how we can add a folded potential in ordered regions, as we do for HP1 in our study.

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
What you will need:
- GROMACS: a version of GROMACS that supports tabulated potentials. These were temporarily removed in 2020, but updated in 2021. (Version 2018.4 was used in published work)
- Structure-based models of biomolecules (SMOG): https://smog-server.org (Version 2.2 was used in published work)
- STRIDE: Use the stide web interface to generate a secondary structure prediction from your input PDB (http://webclu.bio.wzw.tum.de/cgi-bin/stride/stridecgi.py)
- PDB: structure of the protein you want to simulate. This can be taken from the PDB or from structure prediction tools such as RaporX (http://raptorx.uchicago.edu) or i-TASSER (https://zhanglab.ccmb.med.umich.edu/I-TASSER/)


---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Steps:
1.		Run STRIDE calculation
		a. Upload PDB to http://webclu.bio.wzw.tum.de/cgi-bin/stride/stridecgi.py and download results
2.		Run SMOG calculation
		a. Run SMOG on your input PDB. Generate an all atom contact map and project a CG model
3.      Setup template and pdb
        a.      Run write_MOFF.py. There are up to 5 possible inputs:
        (python write_MOFF.py jobid template_file smog_file stride_file)
        		i.		jobid is the name of your pdb file. Output topology and starting structure will be similarly named
        		ii.		template_file is the force field template. Unless doing further method development, this should be template_MOFF.top
        		iii.	smog_file is the SMOG topology generated in step 2. Called smog.top by default
        		iv.		stride_file is the STRIDE file generated in step 1. Called $jobid_stide.dat by default
			v.	pair_eps is the strength of the pair potentials that stabilize alpha-helices. We used 3 kJ/mol for disordered proteins and 6 kJ/mol for ordered proteins. The default value is 3 kJ/mol.
        b.      The output should be a CA only pdb file ($jobid_CA.pdb) and a topology file ($jobid.top). These will be input for GROMACS.
4.      Generate a table potential
        a.      Use write_table.py to make a tabulated potential for MOFF. There are up to 7 input possibilities. Defaults for those beyond the 1st three are recommended:
        (python write_table.py ionic_strength MOFFout SMOGout cut cut2 table_length)
                i.      ionic_strength is used for the debye-huckel electrostatic interactions with a distance dependent dielectric constant, in mM (default=150 mM).
                ii.     MOFFout is the name of the output table for the MOFF portion of the potential (default=table_MOFF.xvg)
                iii.    SMOGout is the name of the output table for the SMOG portion of the potential (default=table_smog.xvg)
                iv.    	cut is the distance at which electrostatic interactions are cut to 0 with a fifth degree polynomial switching function (default=1.5 nm)
                v.     	cut2 is the distance at which electrostatic interactions are dampened with a fifth degree polynomial switching function (default=1.2 nm)
                vi.     table_length is the length of the generated table (default=15 nm)
                vii.    dr is the minimal distance at which the table is calculated (default=0.002 nm)
5.      Optional: For folded proteins, additional potentials can be added to stabilize folded structures (U_fold). See the README in MOFF/Example/HP1alpha_monomer for full details on using write_Ufold.py.

6.      Run simulation
        a.      Place in an empty simulation box
        b.      Minimize energy
        c.      Run your simulation
        
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Check out the Example folder for an example on how to do a simple simulation on IBB from our training set, or, for a more complicated example with U_fold, you can try the HP1_alpha_monomer. The optimization folder walks you through how to use our code to update a force field iteratively. The data to run iteration0 and iteration15 are included, along with the code to generate a new force field. A separate README will walk you through each example.

Email aplatham@mit.edu with any questions
