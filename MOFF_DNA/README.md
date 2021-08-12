# MOFF - MRG-CG-DNA Guide
By Andrew Latham and Bin Zhang
Please cite: Latham, A; Zhang, B BioRxiv, 2021. ()
Latham, A; Zhang, B J. Chem. Theory Comput. 2021, 17, 3134. (https://pubs.acs.org/doi/abs/10.1021/acs.jctc.0c01220)
Savelyev, A; Papoian, GA Proc. Natl. Acad. Sci. U.S.A. 20210, 107, 20340.

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
In this folder, we have the scripts necessary to run MOFF simulations with MRG-CG-DNA (Scripts), an example of a coarse grained DNA simulation (DNA_Example), and an example of a coarse grained protein-DNA simulation (Protein_DNA_Example).

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
What you will need:
- GROMACS: a version of GROMACS that supports tabulated potentials. These were temporarily removed in 2020, but updated in 2021. (Version 2018.4 was used in published work)
- Structure-based models of biomolecules (SMOG): https://smog-server.org (Version 2.2 was used in published work)
- STRIDE: Use the stide web interface to generate a secondary structure prediction from your input PDB (http://webclu.bio.wzw.tum.de/cgi-bin/stride/stridecgi.py)
- PDB: structure of the protein or DNA you want to simulate. This can be taken from the PDB or from structure prediction tools such as RaporX (http://raptorx.uchicago.edu) or i-TASSER (https://zhanglab.ccmb.med.umich.edu/I-TASSER/)

        
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Check out the DNA_Example folder for an example on how to do a simple simulation on 200 bp of DNA, or the Protein_DNA_Example for how to do a simulation with both an HP1 monomer and DNA. A separate README will walk you through each example.

Email aplatham@mit.edu with any questions
