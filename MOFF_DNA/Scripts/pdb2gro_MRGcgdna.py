##########################################################################
# This script will prepare the coarse-grained Gromacs .gro file
#
# Written by Xingcheng Lin, 05/03/2017;
###########################################################################

import time
import subprocess
import os
import math
import sys
import numpy as np

################################################


def my_lt_range(start, end, step):
    while start < end:
        yield start
        start += step


def my_le_range(start, end, step):
    while start <= end:
        yield start
        start += step
#############################################

def switch_resname(argument):
    switcher = {
        "DA": "ADE",
        "DT": "THY",
        "DC": "CYT",
        "DG": "GUA",
    }
    return switcher.get(argument, "Invalid resname")

# Add line to the top of a file;
def line_prepender(filename, line):
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)

def pdb2gro_single_residue_dna(inputFile, outputFile):

    # Get current working directory
    pwd = os.getcwd()

    infile = open(inputFile, 'r')
    outfile = open(outputFile, 'w')

    output_coord_precision = 3;

    # Read in lines from the file;

    lines = [line.strip() for line in infile]

    infile.close()

    length = len(lines)

    aa_coords = []
    aa_resid = 0

    for i in my_lt_range(0, length, 1):

        line = lines[i].split()
        # Find the relevant lines;
        if (line[0] == "ATOM"):
            if (int(line[5]) == (aa_resid+1) and aa_resid == 0):
                # When it is in the first line of the first residue, we just add 
                # the atomistic coordinates in
                aa_coords = np.array([float(line[6]), float(line[7]), float(line[8])])

                aa_resname = line[3]
                aa_resid = int(line[5])

            elif (int(line[5]) == (aa_resid+1) and aa_resid != 0):
                # When it is in the first line of the one (not the first) residue,
                # we perform the average of all the coordinates for the coarse-grained
                # level coordinates and add in the new line of atomistic coordinates;
                
                # Do average and convert to nm unit
                cg_coords = np.average(aa_coords, axis=0)/10.0
                cg_coords = np.round(cg_coords, decimals=output_coord_precision)

                # Output to file
                # Decide the resname
                cg_resname = switch_resname(aa_resname)
                cg_resid = aa_resid

                outfile.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" % (cg_resid, str(cg_resname), str("NUC"), cg_resid, cg_coords[0], cg_coords[1], cg_coords[2]))

                # Reset and record the first line of aa coordinates;
                aa_coords = np.array([float(line[6]), float(line[7]), float(line[8])])

                aa_resname = line[3]
                aa_resid = int(line[5])

            else:
                # When it is in the rest of the lines of that residue, we just add 
                # the atomistic coordinates in
                to_be_added =  np.array([float(line[6]), float(line[7]), float(line[8])])
                aa_coords = np.vstack((aa_coords, to_be_added))

        # Reaching the end of the ATOM section, record the last residue
        elif (line[0] == "END"):
            # When it is in the first line of the one (not the first) residue,
            # we perform the average of all the coordinates for the coarse-grained
            # level coordinates and add in the new line of atomistic coordinates;
            
            # Do average and convert to nm unit
            cg_coords = np.average(aa_coords, axis=0)/10.0
            cg_coords = np.round(cg_coords, decimals=output_coord_precision)

            # Output to file
            # Decide the resname
            cg_resname = switch_resname(aa_resname)
            cg_resid = aa_resid

            outfile.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" % (cg_resid, str(cg_resname), str("NUC"), cg_resid, cg_coords[0], cg_coords[1], cg_coords[2]))

        # Pass the information of PBC into the file;
        elif (line[0] == "CRYST1"):
            pbc_info = [line[1], line[2], line[3]]

    #outfile.write(str(pbc_info[0]) + "\t" + str(pbc_info[1]) + "\t" + str(pbc_info[2]) + "\n")
    outfile.write(str(1000) + "\t" + str(1000) + "\t" + str(1000) + "\n")

    outfile.close()

    # Add the number of residues to the top of the file
    line_prepender(outputFile, str(cg_resid))
    # Add a verbal line to the top of the file (required by the format of .gro file)
    line_prepender(outputFile, "Gro file for the CG DNA")

    return


############################################################################

if __name__ == "__main__":
    # Input all-atom DNA structure;
    inputFile = sys.argv[1]
    # Output CG DNA structure;
    outputFile = sys.argv[2]


    pdb2gro_single_residue_dna(inputFile, outputFile)

    print("When the voice of the Silent touches my words,")
    print("I know him and therefore know myself.")
