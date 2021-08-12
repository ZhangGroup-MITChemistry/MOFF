###########################################################################
# This script will prepare the corresponding Gromacs .top file
#
# Written by Xingcheng Lin, and Andrew Latham, 07/15/2020
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


def pdb2top_single_residue_dna(inputFile, outputFile):
    damp=0.9 # reduce strength of interactions by a scalar

    # Conversion factor from kcal/mol to kj/mol
    cal2j = 4.184 * damp
    deg2rad = 0.0174533

    # Get current working directory
    pwd = os.getcwd()

    infile = open(inputFile, 'r')
    outfile = open(outputFile, 'w')

    output_coord_precision = 3

    # Output the common sections:
    outfile.write("; CG DNA from Garegin A. Papoian\n")
    outfile.write("; Implemented by Xingcheng Lin and Andrew Latham\n")
    outfile.write("; Bin Zhang Group\n")
    outfile.write("\n")
    outfile.write("\n")

    outfile.write("[ defaults ]\n")
    outfile.write("; nbfunc comb-rule gen-pairs\n")
    outfile.write("  1      1         no\n")
    outfile.write("\n")

    outfile.write(" [ atomtypes ] \n")
    outfile.write("; name mass   charge ptype c6            c12\n")
    outfile.write("  NB_2 325.000  -1.000  A     0.00000e+00   1.67772e-05\n")
    outfile.write("\n")

    outfile.write("[ moleculetype ]\n")
    outfile.write("; name            nrexcl\n")
    # nrexcl = 0 here because we don't want to exclude the atoms pairs involved in the Fan interactions;
    outfile.write("  Macromolecule   0\n")
    outfile.write("\n")

    # The [ Atom ] Section;
    outfile.write("[ atoms ]\n")
    outfile.write(";nr  type  resnr residue atom   cgnr\n")
    # Read in lines from the file;

    lines = [line.strip() for line in infile]

    infile.close()

    length = len(lines)

    aa_resid = 0

    for i in my_lt_range(0, length, 1):

        line = lines[i].split()

        # Find the relevant lines;
        if (line[0] == "ATOM"):
            if (int(line[5]) == (aa_resid+1)):
                # When it is in the first line of the one residue,
                # update aa_resid
                aa_resid = int(line[5])
                # start recording

                # CA model, the atom index equals the residue index
                output_nr = line[5]
                output_type = "NB_2"
                output_resnr = line[5]
                aa_resname = line[3]
                output_residue = switch_resname(aa_resname)
                output_atom = "NUC"
                # CA model, cgnr index equals the residue index
                output_cgnr = output_resnr
                outfile.write(str(output_nr) + " " + str(output_type) + " " + str(output_resnr) + " "
                              + str(output_residue) + " " + str(output_atom) + " " + str(output_cgnr) + "\n")

    # At this moment, the total number of residues have been recorded
    tot_num_residues = aa_resid
    # number of residues for each chain
    tot_num_residues_per_chain = tot_num_residues / 2

    # Section separator;
    outfile.write("\n")

    # The [ bonds ] Section;
    outfile.write("[ bonds ]\n")
    outfile.write(";ai     aj      func    r0(nm)  Kb\n")

    l0 = 4.96 * 0.1
    K1 = 2.625 * 100 * cal2j
    K2 = -0.226 * 1000 * cal2j
    K3 = 0.0149 * 10000 * cal2j

    # Reset aa_resid
    aa_resid = 0
# Bonded interactions (single) ---------------------------------------------------------------------------------------
    for i in my_lt_range(0, length, 1):

        # Flag for checking whether to input the bonding / fan interactions;
        input_flag = 0

        line = lines[i].split()
        # Find the relevant lines;
        if (line[0] == "ATOM"):
            if (int(line[5]) == (aa_resid+1)):
                # When it is in the first line of the one residue,
                # update aa_resid
                aa_resid = int(line[5])
                # start recording

                # find the resid for the intrachain partner
                bonding_partner_resid = aa_resid + 1
                # Check if the pair is within one chain
                if (aa_resid <= tot_num_residues_per_chain and bonding_partner_resid <= tot_num_residues_per_chain):
                    # In Chain A
                    input_flag = 1
                elif (aa_resid > tot_num_residues_per_chain and bonding_partner_resid > tot_num_residues_per_chain and bonding_partner_resid <= tot_num_residues):
                    # In Chain B
                    input_flag = 1
                else:
                    input_flag = 0

                if (input_flag == 1):
                    outfile.write(str(aa_resid) + " " + str(bonding_partner_resid) +
                                  " " + str(4) + " " + str(l0) + " " + str(K1) +" "+str(K2/K1) + "\n")
                    outfile.write(str(aa_resid) + " " + str(bonding_partner_resid) +
                                  " " + str(8) + " " + str(0) + " " + str(K3) + "\n")

    # The inter-chain (fan) interactions;

    l0_fan = [17.1, 16.35, 14.7, 13.45, 12.3, 11.3, 9.9, 9.2, 10.2, 12.5, 16.9]
    l0_fan = [float(x) * 0.1 for x in l0_fan]
    l0_fan = ['%.3f' % elem for elem in l0_fan]
    K1_fan = [0.0467, 1.324e-06, 0.085, 0.123, 0.04,
              2.92, 0.115, 0.0955, 0.1378, 0.1386, 0.3626]
    K1_fan = [float(x) * 100 * cal2j for x in K1_fan]
    #K1_fan = ['%.3f' % elem for elem in K1_fan]
    K2_fan = [0.0021, -0.0122, -0.0444, -0.04, -0.01,
              0.41, -0.041, -0.0459, -0.0527, -0.0568, -0.077]
    K2_fan = [float(x) * 1000 * cal2j for x in K2_fan]
    #K2_fan = ['%.3f' % elem for elem in K2_fan]
    K3_fan = [0.000146, 0.00185, 0.005, 0.0037, 0.0008,
              0.072, 0.0058, 0.00502, 0.005, 0.005, 0.005]
    K3_fan = [float(x) * 10000 * cal2j for x in K3_fan]
    K3_fan = ['%.3f' % elem for elem in K3_fan]

    # Reset aa_resid
    aa_resid = 0
# Bonded fan interactions ---------------------------------------------------------------------------------------

    for i in my_lt_range(0, length, 1):

        # Flag for checking whether to input the bonding / fan interactions;
        input_flag = 0

        line = lines[i].split()
        # Find the relevant lines;
        # here Fan is a symmetric interaction, so we only need to consider one DNA chain;
        if (line[0] == "ATOM" and line[4] == "A"):
            if (int(line[5]) == (aa_resid+1)):
                # When it is in the first line of the one residue,
                # update aa_resid
                aa_resid = int(line[5])
                # start recording

                # find the resid for the interchain (fan) partner

                for j in my_le_range(-4, 6, 1):

                    # Offset the index for getting parameters
                    offset = 4
                    fan_idx = j + offset

                    bonding_partner_resid = tot_num_residues + j - aa_resid
                    # Check if the pair is within range of the second chain

                    if (bonding_partner_resid > tot_num_residues_per_chain and bonding_partner_resid <= tot_num_residues):
                        input_flag = 1
                    else:
                        input_flag = 0

                    if (input_flag == 1):
                        outfile.write(str(aa_resid) + " " + str(bonding_partner_resid) + " " + str(4) + " " + str(l0_fan[fan_idx]) + " " + str(K1_fan[fan_idx])+ " "+ str(K2_fan[fan_idx]/K1_fan[fan_idx]) + "\n")
                        outfile.write(str(aa_resid) + " " + str(bonding_partner_resid) + " " + str(8) + " " + str(fan_idx+1) + " " + str(K3_fan[fan_idx]) + "\n")

    # Section separator;
    outfile.write("\n")

# angles  ---------------------------------------------------------------------------------------

    # The [ angles ] Section;
    outfile.write("[ angles ]\n")
    outfile.write(";ai  aj   ak  func  th0(deg)  C0 C1 C2 C3 C4 \n")

    # Gromacs uses the unit of deg for angles;
    theta0 = 156
    K1 = 9.22 * cal2j
    K2 = 4.16 * cal2j
    K3 = 1.078 * cal2j
    C0 = 0.0
    C1 = 0.0
    C2 = K1
    C3 = K2
    C4 = K3

    # Reset aa_resid
    aa_resid = 0

    for i in my_lt_range(0, length, 1):

        # Flag for checking whether to input the angle interactions;
        input_flag = 0

        line = lines[i].split()
        # Find the relevant lines;
        if (line[0] == "ATOM"):
            if (int(line[5]) == (aa_resid+1)):
                # When it is in the first line of the one residue,
                # update aa_resid
                aa_resid = int(line[5])
                # start recording

                # find the resid for the intrachain partner
                bonding_partner_resid_1 = aa_resid + 1
                bonding_partner_resid_2 = aa_resid + 2
                # Check if the triplet is within one chain
                if (aa_resid <= tot_num_residues_per_chain and bonding_partner_resid_1 <= tot_num_residues_per_chain and bonding_partner_resid_2 <= tot_num_residues_per_chain):
                    # In Chain A
                    input_flag = 1
                elif (aa_resid > tot_num_residues_per_chain and bonding_partner_resid_1 > tot_num_residues_per_chain and bonding_partner_resid_1 <= tot_num_residues and bonding_partner_resid_2 > tot_num_residues_per_chain and bonding_partner_resid_2 <= tot_num_residues):
                    # In Chain B
                    input_flag = 1
                else:
                    input_flag = 0

                if (input_flag == 1):
                    outfile.write(str(aa_resid) + " " + str(bonding_partner_resid_1) + " " + str(bonding_partner_resid_2) + " " + str(
                        6) + " " + str(theta0) + " " + str(C0) + " " + str(C1) + " " + str(C2) + " " + str(C3) + " " + str(C4) + "\n")

    # Section separator;
    outfile.write("\n")
# exclusions  ---------------------------------------------------------------------------------------

    # The [ exclusions ] Section;
    # Because nrexcl is 3 here, we only need to exclude contacts from those bonding pairs that are more than 3 bonds away from each other;
    outfile.write("[ exclusions ]\n")
    outfile.write("; ai	aj\n")

    # The bonding interactions, because we are going to set the nrexcl to be 0, we need to explicitly exclude those atoms out from non-bonded 
    # interactions;
    # Reset aa_resid
    aa_resid = 0

    for i in my_lt_range(0, length, 1):

        # Flag for checking whether to input the bonding / fan interactions;
        input_flag = 0

        line = lines[i].split()
        # Find the relevant lines;
        if (line[0] == "ATOM"):
            if (int(line[5]) == (aa_resid+1)):
                # When it is in the first line of the one residue,
                # update aa_resid
                aa_resid = int(line[5])
                # start recording

                # find the resid for the intrachain partner
                bonding_partner_resid = aa_resid + 1
                # Check if the pair is within one chain
                if (aa_resid <= tot_num_residues_per_chain and bonding_partner_resid <= tot_num_residues_per_chain):
                    # In Chain A
                    input_flag = 1
                elif (aa_resid > tot_num_residues_per_chain and bonding_partner_resid > tot_num_residues_per_chain and bonding_partner_resid <= tot_num_residues):
                    # In Chain B
                    input_flag = 1
                else:
                    input_flag = 0

                if (input_flag == 1):
                    outfile.write(str(aa_resid) + " " + str(bonding_partner_resid) + "\n")

    # The angle interactions, because we are going to set the nrexcl to be 0, we need to explicitly exclude those atoms out from non-bonded 
    # interactions;
    # Reset aa_resid
    aa_resid = 0

    for i in my_lt_range(0, length, 1):

        # Flag for checking whether to input the angle interactions;
        input_flag = 0

        line = lines[i].split()
        # Find the relevant lines;
        if (line[0] == "ATOM"):
            if (int(line[5]) == (aa_resid+1)):
                # When it is in the first line of the one residue,
                # update aa_resid
                aa_resid = int(line[5])
                # start recording

                # find the resid for the intrachain partner
                bonding_partner_resid_1 = aa_resid + 1
                bonding_partner_resid_2 = aa_resid + 2
                # Check if the triplet is within one chain
                if (aa_resid <= tot_num_residues_per_chain and bonding_partner_resid_1 <= tot_num_residues_per_chain and bonding_partner_resid_2 <= tot_num_residues_per_chain):
                    # In Chain A
                    input_flag = 1
                elif (aa_resid > tot_num_residues_per_chain and bonding_partner_resid_1 > tot_num_residues_per_chain and bonding_partner_resid_1 <= tot_num_residues and bonding_partner_resid_2 > tot_num_residues_per_chain and bonding_partner_resid_2 <= tot_num_residues):
                    # In Chain B
                    input_flag = 1
                else:
                    input_flag = 0

                if (input_flag == 1):
                    # pairs between aa_resid & bonding_partner_resid_1, between bonding_partner_resid_1 & bonding_partner_resid_2 have been 
                    # included by the excluded_bonding section
                    
                    outfile.write(str(aa_resid) + " " + str(bonding_partner_resid_2) + "\n")

#    # The 1-4 interactions, because we are going to set the nrexcl to be 0, we need to explicitly exclude those atoms out from non-bonded 
#    # interactions;
#    # Reset aa_resid
#    aa_resid = 0
#
#    for i in my_lt_range(0, length, 1):
#
#        # Flag for checking whether to input the angle interactions;
#        input_flag = 0
#
#        line = lines[i].split()
#        # Find the relevant lines;
#        if (line[0] == "ATOM"):
#            if (int(line[5]) == (aa_resid+1)):
#                # When it is in the first line of the one residue,
#                # update aa_resid
#                aa_resid = int(line[5])
#                # start recording
#
#                # find the resid for the intrachain partner
#                bonding_partner_resid_1 = aa_resid + 1
#                bonding_partner_resid_2 = aa_resid + 2
#                bonding_partner_resid_3 = aa_resid + 3
#                # Check if the triplet is within one chain
#                if (aa_resid <= tot_num_residues_per_chain and bonding_partner_resid_1 <= tot_num_residues_per_chain and bonding_partner_resid_2 <= tot_num_residues_per_chain and bonding_partner_resid_3 <= tot_num_residues_per_chain ):
#                    # In Chain A
#                    input_flag = 1
#                elif (aa_resid > tot_num_residues_per_chain and bonding_partner_resid_1 > tot_num_residues_per_chain and bonding_partner_resid_1 <= tot_num_residues and bonding_partner_resid_2 > tot_num_residues_per_chain and bonding_partner_resid_2 <= tot_num_residues and bonding_partner_resid_3 > tot_num_residues_per_chain and bonding_partner_resid_3 <= tot_num_residues):
#                    # In Chain B
#                    input_flag = 1
#                else:
#                    input_flag = 0
#
#                if (input_flag == 1):
#                    # pairs between aa_resid & bonding_partner_resid_1, bonding_partner_resid_1 & bonding_partner_resid_2, bonding_partner_resid2 & bonding_partner_resid_3 have been included by the excluded_bonding section
#                    # pairs between aa_resid & bonding_partner_resid_2, bonding_partner_resid_1 & bonding_partner_resid_3 have been included by the excluded_angle section
#                    outfile.write(str(aa_resid) + " " + str(bonding_partner_resid_3) + "\n")
#
    # The inter-chain (fan) interactions; We cannot exclude the fan interaction here, because we need electrostatics;

#    # Reset aa_resid
#    aa_resid = 0
#
#    for i in my_lt_range(0, length, 1):
#
#        # Flag for checking whether to input the bonding / fan interactions;
#        input_flag = 0
#
#        line = lines[i].split()
#        # Find the relevant lines;
#        # here Fan is a symmetric interaction, so we only need to consider one DNA chain;
#        if (line[0] == "ATOM" and line[4] == "A"):
#            if (int(line[5]) == (aa_resid+1)):
#                # When it is in the first line of the one residue,
#                # update aa_resid
#                aa_resid = int(line[5])
#                # start recording
#
#                # find the resid for the interchain (fan) partner
#
#                for j in my_le_range(-5, 5, 1):
#
#                    bonding_partner_resid = tot_num_residues + j - aa_resid
#
#                    # Check if the bonding partner is 3 bonds away from the first residue;
#                    if (abs(bonding_partner_resid - aa_resid) > 3):
#                        # Check if the pair is within range of the second chain
#
#                        if (bonding_partner_resid > tot_num_residues_per_chain and bonding_partner_resid <= tot_num_residues):
#                            input_flag = 1
#                        else:
#                            input_flag = 0
#
#                    if (input_flag == 1):
#                        outfile.write(str(aa_resid) + " " +
#                                      str(bonding_partner_resid) + "\n")
#
    # Section separator;
    outfile.write("\n")

    # Output the the rest of the common sections:
    outfile.write("[ system ]\n")
    outfile.write("; name\n")
    outfile.write("  Macromolecule\n")
    outfile.write("\n")

    outfile.write(" [ molecules ] \n")
    outfile.write("; name            #molec\n")
    outfile.write("  Macromolecule   1\n")

    return

############################################################################


if __name__ == "__main__":
    # Input all-atom DNA structure;
    inputFile = sys.argv[1]
    # Output CG DNA structure;
    outputFile = sys.argv[2]

    pdb2top_single_residue_dna(inputFile, outputFile)

    print("Life is given to us,")
    print("we earn it by giving it.")
