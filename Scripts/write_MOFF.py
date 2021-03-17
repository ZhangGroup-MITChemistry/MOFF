import sys
import os
import math
import numpy
import MDAnalysis as mda

# Instructions; quit and print instructions if code fails--------------------------------------------------------------
if len(sys.argv)>6 or len(sys.argv)==1:
    print("Incorrect usage. Please input proper files:")
    print("python write_MOFF.py jobid template_file smog_file stride_file pair_eps")
    print("note that the pdb should be the same as the jobid")
    print("Quiting")
    exit()



# Functions to read different parts of smog topology -----------------------------------------------------------------

# returns list (str) of index1, index2, index3, angle, spring constant
def read_angles(top_file2):
    angles_temp = []
    old = open(top_file2, 'r')
    line = old.readline()
    AA_flag = 0
    while line:
        if len(line) > 2:
            line_split = line.split()
            if AA_flag == 1:
                if line_split[0] == '[':
                    AA_flag = 0
                    pass
                elif line_split[0][0] == ';':
                    pass
                else:
                    index1 = line_split[0]
                    index2 = line_split[1]
                    index3 = line_split[2]
                    R0 = line_split[4]
                    K = line_split[5]
                    angles_temp.append([index1, index2, index3, R0, K])
            if len(line_split) > 2:
                if line_split[1] == 'angles':
                    AA_flag = 1
        line = old.readline()
    old.close()

    return angles_temp

# returns list (str) of index1, index2, index3, index4, angle, spring constant
# fix for proper dihedral
def read_dihedrals(top_file2):
    dihedrals_temp = []
    old = open(top_file2, 'r')
    line = old.readline()
    AA_flag = 0
    while line:
        if len(line) > 2:
            line_split = line.split()
            if AA_flag == 1:
                if line_split[0] == '[':
                    AA_flag = 0
                    pass
                elif line_split[0][0] == ';':
                    pass
                elif line_split[7]=='3':
                    pass
                else:
                    index1 = line_split[0]
                    index2 = line_split[1]
                    index3 = line_split[2]
                    index4 = line_split[3]
                    R0 = line_split[5]
                    K = line_split[6]
                    dihedrals_temp.append([index1, index2, index3, index4, R0, K])
            if len(line_split) > 2:
                if line_split[1] == 'dihedrals':
                    AA_flag = 1
        line = old.readline()
    old.close()

    return dihedrals_temp

# returns list (str) of index1, index2, A, B
def read_pairs(top_file2):
    pairs_temp = []
    old = open(top_file2, 'r')
    line = old.readline()
    AA_flag = 0
    while line:
        if len(line) > 2:
            line_split = line.split()
            if AA_flag == 1:
                if line_split[0] == '[':
                    AA_flag = 0
                    pass
                elif line_split[0][0] == ';':
                    pass
                else:
                    index1 = line_split[0]
                    index2 = line_split[1]
                    A = line_split[3]
                    B = line_split[4]
                    pairs_temp.append([int(index1), int(index2), float(A), float(B)])
            if len(line_split) > 2:
                if line_split[1] == 'pairs':
                    AA_flag = 1
        line = old.readline()
    old.close()
    pairs_temp=numpy.asarray(pairs_temp)
    return pairs_temp

# Functions to remove unnecessary pairs form the potential -----------------------------------------------------------
def strip_pairs(pairs_start,eps0,stridef):
    secondary = []

    stride = open(stridef, 'r')
    line = stride.readline()
    while line:
        line_split = line.split()
        if len(line_split) > 0:
            if line_split[0] == 'ASG':
                temp = line_split[6]
                index = int(line_split[4])
                temp2 = [index, temp]
                secondary.append(temp2)
        line = stride.readline()
    stride.close()

    good_pairs = []
    pairs_start[:, 2] = (pairs_start[:, 2] / (6 * eps0)) ** (1 / 10)
    pairs_start[:, 3] = (pairs_start[:, 3] / (5 * eps0)) ** (1 / 12)
    for i in range(0, len(pairs_start)):
        index1 = int(pairs_start[i][0])
        index2 = int(pairs_start[i][1])
        temp = []
        for j in range(index1 - 1, index2):
            temp.append(secondary[j][1])
        # flag detects if secondary structure changes
        flag = 0
        for j in range(0, len(temp)):
            if temp[j] != 'Strand' and temp[j] != 'AlphaHelix':
                flag = 1
            if j < len(temp) - 1:
                if temp[j] != temp[j + 1]:
                    flag = 1
        if flag == 0:
            good_pairs.append(pairs_start[i])
            if round(pairs_start[i, 2], 3) == round(pairs_start[i, 3], 3):
                pass
            else:
                print(
                    "Error!!!! Check pairs. Code is intended for the default (10-12) smog potential, with eps0=1 kj/mol. Check smog file for strength of pair potential or functional form")
    return good_pairs

# write final topolgy ----------------------------------------------------------------------------------------------
# writes topology after reading in angles, dihedrals, pairs from smog
def write_topology(calphas,angles,dihedrals,pairs,pair_eps):
    # Constants
    N=len(calphas)
    bond_l=0.38 # nm
    bond_k=1000 # kJ/mol /nm^2
    angle_k=120 # kJ/mol /deg^2
    dihedral_k=3.0 # kj/mol


    # write topology file
    new2=open(top_file,'w')
    template=open(template_file,'r')
    line=template.readline()
    while line:
        line_split=line.split()
        if len(line_split)>2:
            if line_split[1] == 'atoms':
                new2.write(line)
                line = template.readline()
                new2.write(line)
                # write all CA atoms
                for i in range(0,N):
                    index=calphas.atoms[i].index+1
                    res=calphas.atoms[i].resname
                    new2.write('\t'+str(index)+'\t'+res+'\t'+str(index)+'\t'+res+'\t'+'CA'+'\t'+str(index)+'\n')

            elif line_split[1] == 'bonds':
                new2.write(line)
                line = template.readline()
                new2.write(line)
                # Write all bonds
                for i in range(0,N-1):
                    if calphas.atoms[i].segid==calphas.atoms[i+1].segid:
                        index=calphas.atoms[i].index+1
                        next_i=calphas.atoms[i+1].index+1
                        new2.write(str(index)+'\t'+str(next_i)+'\t'+'1'+'\t'+str(bond_l)+'\t\t'+str(bond_k)+'\n')
            elif line_split[1] == 'angles':
                new2.write(line)
                line = template.readline()
                new2.write(line)
                # Write all bonds
                for i in range(0,len(angles)):
                    index1=angles[i][0]
                    index2=angles[i][1]
                    index3 = angles[i][2]
                    theta0=angles[i][3]
                    new2.write(index1+'\t'+index2+'\t'+index3+'\t1\t'+theta0+'\t'+str(angle_k)+'\n')
            elif line_split[1] == 'dihedrals':
                # returns list (str) of index1, index2, index3, index4, angle, spring constant
                new2.write(line)
                line = template.readline()
                new2.write(line)
                # Write all bonds
                for i in range(0,len(dihedrals)):
                    index1 = dihedrals[i][0]
                    index2 = dihedrals[i][1]
                    index3 = dihedrals[i][2]
                    index4 = dihedrals[i][3]
                    phi0 = float(dihedrals[i][4])
                    phi1= float(dihedrals[i][4])*3
                    new2.write(index1 + '\t' + index2 + '\t' + index3+ '\t' + index4 + '\t1\t' + str(phi0) + ' \t' + str(dihedral_k) + ' \t1\n')
                    new2.write(index1 + '\t' + index2 + '\t' + index3 + '\t' + index4 + '\t1\t' + str(phi1) + ' \t' + str(dihedral_k/2) + ' \t3\n')
            elif line_split[1] == 'pairs':
                new2.write(line)
                line = template.readline()
                new2.write(line)
                # Write all bonds
                for i in range(0,len(pairs)):
                    index1 = str(int(pairs[i][0]))
                    index2 = str(int(pairs[i][1]))
                    R0=pairs[i][2]
                    A=str( 6*pair_eps*(R0)**10 )
                    B= str(5 * pair_eps * (R0) ** 12)
                    new2.write(index1+'\t'+index2+'\t1\t'+A+'\t'+B+'\n')
            else:
                new2.write(line)
        else:
            new2.write(line)
        line=template.readline()

    new2.close()
    template.close()
    return


 # Main ---------------------------------------------------------------------------------------------------------------

 # Read in pdb name
struct_id = sys.argv[1]
if struct_id[-4:].lower()==".pdb":
        pdb = struct_id[-4:]
else:
        pdb = struct_id

# new file types
pdb_file=pdb+'.pdb'
out_pdb=pdb+'_CA.pdb'
top_file=pdb+'.top'


# Check inputs for other inputs
if len(sys.argv)>2:
    template_file=sys.argv[2]
    if template_file[-4:].lower()==".top":
        pass
    else:
        template_file = template_file+".top"
else:
    template_file = 'template_MOFF.top'

if len(sys.argv)>3:
    smog_file=sys.argv[3]
    if smog_file[-4:].lower()==".top":
        pass
    else:
        smog_file = smog_file+".top"
else:
    smog_file = 'smog.top'

if len(sys.argv)>4:
    stride_file=sys.argv[4]
else:
    stride_file = pdb+'_stride'+'.dat'

if len(sys.argv)>5:
    pair_eps2=int(sys.argv[5])
else:
    pair_eps2 = 3.0

# write pdb with only CA
u = mda.Universe(pdb_file,multiframe='no')
calphas1 = u.select_atoms("name CA")
calphas1.write(out_pdb)
# Load new pdb to fix atom indexing
new_universe = mda.Universe(out_pdb,multiframe='no')
calphas1= new_universe.select_atoms("name CA")


# read smog
angles1=read_angles(smog_file)
dihedrals1=read_dihedrals(smog_file)
pairs1=read_pairs(smog_file)
# strip pairs
pairs2=strip_pairs(pairs1,1.0,stride_file)

# Check length of angles and dihedrals:
N1=len(calphas1)
print(len(calphas1.atoms.segments))
if len(angles1)==N1-2*len(calphas1.atoms.segments):
    pass
else:
    print('WARNING: number of dihedral potentials is '+str(len(angles1))+', but there should be '+str(N1-2*len(calphas1.atoms.segments-1)))
    print('Check smog for missing or extra angles')
if len(dihedrals1)==N1-3*len(calphas1.atoms.segments):
    pass
else:
    print('WARNING: number of dihedral potentials is '+str(len(dihedrals1))+', but there should be '+str(N1-3*len(calphas1.atoms.segments)))
    print('Check smog for missing or extra dihedrals')


# write topology
write_topology(calphas1,angles1,dihedrals1,pairs2,pair_eps2)

print('Successfully wrote topology and pdb!')
print('Place in an empty box and generate tables before running a simulation')
