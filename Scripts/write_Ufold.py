import sys
import os
import math
import numpy
import MDAnalysis as mda

# Instructions; quit and print instructions if code fails--------------------------------------------------------------
if len(sys.argv)>7 or len(sys.argv)<4:
    print("Incorrect usage. Please input proper files:")
    print("python write_Ufold.py jobid start1 end1 smog_file pair_eps")
    print("jobid - name of the pdb / topology you want to run")
    print("start1 - starting index for including the folded potential")
    print("end1 - ending index for including the folded potential")
    print("smog_file - name of smog topology file (default = smog.top)")
    print("pair_eps - strength of folded potential (default = 3 kJ/mol, this should be fit to all atom simulations)")
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
def strip_pairs_fold(pairs_start,eps0,start,end):

    good_pairs = []
    pairs_start[:, 2] = (pairs_start[:, 2] / (6 * eps0)) ** (1 / 10)
    pairs_start[:, 3] = (pairs_start[:, 3] / (5 * eps0)) ** (1 / 12)
    for i in range(0, len(pairs_start)):
        index1 = int(pairs_start[i][0])
        index2 = int(pairs_start[i][1])
        if index1>start and index1<end and index2>start and index2<end:
            good_pairs.append(pairs_start[i])

    return good_pairs

# write final topolgy ----------------------------------------------------------------------------------------------
# writes topology after reading in angles, dihedrals, pairs from smog
def write_ufold(top_file,pairs,pair_eps):
    # write topology file
    os.system('cp ' + top_file + ' temp.top')
    new2=open(top_file,'w')
    template=open('temp.top','r')
    line=template.readline()
    pair_list=[]
    exclusion_flag = 1
    while line:
        line_split=line.split()
        if len(line_split)>2:
            if line_split[1] == 'pairs':
                # pair flag indicates the pairs location has been found
                pair_flag=1
                new2.write(line)
                line = template.readline()
                new2.write(line)
                # Write all good pairs from smog
                for i in range(0,len(pairs)):
                    index1 = str(int(pairs[i][0]))
                    index2 = str(int(pairs[i][1]))
                    R0=pairs[i][2]
                    A=str( 6*pair_eps*(R0)**10 )
                    B= str(5 * pair_eps * (R0) ** 12)
                    new2.write(index1+'\t'+index2+'\t1\t'+A+'\t'+B+'\n')
                # Check if pairs repeat. Remove duplicate pairs, but keep new pairs
                while pair_flag==1:
                    line = template.readline()
                    line_split=line.split()
                    if len(line_split)==5:
                        # flag=1 if different, =0 if the pair repeats
                        flag=1
                        index1=int(line_split[0])
                        index2=int(line_split[1])
                        A=float(line_split[3])
                        B = float(line_split[4])
                        pair=[index1, index2, A, B]
                        pair_list.append(pair)
                        for i in range(0,len(pairs)):
                            indexA = int(pairs[i][0])
                            indexB = int(pairs[i][1])
                            if index1==indexA and index2==indexB:
                                flag=0
                        if flag==0:
                            pass
                        else:
                            new2.write(str(index1)+'\t'+str(index2)+'\t1\t'+str(A)+'\t'+str(B)+'\n')
                    else:
                        pair_flag=0
            # write exclusions for added pairs exclusion_flag prevents exclusions section from being added twice
            elif line_split[1] == 'exclusions' or line_split[1] == 'system' and exclusion_flag==1:
                if line_split[1] == 'exclusions':
                    new2.write(line)
                    for i in range(0, len(pairs)):
                        index1 = str(int(pairs[i][0]))
                        index2 = str(int(pairs[i][1]))
                        new2.write(index1 + '\t' + index2 + '\n')
                    exclusion_flag = 0
                elif line_split[1] == 'system':
                    new2.write('[ exclusions ]\n')

                    for i in range(0, len(pairs)):
                        index1 = str(int(pairs[i][0]))
                        index2 = str(int(pairs[i][1]))
                        new2.write(index1 + '\t' + index2 + '\n')
                    new2.write('\n')
                    new2.write(line)
                    exclusion_flag =0
            # Write other files
            else:
                new2.write(line)
        else:
            new2.write(line)
        line=template.readline()

    new2.close()
    template.close()
    os.system('rm temp.top')
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
    start1=int(sys.argv[2])-1


if len(sys.argv)>3:
    end1=int(sys.argv[3])+1

if len(sys.argv)>4:
    smog_file=sys.argv[4]
    if smog_file[-4:].lower()==".top":
        pass
    else:
        smog_file = smog_file+".top"
else:
    smog_file = 'smog.top'

if len(sys.argv)>5:
    pair_eps2=float(sys.argv[5])
else:
    pair_eps2 = 3.0


# Load new pdb to fix atom indexing
new_universe = mda.Universe(out_pdb,multiframe='no')
calphas1= new_universe.select_atoms("name CA")


# read smog
pairs1=read_pairs(smog_file)
# strip pairs
pairs2=strip_pairs_fold(pairs1,1.0,start1,end1)


# write topology
write_ufold(top_file,pairs2,pair_eps2)

print('Successfully rewrote topology!')
print('Place in an empty box and generate tables before running a simulation')
