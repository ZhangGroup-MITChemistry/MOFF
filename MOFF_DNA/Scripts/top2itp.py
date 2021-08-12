import sys
import os
import math
import numpy
import MDAnalysis as mda

# Instructions; quit and print instructions if code fails--------------------------------------------------------------
if len(sys.argv)>2 or len(sys.argv)==1:
    print("Incorrect usage. Please input proper files:")
    print("python top2itp.py jobid")
    print("note that the top, itp, and molecule will have the same name as the jobid")
    print("Quiting")
    exit()

 # Read in pdb name
struct_id = sys.argv[1]
if struct_id[-4:].lower()==".pdb":
        pdb = struct_id[-4:]
else:
        pdb = struct_id
# new file types
top_file=pdb+'.top'
itp_file=pdb+'.itp'


old = open(top_file, 'r')
new = open(itp_file, 'w')
AA_flag=1
line=old.readline()
while line:
    if len(line) > 2:
        line_split = line.split()
        if len(line_split) > 2:
            if line_split[1] == 'moleculetype':
                AA_flag = 0
                new.write(line)
                line=old.readline()
                line_split=line.split()
                while line_split[0]==';':
                    new.write(line)
                    line=old.readline()
                    line_split = line.split()
                new.write(pdb+'\t'+line_split[1]+'\n\n')

            elif line_split[1] == 'system':
                AA_flag = 0
                new.write(line)
                line = old.readline()
                line_split = line.split()
                while line_split[0] == ';':
                    new.write(line)
                    line = old.readline()
                    line_split = line.split()
                new.write(pdb + '\n\n')

            elif line_split[1] == 'atoms' or line_split[1] == 'bonds' or line_split[1] == 'angles'  or line_split[1] == 'exclusions' or line_split[1] == 'pairs' or line_split[1] == 'dihedrals':
                AA_flag = 1
            elif line_split[1] == 'defaults' or line_split[1]=='atomtypes' or line_split[1]=='molecules':
                AA_flag = 0
        if AA_flag == 1:
            new.write(line)
        else:
            pass
    else:
        new.write(line)
    line=old.readline()
new.close()
old.close()




