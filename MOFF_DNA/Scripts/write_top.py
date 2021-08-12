import sys
import os
import math
import numpy

# Instructions; quit and print instructions if code fails--------------------------------------------------------------
if len(sys.argv)%2==1:
    print("Incorrect usage. Please input proper files:")
    print("python write_top.py topol.top mol1 num1 mol2 num2 ...")
    print("note that topol.top is the name of the topology file")
    print("mol1 is the name of the first molecule. Must be the same name as the .itp file")
    print("num1 is the number of mol1 in the simulation (must be integer). Remember to add the same number of mol1 to the structure using GROMAC's insert-molecule command")
    print("mol2 is the name of the second molecule. Must be the same name as the .itp file")
    print("num2 is the number of mol2 in the simulation (must be integer). Remember to add the same number of mol2 to the structure using GROMAC's insert-molecule command")
    print("More than 2 molecules can be added following the same formula")
    print("Quiting")
    exit()

 # Read in pdb name
struct_id = sys.argv[1]
if struct_id[-4:].lower()==".top":
        top = struct_id[:-4]
else:
        top = struct_id
top2=top+'.top'

mol=[]
nmol=[]
for i in range(2,len(sys.argv),2):
    mol_temp = sys.argv[i]
    if struct_id[-4:].lower() == ".itp":
        mol.append(mol_temp[-4:])
    else:
        mol.append(mol_temp)
    nmol_temp=str(int(sys.argv[i+1]))
    nmol.append(nmol_temp)

new=open(top2,'w')
print(top2)
new.write('#include \"prot_dna.itp\"\n')
for i in range(0,len(mol)):
    new.write('#include \"'+mol[i]+'.itp\"\n')
new.write('\n[ molecules ]\n')
new.write('; name            #molec\n')
for i in range(0,len(mol)):
    new.write(mol[i]+'\t'+nmol[i]+'\n')
new.close()

