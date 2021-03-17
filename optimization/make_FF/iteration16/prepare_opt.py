import sys
import os
import math
import numpy
from scipy import linalg

# Function to calculate distance
def distance(pos1,pos2):
    diff=[0,0,0]
    diff[0]=pos1[0]-pos2[0]
    diff[1] = pos1[1] - pos2[1]
    diff[2]=pos1[2]-pos2[2]
    dist=math.sqrt(diff[0]**2+diff[1]**2+diff[2]**2)
    return dist

# Function reads in PDB. Stores 3 elements for each atom: type[0], index[1], position[2]. Atoms are stored in array by timestep and index
def read_pdb(pdb_file):
    pdb=open(pdb_file,'r')
    line=pdb.readline()
    atom_list=[]
    while line:
        line_split=line.split()
        if len(line_split) > 0:
            if line_split[0]=='MODEL':
                atom_flag = 1
                atom_list_timestep = []
                while atom_flag == 1:
                    line = pdb.readline()
                    line_split=line.split()
                    if line_split[0]=='ENDMDL':
                        atom_list.append(atom_list_timestep)
                        atom_flag = 0
                        pass
                    elif line_split[0]=='ATOM':
                        type=line_split[3]
                        index=int(line_split[1])
                        # Get postion
                        pos = [0, 0, 0]
                        pos[0] = float(line[30:37])
                        pos[1] = float(line[38:45])
                        pos[2] = float(line[46:53])
                        #pos[0] = float(line_split[5])
                        #pos[1] = float(line_split[6])
                        #pos[2] = float(line_split[7])
                        # Always add to toal atom list
                        atom=[type,index,pos]
                        atom_list_timestep.append(atom)
                    else:
                        pass
        line = pdb.readline()
    pdb.close()
    return atom_list


# Calculates contacts at for a trajectory and stores it in an array according to amino acids. Cutoff is distance cutoff (in Angrstoms) to be considered a contact
def calc_contacts(trj,cutoff):
    skip=4 # number of neighbors to skip along the chain and exclude from contact calculation
    eta=0.7
    r_cut=8
    timesteps=len(trj)
    N=len(trj[0])
    types=20
    interactions=int(types*(types+1)/2)
    tot_contacts=numpy.zeros((timesteps,interactions))
    AA_conv={'CYS':0,'MET':1,'PHE':2,'ILE':3,'LEU':4,'VAL':5,'TRP':6,'TYR':7,'ALA':8,'GLY':9,'THR':10,'SER':11,'ASN':12,'GLN':13,'ASP':14,'GLU':15,'HIS':16,'ARG':17,'LYS':18,'PRO':19}
    for i in range(0,timesteps):
        contacts = numpy.zeros((types, types))
        for j in range(0,N-skip):
            for k in range(j+skip,N):
                posA=trj[i][j][2]
                posB=trj[i][k][2]
                r_abs=distance(posA, posB)
                if r_abs < cutoff:
                    typeA=trj[i][j][0]
                    typeB=trj[i][k][0]
                    indexA=AA_conv[typeA]
                    indexB=AA_conv[typeB]
                    # Debug
                    # if typeA==typeB and typeA=='PRO':
                    #    print(j)
                    #    print(posA)
                    #    print(k)
                    #    print(posB)
                    if indexA<indexB:
                        contacts[indexA][indexB]=contacts[indexA][indexB]+ 0.5 * (1.0 + numpy.tanh( eta*(r_cut-r_abs) ) )
                    else:
                        contacts[indexB][indexA] = contacts[indexB][indexA] + 0.5 * (1.0 + numpy.tanh( eta*(r_cut-r_abs) ) )
        counter=0
        #numpy.savetxt('contacts.txt',contacts,'%.1e')
        for a in range(0,types*types):
            indexAA=int(a/types)
            indexBB=(a)%types
            if indexAA <= indexBB:
                tot_contacts[i][counter]=contacts[indexAA][indexBB]
                counter=counter+1
    #numpy.savetxt('tot_contacts.txt',tot_contacts.transpose(),'%.1e')
    return tot_contacts

# Presets and parameters
AA=['CYS',    'MET',  'PHE',  'ILE',  'LEU',  'VAL',  'TRP',  'TYR',  'ALA',  'GLY',  'THR',  'SER',  'ASN',  'GLN',  'ASP',  'GLU',  'HIS',  'ARG',  'LYS',  'PRO']
type=20
n=int(type*(type+1)/2)
N=len(sys.argv)-1
timesteps=1501
ordered=7
M=timesteps*N
m=timesteps
bias_tot=numpy.zeros((2*M,1))
contact_list=numpy.zeros((2*M,n))
M_ordered=timesteps*2*ordered
pdb_list_contacts=numpy.zeros((M_ordered,n))
pdb_list=numpy.zeros((M_ordered,n))

# 1soy  1ubq  1wla  3mzq  5tvz  6eez  6h8m  ACTR  An16 asynuclein  ERMTADn hNHE1cdt IBB N49 N98 NLS NSP NUL NUS P53 ProTa sh4ud Sic
rg_exp=[1.530,1.311,1.650,1.610,1.824,1.884,1.535,2.51,4.44,3.31,3.96,3.63,3.20,1.59,2.86,2.40,4.10,3.00,2.49,2.87,3.79,2.90,3.21]
alpha=[3.5,18,7,6.5,6.5,2.5,10,-5.5,-1.5,-1.5,-6.4,-1.5,-2.3,3,1,-13,0.3,-1.5,-0.6,-2.7,-3,-3,-3.6]

for i in range(1,N+1):
    # timesteps for this range of values
    start1=(2*i-2)*timesteps
    finish1=(2*i-1)*timesteps
    start2 = (2 * i - 1) * timesteps
    finish2 = (2 * i ) * timesteps

    # import PDB, get contacts
    pdb_id = sys.argv[i]
    print(pdb_id)
    new_pdb = 'bias_dat/'+pdb_id + '_bias.pdb'
    new_trj = read_pdb(new_pdb)
    new_contacts = calc_contacts(new_trj, 12.0)
    contact_list[start1:finish1,:]=new_contacts
    # Get the strength of the bias (E), and Rg for the AA used
    rg_file='bias_dat/'+pdb_id+'_bias.xvg'
    rg_dat=numpy.loadtxt(rg_file)
    rg=numpy.zeros((timesteps,1))
    bias=numpy.zeros((timesteps,1))
    # deal with different beginning of simulation
    skip=1
    begin=0
    for j in range(0,timesteps):
        #rg[j]=rg_dat[5*skip*j+begin][1]
        #bias[j]=alpha[i-1] * (rg_dat[5*skip*j+begin][1]-rg_exp[i-1])
        rg[j]=rg_dat[j+begin][1]
        bias[j]=alpha[i-1] * (rg_dat[j+begin][1]-rg_exp[i-1])
    #bias[j] = alpha[i - 1] * (rg_dat[5 * skip * j + begin][1] - rg_max[i - 1])
    bias_tot[start1:finish1]=bias

    # Normalize number of contacts by average number of contacts with similar Rg
    Rg_min = rg_exp[i - 1] - 0.05
    Rg_max = rg_exp[i - 1] + 0.05
    indeces = []
    for j in range(0, len(rg)):
        if rg[j] > Rg_min and rg[j] < Rg_max:
            indeces.append(j)
    N_rg = len(indeces)
    if N_rg > 0:
        pass
    else:
        print('Error!!! No Rg within 0.05 of the experimental Rg!!!')
    print(N_rg)
    Rg_contacts = numpy.zeros((N_rg, n))
    for j in range(0, len(indeces)):
        index_rg = indeces[j]
        Rg_contacts[j, :] = new_contacts[index_rg, :]
    ave_contacts = numpy.mean(Rg_contacts, axis=0)


    # Get values for unbiased simulation
    old_pdb='nb_dat/'+pdb_id+ '_unbiased.pdb'
    old_trj=read_pdb(old_pdb)
    old_contacts=calc_contacts(old_trj,12.0)
    contact_list[start2:finish2,:]=old_contacts
    # Get strength of theoretical bias
    rg_file='nb_dat/'+pdb_id+'_unbiased.xvg'
    rg_dat = numpy.loadtxt(rg_file)
    begin=0
    rg = numpy.zeros((timesteps, 1))
    bias = numpy.zeros((timesteps, 1))
    for j in range(0,timesteps):
        rg[j]=rg_dat[j+begin][1]
        bias[j]=alpha[i-1] * (rg_dat[j+begin][1]-rg_exp[i-1])
        #bias[j] = alpha[i - 1] * (rg_dat[j + begin][1] - rg_max[i - 1])
    bias_tot[start2:finish2] = bias

    if i<ordered+1:
        pdb_pdb = 'pdb_dat/' + pdb_id + '_pdb.pdb'
        pdb_trj = read_pdb(pdb_pdb)
        pdb_contacts = calc_contacts(pdb_trj, 12.0)
        for j in range(start1, finish2):
            pdb_list[j, :] = pdb_contacts[0, :]
            pdb_list_contacts[j,:]=contact_list[j, :]
            # subtract pdb contacts from ordered proteins
            contact_list[j, :] = contact_list[j, :] - ave_contacts
    # Subtract average contacts from IDPs
    else:
        for j in range(start1, finish2):
            contact_list[j, :] = contact_list[j, :] - ave_contacts


# Calculate SVD
U,D1,V=linalg.svd(contact_list)
numpy.savetxt('D1.txt',D1)

# parameters to remove noise from SVD
tot=numpy.sum(numpy.square(D1))
cutEig=0.95
sum=0
count=0
# cutEig based on sum variance of vector
for i in range(0,n):
    if sum<cutEig:
        sum=sum+(D1[i]**2/tot)
    else:
        D1[i]=0
        print('Cutting!')
        count=count+1
print(count)

# remake matrix and save
D=numpy.zeros((2*M,n))
numpy.fill_diagonal(D,D1)
numpy.savetxt('D1.txt',D1)
numpy.savetxt('D.txt',D)
contact_list2=numpy.dot(U,numpy.dot(D,V))


# Needed for next step: C - contact list, d - bias_tot, A - pdb_list - pdb_list_contacts, b - 0
# Solving C*x=d where Ax<b
numpy.savetxt('contact_list.txt',contact_list2)
numpy.savetxt('bias_tot.txt',bias_tot)
numpy.savetxt('pdb_list.txt',pdb_list)
numpy.savetxt('pdb_list_contacts.txt',pdb_list_contacts)



