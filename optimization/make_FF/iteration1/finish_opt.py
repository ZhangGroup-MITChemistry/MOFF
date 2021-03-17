import sys
import os
import math
import numpy

AA=['CYS',    'MET',  'PHE',  'ILE',  'LEU',  'VAL',  'TRP',  'TYR',  'ALA',  'GLY',  'THR',  'SER',  'ASN',  'GLN',  'ASP',  'GLU',  'HIS',  'ARG',  'LYS',  'PRO']
type=20

eps=numpy.loadtxt('delta_eps.txt')

# default assumes positive epsilon correspond lower energy
eps=-1*eps



# store back into 20*20 matrix
d_eps=numpy.zeros((type,type))
counter=0
for i in range(0,type*type):
    indexAA = int(i / type)
    indexBB = (i) % type
    if indexAA <= indexBB:
        d_eps[indexAA][indexBB]=eps[counter]
        counter = counter + 1
# add new potential to old potential, remove nonsensical values of eps
old_eps=numpy.loadtxt('iter15.dat')

# debug
#numpy.savetxt('eps.txt',eps)
#numpy.savetxt('d_eps1.txt',d_eps)

new_eps=old_eps+d_eps


numpy.savetxt('new_eps.dat',new_eps)




sigma=numpy.loadtxt('sigma.dat')

new_file='new_ff.dat'
new=open(new_file,'w')
new.write('[ nonbond_params ]\n')
new.write(';i\tj\tfunc\tC\t\tA\n')


dr=0.0001
ndr=1/dr
r=numpy.linspace(dr,1,int(1/dr))


k=len(AA)
for i in range(0,k):
    for j in range(i,k):
        if new_eps[i][j]>0:
            s=sigma[i][j]
            e=new_eps[i][j]
            eta=7
            r_0=0.8
            f=e*(s/r)**12-(e/2)*(1+numpy.tanh(eta*(r_0-r)))
            minimum = min(f)
            eps_2=e*(e/(abs(minimum)))
            f2=eps_2*(s/r)**12-(eps_2/2)*(1+numpy.tanh(eta*(r_0-r)))
            print(min(f2))
            repulsive=eps_2*s**12
            to_write=AA[i]+'\t'+AA[j]+'\t1\t'+str(eps_2)+'\t'+str(repulsive)+'\n'
            new.write(to_write)
        if new_eps[i][j]<0:
            s = sigma[i][j]
            e = abs(new_eps[i][j])
            eta = 7
            r_0 = 0.8
            n2 = int(ndr * s) - 1
            f = e * (s / r) ** 12 + (e / 2) * (1 + numpy.tanh(eta * (r_0 - r)))
            eps_2 = e * e / f[n2]
            repulsive = eps_2 * s ** 12
            to_write = AA[i] + '\t' + AA[j] + '\t1\t' + str(-1*eps_2) + '\t' + str(repulsive) + '\n'
            new.write(to_write)
new.close()
