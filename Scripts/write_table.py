import sys
import os
import math
import numpy

# Potential of the form U=q_i*q_j*exp(-kappa*r)/r+A / r^12-B 0.5*(1+tanh(eta * (rmax-x) ) )
# Expects f(r), -f'(r), g(r), -g'(r), h(r), -h'(r).
# f(r) is charge, g(r) is attractive term, h(r) is repulsive term
# dr=0.002 nm in most tables


# Instructions
if len(sys.argv)>7:
    print("Error: Incorrect usage. Please input proper files:")
    print("python write_table.py ionic_strength MOFFout SMOGout cut cut2 table_length dr\n")
    print("ionic_strength is for the debye-huckel electrostatic interactions, in mM (default=150 mM)")
    print("outputfile is the name of the table (default=table.xvg)")
    print("cut is the distance at which electrostatic interactions are cut to 0 with a fifth degree polynomial switching function (default=1.5 nm)")
    print("cut2 is the distance at which electrostatic interactions are dampened with a fifth degree polynomial switching function (default=1.2 nm)")
    print("table_length is the length of the generated table (default=15 nm)")
    print("dr is the minimal distance at which the table is calculated (default=0.002 nm)\n")
    print("Quiting")
    exit()

# deal with optional inputs
if len(sys.argv)>1:
    ionic_strength=int(sys.argv[1])
else:
    ionic_strength=150
#Convert ionic strength from mM to M
ionic_strength=ionic_strength*0.001
# Constants
k=1.380649*10**-23 # boltzmann constant
T=300 # temperature
#eps=80 # dielectric of water
eps0= 8.8541878128*10**-12 # Permativity of free space
charge=1.60217662*10**-19 # charge of an electron
No=6.023*10**23 # avagagro's number
conv=2000 # converts mol / L to mol / m^3 and takes sum over positive and negative ions
# calculate debye length




if len(sys.argv)>2:
    MOFFout=sys.argv[2]
else:
    MOFFout='table_MOFF.xvg'

if len(sys.argv) > 2:
    SMOGout = sys.argv[2]
else:
    SMOGout = 'table_smog.xvg'

if len(sys.argv)>4:
    cut=float(sys.argv[4])
else:
    cut=1.5

if len(sys.argv)>5:
    cut2=float(sys.argv[5])
else:
    cut2=1.2

if len(sys.argv)>6:
    table_length=float(sys.argv[6])
else:
    table_length=15.0

if len(sys.argv)>7:
    dr=float(sys.argv[7])
else:
    dr=0.002



N=int(table_length/dr)+1
new_mat=numpy.zeros((N,7))
new_mat2=numpy.zeros((N,7))


# parameters for attractive well in nm
eta=7 # eta is smoothness parameter
r_cut=0.8 # r_cut is the switching distance

# switching function parameters
c0=1
c1=0
c2=0
c3=-10
c4=15
c5=-6
r1=cut2
ru=cut


# electrostatic parameters
eps_0=78.4
A=-8.5525
kappa=7.7839
lam=0.003627*10
B=eps_0-A


for i in range(0,N):
    r=i*dr

    # Prevent infinities at small distances
    if r < 0.015:
        f=0
        f2=0
        g=0
        g2=0
        h=0
        h2=0


    else:
        # distance dependent dielectric constant calculation
        eps = A + B / (1 + kappa * numpy.exp(-1*lam *B*r))
        # recalculate lambda_D
        lambda_D = numpy.sqrt((eps0 * eps * k * T) / (conv * No * ionic_strength * charge ** 2))
        # convert units to nm. 1 m = 10**9 nm
        lambda_D = lambda_D * 10 ** 9
        # calculate electrostatic potential. Apply cuttoff if in range
        f_i = (1 / (eps * r)) * numpy.exp((-r / lambda_D))
        df=( numpy.exp((-r / lambda_D)) /(eps*r**2))+( numpy.exp((-r / lambda_D)) /(eps*r*lambda_D))

        if r>=r1 and r<=ru:
            S_r = c0 + c1 * ((r - r1) / (ru - r1)) + c2 * ((r - r1) / (ru - r1)) ** 2 + c3 * ((r - r1) / (ru - r1)) ** 3 + c4 * ((r - r1) / (ru - r1)) ** 4 + c5 * ((r - r1) / (ru - r1)) ** 5
            dS=c1*(1/(ru-r1))+2*c2*((r-r1) / (ru-r1))+3*c3*((r-r1) / (ru-r1))**2+4*c4*((r-r1) / (ru-r1))**3+5*c5*((r-r1) / (ru-r1))**4
            f=f_i*S_r
            f2=S_r*df+dS*f_i
        elif r<r1: # Before cutoff use normal value
            f=f_i
            f2=df
        else: # After cutoff use 0
            f=0
            f2=0



        # calculate nonbonded potential
        g= -0.5 * (1.0 + numpy.tanh( eta*(r_cut-r) ) )
        g2= -0.5 * eta / numpy.power((numpy.cosh( eta*(r_cut-r)  )),2)



        h= numpy.power((1.0/r),12)
        h2 = 12 * numpy.power((1.0 / r), 13)






    # Store values in table
    new_mat[i, 0] = r
    new_mat[i, 1] = f
    new_mat[i, 2] = f2
    new_mat[i, 3]=g
    new_mat[i,4]=g2
    new_mat[i,5]=h
    new_mat[i,6]=h2

for i in range(0, N):
    r = i * dr

    # Prevent infinities at small distances
    if r < 0.015:
        f = 0
        f2 = 0
        g = 0
        g2 = 0
        h = 0
        h2 = 0


    else:
        # distance dependent dielectric constant calculation
        eps = A + B / (1 + kappa * numpy.exp(-1 * lam * B * r))
        # recalculate lambda_D
        lambda_D = numpy.sqrt((eps0 * eps * k * T) / (conv * No * ionic_strength * charge ** 2))
        # convert units to nm. 1 m = 10**9 nm
        lambda_D = lambda_D * 10 ** 9
        # calculate electrostatic potential. Apply cuttoff if in range
        f_i = (1 / (eps * r)) * numpy.exp((-r / lambda_D))
        df = (numpy.exp((-r / lambda_D)) / (eps * r ** 2)) + (numpy.exp((-r / lambda_D)) / (eps * r * lambda_D))

        if r >= r1 and r <= ru:
            S_r = c0 + c1 * ((r - r1) / (ru - r1)) + c2 * ((r - r1) / (ru - r1)) ** 2 + c3 * ((r - r1) / (ru - r1)) ** 3 + c4 * ((r - r1) / (ru - r1)) ** 4 + c5 * ((r - r1) / (ru - r1)) ** 5
            dS = c1 * (1 / (ru - r1)) + 2 * c2 * ((r - r1) / (ru - r1)) + 3 * c3 * ((r - r1) / (ru - r1)) ** 2 + 4 * c4 * ((r - r1) / (ru - r1)) ** 3 + 5 * c5 * ((r - r1) / (ru - r1)) ** 4
            f = f_i * S_r
            f2 = S_r * df + dS * f_i
        elif r < r1:  # Before cutoff use normal value
            f = f_i
            f2 = df
        else:  # After cutoff use 0
            f = 0
            f2 = 0

        # calculate nonbonded potential
        g = -numpy.power((1.0 / r), 10)
        g2 = -10 * numpy.power((1.0 / r), 11)

        h = numpy.power((1.0 / r), 12)
        h2 = 12 * numpy.power((1.0 / r), 13)

    # Store values in table
    new_mat2[i, 0] = r
    new_mat2[i, 1] = f
    new_mat2[i, 2] = f2
    new_mat2[i, 3] = g
    new_mat2[i, 4] = g2
    new_mat2[i, 5] = h
    new_mat2[i, 6] = h2



# save new matrix to file
numpy.savetxt(MOFFout,new_mat)
numpy.savetxt(SMOGout,new_mat2)
print('Successfully wrote tables!')
