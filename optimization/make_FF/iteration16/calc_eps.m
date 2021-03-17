clear all;
close all;

% N- number of structures for each protein
N=3002;
% n- number of energies to calculate (210=upper diagonal of 20*20)
n=210;

% import old data
old=importdata('iter15.dat');
old2=zeros(n,1);

% fill out full matrix
count=1;
for i=1:20
    for j=1:20
        if i<=j
            old(j,i)=old(i,j);
            old2(count)=old(i,j);
            count=count+1;
        end
    end
end

% old is positive is attractive. Make negative attractive
old2=-old2;


% Import data from prepare_opt.py
contact_list=importdata('contact_list.txt');
bias_tot=importdata('bias_tot.txt');
pdb_list=importdata('pdb_list.txt');
pdb_list_contacts=importdata('pdb_list_contacts.txt');


E_pdb2=pdb_list*old2;
E_sim2=pdb_list_contacts*old2;

save('E_pdb_old.txt','E_pdb2','-ascii');
save('E_sim_old.txt','E_sim2','-ascii');

% Subtract contacts from pdb contacts
pdb_diff=pdb_list-pdb_list_contacts;
% energy of pdb sturcture
min_E=-pdb_diff*old2;




% std_dev of MJ matrix, for lower and upper bounds of tolerance
std_dev=0.2439;
lb=ones(n,1)*2*-std_dev;
ub=ones(n,1)*2*std_dev;

% make matrices sparse
contact_list=sparse(contact_list);
pdb_diff=sparse(pdb_diff);


% Calculate std dev of pdb E for tolerance
sigma_pdb=zeros(7,1);

% sigma_pdb has 7 values. Each is the std deviation of contact energies of the simulated structures
pdb_E=pdb_list_contacts*old2;
for i=1:7
    sigma_pdb(i)=std(pdb_E(((i-1)*N+1) : (i*N)));    
end


% Add gamma*sigma_pdb to add in the error tolerance
gamma=1.5;
for i=1:7
    min_E(((i-1)*N+1) : (i*N))=min_E(((i-1)*N+1) : (i*N))+gamma*sigma_pdb(i);
end

% fit with constrained least squares
deps=lsqlin(contact_list,bias_tot,pdb_diff,min_E,[],[],lb,ub,[]);

save('delta_eps.txt','deps','-ascii');

% add new energy to old energy
eps=old2+deps;

% Calculate pdb and contact energies to check convergence
E_pdb=pdb_list*eps;
E_sim=pdb_list_contacts*eps;

save('E_pdb.txt','E_pdb','-ascii');
save('E_sim.txt','E_sim','-ascii');

