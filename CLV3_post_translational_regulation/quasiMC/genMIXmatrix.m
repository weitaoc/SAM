% generate mixed parameter sets for Quasi-Monte Carlo 
clear all;

load('sampling1_mat.mat');
m1 = sampling_mat;
load('sampling2_mat.mat');
m2 = sampling_mat;

for i = 1:size(m1,2)
    eval(['n' num2str(i) ' = m2;']);
    eval(['n' num2str(i) '(:,' num2str(i) ') = m1(:,' num2str(i) ');']);
    filename = sprintf('N%d.mat',i);
    varname = sprintf('n%d',i);
    save(filename,varname);
end


for i = 1:size(m1,2)
    eval(['nn' num2str(i) ' = m1;']);
    eval(['nn' num2str(i) '(:,' num2str(i) ') = m2(:,' num2str(i) ');']);
    filename = sprintf('NN%d.mat',i);
    varname = sprintf('nn%d',i);
    save(filename,varname);
end