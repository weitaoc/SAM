clear all

% Quasi-Monte Carlo method: Sobol sequences for M parameters
M = 5;
sampling_mat = sobolset(M); % each component is between 0 and 1
sampling_mat = scramble(sampling_mat,'MatousekAffineOwen'); % apply a linear scramble to shuffle the points
N = 10000; % sample size

% generate samples in a large range:
% rw: [10^(-1) 10^2]
% Dw: [10^(-2) 10^1]
% dwn: [10^(-2) 10^1]
% dwc: [10^(-2) 10^1]
% rex: [10^(-2) 10^1]

% project the sampling from [0,1] to the ranges of parameters
sampling_mat = sampling_mat(1:N,:);
sampling_mat(:,1)   = 3*sampling_mat(:,1)-1;
sampling_mat(:,2:M) = 3*sampling_mat(:,2:M)-2;
sampling_mat = 10.^(sampling_mat); 
save('sampling_mat.mat')







