% source functions
%wus mRNA
A1  = 70;
xc1 = 0;
yc1 = 0;
zc1 = 5.5;
sigmax1 = 2;
sigmay1 = 2;
sigmaz1 = 4;
L_wus = 8.5;
% clv3
A2  = 0;%25;
xc2 = 0;
yc2 = 0;
zc2 = 9.5;
sigmax2 = 5;
sigmay2 = 5;
sigmaz2 = 5;
L_clv = 3;
%cytokinin
A3  = 1;
xc3 = 0;
yc3 = 0;
zc3 = 0;
sigmax3 = 3;
sigmay3 = 3;
sigmaz3 = 6;
%ck receptor
A4  = 1;
L_ckR = 7.5;



% diffusion coefficients
Dc = 0.1;
Dw = 3;
Dck = 1;

% protein production
rn = 0;
rc = 0.5;

ckR0 = 1;

kon = 0.5;

wthr1 = 4;
wthr2 = 6;
klow  = 2;

% degradation
dc = 0.1;
dw = 0.2;
dwn = 1;
dwc = 2.5;
dck = 0.1;
dckR = 0.1;
dCK = 0.1;

% EC50
kwc = 0.01;
kcw1 = 80;
kcw2 = 151;
kcw3 = 151;
kckw1 = 1;
kckw2 = 1;
kww = 14;

% Hill power
n = 2;

% transport between nuclei and cytoplasm
rmin = 0.8; % export
rmax = 0.8; % export
rim  = 2.4; % import

% flux rate in boundary conditions
falpha = 0.02;
