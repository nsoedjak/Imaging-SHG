MinVar=["Ref","Sigma","gamma"]; % Subset of ["Ref","Sigma","gamma", "Gamma"]

noiselevel=0.00; % set noise level
betan=0e-9; betaS=1*betan; betaG=1*betan; betag=1*betan; % regularization parameters

MaxIT=10000;

iter_plot = 15; %plots intermediate reconstructions after every iter_plot iterations

Ns=36; %number of sources


% Load geometrical information on the domain
load geo-2b2

SHG(MinVar, Ns, noiselevel, MaxIT, iter_plot, betan, betaS, betaG, betag, geo)