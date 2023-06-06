MinVar=["gamma"]; % Subset of ["Ref","Sigma","gamma"]

noiselevel=0.00; % set noise level
betan=0e-9; betaS=1*betan; betaG=1*betan; betag=1*betan; % regularization parameters

MaxIT=50;

Ns=36; %number of sources


% Load geometrical information on the domain
load geo-2b2

SHG(MinVar, Ns, noiselevel, MaxIT, betan, betaS, betaG, betag, geo)