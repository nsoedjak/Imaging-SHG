MinVar=["Gamma"]; % Subset of ["Ref","Sigma","Gamma","gamma"]

noiselevel=0.00; % set noise level
betan=0e-9; betaS=1*betan; betaG=1*betan; betag=1*betan; % regularization parameters

MaxIT=7;

Ns=1; %number of sources

SHG(MinVar, Ns, noiselevel, MaxIT, betan, betaS, betaG, betag)