% Coefficients to reconstruct; subset of ["Ref","Sigma","gamma", "Gamma"]
MinVar=["Ref","Sigma","gamma","Gamma"];

noiselevel=0.00; % set noise level
betan=0e-7; betaS=0e-7; betag=0e-7; % regularization parameters

MaxIT=5000; % max number of iterations in optimization algorithm

iter_plot=100; %plots intermediate reconstructions after every iter_plot iterations

Ns=36; %number of sources (maximum is 36)

%Set true coefficients
%True Ref
reft = Profile(); 
reft.background = 1;
reft.rectangles = [Rectangle([0.4 1.3; 1 1.5],0.5)];

%True Sigma
sigmat = Profile(); 
sigmat.background = 1;
sigmat.rectangles = [Rectangle([0.7 1.3; 0.7 1.3],0.5)];

%True gamma
gammat = Profile();
gammat.background = 2;
gammat.rectangles = [Rectangle([1.2 1.8; 0.4 1.6],1)];

%True Gamma
Gammat = Profile();
Gammat.background = 1;
Gammat.circles = [Circle([1.5 1.0],0.3,1)];
Gammat.rectangles = [Rectangle([0.5 0.8; 0.5 0.8],2)];


% Load geometrical information on the domain
load geo-2b2

SHG(MinVar, reft, sigmat, gammat, Gammat, Ns, noiselevel, MaxIT, ...
    iter_plot, betan, betaS, betaG, betag, geo)