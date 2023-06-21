%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SHG: Second harmonic generation image reconstruction with least-squares 
%       minimization method
%
% Authors: Kui Ren  (kr2002@columbia.edu), Nathan Soedjak (ns3572@columbia.edu)
% Institution: Department of Applied Physics and Applied Mathematics, Columbia University
% Last update: 2023-2-18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mathematical model:
% \Delta u + k^2 (1+n) u +i k \sigma u =0
% u = g
% 
% \Delta v + (2k)^2 (1+n) v +i 2k \sigma v = -(2k)^2 \gamma u^2
% v + i2k n*grad v = 0
%
% Internal data: H = \Gamma\sigma (|u|^2 + |v|^2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear all;  close all;
% 

function SHG(MinVar, reft, sigmat, gammat, Gammat, Ns, noiselevel, ...
    MaxIT, iter_plot, betan, betaS, betag, geo)

tic; tb=toc;

% Create the discretized domain (0,2)x(0,2)
dx=0.025; x=0:dx:2;
dy=0.025; y=0:dy:2;
Nx=length(x);
Ny=length(y);
% [X,Y]=meshgrid(x,y);

% Generate regular finite element mesh on rectangular geometry
[P,E,T]=poimesh(geo,Nx-1,Ny-1); 
% figure;
% pdemesh(P,E,T);
% axis tight; axis square; box on; axis off;
% title('Finite element mesh');

M=Nx*Ny; % total number of nodes in the mesh

wnum=1; % wave number k



% Set up initial guesses
ref0=reft.background*ones(M,1);
sigma0=sigmat.background*ones(M,1);
gamma0=gammat.background*ones(M,1);


%Evaluate the true coefficients at the grid points
reft = reft.evaluate(P);
sigmat = sigmat.evaluate(P);
gammat = gammat.evaluate(P);
Gammat = Gammat.evaluate(P);


% Set up initial guesses only for coefficients we want to reconstruct
if ~ismember("Ref",MinVar)
    ref0=reft;
end

if ~ismember("gamma",MinVar)
    gamma0=gammat;
end

if ~ismember("Sigma",MinVar)
    sigma0=sigmat;
end

X0=[ref0' sigma0' gamma0']';



reft_min = min(reft);
reft_max = max(reft);
sigmat_min = min(sigmat);
sigmat_max = max(sigmat);
gammat_min = min(gammat);
gammat_max = max(gammat);
Gammat_min = min(Gammat);
Gammat_max = max(Gammat);


%Plot initial guesses
if ismember("Ref",MinVar)
    ref0g=tri2grid(P,T,ref0,x,y);
    figure;
    ph = pcolor(x,y,ref0g); axis tight; colorbar('SouthOutside');
    axis square; axis off; shading interp;
    ph.ZData = ph.CData;
    caxis([reft_min reft_max]);
    title('initial guess of \eta');
    drawnow;
end
if ismember("Sigma",MinVar)
    sigma0g=tri2grid(P,T,sigma0,x,y);
    figure;
    ph = pcolor(x,y,sigma0g); axis tight; colorbar('SouthOutside');
    axis square; axis off; shading interp;
    ph.ZData = ph.CData;
    caxis([sigmat_min sigmat_max]);
    title('initial guess of \sigma');
    drawnow;
end
if ismember("gamma",MinVar)
    gamma0g=tri2grid(P,T,gamma0,x,y);
    figure;
    ph = pcolor(x,y,gamma0g); axis tight; colorbar('SouthOutside');
    axis square; axis off; shading interp;
    ph.ZData = ph.CData;
    caxis([gammat_min gammat_max]);
    title('initial guess of \gamma');
    drawnow;
end

% Setup information on sources and boundary edges
SrcInfo=SetSources(Ns);
BdaryInfo=SetBdaryInfo(P,E);


%Plot true coefficients
if ismember("Ref", MinVar)
    reftg=tri2grid(P,T,reft,x,y);
    figure;
    ph = pcolor(x,y,reftg); axis tight; colorbar('SouthOutside');
    axis square; axis off; shading interp;
    ph.ZData = ph.CData;
    caxis([reft_min reft_max]);
    title('true \eta');
    drawnow;
end
if ismember("Sigma", MinVar)
    sigmatg=tri2grid(P,T,sigmat,x,y);
    figure;
    ph = pcolor(x,y,sigmatg); axis tight; colorbar('SouthOutside');
    axis square; axis off; shading interp;
    ph.ZData = ph.CData;
    caxis([sigmat_min sigmat_max]);
    title('true \sigma');
    drawnow;
end
if ismember("gamma", MinVar)
    gammatg=tri2grid(P,T,gammat,x,y);
    figure;
    ph = pcolor(x,y,gammatg); axis tight; colorbar('SouthOutside');
    axis square; axis off; shading interp;
    ph.ZData = ph.CData;
    caxis([gammat_min gammat_max]);
    title('true \gamma');
    drawnow;
end
if ismember("Gamma", MinVar)
    Gammatg=tri2grid(P,T,Gammat,x,y);
    figure;
    ph = pcolor(x,y,Gammatg); axis tight; colorbar('SouthOutside');
    axis square; axis off; shading interp;
    ph.ZData = ph.CData;
    caxis([Gammat_min Gammat_max]);
    title('true \Gamma');
    drawnow;
end


% Generating synthetic data
disp(' ');
disp(' ');
disp('Generating synthetic data .......');
disp(' ');

srczero=zeros(M,1);
Hm=zeros(M,Ns);
for ks=1:Ns
    
    % Solve the Helmholtz equations
    ut=HelmholtzSolve('u_Forward',SrcInfo,BdaryInfo,ks,P,E,T,wnum,reft,sigmat,srczero);

    srcv = -(2*wnum)^2 * gammat .* ut.^2;
    vt=HelmholtzSolve('Homogeneous_Robin',SrcInfo,BdaryInfo,ks,P,E,T,2*wnum,reft,sigmat,srcv);

    Ht=Gammat.*sigmat.*(abs(ut).^2 + abs(vt).^2);


    %Plot data/solutions
    %utg=tri2grid(P,T,vt,x,y);
    %figure;
    %ph = contourf(x,y,angle(utg)); axis tight; colorbar('SouthOutside');
    %ph = pcolor(x,y,real(utg)); axis tight; colorbar('SouthOutside');
    %axis square; axis off; shading interp;
    %ph.ZData = ph.CData;
    %drawnow;
    %pause
   
    % Add noise to data
	Hm(:,ks)=Ht.*(1+noiselevel*2*(rand(M,1)-0.5));
    
%   disp(['Synthetic data generated for source #: ' num2str(ks)]);
%   disp('  ');

    clear ut Ht;
end
disp('Finished generating synthetic data .......');


% This short part is only for debugging
%[f0 g0]=qTATObj(X0,MinVar,x,y,dx,dy,Nx,Ny,P,E,T,Ns,Hm,SrcInfo,BdaryInfo,wnum);
%g0g=tri2grid(P,T,g0(1:M),x,y);
%g0g=tri2grid(P,T,g0(M+1:2*M),x,y);
%figure;
%pcolor(x,y,g0g); axis tight; colorbar('SouthOutside');
%axis square; shading interp;
%title('Gradient');
%drawnow;

OptimMethod='UNCON';

% Setup the minimization algorithm
disp(' ');
disp(' ');
disp('Minimizing objective function .......');
disp(' ');


% Set up shared variables with outfun
% history_x = [];
% history_fval = [];
% history_iter = [];

% This function is used to plot the intermediate reconstruction results
function stop = outfun(x_,optimValues,state)
     stop = false;
 
     switch state
         case 'iter'
             if rem(optimValues.iteration, iter_plot) == 0 && optimValues.iteration ~= 0
                 % Concatenate current point and objective function
                 % value with history. x must be a row vector.
%                  history_iter = [history_iter; optimValues.iteration];
%                  history_fval = [history_fval; optimValues.fval];
%                  history_x = [history_x; x];
                 refr=x_(1:M);
                 sigmar=x_(M+1:2*M);
                 gammar=x_(2*M+1:3*M);
                 % Plot intermediate reconstruction results
                 if ismember("Ref",MinVar)
                     refrg=tri2grid(P,T,refr,x,y);
                     figure;
                     ph = pcolor(x,y,refrg); axis tight; colorbar('SouthOutside');
                     axis square; axis off; shading interp;
                     ph.ZData = ph.CData;
                     caxis([reft_min reft_max]);
                     title(sprintf('recovered \\eta after %.0f iterations', optimValues.iteration));
                     drawnow;
                end
                if ismember("Sigma",MinVar)
                     sigmarg=tri2grid(P,T,sigmar,x,y);
                     figure;
                     ph = pcolor(x,y,sigmarg); axis tight; colorbar('SouthOutside');
                     axis square; axis off; shading interp;
                     ph.ZData = ph.CData;
                     caxis([sigmat_min sigmat_max]);
                     title(sprintf('recovered \\sigma after %.0f iterations', optimValues.iteration));
                     drawnow;
                end
                if ismember("gamma",MinVar)
                     gammarg=tri2grid(P,T,gammar,x,y);
                     figure;
                     ph = pcolor(x,y,gammarg); axis tight; colorbar('SouthOutside');
                     axis square; axis off; shading interp;
                     ph.ZData = ph.CData;
                     caxis([gammat_min gammat_max]);
                     title(sprintf('recovered \\gamma after %.0f iterations', optimValues.iteration));
                     drawnow;
                end
                %Reconstruct and plot Gamma
                if ismember("Gamma",MinVar)
                    Gammar = zeros(M,1);
                    for ks=1:Ns
                        % Solve the Helmholtz equations
                        ur=HelmholtzSolve('u_Forward',SrcInfo,BdaryInfo,ks,P,E,T,wnum,refr,sigmar,srczero);
                    
                        srcv = -(2*wnum)^2 * gammar .* ur.^2;
                        vr=HelmholtzSolve('Homogeneous_Robin',SrcInfo,BdaryInfo,ks,P,E,T,2*wnum,refr,sigmar,srcv);
                    
                        Gammar = Gammar + Hm(:,ks)./(sigmar.*(abs(ur).^2 + abs(vr).^2));
                    end
                    Gammar = Gammar./Ns;
                
                    %Plot Gamma
                    Gammarg=tri2grid(P,T,Gammar,x,y);
                    figure;
                    ph = pcolor(x,y,Gammarg); axis tight; colorbar('SouthOutside');
                    axis square; axis off; shading interp;
                    ph.ZData = ph.CData;
                    caxis([Gammat_min Gammat_max]);
                    title(sprintf('recovered \\Gamma after %.0f iterations', optimValues.iteration));
                    drawnow;
                end
             end
         otherwise
     end
end

f=@(X) SHGObj(X,Gammat,MinVar,x,y,dx,dy,Nx,Ny,P,E,T,Ns,Hm,...
    SrcInfo,BdaryInfo,wnum,betan,betaS,betag);

if strcmp(OptimMethod,'UNCON')
    options=optimoptions(@fminunc,'OutputFcn',@outfun,'Algorithm','quasi-newton', ...
    'Display','iter-detailed','GradObj','on','TolFun',1e-12,...
    'MaxIter',MaxIT);
    [X,fval,exitflag,output,grad]=fminunc(f,X0,options);
else
    % Set inequality constraint
    Aieq=zeros(1,2*M);
    Bieq=0;
    % Set equality constraint
    Aeq=zeros(1,2*M);
    Beq=0;
    % Set upper and lower bounds
    LB=[0.1*ones(1,M) 0.1*ones(1,M)]';
    UB=[0.9*ones(1,M) 0.4*ones(1,M)]';

    options=optimoptions(@fmincon,'Algorithm','trust-region', ...
    'Display','iter-detailed','GradObj','on','TolFun',1e-12,...
    'MaxIter',MaxIT);
    %options=optimset('Display','iter-detailed','GradObj','on','TolFun',1e-12,'MaxIter',MaxIT);
    %options = optimset('algorithm','sqp','maxfunevals',5000,'maxiter',100);
    %options = optimset(options,'tolx',1e-9,'tolcon',1e-9,'tolfun',1e-6);
    %options = optimset(options,'GradObj','on','GradConstr','off');
    
    [X,fval,exitflag,output,lambda]=fmincon(f,X0,Aieq,Bieq,Aeq,Beq,LB,UB,[],options);
    %[X,fval,exitflag,output]=fmincon(f,X0,zeros(M,M),zeros(M,1),[],[],LB,UB);
end


disp(' ');
disp(' ');
disp('Finished minimizing objective function .......');

disp(' ');
disp(' ');
disp('Plotting final results .......');
disp(' ');

refr=X(1:M);
sigmar=X(M+1:2*M);
gammar=X(2*M+1:3*M);
Gammar=Gammat;
% Plot reconstruction results
if ismember("Ref",MinVar)
    refrg=tri2grid(P,T,refr,x,y);
    figure;
    ph = pcolor(x,y,refrg); axis tight; colorbar('SouthOutside');
    axis square; axis off; shading interp;
    ph.ZData = ph.CData;
    caxis([reft_min reft_max]);
    title(sprintf('recovered \\eta after %.0f iterations', output.iterations));
    drawnow;
end
if ismember("Sigma",MinVar)
    sigmarg=tri2grid(P,T,sigmar,x,y);
    figure;
    ph = pcolor(x,y,sigmarg); axis tight; colorbar('SouthOutside');
    axis square; axis off; shading interp;
    ph.ZData = ph.CData;
    caxis([sigmat_min sigmat_max]);
    title(sprintf('recovered \\sigma after %.0f iterations', output.iterations));
    drawnow;
end
if ismember("gamma",MinVar)
    gammarg=tri2grid(P,T,gammar,x,y);
    figure;
    ph = pcolor(x,y,gammarg); axis tight; colorbar('SouthOutside');
    axis square; axis off; shading interp;
    ph.ZData = ph.CData;
    caxis([gammat_min gammat_max]);
    title(sprintf('recovered \\gamma after %.0f iterations', output.iterations));
    drawnow;
end

%Reconstruct and plot Gamma
if ismember("Gamma",MinVar)
    Gammar = zeros(M,1);
    for ks=1:Ns
        % Solve the Helmholtz equations
        ur=HelmholtzSolve('u_Forward',SrcInfo,BdaryInfo,ks,P,E,T,wnum,refr,sigmar,srczero);
    
        srcv = -(2*wnum)^2 * gammar .* ur.^2;
        vr=HelmholtzSolve('Homogeneous_Robin',SrcInfo,BdaryInfo,ks,P,E,T,2*wnum,refr,sigmar,srcv);
    
        Gammar = Gammar + Hm(:,ks)./(sigmar.*(abs(ur).^2 + abs(vr).^2));
    end
    Gammar = Gammar./Ns;

    %Plot Gamma
    Gammarg=tri2grid(P,T,Gammar,x,y);
    figure;
    ph = pcolor(x,y,Gammarg); axis tight; colorbar('SouthOutside');
    axis square; axis off; shading interp;
    ph.ZData = ph.CData;
    caxis([Gammat_min Gammat_max]);
    title(sprintf('recovered \\Gamma after %.0f iterations', output.iterations));
    drawnow;
end



disp('Finished plotting final results .......');

save Exp01-Info geo P E T SrcInfo BdaryInfo wnum Ns MaxIT ...
                  OptimMethod noiselevel dx dy -ASCII
save Exp01-Results reft ref0 refr sigmat sigma0 sigmar ...
                  gammat gamma0 gammar Gammat Gammar -ASCII

te=toc;
disp(' ');
disp(' ');
disp(['The code ran for: ' num2str(te-tb) ' seconds']);
disp(' ');
disp(' ');

% This last line is used to close MATLAB after the computation. It is 
% only used when runing the code in background.

%exit; % to exit MATLAB
end