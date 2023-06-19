%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function setup boundary conditions for the forward and the adjoint 
% Helmholtz problems.
%
% Type="u_Forward" solves the forward Helmholtz model for u:
%
% \Delta u + k^2(1+n)u+ ik\sigma u = S  in \Omega
% u = g, on \partial\Omega
%
% with S=0, g given by the boundary source
%
% Type="Homogeneous_Dirichlet" solves the Helmholtz model:
%
% \Delta w + k^2(1+n)u + ik\sigma w = S  in \Omega
% w=0, on \partial\Omega
%
%
% Type="Homogeneous_Robin" solves the Helmholtz model:
%
% \Delta w + k^2(1+n)u + ik\sigma w = S  in \Omega
% w + ik n*grad w=0, on \partial\Omega
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [qmatrix,gmatrix,hmatrix,rmatrix] = ...
         HelmholtzBC(Type,SrcInfo,BdaryInfo,ks,wnum,p,e,u,time)

ne = size(e,2); % number of edges on the domain boundary

if strcmp(Type,'Homogeneous_Robin')
    % Robin BCs
    %w + ik n*grad w=0 is equivalent to n*(-grad w) + i/k w = 0
    qmatrix = 1i/wnum*ones(1,ne);
    gmatrix = zeros(1,ne);

    %No Dirichlet BCs
    hmatrix = zeros(1,2*ne);
    rmatrix = zeros(1,2*ne);
else
    % This two lines indicate that the BC is NOT Robin
    qmatrix = zeros(1,ne);
    gmatrix = zeros(1,ne);

    % The following two lines set the BC to homogeneous Dirichlet type
    hmatrix = ones(1,2*ne);
    rmatrix = zeros(1,2*ne);
    
    % Set the rmatrix
    if strcmp(Type,'Homogeneous_Dirichlet')
        % boundary source: zero
        for k = 1:ne
            rmatrix(k)=0.0;
            rmatrix(ne+k)=0.0;
        end
    elseif strcmp(Type,'u_Forward') % boundary sources for forward: Gaussians
        xs=SrcInfo(1,ks);
        ys=SrcInfo(2,ks);
        srcseg=SrcInfo(3,ks);
        for k = 1:ne
            x1 = p(1,e(1,k)); % x at first point in segment
            y1 = p(2,e(1,k)); % y at first point in segment
            x2 = p(1,e(2,k)); % x at second point in segment
            y2 = p(2,e(2,k)); % y at second point in segment
            rmatrix(k)=1.0;
            if BdaryInfo(2,k)==srcseg % if the edge lives on the same side with the source
                rmatrix(k) = 2*exp(-((x1-xs)^2+(y1-ys)^2)/0.1);
                rmatrix(k+ne) = 2*exp(-((x2-xs)^2+(y2-ys)^2)/0.1);
            end        
        end
    else
        disp('Must specify problem type to fix BC!');
    end
end