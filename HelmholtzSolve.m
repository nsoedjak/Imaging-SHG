%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function solves the the forward and adjoint Helmholtz 
% equations with a P_1 finite element method
%
% Type="u_Forward" solves the forward Helmholtz model for u:
%
% \Delta u + k^2(1+n)u+ ik\sigma u = S  in \Omega
% u = f, on \partial\Omega
%
% with S=0, f given by the boundary source
%
% Type="Homogeneous_Dirichlet" solves the Helmholtz model:
%
% \Delta w + k^2(1+n)u + ik\sigma w = S  in \Omega
% w=f, on \partial\Omega
%
% with f=0
%
% Type="Homogeneous_Robin" solves the Helmholtz model:
%
% \Delta w + k^2(1+n)u + ik\sigma w = S  in \Omega
% w + ik n*grad w=f, on \partial\Omega
%
% with f=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u=HelmholtzSolve(Type,SrcInfo,BdaryInfo,ks,P,E,T,wnum,ref,sigma,S)

% interpolation to triangle middle point
refm=pdeintrp(P,T,ref);
sigmam=pdeintrp(P,T,sigma);
Sm=pdeintrp(P,T,S);

% construct mass matrices
[K,M,F]=assema(P,T,-1,wnum^2*(1+refm)+1i*wnum*sigmam,Sm);

% construct boundary conditions
pdebound =@(p,e,u,time)HelmholtzBC(Type,SrcInfo,BdaryInfo,ks,wnum,p,e,[],[]);
[Q,G,H,R] = assemb(pdebound,P,E);

% solve the PDE
u = assempde(K,M,F,Q,G,H,R);