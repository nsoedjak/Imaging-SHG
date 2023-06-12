%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function evaluate the objective function and its gradients with 
% respect to the optimization variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f g]=SHGObj(X,GammaFlag,Gammat,MinVar,x,y,dx,dy,Nx,Ny,P,E,T,...
    Ns,Hm,SrcInfo,BdaryInfo,wnum,betan,betaS,betaG,betag)

M=Nx*Ny; % total number of nodes in the mesh
ne = size(SrcInfo,2); % number of edges/nodes on the domain boundary

refc=X(1:M);% current value of n
sigmac=X(M+1:2*M); % current value of sigma
gammac=X(2*M+1:3*M); %current value of gamma

f=0.0;
g=zeros(3*M,1);

if ~GammaFlag %Gamma is known
    for ks=1:Ns
        
        Hc=zeros(M,1); % predicted data
        rz=zeros(M,1); % residual on measurement locations
        srczero=zeros(M,1); % zero volume source for forward problems
     
        uc=HelmholtzSolve('u_Forward',SrcInfo,BdaryInfo,ks,P,E,T,wnum,refc,sigmac,srczero);
    
        srcv = -(2*wnum)^2 * gammac .* uc.^2;
        vc=HelmholtzSolve('Homogeneous_Robin',SrcInfo,BdaryInfo,ks,P,E,T,2*wnum,refc,sigmac,srcv);
    
        Hc=Gammat.*sigmac.*(abs(uc).^2 + abs(vc).^2);
        
        %Hcg=tri2grid(P,T,Hc,x,y);
        %figure;
        %pcolor(x,y,Hcg); axis tight; colorbar('SouthOutside');
        %axis square; axis off; shading interp;
        %drawnow;
        
        HmL=Hm(:,ks);
        rz=(Hc-HmL); % for unnormalized objective function
        %rz=(Hc-HmL)./HmL; % for normalized objective function
        
        % the contribution to the objective function from source ks
        f=f+0.5*sum(rz.^2)*dx*dy;
        
        % the contribution to the gradient from source ks
        if nargout > 1         
    
            % solve the adjoint equations
            src_uc_2=-Gammat.*sigmac.*rz.*conj(uc);        
            uc_2=HelmholtzSolve('Homogeneous_Dirichlet',SrcInfo,BdaryInfo,ks,P,E,T,wnum,refc,sigmac,src_uc_2);
    
            src_vc_2=-Gammat.*sigmac.*rz.*conj(vc);        
            vc_2=HelmholtzSolve('Homogeneous_Robin',SrcInfo,BdaryInfo,ks,P,E,T,2*wnum,refc,sigmac,src_vc_2);
    
            src_uc_3=-2*(2*wnum)^2*gammac.*uc.*vc_2;        
            uc_3=HelmholtzSolve('Homogeneous_Dirichlet',SrcInfo,BdaryInfo,ks,P,E,T,wnum,refc,sigmac,src_uc_3);
            
            %wcg=tri2grid(P,T,wc,x,y);
            %figure;
            %pcolor(x,y,real(wcg)); axis tight; colorbar('SouthOutside');
            %axis square; axis off; shading interp;
            %drawnow;
            %pause;
        
            % the gradient w.r.t n            
            if ismember("Ref",MinVar)
                g(1:M)=g(1:M)+2*wnum^2*real(uc.*uc_2 + 2^2*vc.*vc_2 + uc.*uc_3)*dx*dy;
            end
            % the gradient w.r.t sigma
            if ismember("Sigma",MinVar)
                g(M+1:2*M)=g(M+1:2*M)+(rz.*Gammat.*(abs(uc).^2 + abs(vc).^2) ...
                    +2*wnum*real(1i*uc.*uc_2 + 2i*vc.*vc_2 + 1i*uc.*uc_3))*dx*dy;
            end
            % the gradient w.r.t. gamma
            if ismember("gamma",MinVar)
                g(2*M+1:3*M)=g(2*M+1:3*M)+2*(2*wnum)^2*real(uc.^2.*vc_2)*dx*dy;
            end
            
        end
        
    end
else %Gamma is unknown
    srczero=zeros(M,1); % zero volume source for forward problems
 
    u1=HelmholtzSolve('u_Forward',SrcInfo,BdaryInfo,1,P,E,T,wnum,refc,sigmac,srczero);

    srcv = -(2*wnum)^2 * gammac .* u1.^2;
    v1=HelmholtzSolve('Homogeneous_Robin',SrcInfo,BdaryInfo,1,P,E,T,2*wnum,refc,sigmac,srcv);

    K1 = abs(u1).^2 + abs(v1).^2;
    for ks=2:Ns
        
        Kc=zeros(M,1);
        rz=zeros(M,1); % residual on measurement locations
        srczero=zeros(M,1); % zero volume source for forward problems
     
        uc=HelmholtzSolve('u_Forward',SrcInfo,BdaryInfo,ks,P,E,T,wnum,refc,sigmac,srczero);
    
        srcv = -(2*wnum)^2 * gammac .* uc.^2;
        vc=HelmholtzSolve('Homogeneous_Robin',SrcInfo,BdaryInfo,ks,P,E,T,2*wnum,refc,sigmac,srcv);
    
        Kc = abs(uc).^2 + abs(vc).^2;
        %Hcg=tri2grid(P,T,Hc,x,y);
        %figure;
        %pcolor(x,y,Hcg); axis tight; colorbar('SouthOutside');
        %axis square; axis off; shading interp;
        %drawnow;
        
        rz = Kc./K1 - Hm(:,ks)./Hm(:,1); % for unnormalized objective function
        
        % the contribution to the objective function from source ks
        f=f+0.5*sum(rz.^2)*dx*dy;
        
        % the contribution to the gradient from source ks
        if nargout > 1         
    
            % solve the adjoint equations
            src_uc_2=-rz.*conj(uc)./K1;        
            uc_2=HelmholtzSolve('Homogeneous_Dirichlet',SrcInfo,BdaryInfo,ks,P,E,T,wnum,refc,sigmac,src_uc_2);

            src_u1_2=rz.*Kc.*conj(u1)./K1.^2;        
            u1_2=HelmholtzSolve('Homogeneous_Dirichlet',SrcInfo,BdaryInfo,1,P,E,T,wnum,refc,sigmac,src_u1_2);
    
            src_vc_2=-rz.*conj(vc)./K1;        
            vc_2=HelmholtzSolve('Homogeneous_Robin',SrcInfo,BdaryInfo,ks,P,E,T,2*wnum,refc,sigmac,src_vc_2);

            src_v1_2=rz.*Kc.*conj(v1)./K1.^2;        
            v1_2=HelmholtzSolve('Homogeneous_Robin',SrcInfo,BdaryInfo,1,P,E,T,2*wnum,refc,sigmac,src_v1_2);
    
            src_uc_3=-2*(2*wnum)^2*gammac.*uc.*vc_2;        
            uc_3=HelmholtzSolve('Homogeneous_Dirichlet',SrcInfo,BdaryInfo,ks,P,E,T,wnum,refc,sigmac,src_uc_3);

            src_u1_3=-2*(2*wnum)^2*gammac.*u1.*v1_2;        
            u1_3=HelmholtzSolve('Homogeneous_Dirichlet',SrcInfo,BdaryInfo,1,P,E,T,wnum,refc,sigmac,src_u1_3);
            
            %wcg=tri2grid(P,T,wc,x,y);
            %figure;
            %pcolor(x,y,real(wcg)); axis tight; colorbar('SouthOutside');
            %axis square; axis off; shading interp;
            %drawnow;
            %pause;
        
            % the gradient w.r.t n            
            if ismember("Ref",MinVar)
                g(1:M)=g(1:M)+2*wnum^2*real(uc.*uc_2 + 2^2*vc.*vc_2 + uc.*uc_3 ...
                    + u1.*u1_2 + 2^2*v1.*v1_2 + u1.*u1_3)*dx*dy;
            end
            % the gradient w.r.t sigma
            if ismember("Sigma",MinVar)
                g(M+1:2*M)=g(M+1:2*M)+2*wnum*real(1i*uc.*uc_2 + 2i*vc.*vc_2 + 1i*uc.*uc_3 ...
                    + 1i*u1.*u1_2 + 2i*v1.*v1_2 + 1i*u1.*u1_3)*dx*dy;
            end
            % the gradient w.r.t. gamma
            if ismember("gamma",MinVar)
                g(2*M+1:3*M)=g(2*M+1:3*M)+2*(2*wnum)^2*real(uc.^2.*vc_2 ...
                    + u1.^2.*v1_2)*dx*dy;
            end
            
        end
        
    end    

end

% Add regularization terms to both the objective function and its gradients

if ismember("Ref", MinVar)
    [Rx,Ry] = pdegrad(P,T,refc);
    Rx1=pdeprtni(P,T,Rx); Ry1=pdeprtni(P,T,Ry);
    f=f+0.5*betan*sum(Rx1.^2+Ry1.^2)*dx*dy;
    if nargout > 1
        [Rxx, Rxy]=pdegrad(P,T,Rx1); [Ryx, Ryy]=pdegrad(P,T,Ry1);
        Rx2=pdeprtni(P,T,Rxx); Ry2=pdeprtni(P,T,Ryy);
        Deltan=Rx2+Ry2;
        g(1:M)=g(1:M)-betan*Deltan*dx*dy;
        for j=1:ne
            nd=BdaryInfo(1,j);
            g(nd)=g(nd)+betan*(BdaryInfo(3,j)*Rx1(nd)+BdaryInfo(4,j)*Ry1(nd))*BdaryInfo(5,j);
        end
    end
end
if ismember("Sigma", MinVar)
    [Sx,Sy] = pdegrad(P,T,sigmac);
    Sx1=pdeprtni(P,T,Sx); Sy1=pdeprtni(P,T,Sy);
    f=f+0.5*betaS*sum(Sx1.^2+Sy1.^2)*dx*dy;
    if nargout > 1
        [Sxx, Sxy]=pdegrad(P,T,Sx1); [Syx, Syy]=pdegrad(P,T,Sy1);
        Sx2=pdeprtni(P,T,Sxx); Sy2=pdeprtni(P,T,Syy);
        DeltaSigma=Sx2+Sy2;
        g(M+1:2*M)=g(M+1:2*M)-betaS*DeltaSigma*dx*dy;
        for j=1:ne
            nd=BdaryInfo(1,j);
            g(M+nd)=g(M+nd)+betaS*(BdaryInfo(3,j)*Sx1(nd)+BdaryInfo(4,j)*Sy1(nd))*BdaryInfo(5,j);
        end
    end
end
if ismember("gamma", MinVar)
    [gx,gy] = pdegrad(P,T,gammac);
    gx1=pdeprtni(P,T,gx); gy1=pdeprtni(P,T,gy);
    f=f+0.5*betag*sum(gx1.^2+gy1.^2)*dx*dy;
    if nargout > 1
        [gxx, gxy]=pdegrad(P,T,gx1); [gyx, gyy]=pdegrad(P,T,gy1);
        gx2=pdeprtni(P,T,gxx); gy2=pdeprtni(P,T,gyy);
        Deltagamma=gx2+gy2;
        g(2*M+1:3*M)=g(2*M+1:3*M)-betag*Deltagamma*dx*dy;
        for j=1:ne
            nd=BdaryInfo(1,j);
            g(2*M+nd)=g(2*M+nd)+betag*(BdaryInfo(3,j)*gx1(nd)+BdaryInfo(4,j)*gy1(nd))*BdaryInfo(5,j);
        end
    end
end