%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function collect information about each node on the boundary of the 
% domain
%
% The matrix BdaryInfo is arranged as follows:
%
% BdaryInfo(1,nd): the boundary node where detector is located
% BdaryInfo(2,nd): which segment of the boundary each node is located
%                  1, 2, 3, 4, 5 for bottom, right, top, left and corner 
% BdaryInfo(3,nd): x component of outer normal vector
% BdaryInfo(4,nd): y component of outer normal vector
% BdaryInfo(5,nd): length of the boundary edge where the node is located
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function BdaryInfo=SetBdaryInfo(P,E)

ne=length(E(1,:)); % total # of boundary nodes, same total # bdary edges
BdaryInfo=zeros(5,ne);

BdaryInfo(1,:)=E(1,:); % all boundary nodes

% Find which segment of boundary a node belongs to
for k=1:ne
    if abs(P(1,BdaryInfo(1,k))-0)<1e-8 % left
        if abs(P(2,BdaryInfo(1,k))-0)<1e-8
            BdaryInfo(2,k)=5; % corner segments are denoted with 5
        elseif abs(P(2,BdaryInfo(1,k))-2)<1e-8
            BdaryInfo(2,k)=5;
        else
            BdaryInfo(2,k)=4;
        end
    elseif abs(P(1,BdaryInfo(1,k))-2)<1e-8 % right
        if abs(P(2,BdaryInfo(1,k))-0)<1e-8
            BdaryInfo(2,k)=5;
        elseif abs(P(2,BdaryInfo(1,k))-2)<1e-8
            BdaryInfo(2,k)=5;
        else
            BdaryInfo(2,k)=2;
        end
    else
        if abs(P(2,BdaryInfo(1,k))-0)<1e-8 % bottom & top
            BdaryInfo(2,k)=1; % bottom
        else
            BdaryInfo(2,k)=3; % top
        end
    end
end  

% Find the outer normal vector and the length of each boundary edge 
% The boundary edges given in E are assumed to be ordered in the usual 
% convention: when traveling on the boundary, the domain is always on 
% the left. The outer normal vectors are found by rotating the tangent 
% vector by 90 degrees clockwise. This is done by applying the rotation 
% matrix:
% \cos\theta -\sin\theta
% \sin\theta \cos\theta
% with \theta=-90 since the rotation matrix rotates in the anticlock direction.
BdryEdgeVec=zeros(2,ne);
for j=1:ne
    ind1=E(1,j); ind2=E(2,j);
    vecx=P(1,ind2)-P(1,ind1);
    vecy=P(2,ind2)-P(2,ind1);
    veclength=sqrt(vecx^2+vecy^2);
    vecx=vecx/veclength;
    vecy=vecy/veclength;
    BdryEdgeVec(1,j)=vecy;
    BdryEdgeVec(2,j)=-vecx;
    BdaryInfo(5,j)=veclength;
end

% Find the outer normal vectors of boundary nodes
% The normal vector at a node is defined here as the average of the outer 
% normal vector of the boundary edges joined together by the node
BdaryInfo(3:4,1)=(BdryEdgeVec(:,ne)+BdryEdgeVec(:,1))/2;
for j=2:ne
    BdaryInfo(3:4,j)=(BdryEdgeVec(:,j-1)+BdryEdgeVec(:,j))/2;
end