%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function set up boundary sources for the simulations
%
% SrcInfo(1,ns): x coordinate of source
% SrcInfo(2,ns): y coordinate of source
% SrcInfo(3,ns): which part of boundary the source is located at: 
%                             bottom 1, right 2, top 3 and left 4
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SrcInfo=SetSources(Ns)

% The following setup is for the domain (0,2)x(0,2). This can
% be changed for other type of domains.

% Source locations: bottom, right, top, left
location=[0.2:0.2:1.8 2*ones(1,9) 1.8:-0.2:0.2 zeros(1,9);...
          zeros(1,9) 0.2:0.2:1.8 2*ones(1,9) 1.8:-0.2:0.2];
% The side of the boundary where the each source is located
segment=[ones(1,9) 2*ones(1,9) 3*ones(1,9) 4*ones(1,9)]; 

if Ns>length(segment)
	disp('Error in setting # of sources!');
	exit;
end

SrcInfo(1:2,:)=location(:,1:Ns);
SrcInfo(3,:)=segment(1:Ns);