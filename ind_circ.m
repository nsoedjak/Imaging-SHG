%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The indicator function of the circle: circ
% 
% Input:
% 
% x: spatial variable
% circ: parameters for the circle located at (x,y)=(circ(1), circ(2)) with 
%       radius r=circ(3)
%
% Output: the value of the characteristic function of circ at x
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=ind_circ(x,circ)
N=length(x(1,:));
f(1:N)=0.0;
for k=1:N
    if norm(circ(1:2)'-x(:,k))<=circ(3)
        f(k)=1;
    end
end
f=f';  