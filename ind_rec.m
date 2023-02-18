%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The indicator function of the rectangle: rec
% 
% Input:
% 
% x: spatial variable
% rec: parameters for the rectangle [a b]x[c d]=[rec(1,1) rec(1,2)]x
%      [rec(2,1) rec(2,2)]
%
% Output: the value of the characteristic function of rec at x
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=ind_rec(x,rec)
N=length(x(1,:));
f(1:N)=0.0;
for k=1:N
    if rec(1,1)<x(1,k) & x(1,k)<=rec(1,2)
        if rec(2,1)<x(2,k) & x(2,k)<=rec(2,2)
            f(k)=1;
        end
    end
end
f=f';  