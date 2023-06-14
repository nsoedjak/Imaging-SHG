classdef Rectangle

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties:
    %
    % rec: 2x2 array
    % height: scalar
    % Represents the rectangle [a b]x[c d]=[rec(1,1) rec(1,2)]x
    %   [rec(2,1) rec(2,2)] with z-height equal to height
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        rec
        height
    end

    methods 
        % constructor
        function obj = Rectangle(rec,height)
            obj.rec = rec;
            obj.height = height;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % The indicator function of the rectangle: rec
        % 
        % Input:
        % 
        % x: spatial variable
        %
        % Output: the value of the characteristic function of rectangle at x
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function f=indicator(obj, x)
            rec = obj.rec;
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
        end
    end
end