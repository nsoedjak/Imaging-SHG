classdef Circle
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties:
    %
    % center: array of length 2
    % radius: scalar
    % height: scalar
    % Represents the circle with center (center(1),center(2)), radius equal
    %  to radius, and z-height equal to height
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        center
        radius
        height
    end

    methods 
        % constructor
        function obj = Circle(center,radius,height)
            obj.center = center;
            obj.radius = radius;
            obj.height = height;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % The indicator function of the circle: circ
        % 
        % Input:
        % 
        % x: spatial variable
        %
        % Output: the value of the characteristic function of circ at x
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function f=indicator(obj,x)
            center = obj.center;
            radius = obj.radius;
            N=length(x(1,:));
            f(1:N)=0.0;
            for k=1:N
                if norm(center'-x(:,k))<=radius
                    f(k)=1;
                end
            end
            f=f';
        end
    end
end