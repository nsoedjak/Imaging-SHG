classdef Profile
    
    % A Profile object is used to describe the profile of each of the 
    % coefficients of the PDE system. The simple profiles we consider here
    % are composed of rectangles and circles on top of a constant
    % background.
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties:
    %
    % background: scalar
    % rectangles: array of Rectangle objects
    % circles: array of Circle objects
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        background
        rectangles
        circles
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Evaluates the profile at the point x
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function val = evaluate(obj, x)
            val = obj.background;
            for i = 1:length(obj.rectangles)
                rectangle = obj.rectangles(i);
                val = val + rectangle.height * rectangle.indicator(x);
            end
            for i = 1:length(obj.circles)
                circle = obj.circles(i);
                val = val + circle.height * circle.indicator(x);
            end
        end
    end
end