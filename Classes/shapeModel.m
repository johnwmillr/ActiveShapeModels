% SHAPEMODEL
%
%
% John W. Miller
% 24-Jul-2017

classdef shapeModel < handle
    properties (SetAccess = private)
                
    end
    properties (SetAccess = protected)
                
    end
    methods
        % Constructor
        function shapeObj = shapeModel(n_images)
            if nargin > 0 
                shapeObj.n_images = n_images;                                                
            end
        end
        function view(obj)                                                                                                                               
            
        end
                                                                
    end
end


















