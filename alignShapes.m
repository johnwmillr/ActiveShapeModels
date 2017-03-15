function [alignedShapes, avgDiff] = alignShapes(allShapes, scaling)
% ALIGNSHAPES uses Procrustes analysis to align a set of shapes (with or without
% scaling).
%
%	INPUT
%       allShapes: [2*n_landmarks x n_shapes], n_shapes is basically n_subjects
%           Collected from the HD scan semi-automated segmentation           
%           for each row: 20 points (40 elements): x1, y1, x2, y2, ..., x20, y20%           %
%       scaling: (bool) Do you want to scale your images or not? (default = 0)
%
%	OUTPUT
%       alignedShapes: The realigned shapes. Same shape as totalShapes
%       avgDiff: The average difference (error?) between each shape and the mean shape
%
%   Shape analysis techniques based on this paper:
%   Cootes, T. F., Taylor, C. J., Cooper, D. H., & Graham, J. (1995).
%       "Active Shape Models-Their Training and Application. Computer Vision and
%       Image Understanding."
%
%   See also PROCRUSTES, PLOTLANDMARKS
%
% John W. Miller
% 15-Mar-2017

% Pre-allocate
n_shapes = size(allShapes,2);
alignedShapes = zeros(size(allShapes));
totalDiff = 0;

% Mean shape across all subjects (assuming each shape in totalShapes is a different subj)
meanShape = mean(allShapes,2); % x1, y1, x2, y2, ..., x20, y20
meanShape = [meanShape(1:2:end) meanShape(2:2:end)]; % Reshape for Procrustes

%% Loop thru each shape in totalShapes and transform via Procrustes analysis
for n_shape = 1:n_shapes
    % Landmarks shape for the current subject (if 1 shape per subject)    
    iShape = [allShapes(1:2:end,n_shape) allShapes(2:2:end,n_shape)];
    
    % Do the Procrustes alignment
    [d, iShapeAligned] = procrustes(meanShape,iShape,'scaling',scaling,'reflection','best');
    totalDiff = totalDiff + d; % Total difference between iShape and the mean
    
    % Store the aligned shape in a similar manner as how totalShapes is passed in            
    alignedShapes(1:2:end,n_shape) = iShapeAligned(:,1)';
    alignedShapes(2:2:end,n_shape) = iShapeAligned(:,2)';            
end

% On average, how far off was each shape from the mean shape?
avgDiff = totalDiff/n_shapes;

end % End of main