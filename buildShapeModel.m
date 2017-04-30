function shapeModel = buildShapeModel(unalignedShapes,pathToTrainingImages)
% BUILDSHAPEMODEL performs PCA on a set of aligned shapes to create an active shape
% model.
%
%	INPUT
%       unalignedLandmarks: Unaligned shapes, placed on training images [2*n_landmarks x n_shapes]
%
%	OUTPUT
%       ShapeModel (struct)
%       xBar: Mean shape [2*n_landmarks x 1]
%       V: Eigenvectors (decreasing energy)
%       D: Eigenvalues  (decreasing energy)
%
% TODO: Do I need to account for resolution somehow with the eigenvectors and values?
%
% John W. Miller
% 25-Apr-2017

% Align shapes from training images using Procrustes
scaling = 0; % Almost always best to set scaling to 0
x = alignShapes(unalignedShapes,scaling);

% Use PCA to create model
xBar = mean(x,2);  % Mean shape
S = cov(x');       % Covariance matrix
[V,D] = eig(S);    % Eigenvectors
D = sort(diag(D),'descend');
V = fliplr(V);

% Store model as a struct
shapeModel = struct();
shapeModel.meanShape = xBar;
shapeModel.eVectors = V;
shapeModel.eValues  = D;
shapeModel.alignedShapes = x;
shapeModel.unalignedShapes = unalignedShapes;
if nargin > 1
    shapeModel.trainingImages = pathToTrainingImages;
end
shapeModel.n_shapes = size(x,2);

end % End of main