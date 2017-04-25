function [xBar, V, D] = buildShapeModel(alignedShapes)
% BUILDSHAPEMODEL performs PCA on a set of aligned shapes to create an active shape
% model.
%
%	INPUT
%       alignedShapes: [2*n_landmarks x n_shapes]
%
%	OUTPUT
%       xBar: Mean shape [2*n_landmarks x 1]
%       V: Eigenvectors (decreasing energy)
%       D: Eigenvalues  (decreasing energy)
%
% John W. Miller
% 25-Apr-2017


%% Create PCA model (on a subset of shapes)
x = alignedShapes(:,1:45);
xBar = mean(x,2);  % Mean shape
S = cov(x');       % Covariance matrix
[V,D] = eig(S);    % Eigenvectors
D = sort(diag(D),'descend');
V = fliplr(V);


end % End of main