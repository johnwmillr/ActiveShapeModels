% SHAPEMODEL_COOTESTAYLOR
%
%   Shape analysis techniques based on this paper:
%   Cootes, T. F., Taylor, C. J., Cooper, D. H., & Graham, J. (1995).
%       "Active Shape Models-Their Training and Application. Computer Vision and
%       Image Understanding."
%
%
% John W. Miller
% 16-Feb-2017

% Here's an annotated faces dataset:
% http://personalpages.manchester.ac.uk/staff/timothy.f.cootes/data/bioid_points.html


% Idea: use ASMs as a classifier of ILM (normal and swollen)
% start by testing the classifier just between ILM and RPE (should be easy)


%% Faces
% Place landmarks
% landmarks = placeLandmarks(pathToImages,20,10);
addpath('./Utilities/')
load('./Landmarks/landmarks_faces_A_1-50')
alignedShapes = alignShapes(allLandmarks,0);
plotLandmarks(alignedShapes)

%% Create PCA model (on a subset of shapes)
x = alignedShapes(:,1:45);
xBar = mean(x,2);  % Mean shape
S = cov(x');       % Covariance matrix
[V,D] = eig(S);    % Eigenvectors
D = sort(diag(D),'descend');
V = fliplr(V);

%% Get weights for first 3 PCs and plot them to check for independence between the PCs
P = V(:,1:3);
weights = P'*(alignedShapes-repmat(xBar,1,size(alignedShapes,2)));
figure
plot3(weights(1,:),weights(2,:),weights(3,:),'o')
xlabel('b1','fontsize',FS), ylabel('b2','fontsize',FS), zlabel('b3','fontsize',FS), grid on

%% Get weights for a new shape
n_pcs = 3;
y = alignedShapes(:,end);

% Solve for the weights that will approximate the new shape using the model
P = V(:,1:n_pcs);
b = P'*(y-xBar);

newShape = xBar + P*b;
figure, hold on
plot(y(1:2:end),y(2:2:end),'ko','linewidth',1)
plot(newShape(1:2:end),newShape(2:2:end),'r.','markersize',10,'linewidth',2)
set(gca,'ydir','reverse')
legend({'Real','Model'},'location','best')
title(sprintf('Approximation of new shape using first %d PCs from the model',n_pcs),'fontsize',FS)

%% GUI variations along the various PCs
guiPrinComps(xBar,V,D,'show_image',0);

%% Examine variations from individual PCs
plotPrinComp(V,D,xBar,1)

%% Edge detection using ASMs

% Probably going to need to do some image processing to enhance edges in the images
imDir = './Faces/faces_B';
imFile = 'B_40_0.jpg';
im = imread(fullfile(imDir,imFile));
imshow(im), hold on
plot(xBar(1:2:end),xBar(2:2:end),'ro','linewidth',2)

hFilt = fspecial('average',4*[1 1]);
im_filt = imfilter(im,hFilt);
im_gMag = imgradient(im_filt); % Image gradient
figure, imshow(im_gMag,[])

%% Multi-resolution
% Roughly align mean shape to face in image using multi-resolution
x_aligned = asm_multiResolution(im,alignedShapes);

%% Loop through each landmark point, calculating the normal vector
asm_findFace(im,x_aligned);

%%
% After looping through each point and looking along the normal, we'll end up with a
% set of suggested points, the same dimension as our mean shape. It's then just a
% simple matter of using Procrustes to translate the current points to the suggested
% points.
suggestedPoints = [];

currentPoints = procrustes(suggestedPoints,currentPoints,'scaling',0);

P = V(:,1:5);
dx; % adjustment vectors w/in the image space
db = P'*dx; % adjustments w/in the model space











