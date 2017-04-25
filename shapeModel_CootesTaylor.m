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

% Another one: http://www.milbo.org/muct/index.html (w/ more landmarks around face)

% ASM PhD dissertation: http://www.milbo.org/stasm-files/phd-milborrow.pdf#page=16

% Idea: use ASMs as a classifier of ILM (normal and swollen)
% start by testing the classifier just between ILM and RPE (should be easy)

addpath(fullfile(go('iowa'),'ECE_7480_AdvancedDigitalImageProcessing','Project','Code'))

%% Place landmarks
% landmarks = placeLandmarks(pathToImages,20,10);
addpath('./Utilities/')
load('./Landmarks/landmarks_faces_A_1-50')
alignedShapes = alignShapes(allLandmarks,0);
plotLandmarks(alignedShapes)

%% Create the shape model
[xBar, eVectors, eVals] = buildShapeModel(alignedShapes);

%% Explore the model

% Effect of PC weights on shape (GUI)
guiPrinComps(xBar,eVectors,eVals,'show_image',0)

% Plot variations of the different PCs
plotPrinComp(eVectors,eVals,xBar,1)

%% Edge detection using ASMs
imDir = './Faces/faces_B';
imFile = 'B_40_0.jpg';
im = imread(fullfile(imDir,imFile));

%% Multi-resolution
% Roughly align mean shape to face in image using multi-resolution
% x_aligned = asm_multiResolution(im,alignedShapes);
x_aligned = placeShape(im,xBar);

%% Loop through each landmark point, calculating the normal vector
asm_findFace(im,x_aligned,eVectors,eVals,xBar);








