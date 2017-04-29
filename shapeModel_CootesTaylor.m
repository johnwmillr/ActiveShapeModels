% SHAPEMODEL_COOTESTAYLOR
%
%   Shape analysis techniques based on this paper:
%   Cootes, T. F., Taylor, C. J., Cooper, D. H., & Graham, J. (1995).
%       "Active Shape Models-Their Training and Application. Computer Vision and
%       Image Understanding."
%
% John W. Miller
% 16-Feb-2017

% Here's an annotated faces dataset:
% http://personalpages.manchester.ac.uk/staff/timothy.f.cootes/data/bioid_points.html

% Another one: http://www.milbo.org/muct/index.html (w/ more landmarks around face)

% ASM PhD dissertation: http://www.milbo.org/stasm-files/phd-milborrow.pdf#page=16

% Idea: use ASMs as a classifier of ILM (normal and swollen)
% start by testing the classifier just between ILM and RPE (should be easy)

%% Add necessary paths
projectDir = fullfile(go('iowa'),'ECE_7480_AdvancedDigitalImageProcessing','Project','Code');
addpath(projectDir,fullfile(projectDir,'Utilities'),fullfile(projectDir,'Visualization')), cd(projectDir)

%% Place landmarks on images (or load them from disk)
place_new_landmarks = 0;
pathToImages = fullfile(projectDir,'Faces','faces_A_50');
if place_new_landmarks
    allLandmarks = placeLandmarks(pathToImages,20,10);
else
    load(fullfile(projectDir,'Landmarks','landmarks_faces_A_1-50'))
end

%% Create the shape model from the unaligned shapes
shapeModel = buildShapeModel(allLandmarks,pathToImages);

%% Explore the shape model
guiPrinComps(shapeModel,'show_image',1) % Effect of PC weights on shape (GUI)
plotPrinComp(shapeModel,1) % Plot variations of the different PCs

%% Create the gray-level 2D profile model
grayModel = buildGrayLevelModel(pathToImages,shapeModel);

%% Edge detection using ASMs
imDir = fullfile(projectDir,'Faces','faces_B');
n_im = 35;
imFile = sprintf('B_%02d_0.jpg',n_im);
im = imread(fullfile(imDir,imFile));
figure(34), hold on, imshow(im,[]), hold on

%% Multi-resolution
% Roughly align mean shape to face in image using multi-resolution
% x_aligned = asm_multiResolution(im,shapeModel.meanShape);
x_aligned = placeShape(im,xBar);

%% Find a face!
findFace(im,shapeModel,grayModel)







