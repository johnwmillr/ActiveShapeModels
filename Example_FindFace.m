% EXAMPLE_FINDFACE steps the user through the process of detecting a face in a single
% example image. The models are built from 50 training images from the link below.
%
%   Run this script for an example of how to use the active shape model (ASM) code.
%   This script needs to be within the main ActiveShapeModels directory in order to
%   load all of the paths properly.
%
%   Link to training images:
%       http://robotics.csie.ncku.edu.tw/Databases/FaceDetect_PoseEstimate.htm#Our_Database_
%
%   Shape analysis techniques based on this paper:
%       Cootes, T. F., Taylor, C. J., Cooper, D. H., & Graham, J. (1995).
%           "Active Shape Models-Their Training and Application. Computer Vision and
%           Image Understanding."
%
%   2D gray-level profiles inspired by:
%       Milborrow, S. (2016). "Multiview Active Shape Models with SIFT Descriptors."
%
%   See also BUILDSHAPEMODEL, BUILDGRAYLEVELMODEL, FINDFACE
%
% John W. Miller
% Electrical & Computer Engineering
% University of Iowa
% 30-Apr-2017

close all; clear all; clc
%% Add necessary paths
% This will work if you run this script. If you're just executing these lines into
% the command line, make sure you set filename = 'Example_FindFace'
filename = mfilename;
project_dir = fileparts(which(filename));
addpath(project_dir,fullfile(project_dir,'Utilities'),fullfile(project_dir,'Visualization')), cd(project_dir)

%% Load image landmarks from disk
% This will load the variable 'allLandmarks' into the workspace.
% allLandmarks is a 40x50 matrix, corresponding to [2*n_landmarks x n_images]
load(fullfile(project_dir,'Landmarks','Example_FindFace_Landmarks'))

% View the landmarks we just loaded (they are not aligned, we will align them)
plotLandmarks(allLandmarks), pause(1), close

%% Create the shape model from the unaligned shapes
shapeModel = buildShapeModel(allLandmarks);

%% Explore the shape model
if ~strcmp(mfilename,'Example_FindFace') % Only view GUI if user didn't call the entire script
    guiPrinComps(shapeModel,'show_image',1) % Effect of PC weights on shape (GUI)
end

%% Create the gray-level 2D profile model (or load it from disk more likely)
if isdir(fullfile(project_dir,'Faces','faces_A_50')) && 0
    pathToImages = fullfile(project_dir,'Faces','faces_A_50');
    grayModel = buildGrayLevelModel(pathToImages,shapeModel); % This takes about 30 seconds
else
    % Load the 'grayModel' struct into the workspace.
    % grayModel is a 4x1 ([n_resolutions x 1]) struct containing the gray-level
    % gradient information for each of the training images at each resolution for
    % each landmark making up the face shape.
    fprintf('\nLoading the gray-level model...')
    load(fullfile(project_dir,'SavedModels','grayModel_Example_FindFace')); fprintf(' Done loading.\n')
end

%% Load the example face
im = imread(fullfile(project_dir,'Faces','Face.jpg')); close all

%% Find a face!
disp('Finding a face...')
findFace(im,shapeModel,grayModel,'visualize',1,'facefinder','cick')

