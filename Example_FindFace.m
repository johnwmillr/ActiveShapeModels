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
%   See also BUILDSHAPEMODEL, BUILDGRAYLEVELMODEL, FINDFACE, EXAMPLE_BUILDMUCTMODEL
%
% John W. Miller
% Electrical & Computer Engineering
% University of Iowa
% 15-Dec-2017

close all; clear all; clc
%% Add necessary paths
% Make sure you run this example script while within its containing folder
project_dir = pwd;
addpath(project_dir,fullfile(project_dir,'Utilities'),fullfile(project_dir,'Visualization'))

%% Load image landmarks from disk
landmark_style = 'MUCT';
load(fullfile(project_dir,'Landmarks','Example_FindFace_Landmarks_MUCT'))    

% View the landmarks we just loaded (they are not aligned, we will align them)
plotLandmarks(allLandmarks), pause(1), close

%% Create the shape model from the unaligned shapes
shapeModel = buildShapeModel(allLandmarks);

%% Explore the shape model
view_gui = 0;
if view_gui
    guiPrinComps(shapeModel,'layout',landmark_style) % Effect of PC weights on shape (GUI)
end

%% Create the gray-level 2D profile model (or load it from disk more likely)
create_new_gray_model = false;
if create_new_gray_model
    faceFiles = dir(fullfile(pathToImages,'*.jpg'));
    faceFiles = faceFiles(~cellfun(@isempty,regexp({faceFiles(:).name}','i\d{3}qa-[fm]n')));
    shapeModel.trainingImages = fullfile(project_dir,'Faces','MUCT','muct_images');
    grayModel = buildGrayLevelModel(faceFiles,shapeModel); % This takes about 30 seconds   
else
    % Load the 'grayModel' struct into the workspace.
    % grayModel is a struct ([n_resolutions x 1]) containing the gray-level
    % gradient information for each of the training images at each resolution for
    % each landmark making up the face shape.
    fprintf('\nLoading the gray-level model...')
    load(fullfile(project_dir,'SavedModels','grayModel_MUCT')); fprintf(' Done loading.\n')
end

%% Explore the gray model
n_landmark = 67;
viewGrayProfile(grayModel,1,n_landmark)

%% Find a face!
disp('Finding a face...')
im = imread(fullfile(project_dir,'Faces','Face.jpg')); close all
findFace(im,shapeModel,grayModel,'visualize',1,'facefinder','click','layout','muct')

