% Build MUCT model from scratch
%
% John W. Miller
% 15-Dec-2017

% Download the faces dataset from the link below
% http://www.milbo.org/muct/index.html
project_dir = pwd;
addpath(project_dir,fullfile(project_dir,'Utilities'),fullfile(project_dir,'Visualization')), cd(project_dir)

% Load landmarks and images
pathToLandmarks = fullfile(project_dir,'Faces','MUCT','muct_landmarks');
pathToImages = fullfile(project_dir,'Faces','MUCT','muct_images');
faceFiles = dir(fullfile(pathToImages,'*.jpg'));

% Load landmarks    
C = importdata(fullfile(pathToLandmarks,'muct76.csv'));
landmarks = C.data(1:end,2:end);
subj_labels = C.textdata(2:end,1);           

% Extract the landmarks and images we want
expression = 'i\d{3}[qrs]a-[fm]n';
mask_landmarks = ~cellfun(@isempty,regexpi(subj_labels,expression));
mask_images = ~cellfun(@isempty,regexp({faceFiles(:).name}',expression));

faceFiles = faceFiles(mask_images);
allLandmarks = landmarks(mask_landmarks,:)';
im = rgb2gray(imread(fullfile(pathToImages,faceFiles(1).name)));
allLandmarks(1:2:end,:) =    allLandmarks(1:2:end,:) + size(im,2)/2; % Adjust for coordinate differences
allLandmarks(2:2:end,:) = -1*allLandmarks(2:2:end,:) + size(im,1)/2;   
% save('Example_FindFace_landmarks_MUCT','allLandmarks')

%% View an image with landmarks
n_im = 2;
im = rgb2gray(imread(fullfile(pathToImages,faceFiles(n_im).name)));
figure, imshow(im,[]), hold on
plot(allLandmarks(1:2:end,n_im),allLandmarks(2:2:end,n_im),'ro'), hold off

%% Create the shape model from the unaligned shapes
shapeModel = buildShapeModel(allLandmarks,pathToImages);
plotLandmarks(shapeModel.alignedShapes,'layout','muct')

%% Explore the shape model
guiPrinComps(shapeModel,'layout','muct') % Effect of PC weights on shape (GUI)
plotPrinComp(shapeModel,1,'muct') % Plot variations of the different PCs

%% Create the gray-level 2D profile model
grayModel = buildGrayLevelModel(faceFiles,shapeModel,'resolutions',[16 6 2 1],'region_size',[15 3]);

%% Explore the gray model
n_resolution = 1;
n_landmark = 67;
viewGrayProfile(grayModel,n_resolution,n_landmark)

%% Find a face!
n_im = 30;
im = rgb2gray(imread(fullfile(pathToImages,faceFiles(n_im).name)));
findFace(im,shapeModel,grayModel,'layout','muct','evolutions',3)
