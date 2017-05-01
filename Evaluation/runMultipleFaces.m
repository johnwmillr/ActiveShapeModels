% RUNMULTIPLEFACES 
%
%	INPUT
%
%
%
%	OUTPUT
%
%
% John W. Miller
% 30-Apr-2017

project_dir = fullfile(go('iowa'),'ECE_7480_AdvancedDigitalImageProcessing','Project','Code');

imageDir = fullfile(project_dir,'Faces','faces_A_50');
load(fullfile(project_dir,'Landmarks','landmarks_faces_A_1-50'))
files = dir([imageDir filesep '*.jpg']);
n_files = length(files);

% Create the shape model from the unaligned shapes
shapeModel = buildShapeModel(allLandmarks,imageDir);

% Gray Model
load('../SavedModels/grayModel_WideCoverage.mat')

% Automatically find face for many images, compare to their manually-labeled "truth" 
all_found_shapes = zeros(size(allLandmarks,1),n_files);
all_true_shapes = zeros(size(all_found_shapes));
t_start = tic;
for n_file = 1:n_files
    fprintf('\nFile #%d out of %d\n',n_file,n_files)
    filename = sprintf('A_%02d_0.jpg',n_file);
    [~,truth_file,ext] = fileparts(landmarkInfo(n_file).path);
    if ~strcmpi(filename,[truth_file ext])
        error('File %d names don''t match.',n_file)
    else        
        im = imread(fullfile(imageDir,filename));
        all_found_shapes(:,n_file) = findFace_NoPlot(im,shapeModel,grayModel);
        all_true_shapes(:,n_file) = allLandmarks(:,n_file);
    end    
end
toc(t_start)

% Compare the found shapes and the true shapes
all_errors = all_found_shapes-all_true_shapes;

























