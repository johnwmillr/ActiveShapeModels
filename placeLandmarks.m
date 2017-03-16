function allLandmarks = placeLandmarks(pathToImages,n_landmarks,n_scans_to_label,varargin)
% PLACELANDMARKS allows user to interactively place landmarks on an image.
%
%	INPUT
%       pathToImages: Directory containing images (string)
%       n_landmarks:  Number of landmarks to assign to each image
%       n_scans_to_label: Number of scans to label from the directory
%       OPTIONAL (key-value pairs)
%           image_idxs:
%           save_landmarks:
%
%	OUTPUT
%       allLandmarks: Matrix containing the assigned landmarks for each image
%                       [2*n_landmarks x n_images]
%
%   See also PLOTLANDMARKS, ALIGNSHAPES
%
% John W. Miller
% 14-Feb-2017

% Key-value pair varargin
keys = {'image_idxs','save_landmarks','file_ext'}; default_values = {[],1,'.jpg'};
[image_idxs,save_landmarks,file_ext] = parseKeyValuePairs(varargin,keys,default_values);

% Get paths for every single patient scan in directory
files = dir(fullfile(pathToImages,['*' file_ext]));
pathsToAllImages = strcat([pathToImages filesep],{files.name}');

% Label all images in directory, unless scan indices were specified
if isempty(image_idxs)
    imagesToLabel = pathsToAllImages;   
else
    imagesToLabel = pathsToAllImages(image_idxs);
end
allLandmarks = zeros(2*n_landmarks,n_scans_to_label);
landmarkInfo = struct();

%% Loop through each scan, user placing landmarks on each
tic
for n_scan = 1:n_scans_to_label
    % Load and view the image
    im = imread(imagesToLabel{n_scan});
    if n_scan==1
        figure(1)
    else
        figure(2)
    end
    imshow(im,[]), hold on
    text(0.01,0.1, sprintf('Place %d landmarks',n_landmarks),'units','normalized','color','r','fontsize',14)
    text(0.01,0.95,sprintf('Scan %d / %d',n_scan,n_scans_to_label),'units','normalized','color','r','fontsize',14)
    
    % Let user click on the image
    [x,y] = deal(zeros(n_landmarks,1));
    for n = 1:n_landmarks
        [x(n),y(n)] = ginput(1);
        plot(x(n),y(n),'rs','markerfacecolor','r','linewidth',3)
        text(x(n),y(n),sprintf('%d',n),'color','g','fontsize',18)
    end
    
    % Store the coordinates in the landmarks vector
    landmarks = zeros(2*n_landmarks,1);
    landmarks(1:2:end,:) = x;
    landmarks(2:2:end,:) = y;
    allLandmarks(:,n_scan) = landmarks;
    
    % Store info about the landmarks and image
    landmarkInfo(n_scan).path = imagesToLabel{n_scan};
    landmarkInfo(n_scan).landmarks = landmarks;
    pause(0.3)
end
close, toc

% Save the landmarks?
if save_landmarks
    save_name = ['landmarks_' input('Save name? landmarks_','s') '.mat'];
    save(save_name,'allLandmarks','image_idxs','landmarkInfo');
end

end % End of main
