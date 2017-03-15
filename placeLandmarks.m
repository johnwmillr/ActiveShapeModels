function [allLandmarks,image_idxs] = placeLandmarks(pathToImages,n_landmarks,n_scans_to_label,image_idxs)
% PLACELANDMARKS allows user to interactively place landmarks on an image.
%
%	INPUT
%
%
%
%	OUTPUT
%
%
% John W. Miller
% 14-Feb-2017

% Varargin
if nargin < 3
    n_scans_to_label = 3;
end
if nargin < 4
    image_idxs = [];
end

% Get paths for every single patient scan in directory
file_extension = '*.jpg';
files = dir(fullfile(pathToImages,file_extension));
pathsToAllImages = strcat([pathToImages filesep],{files.name}');

% Randomly select scans to label, unless scan indices were specified
if isempty(image_idxs)
%     [imagesToLabel,image_idxs] = datasample(pathsToAllImages,n_scans_to_label,'replace',false);
    imagesToLabel = pathsToAllImages;
    image_idxs = 1:length(imagesToLabel);
else
    imagesToLabel = pathsToAllImages(image_idxs);
end
allLandmarks = zeros(2*n_landmarks,n_scans_to_label);

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
    
    % Connect the landmarks w/ a spline
    %     pp = spline(x,y);
    %     xx = linspace(min(x),max(x));
    %     plot(xx,ppval(pp,xx),'-')
    pause(0.3)
end
close, toc

% Save the landmarks?
save_landmarks = 1;
if save_landmarks
    save_name = ['landmarks_' input('Save name? landmarks_','s') '.mat'];
    save(save_name,'allLandmarks','image_idxs');
end


end % End of main
