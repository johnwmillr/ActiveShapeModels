function grayModel = buildGrayLevelModel(pathToImages,shapeModel,varargin)
% BUILDGRAYLEVELMODEL
%
%	INPUT
%       pathToImages: Full path to directory containing face images.
%       shapeModel: Comes from buildShapeModel();
%
%	OUTPUT
%       grayProfileModel: (struct) Contains all ya need for the model.
%
%
%   TODO: Make the process for getting gray-level info its own function, so it can be
%   called during the search process for new images (right now I'm repeating code)
%
%   TODO: Add varargin to specify model parameters (e.g. resolution, filter, etc.)
%   TODO: Add ability to save a model to disk once it's trained
%
% John W. Miller
% 25-Apr-2017
tic

% Key-value pair varargin
keys = {'save_model','resolutions'}; default_values = {0,6:-1:3};
[save_model, downsample_factors] = parseKeyValuePairs(varargin,keys,default_values);

if save_model
    fprintf('\n\n')
    save_name = ['grayModel_' input('Save name? grayModel_','s') '.mat'];
    save_dir = fullfile(fileparts(mfilename),'SavedModels');
    disp('Saving file...')
end

% Multi-resolution
n_resolutions = numel(downsample_factors);

% Change parameters based on resolution
interp_step_sizes = 1*ones(n_resolutions,1); % Probably just leave this as 1

% These model parameters are very important.
filter_sigs = linspace(0.8,0.35,n_resolutions); % Should go from about 0.8 to 0.3
region_size = linspace(10,4,n_resolutions);     % Should go from about 20 to 5

% Build a 2D gray-level profile for each landmark at each resolution
for n_resolution = 1:n_resolutions
    downsample_factor = downsample_factors(n_resolution);
    % Pre-allocation and stuff
    % Downsample the image and shape coordinates
    x_bar = shapeModel.meanShape./downsample_factor;
    imageLandmarks = shapeModel.unalignedShapes./downsample_factor;
    
    % Image filenames
    imFiles = dir([pathToImages filesep '*.jpg']);
    n_images = length(imFiles);
    
    % Sample a square region of pixels around each landmark
    rc = region_size(n_resolution); % Size of square (rc+1)
    if mod(rc,2)~=0;rc=rc+1;end
    
    % Convolution mask (for gradient of the square regions)
    kernel = [0 -1 0; -1 2 0; 0 0 0];
    C = 5; % Sigmoid equalization
    
    % Pre-allocate to store square profiles at each landmark
    n_landmarks = length(x_bar)/2;
    im_profiles = cell(n_landmarks,1); % One profile per landmark
    
    % Filter for smoothing before measuring profile
    [A,mew,sig] = deal(1,0,filter_sigs(n_resolution)); % Increase sig for more smoothing
    x_vals = (-3*sig+mew):0.1:(3*sig+mew);
    h_filt = A*exp(-((x_vals-mew).^2)./(2*sig^2));
    
    %% Start loopin'
    % For each image in the training set, calculate the image gradient (kinda) in a
    % square region of pixels around each landmark
    for n_image = 1:n_images
        im = imread(fullfile(pathToImages,imFiles(n_image).name));
        
        % Smooth and downsample the image
        im = conv2(h_filt,h_filt,im,'same');
        im = imresize(im,1./downsample_factor);
        %     figure(1), hold off, imshow(im,[]), hold on
        %     plotLandmarks(imageLandmarks(:,n_image),'hold',1)
        for n_landmark = 1:n_landmarks
            iPixel = imageLandmarks(2*n_landmark+[-1 0],n_image); % Coordinates of current landmark
            
            % Interpolate the square around the pixel
            step_size = interp_step_sizes(n_resolution);
            [Xq,Yq] = meshgrid((iPixel(1)-rc/2):step_size:(iPixel(1)+rc/2),(iPixel(2)-rc/2):step_size:(iPixel(2)+rc/2));
%             im_region = interp2(double(im),Xq,Yq,'cubic'); % You could also put a scalar after 'cubic' to fill NaNs
            im_region = interp2(im2double(im),Xq,Yq,'cubic'); % You could also put a scalar after 'cubic' to fill NaNs
            if any(any(isnan(im_region)))
                imr = reshape(im_region,numel(im_region),1);
                imr(isnan(imr)) = nanmedian(imr);
                im_region = reshape(imr,size(im_region,1),size(im_region,2));
            end
            
            % Calculate the gradient for this region
            im_region_filt = conv2(im_region,kernel,'valid');
            abs_sum = sum(abs(im_region_filt(:)));
            if abs_sum ~= 0
                im_region_filt = im_region_filt./abs_sum;
            end
            
            % Sigmoid equalization
            im_region_filt = im_region_filt./(abs(im_region_filt)+C);
            %         figure(2), hold on, imshow(im_region_filt,[]), hold on
            %         figure(1), hold on, plot(iPixel(1),iPixel(2),'ys'), hold on
            
            % Store as 1D vector
            im_profile = reshape(im_region_filt,size(im_region_filt,1).^2,1); % minus 1 b/c of 'valid' above
            im_profiles{n_landmark}(:,n_image) = im_profile;
        end % Looping through landmarks
    end % Looping through images
    
    %% Build gray-level model for each landmark (using PCA) from all images
    % Covariance matrices (one per landmark)
    [gBar,S,eVectors,eValues] = deal(cell(n_landmarks,1)); % This deal might be slow
    for n_landmark = 1:n_landmarks
        gBar{n_landmark} = mean(im_profiles{n_landmark},2);
        S{n_landmark} = cov(im_profiles{n_landmark}'); % Must transpose here
        [V,D] = eig(S{n_landmark});
        D = sort(diag(D),'descend'); V = fliplr(V);
        eVectors{n_landmark} = V;
        eValues{n_landmark} = D;
    end
    
    %% Visualize the model (optional)
    view_model = 0;
    if view_model
        n_landmark = 14;
        % Mean profile at a specific landmark
        figure,imshow(reshape(gBar{n_landmark},rc-1,rc-1),[])
    end
    
    %% Store in a struct
    if n_resolution == 1
        grayModel(n_resolution,1) = struct(); %#ok<*AGROW>
    else
        grayModel(n_resolution,1) = cell2struct(cell(length(fields(grayModel)),1),fields(grayModel),1);
    end
    
    % Populate the struct
    grayModel(n_resolution).meanProfile = gBar;
    grayModel(n_resolution).covMatrix = S;
    grayModel(n_resolution).eVectors  = eVectors;
    grayModel(n_resolution).eValues   = eValues;
    grayModel(n_resolution).Info.n_images = n_images;
    grayModel(n_resolution).Info.imageDir = pathToImages;
    grayModel(n_resolution).Info.downsampleFactor = downsample_factor;
    grayModel(n_resolution).Info.rc_squaresize = rc;
    grayModel(n_resolution).Info.SigmoidEQ = C;
    grayModel(n_resolution).Info.SmoothingFilter = h_filt;
    grayModel(n_resolution).Info.EdgeKernel = kernel;
    grayModel(n_resolution).Info.interp_step_size = interp_step_sizes(n_resolution);
    
    fprintf('\nResolution scale: 1/%d. %d remaining.',downsample_factor,numel(downsample_factors)-n_resolution)
end, toc % End looping through downsampling factors

% Save the model (optional)
if save_model
    save(fullfile(save_dir,save_name),'grayModel');
end

fprintf('\nAll done. Have a nice day!\n')
end % End of main