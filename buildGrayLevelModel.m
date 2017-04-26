function grayProfileModel = buildGrayLevelModel(pathToImages,shapeModel)
% BUILDGRAYLEVELMODEL
%
%	INPUT
%       pathToImages: Full path to directory containing face images.
%       shapeModel: Comes from buildShapeModel();
%
%	OUTPUT
%       grayProfileModel: (struct) Contains all ya need for the model.
%
% John W. Miller
% 25-Apr-2017

% Multi-resolution
downsample_factors = [6 4 3 2];

for n_resolution = 1:numel(downsample_factors)    
    downsample_factor = downsample_factors(n_resolution);    
    % Pre-allocation and stuff
    % Downsample the image and shape coordinates
    x_bar = shapeModel.meanShape./downsample_factor;
    imageLandmarks = shapeModel.unalignedShapes./downsample_factor;
    
    % Image filenames
    imFiles = dir([pathToImages filesep '*.jpg']);
    n_images = length(imFiles);
    
    % Sample a square region of pixels around each landmark
    rc = 8; % Size of square (rc+1)
    if mod(rc,2)~=0;rc=rc+1;end
    
    % Convolution mask (for gradient of the square regions)
    kernel = [0 -1 0; -1 2 0; 0 0 0];
    C = 5; % Sigmoid equalization
    
    % Pre-allocate to store square profiles at each landmark
    n_landmarks = length(x_bar)/2;
    im_profiles = cell(n_landmarks,1); % One profile per landmark
    for n_landmark = 1:n_landmarks
        im_profiles{n_landmark} = zeros((rc-1)^2,n_images);
    end
    
    % Filter for smoothing before measuring profile
    [A,mew,sig] = deal(1,0,0.1);
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
            step_size=1;
            [Xq,Yq] = meshgrid((iPixel(1)-rc/2):step_size:(iPixel(1)+rc/2),(iPixel(2)-rc/2):step_size:(iPixel(2)+rc/2));
            im_region = interp2(double(im),Xq,Yq,'cubic');
            
            % Calculate the gradient for this region
            im_region_filt = conv2(im_region,kernel,'valid');
            abs_sum = sum(abs(im_region_filt(:)));
            if abs_sum > 0.1 % (Don't divide by zero)
                im_region_filt = im_region_filt./abs_sum;
            end
            
            % Sigmoid equalization
            im_region_filt = im_region_filt./(abs(im_region_filt)+C);
            %         figure(2), hold on, imshow(im_region_filt,[]), hold on
            %         figure(1), hold on, plot(iPixel(1),iPixel(2),'ys'), hold on
            
            % Store as 1D vector
            im_profile = reshape(im_region_filt,(rc-1)^2,1); % minus 1 b/c of 'valid' above
            im_profiles{n_landmark}(:,n_image) = im_profile;
        end % Looping through landmarks
    end % Looping through images
    
    %% Build gray-level model for each landmark (using PCA) from all images
    % Covariance matrices (one per landmark)
    [gBar,S,eVectors,eValues] = deal(cell(n_landmarks,1));
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
        n_landmark = 10;
        % Mean profile at a specific landmark
        figure,imshow(reshape(gBar{n_landmark},rc-1,rc-1),[])
        
        % Generate new profile
        n_pcs = 5;
        prof_mew = gBar{n_landmark};
        P = eVectors{n_landmark}(:,1:n_pcs);
        b = 0.1*ones(size(P,2),1);
        prof_new = prof_mew + P*b;
        figure,imshow(reshape(prof_new,rc-1,rc-1),[])
        
    end
    
    %% Store in a struct
    if n_resolution == 1
        grayProfileModel(n_resolution,1) = struct(); %#ok<*AGROW>
    else
        grayProfileModel(n_resolution,1) = cell2struct(cell(length(fields(grayProfileModel)),1),fields(grayProfileModel),1);
    end    
    
    % Populate the struct
    grayProfileModel(n_resolution).meanProfile = gBar;
    grayProfileModel(n_resolution).covMatrix = S;
    grayProfileModel(n_resolution).eVectors  = V;
    grayProfileModel(n_resolution).eValues   = D;
    grayProfileModel(n_resolution).Info.n_images = n_images;
    grayProfileModel(n_resolution).Info.imageDir = pathToImages;
    grayProfileModel(n_resolution).Info.n_downsample = downsample_factor;
    grayProfileModel(n_resolution).Info.rc_squaresize = rc;
    
    fprintf('\nResolution: 1/%d. %d remaining.',downsample_factor,numel(downsample_factors)-n_resolution)
end % End looping through downsampling factors
end % End of main