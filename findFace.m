function varargout = findFace(im_original,shapeModel,grayModel,varargin)
% FINDFACE uses an active shape model to locate the facial features in the supplied
% image.
%
%	INPUT
%       im_original: Image file, roughly 480x640. Grayscale.
%       shapeModel:  Active shape model, comes from BUILDSHAPEMODEL
%       grayModel:   Gray-level 2D profiles, comes from BUILDGRAYMODEL
%       OPTIONAL (key-value pairs)
%           save_video: (bool) Do you want to save a video?
%           visualize:  (bool) Do you want to watch the iterations?
%           facefinder:  'click' or 'auto', method for initial face localization, default: 'click'
%           dist_metric: 'pca'   or 'maha', method for measuring accuracy of shifted profiles, default: 'pca'
%           evolutions:  Number of shape evolutions at each resolution level, default: 4
%
%	OUTPUT
%       x_final: Final estimate for face shape position [2*n_landmarks x 1]
%       b: Final model weights for the PCA shape model
%
% John W. Miller
% 25-Apr-2017

% Key-value pair varargin
keys = {'save_video','visualize','facefinder','dist_metric','evolutions','vis_grid'};
default_values = {0,1,'click','pca',4,0};
[save_video, vis, face_find_method, dist_metric, n_evolutions, vis_grid] =...
    parseKeyValuePairs(varargin,keys,default_values);

% Save a video?
if save_video, close all, vis=1;
    videoFileName = [pwd filesep sprintf('ASM_FaceDetection_%s',date)];
    vidObj = VideoWriter(videoFileName,'MPEG-4'); % Video File
    vidObj.FrameRate = 50;
    open(vidObj); disp('Recording video...')
end

% Visualization stuff
if vis_grid, vis=1; end;
if ~vis, face_find_method = 'auto'; end

% Inverse of covariance matrix sometimes badly scaled
warning('off','MATLAB:nearlySingularMatrix');

% Shape model and constraints on new shapes
[~,V,D] = deal(shapeModel.meanShape,shapeModel.eVectors,shapeModel.eValues);

% How many resolution levels are in the gray model?
n_resolutions = length(grayModel);
search_sizes = round(linspace(8,3,n_resolutions));
for n_resolution = 1:n_resolutions
    GM = grayModel(n_resolution);
    step_size=GM.Info.interp_step_size; % Interpolation step size (for im_region)
    
    % Smooth the image (Anti-aliasing?)
    h_filt = GM.Info.SmoothingFilter;
    im = conv2(h_filt,h_filt,im_original,'same');
    
    % Downsample everything that needs to be downsampled
    downsampleFactor = GM.Info.downsampleFactor;
    im = imresize(im,1/downsampleFactor);            % Resize image
    x_mean = shapeModel.meanShape./downsampleFactor; % Resize mean shape
    
    % Gray-level profiles model
    [g_mean, g_S, g_eVecs, g_eVals] = deal(GM.meanProfile, GM.covMatrix, GM.eVectors, GM.eValues);
    n_pcs = 15;
    if n_pcs >= size(g_eVecs{1},1)
        n_pcs = round(0.8*size(g_eVecs{1},1));
    end
    
    % Shape model constraints
    P = V(:,1:n_pcs);
    maxb=3*sqrt(D(1:n_pcs));
    
    % Sample a square region of pixels around each landmark
    rc = GM.Info.rc_squaresize; % Size of square (rc+1)
    n_landmarks = length(x_mean)/2;
    
    % Convolution mask (for gradient of the square regions)
    kernel = GM.Info.EdgeKernel; % Basically the image gradient
    C = GM.Info.SigmoidEQ; % Sigmoid equalization
    
    % Place mean shape over face (Can be automatic or user clicks on image)
    if n_resolution == 1
        if strcmpi(face_find_method,'auto')
            [x_original_estimate, h_im] = estimateFaceLocation(im,x_mean,GM.Info.SmoothingFilter,vis);
        else
            [x_original_estimate, h_im] = placeShape(im,x_mean);
        end, tic
        [x_aligned,x_current] = deal(x_original_estimate);
        if save_video && ishandle(h_im)
            try vidFrame = getframe(h_im);
                writeVideo(vidObj,vidFrame);
            catch
                close(vidObj)
            end
        end
    else
        resolutionScale = grayModel(n_resolution-1).Info.downsampleFactor/downsampleFactor;
        x_original_estimate = x_original_estimate*resolutionScale; % For final display
        x_current = x_current*resolutionScale;
        x_aligned = x_current; % Use the position determined from the lower resolution
        if vis, figure(h_im), hold off
            imshow(im,[]),plotLandmarks(x_aligned,'hold',1), end
    end
    
    % Evolve estimate of face location, adjusting landmarks w/in model space
    for n_evolution = 1:n_evolutions;
        if vis, title(sprintf('Downsample: 1/%d. %d remaining.\nEvolution %d/%d',...
                downsampleFactor,n_resolutions-n_resolution,n_evolution,n_evolutions),'fontsize',FS), drawnow('expose'), end
        x_suggested = zeros(size(x_aligned));
        
        for n_landmark = 1:n_landmarks
            %% Calculate 2D profiles near iPixel
            % Think of this as moving a square grid around iPixel and calculating a profile
            % for the grid in each shifted location
            search_size = search_sizes(n_resolution); % Shift amount in x,y directions (xy^2 locations in total)
            n_shift = 0;
            dist_min = [];
            iPixel_startingPosition = x_current(2*n_landmark+[-1 0]); % Starting position of current pixel for current evolution
            
            % Shift the 2D profile around the current landmark
            for c = -(search_size/2):(search_size/2)
                for r = -(search_size/2):(search_size/2)
                    n_shift = n_shift+1;
                    iPixel = iPixel_startingPosition+[r c]'; % Coordinates of current pixel
                    if vis_grid, plot(iPixel(1),iPixel(2),'bs'), hold on, end
                    
                    % Interpolate the square around the pixel
                    [Xq,Yq] = meshgrid((iPixel(1)-rc/2):step_size:(iPixel(1)+rc/2),(iPixel(2)-rc/2):step_size:(iPixel(2)+rc/2));
                    im_region = interp2(im2double(im),Xq,Yq,'cubic'); % You could also put a scalar after 'cubic' to fill NaNs
                    if any(any(isnan(im_region))) % Fill nans that show up when region extends out of image boundaries
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
                    
                    % Store the current 2D profile as a 1D vector
                    g_new = reshape(im_region_filt,size(im_region_filt,1).^2,1);
                    
                    % Compute the 'distance' from the current profile to the mean profile
                    g_bar = g_mean{n_landmark};
                    if strcmpi(dist_metric,'pca')
                        % Approximate current gray profile in the gray model space
                        g_P = g_eVecs{n_landmark}(:,1:n_pcs); eVals = g_eVals{n_landmark}(1:n_pcs);
                        g_new_b = g_P'*(g_new-g_bar);
                        
                        % How well does the model approximation fit the mean profile?
                        R = (g_new-g_bar)'*(g_new-g_bar) - g_new_b'*g_new_b; % This is actually R^2
                        F = sum((g_new_b.^2)./eVals) + 2*(R./eVals(end)); % Not sure if the second term should be in the sum
                        dist = F;
                    elseif strcmpi(dist_metric,'maha')
                        % Measure Mahalanobis distance for this shift
                        % (This is the distance between this image's profile at the current
                        % shift, compared to the mean profile for this landmark from the
                        % gray-level model that has been passed into this function)
                        md = (g_new-g_bar)'*inv(g_S{n_landmark})*(g_new-g_bar);
                        
                        % Should I be doing the abs(md) ?
                        md = abs(md);
                        dist = md;
                    end
                    
                    % Keep track of best shift
                    if isempty(dist_min)
                        dist_min = dist;
                        best_pixel = iPixel';
                    elseif dist < dist_min
                        dist_min = dist;
                        best_pixel = iPixel';
                    end
                end
            end % End shifting square grid around iPixel
            
            % TODO: Move the distance calculations outside of the loops
            x_suggested(2.*n_landmark+[-1 0]) = best_pixel;
            
            % Visualize & save video (optional)
            if vis_grid, plot(best_pixel(1),best_pixel(2),'y.'), drawnow(), end
            if save_video && ishandle(h_im)
                try vidFrame = getframe(h_im);
                    writeVideo(vidObj,vidFrame);
                catch
                    close(vidObj)
                end
            end
            
        end % Looping through landmarks
        
        % Update pose parameters towards suggested shape
        [~,x_posed] = procrustes(x_suggested,x_current);
        
        % Deform shape towards suggested points
        b_suggested = P'*(x_posed-x_mean);
        b=max(min(b_suggested,maxb),-maxb); % Keep adjustments within model limits
        
        % Generate new shape (within model space)
        x_new = x_mean + P*b;
        
        % Transfer x_new to image space (for some reason we need to change the array shape)
        xn = [x_new(1:2:end) x_new(2:2:end)]; xs = [x_suggested(1:2:end) x_suggested(2:2:end)];
        [~,xn] = procrustes(xs,xn);
        x_new = zeros(size(x_mean));
        x_new(1:2:end) = xn(:,1); x_new(2:2:end) = xn(:,2);
        x_current = x_new;
        
        if vis % View the current evolution
            imshow(im,[]), plotLandmarks(x_current,'hold',1,'linewidth',3), hold on
            plot(x_suggested(1:2:end),x_suggested(2:2:end),'bo','linewidth',4)
            
            if save_video && ishandle(h_im)
                try vidFrame = getframe(h_im);
                    writeVideo(vidObj,vidFrame);
                catch
                    close(vidObj)
                end
            end
        end
    end % End looping through evolutions
end % End looping through resolution levels

% Scale the final shape estimation to original image size
x_final = x_new*downsampleFactor;
x_original_estimate = x_original_estimate*downsampleFactor; % For final display

% Compare the original estimate with the final evolution
if vis || 1
    figure(gcf), hold off, imshow(im_original,[])
    h_orig  = plotLandmarks(x_original_estimate,'hold',1,'linestyle','--');
    h_final = plotLandmarks(x_final,'hold',1,'color','g','linewidth',4);
    legend([h_orig,h_final],{'Original','Final'},'location','nw','fontsize',FS)
    title('Final shape','fontsize',FS)
    if save_video
        try vidFrame = getframe(h_im);
            writeVideo(vidObj,vidFrame);
        catch
            close(vidObj)
        end, close(vidObj), disp('Video closed.') % Close video
    end
end

% Output final shape and model parameters
if nargout == 1
    varargout{1} = x_final;
elseif nargout == 2
    varargout{1} = x_final;
    varargout{2} = b;
end, toc

end % End of main