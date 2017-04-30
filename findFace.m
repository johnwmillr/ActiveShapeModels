function varargout = findFace(im_original,shapeModel,grayModel)
% FINDFACE
%
%	INPUT
%
%
%
%	OUTPUT
%
%   TODO: Add varargin for selecting parameters (e.g. dist_metric, n_evolutions, etc.)
%
% John W. Miller
% 25-Apr-2017

% Use Mahalanobis distance or gray-level PCA space for search?
dist_metric = 'pca'; % Options: 'maha' or 'pca'

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
    
    % Place mean shape over face (this should be done automatically)
    if n_resolution == 1
        [x_original_estimate, h_im] = placeShape(im,x_mean);
%         [x_original_estimate, h_im] = asm_multiResolution(im,x_mean,GM.Info.SmoothingFilter);
        tic
        x_aligned = x_original_estimate;
        x_current = x_aligned;
    else
        resolutionScale = grayModel(n_resolution-1).Info.downsampleFactor/downsampleFactor;
        x_original_estimate = x_original_estimate*resolutionScale; % For final display
        x_current = x_current*resolutionScale;
        x_aligned = x_current; % Use the position determined from the lower resolution
        figure(h_im), hold off
        imshow(im,[]),plotLandmarks(x_aligned,'hold',1)
    end
    
    % Evolve estimate of face location, adjusting landmarks w/in model space
    n_evolutions = 4;
    for n_evolution = 1:n_evolutions;
        title(sprintf('Downsample: 1/%d. %d remaining.\nEvolution %d/%d',...
            downsampleFactor,n_resolutions-n_resolution,n_evolution,n_evolutions),'fontsize',FS), drawnow('expose')
        x_suggested = zeros(size(x_aligned));
        
        for n_landmark = 1:n_landmarks
            %% Calculate 2D profiles near iPixel
            % Think of this as moving a square grid around iPixel and calculating a profile
            % for the grid in each shifted location
            search_size = search_sizes(n_resolution); % Shift amount in x,y directions (xy^2 locations in total)
            %             gShifting = zeros((rc-1)^2,(xy+1)^2,1);
            n_shift = 0;
            dist_min = [];
            iPixel_startingPosition = x_current(2*n_landmark+[-1 0]); % Starting position of current pixel for current evolution
            
            % Shift the 2D profile around the current landmark
            for c = -(search_size/2):(search_size/2)
                for r = -(search_size/2):(search_size/2)
                    n_shift = n_shift+1;
                    iPixel = iPixel_startingPosition+[r c]'; % Coordinates of current pixel
                    plot(iPixel(1),iPixel(2),'bs'), hold on
                    
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
%             plot(best_pixel(1),best_pixel(2),'y.')
            
        end % Looping through landmarks
        
        % Update pose parameters towards suggested shape
        xs = [x_suggested(1:2:end) x_suggested(2:2:end)]; xc = [x_current(1:2:end) x_current(2:2:end)];
%         [~,xp] = procrustes(xs,xc);
%         x_posed = zeros(size(x_mean));
%         x_posed(1:2:end) = xp(:,1);
%         x_posed(2:2:end) = xp(:,2);
        [~,x_posed] = procrustes(x_suggested,x_current);
        
        % Deform shape towards suggested points
        b_suggested = P'*(x_posed-x_mean);
        b=max(min(b_suggested,maxb),-maxb); % Keep adjustments within model limits
        
        % Generate new shape (within model space)
        x_new = x_mean + P*b;
        
        % Transfer x_new to image space (for some reason we need to change the array size)
        xn = [x_new(1:2:end) x_new(2:2:end)];
        [~,xn] = procrustes(xs,xn);
        x_new = zeros(size(x_mean));
        x_new(1:2:end) = xn(:,1); x_new(2:2:end) = xn(:,2);
        x_current = x_new;
        
        % View the current evolution
        imshow(im,[]), plotLandmarks(x_current,'hold',1,'linewidth',3), hold on
        plot(x_suggested(1:2:end),x_suggested(2:2:end),'bo','linewidth',4)
    end % End looping through evolutions
end % End looping through resolution levels

% Scale the final shape estimation to original image size
x_final = x_new*downsampleFactor;
x_original_estimate = x_original_estimate*downsampleFactor; % For final display

% Compare the original estimate with the final evolution
figure(gcf), hold off, imshow(im_original,[])
h_orig  = plotLandmarks(x_original_estimate,'hold',1,'linestyle','--');
h_final = plotLandmarks(x_final,'hold',1,'color','g','linewidth',4);
legend([h_orig,h_final],{'Original','Final'},'location','nw','fontsize',FS)
title('Final shape','fontsize',FS)

% Output final shape and model parameters
if nargout == 1
    varargout{1} = x_final;
elseif nargout == 2
    varargout{1} = x_final;
    varargout{2} = b;
end, toc

end % End of main