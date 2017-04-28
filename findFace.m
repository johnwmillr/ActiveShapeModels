function [x_final, b] = findFace(im_original,shapeModel,grayModel)
% FINDFACE
%
%	INPUT
%
%
%
%	OUTPUT
%
%
% John W. Miller
% 25-Apr-2017

% The initial placement has to be very close to the center of the face...

% Shape model and constraints on new shapes
[x_bar,V,D] = deal(shapeModel.meanShape,shapeModel.eVectors,shapeModel.eValues);

% How many resolution levels are in the gray model?
n_resolutions = length(grayModel);
for n_resolution = 1:n_resolutions
    GM = grayModel(n_resolution);
    step_size=GM.Info.interp_step_size; % Interpolation step size (for im_region)
    
    % Gray-level profiles model
    [gBar,S] = deal(GM.meanProfile, GM.covMatrix);
    
    % Smooth the image (Anti-aliasing?)
    h_filt = GM.Info.SmoothingFilter;
    im = conv2(h_filt,h_filt,im_original,'same');
    
    % Downsample everything that needs to be downsampled
    downsampleFactor = GM.Info.downsampleFactor;
    im = imresize(im,1/downsampleFactor);           % Resize image
    x_bar = x_bar./downsampleFactor;                % Resize mean shape
    %     V = V./downsampleFactor;    % Is this a legit method for downsampling V and D?
    %     D = D./downsampleFactor;
    
    % Shape model constraints
    n_pcs = 5;
    P = V(:,1:n_pcs);
    maxb=3*sqrt(D(1:n_pcs));
    
    % Sample a square region of pixels around each landmark
    rc = GM.Info.rc_squaresize; % Size of square (rc+1)
    n_landmarks = length(x_bar)/2;
    
    % Convolution mask (for gradient of the square regions)
    kernel = GM.Info.EdgeKernel; % Basically the image gradient
    C = GM.Info.SigmoidEQ; % Sigmoid equalization
    
    % Inverse of covariance matrix sometimes badly scaled
    warning('off','MATLAB:nearlySingularMatrix');
    
    % Place mean shape over face (this should be done automatically)
    if n_resolution == 1
        [x_original_estimate,h_im] = placeShape(im,x_bar);
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
    n_evolutions = 5;
    for n_evolution = 1:n_evolutions;
        title(sprintf('Downsample: 1/%d.\n%d remaining. Evolution %d/%d',...
            downsampleFactor,n_resolutions-n_resolution,n_evolution,n_evolutions),'fontsize',FS)
        x_suggested = zeros(size(x_aligned));
        
        for n_landmark = 1:n_landmarks
            %% Calculate 2D profiles near iPixel
            % Think of this as moving a square grid around iPixel and calculating a profile
            % for the grid in each shifted location
            xy = 5; % Shift amount in x,y directions (xy^2 locations in total)
            %             gShifting = zeros((rc-1)^2,(xy+1)^2,1);
            n_shift = 0;
            md_min = [];
            iPixel_startingPosition = x_current(2*n_landmark+[-1 0]); % Starting position of current pixel for current evolution
            for c = -(xy/2):(xy/2)
                for r = -(xy/2):(xy/2)
                    n_shift = n_shift+1;
                    iPixel = iPixel_startingPosition+[r c]'; % Coordinates of current pixel
                    %                     plot(iPixel(1),iPixel(2),'bs'), hold on
                    
                    % Interpolate the square around the pixel
                    [Xq,Yq] = meshgrid((iPixel(1)-rc/2):step_size:(iPixel(1)+rc/2),(iPixel(2)-rc/2):step_size:(iPixel(2)+rc/2));
                    im_region = interp2(double(im),Xq,Yq,'cubic'); % You could also put a scalar after 'cubic' to fill NaNs
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
                    
                    % Store as 1D vector
                    profileAtCurrentShift = reshape(im_region_filt,size(im_region_filt,1).^2,1);
                    gShifting(:,n_shift) = profileAtCurrentShift;
                    
                    % Measure Mahalanobis distance for this shift
                    % (This is the distance between this image's profile at the current
                    % shift, compared to the mean profile for this landmark from the
                    % gray-level model that has been passed into this function)
                    md = (gShifting(:,n_shift)-gBar{n_landmark})'*inv(S{n_landmark})*(gShifting(:,n_shift)-gBar{n_landmark});
                    
                    % Should I be doing the abs(md) ?
                    md = abs(md);
                    
                    if isempty(md_min)
                        md_min = md;
                        best_pixel = iPixel';
                    elseif md < md_min
                        md_min = md;
                        best_pixel = iPixel';
                        %                         plot(iPixel(1),iPixel(2),'g.')
                    end
                end
            end % End shifting square grid around iPixel
            
            %% Mahalanobis Distance
            %             MDs = bsxfun(@minus,gShifting,gBar{n_landmark})'*inv(S{n_landmark})*bsxfun(@minus,gShifting,gBar{n_landmark});
            
            %
            %% MAHADIST is symmetric, multiple options for shift direciton, figure this out!
            %figure(4),imshow(mahaDist,[]), hold on, plot(I,J,'ro','linewidth',3), hold off
            
            % x,y shift from iPixel to the smallest maha distance
            %         iShift = (((size(gShifting,2)-1)/2+1)*[1 1])-[I J];
            %         suggested_shifts(2.*n_landmark+[-1 0]) = iShift;
            %         x_suggested(2.*n_landmark+[-1 0]) = iPixel'+iShift;
            x_suggested(2.*n_landmark+[-1 0]) = best_pixel;
            %             plot(best_pixel(1),best_pixel(2),'y.')
            
        end % Looping through landmarks
        
        % First you do Procrustes (Pose) then you deform (Shape)
        % Procrustes
        [~, x_posed] = procrustes(x_current,x_suggested);
        
        % Determine weights to adjust shape w/in model space
        % TODO: Do I somehow need to account for resolution scaling w/ the
        % eigenvectors and values?
        b = P'*(x_current-x_aligned); % What should the subtraction be between?
        dx = x_posed-x_current;
        db = P'*dx;
        b = b+db;
        
        % Limit the model parameters based on what is considered a nomal
        % contour, using the eigenvalues of the PCA-model
        
        % FROM THE MATHEXCHANGE ASM CODE
        b=max(min(b,maxb),-maxb);
        
        % Generate new shape within constraints of model
        %Wb = diag(D); % Eq. 26 (can also just be identity)
        x_new = x_aligned + P*(b);
        x_current = x_new; % Update the current shape position
        
        imshow(im,[]), plotLandmarks(x_new,'hold',1)
    end % End looping through evolutions
    
end % End looping through resolution levels

% Scale the final shape estimation to original image size
x_final = x_new*downsampleFactor;
x_original_estimate = x_original_estimate*downsampleFactor; % For final display

% Compare the original estimate with the final evolution
figure(gcf), hold off, imshow(im_original,[])
h_orig  = plotLandmarks(x_original_estimate,'hold',1);
h_final = plotLandmarks(x_final,'hold',1,'color','g');
legend([h_orig,h_final],{'Original','Final'},'location','nw','fontsize',FS)
title('Final shape','fontsize',FS)

end % End of main