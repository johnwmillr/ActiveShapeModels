function [] = asm_findFace(im_original,x_bar_aligned,V,D)
% ASM_FINDFACE
%
%	INPUT
%       V: eigenvectors
%       D: eigenvalues
%       x_bar_aligned: Mean shape, roughly aligned over face via multi resolution
%
%	OUTPUT
%
%
% John W. Miller
% 16-Apr-2017

% Model stuff
n_pcs = 4;
P = V(:,1:n_pcs);
D = D(1:n_pcs);
b = zeros(n_pcs,length(-3));
for n_pc = 1:n_pcs
    b(n_pc,:) = sqrt(D(n_pc))*(-3);
end

%% Downsample
n_downsample = 2;
im = double(im_original(1:n_downsample:end,1:n_downsample:end));
x_bar_aligned = x_bar_aligned./n_downsample;

save_video = 0;
if save_video
    close all
    videoFileName = [go('down') filesep sprintf('ASM_FaceDetection_%s_DownampleBy_%02d',date,n_downsample)];
    vidObj = VideoWriter(videoFileName,'MPEG-4'); % Video File
    vidObj.FrameRate = 5;
    open(vidObj);
    disp('Recording video...')
end

% Filter (smooth image)
A = 1;
sig = .4;
mew = 1;
x = (-3*sig+mew):0.1:(3*sig+mew);
h = A*exp(-((x-mew).^2)./(2*sig^2));
im_filt = conv2(h,h,im,'same');
im_gMag = imgradient(im_filt); % Image gradient

% Crop the image (Kludge, dealing with bright ring around edge due to (?)
% convolution)
rc = round(0.05*length(im_gMag)); % rows and columns to cut off the edge off the image (Kludge)
im_gMag(:,1:rc) = 0;
im_gMag(:,end-rc:end) = 0;
im_gMag(1:rc,:) = 0;
im_gMag(end-rc:end,:) = 0;

%% Iterate through landmarks
n_landmarks = length(x_bar_aligned)/2;
[x_new, x_suggested] = deal(x_bar_aligned);

% Point slope form of normal line through the current point
y = @(p1,m,x) round((x-p1(1))*m+p1(2))'; % The output of this will be pixels (right?)

n_evolutions = 50;
for n_evolution = 1:n_evolutions
    xy = [x_new(1:2:end) x_new(2:2:end)];
    x_prev = x_new;
    
    % Plot the mean shape on image
    h = figure(3); imshow(im_gMag,[]), hold on
    %     for n = 1:n_landmarks,plot(xy(n,1),xy(n,2),'go','linewidth',2);hold on;end;
    % Connect the dots
    plotLandmarks(x_prev,'show_lines',1,'hold',1)
    text(0.2,0.80,sprintf('Iteration #: %d',n_evolution),'units','normalized','color','r','fontsize',FS)
    
    for n_landmark = 1:n_landmarks
        [m,p1] = calcNormalSlope(n_landmark,xy); % Internal function
        
        c = p1(1)+(-.6:0.05:.6);
        r = y(p1,m,c);
        if range(r)/size(im,1) > 0.08 % Subtract pixels
            if range(r)/size(im,1) > 0.2
                p = 3;
            elseif range(r)/size(im,1) > 0.1
                p = 4;
            else
                p = 5;
            end
            t = median(1:length(r));
            r = r(t-p:t+p);
            c = c(t-p:t+p);
        elseif range(r)/size(im,1) < 0.01 % Add pixels
            p = 3;
            c = [min(round(c))+(-p:1:-1) c max(round(c))+(1:p)];
            r = y(p1,m,c);
        end
        
        if any(isinf(r)) % Vertical line
            r = p1(2)+(-3:3)';
            c = p1(1)*ones(size(r));
        else
            c = round(c)';
        end
        
        %         figure(3), plot(c,r,'r.','linewidth',2)
        
        % Make sure r and c are within image
        r = max(r,1);
        c = max(c,1);
        r = min(r,size(im_gMag,1));
        c = min(c,size(im_gMag,2));
        
        % Find max value along normal
        idxs_along_norm = sub2ind(size(im_gMag),r,c);
        [~,idx_max_val] = max(im_gMag(idxs_along_norm));
        
        [I,J] = ind2sub(size(im),idxs_along_norm(idx_max_val));
        x_suggested((2*n_landmark-1):(2*n_landmark)) = [J I];
        %         figure(3), plot(J,I,'co','linewidth',2,'markerfacecolor','c'), hold on
    end
    
    % First you do Procrustes (Pose) then you deform (Shape)
    % Procrustes
    [~, x_posed] = procrustes(x_prev,x_suggested);
    
    % Determine weights to adjust shape w/in model space
    b = P'*(x_prev-x_bar_aligned);
    dx = x_posed-x_prev;
    db = P'*dx;
    
    % Generate new shape
    %Wb = diag(D); % Eq. 26 (can also just be identity)
    x_new = x_bar_aligned + P*(b+db);
    
    if save_video && ishandle(h)
        try vidFrame = getframe(h);
            writeVideo(vidObj,vidFrame);
        catch
            close(vidObj)
        end        
    end
end % End looping through evolutions of ASM


if save_video    
    close(vidObj)
    disp('Video closed.')
end

end % End of main

%% ----------------------------------------- %
%           INTERNAL FUNCTIONS               %
%  ----------------------------------------- %

function [slope,p1] = calcNormalSlope(n,xy)
faceRegions = getFaceRegions();
if any(n==faceRegions{1})     % Left eye
    [p0,p1,p2] = deal(xy(1,:),xy(n,:),xy(3,:));
elseif any(n==faceRegions{2}) % Right eye
    [p0,p1,p2] = deal(xy(4,:),xy(n,:),xy(6,:));
elseif any(n==faceRegions{3}) % Left eyebrow
    if n==7
        [p0,p1,p2] = deal(xy(1,:),xy(n,:),xy(8,:));
    else
        [p0,p1,p2] = deal(xy(7,:),xy(n,:),xy(9,:));
    end
elseif any(n==faceRegions{4}) % Right eyebrow
    if n==12
        [p0,p1,p2] = deal(xy(11,:),xy(n,:),xy(6,:));
    else
        [p0,p1,p2] = deal(xy(10,:),xy(n,:),xy(12,:));
    end
elseif any(n==faceRegions{5}) % Nose
    [p0,p1,p2] = deal(xy(13,:),xy(n,:),xy(15,:));
elseif any(n==faceRegions{6}) % Mouth
    switch n
        case 16
            [p0,p1,p2] = deal(xy(19,:), xy(n,:),xy(n+1,:));
        case 17
            [p0,p1,p2] = deal(xy(n-1,:),xy(n,:),xy(n+1,:));
        case 18
            [p0,p1,p2] = deal(xy(n-1,:),xy(n,:),xy(n+1,:));
        case 19
            [p0,p1,p2] = deal(xy(n-1,:),xy(n,:),xy(16,:));
    end
elseif any(n==faceRegions{7}) % Chin
    [p0,p1,p2] = deal(xy(16,:),xy(n,:),xy(18,:));
end
% Normal vector
R = [0 1; -1 0]; % Rotate 90 degrees
vNorm = (p2-p0)*R;
slope = vNorm(2)/vNorm(1);
end