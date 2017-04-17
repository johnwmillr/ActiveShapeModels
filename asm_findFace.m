function [] = asm_findFace(im_original,x_aligned)
% ASM_FINDFACE
%
%	INPUT
%
%
%
%	OUTPUT
%
%
% John W. Miller
% 16-Apr-2017

%% Downsample
n_downsample = 3;
im = double(im_original(1:n_downsample:end,1:n_downsample:end));

% Filter (smooth image)
A = 1;
sig = .4;
mew = 1;
x = (-3*sig+mew):0.1:(3*sig+mew);
h = A*exp(-((x-mew).^2)./(2*sig^2));
im_filt = conv2(h,h,im,'same');
im_gMag = imgradient(im_filt); % Image gradient


% figure(1), imshowpair(im_filt,im_gMag)
% figure(2),imshow(im_gMag,[])
% plotLandmarks(x_aligned./n_downsample,'show_lines',1,'hold_on',1)

%% Iterate through landmarks
n_landmarks = length(x_aligned)/2;


x_evolving = x_aligned;

for n_evolutions = 1:3
    xy = [x_evolving(1:2:end) x_evolving(2:2:end)]./n_downsample;
    
    % Plot the mean shape on image
    figure, imshow(im_gMag,[]), hold on
    for n = 1:n_landmarks,plot(xy(n,1),xy(n,2),'go','linewidth',2);end;
    for n_landmark = 1:n_landmarks
        [m,p1] = calcNormalSlope(n_landmark,xy);
        
        % Point slope form of normal line through the current point
        y = @(p1,m,x) round((x-p1(1))*m+p1(2))'; % The output of this will be pixels (right?)
        
        c = -0.6:0.05:0.6;
        r = y(p1,m,c+p1(1));
        if range(r)/size(im,1) > 0.08
            if range(r)/size(im,1) > 0.2
                p = 1;
            elseif range(r)/size(im,1) > 0.1
                p = 2;
            else
                p = 3;
            end
            t = median(1:length(r));
            r = r(t-p:t+p);
            c = c(t-p:t+p);
        end
        c = round(p1(1)+c)';
        plot(c,r,'ro-','linewidth',2)
        
        % Find max value along normal
        idxs = sub2ind(size(im_gMag),r,c);
        [~,idx_max_val] = max(im_gMag(idxs));
        
        [I,J] = ind2sub(size(im),idxs(idx_max_val));
        x_evolving(n_landmark:n_landmark+1) = [J I];
    end
    pause
end














end % End of main




function [slope,p1] = calcNormalSlope(n,xy)
faceRegions = getFaceRegions();
if any(n==faceRegions{1})     % Left eye
    [p0,p1,p2] = deal(xy(1,:),xy(n,:),xy(3,:));
elseif any(n==faceRegions{2}) % Right eye
    [p0,p1,p2] = deal(xy(4,:),xy(n,:),xy(6,:));
elseif any(n==faceRegions{3}) % Left eyebrow
    [p0,p1,p2] = deal(xy(7,:),xy(n,:),xy(9,:));
elseif any(n==faceRegions{4}) % Right eyebrow
    [p0,p1,p2] = deal(xy(10,:),xy(n,:),xy(12,:));
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














