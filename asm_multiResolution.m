function [x_aligned, f] = asm_multiResolution(im_original,xBar,h_filt)
% ASM_MULTIRESOLUTION 
%
%	INPUT
%
%
%
%	OUTPUT
%
%
% John W. Miller
% 12-Apr-2017

% Filter (smooth image)
if nargin < 3
    A = 1;
    sig = .3; % Need to play w/ these parameters
    mew = 1;
    x = (-3*sig+mew):0.1:(3*sig+mew);
    h_filt = A*exp(-((x-mew).^2)./(2*sig^2));
end
im_filt = conv2(h_filt,h_filt,im_original,'same');

% Downsample
n = 1;
im = double(im_filt(1:n:end,1:n:end));

% Convolve with head shape
hair_shape = ones(round(size(im).*[1/5 1/3]));
hair_response = conv2(1./im,hair_shape,'same');
rc = 8; % rows and columns to cut off the edge off the image (Kludge)
hair_response(:,1:rc) = 0;
hair_response(:,end-rc:end) = 0;
hair_response(1:rc,:) = 0;
hair_response(end-rc:end,:) = 0;
[~,idx] = max(hair_response(:));
[I,J] = ind2sub(size(im),idx);

% Determine where to place mean shape
I = I + round(size(hair_shape,2)/2.5); % KLUDGE (try to move from hair to eyes)

[x1,y1] = deal(xBar(17),xBar(18));
[x2,y2] = deal(xBar(19),xBar(20));
[xM, yM] = deal(mean([x1 x2]),mean([y1 y2])); % Top middle of face
subMean = zeros(size(xBar));
subMean(1:2:end) = xM;
subMean(2:2:end) = yM;
x_aligned = xBar-subMean;
x_aligned(1:2:end) = x_aligned(1:2:end) + J*n; % Scale back up the detected region and shift
x_aligned(2:2:end) = x_aligned(2:2:end) + I*n;

% Add mean shape to image
f = figure('units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
hold on, imshow(im_original,[],'InitialMagnification','fit')
plotLandmarks(x_aligned,'show_lines',1,'hold',1);

end % End of main