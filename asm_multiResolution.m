function x_aligned = asm_multiResolution(im_original)
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


%%
% imDir = './Faces/faces_B';
% imFile = 'B_40_0.jpg';
% im_original = imread(fullfile(imDir,imFile));

% Downsample
n = 8;
im = double(im_original(1:n:end,1:n:end));

% Filter (smooth image)
A = 1;
sig = .5;
mew = 1;
x = (-3*sig+mew):0.1:(3*sig+mew);
h = A*exp(-((x-mew).^2)./(2*sig^2));
im_filt = conv2(h,h,im,'same');
% figure,imshow(im_filt,[])

% Convolve with head shape
hair_shape = ones(round(size(im).*[1/3 1/5]));
hair_response = conv2(1./im_filt,hair_shape,'same');
rc = 8;
hair_response(:,1:rc) = 0;
hair_response(:,end-rc:end) = 0;
hair_response(1:rc,:) = 0;
hair_response(end-rc:end,:) = 0;
[~,idx] = max(hair_response(:));
[I,J] = ind2sub(size(im),idx);
% figure(1), imshow(im,[]); hold on
% plot(J,I,'ro'), hold off


% Determine where to place mean shape
%
figure(1), imshow(im_original,[]), hold on
plot(J*n,I*n,'ro')
I = I + round(size(hair_shape,2)/2); % KLUDGE (try to move from hair to eyes)
plot(J*n,I*n,'ro')
x = mean(alignedShapes,2);

[x1,y1] = deal(x(17),x(18));
[x2,y2] = deal(x(19),x(20));
[xM, yM] = deal(mean([x1 x2]),mean([y1 y2])); % Top middle of face
subMean = zeros(size(x));
subMean(1:2:end) = xM;
subMean(2:2:end) = yM;
x_aligned = x-subMean;
x_aligned(1:2:end) = x_aligned(1:2:end) + J*n; % Scale back up the detected region and shift
x_aligned(2:2:end) = x_aligned(2:2:end) + I*n;

% Add mean shape to image
plotLandmarks(x_aligned,'show_lines',1,'hold_on',1)













end % End of main