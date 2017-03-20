% SHAPEMODEL_COOTESTAYLOR
%
%
%
%
%   Shape analysis techniques based on this paper:
%   Cootes, T. F., Taylor, C. J., Cooper, D. H., & Graham, J. (1995).
%       "Active Shape Models-Their Training and Application. Computer Vision and
%       Image Understanding."
%
%
%
% John W. Miller
% 16-Feb-2017


%% Faces
% Place landmarks
% landmarks = placeLandmarks(pathToImages,20,10);
addpath('./Utilities/')
load('./Landmarks/landmarks_faces_A_1-50')
alignedShapes = alignShapes(allLandmarks,0);
plotLandmarks(alignedShapes)

%% Create PCA model (on a subset of shapes)
x = alignedShapes(:,1:45);
xBar = mean(x,2);  % Mean shape
S = cov(x');       % Covariance matrix
[V,D] = eig(S);    % Eigenvectors
D = sort(diag(D),'descend');
V = fliplr(V);

%% GUI variations along the various PCs
guiPrinComps(xBar,V,D);

%% Get weights for first 3 PCs and plot them to check for independence between the PCs
P = V(:,1:3);
weights = P'*(alignedShapes-repmat(xBar,1,size(alignedShapes,2)));
figure
plot3(weights(1,:),weights(2,:),weights(3,:),'o')
xlabel('b1','fontsize',FS), ylabel('b2','fontsize',FS), zlabel('b3','fontsize',FS), grid on

%% Get weights for a new shape
n_pcs = 3;
y = alignedShapes(:,end);

% Solve for the weights that will approximate the new shape using the model
P = V(:,1:n_pcs);
b = P'*(y-xBar);

newShape = xBar + P*b;
figure, hold on
plot(y(1:2:end),y(2:2:end),'ko','linewidth',1)
plot(newShape(1:2:end),newShape(2:2:end),'r.','markersize',10,'linewidth',2)
set(gca,'ydir','reverse')
legend({'Real','Model'},'location','best')
title(sprintf('Approximation of new shape using first %d PCs from the model',n_pcs),'fontsize',FS)

%% Examine variations from individual PCs
plotPrinComp(V,D,xBar,1)

%% Edge detection using ASMs

% Probably going to need to do some image processing to enhance edges in the images
imDir = './Faces/faces_B';
imFile = 'B_49_0.jpg';
im = imread(fullfile(imDir,imFile));
imshow(im), hold on
plot(xBar(1:2:end),xBar(2:2:end),'ro','linewidth',2)

hFilt = fspecial('average',4*[1 1]);
im_filt = imfilter(im,hFilt);
im_gMag = imgradient(im_filt); % Image gradient
figure, imshow(im_gMag,[])


%% Connect dots around the face
faceRegions = getFaceRegions();

%% Loop through each landmark point, calculating the normal vector

n_landmarks = length(xBar)/2;
xy = [xBar(1:2:end) xBar(2:2:end)];
R = [0 1; -1 0]; % Rotate 90 degrees

% Adjustments along the normal (in image space)
dX = zeros(2*n_landmarks); % {dX0, dY0, dX1 dY1, ..., dXn-1 dYn-1)

% Plot the mean shape on image
figure, imshow(im_gMag,[]), hold on
for n = 1:n_landmarks,plot(xy(n,1),xy(n,2),'mo');end;

for n = 1:n_landmarks        
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
    vNorm = (p2-p0)*R;
    
    % Point slope form of normal line through the current point
    m = vNorm(2)/vNorm(1);
    y = @(p1,m,x) (x-p1(1))*m+p1(2); % The output of this will be pixels (right?)

    % Put the normal vector on the image
    cols = 1:size(im,2);
    rows = y(p1,m,cols);               
    
    mask_bad_cols = ~and(rows>1,rows<size(im,1));
    rows(mask_bad_cols) = [];
    cols(mask_bad_cols) = [];                            
    
    plot(cols,rows,'go')
    
    
    % Identify max edge value along line            
    pixels = im(round(rows),round(cols));
    [val,idx] = max(pixels(:));    
    pause
    
    
    %     dX(n,:) = [];
    
    
end


%%
% After looping through each point and looking along the normal, we'll end up with a
% set of suggested points, the same dimension as our mean shape. It's then just a
% simple matter of using Procrustes to translate the current points to the suggested
% points.
suggestedPoints = [];

currentPoints = procrustes(suggestedPoints,currentPoints,'scaling',0);

P = V(:,1:5);
dx; % adjustment vectors w/in the image space
db = P'*dx; % adjustments w/in the model space











