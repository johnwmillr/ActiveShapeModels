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
% [landmarks,idxs] = placeLandmarks(pathToImages,20,10);
% load('landmarks_faces_set01')
load('./Landmarks/landmarks_faces_A_23')
landmarksTraining = allLandmarks;
alignedShapes = alignShapes(landmarksTraining,0);
plotLandmarks(alignedShapes)

% Create PCA model (on a subset of shapes)
x = alignedShapes(:,1:20);
xBar = mean(x,2);  % Mean shape
S = cov(x');       % Covariance matrix
[V,D] = eig(S);    % Eigenvectors
D = sort(diag(D),'descend');
V = fliplr(V);

%% Get weights for first 3 PCs and plot them to check for independence between the PCs
P = V(:,1:3);
weights = P'*(alignedShapes-repmat(xBar,1,size(alignedShapes,2)));
figure
plot3(weights(1,:),weights(2,:),weights(3,:),'o')
xlabel('b1','fontsize',FS), ylabel('b2','fontsize',FS), zlabel('b3','fontsize',FS), grid on


%% Get weights for a new shape
n_pcs = 3;
y = alignedShapes(:,21);

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


%% Connect dots around the face
faceLabels = cell(7,1);
faceLabels{1} = 1:3;
faceLabels{2} = 4:6;
faceLabels{3} = 7:9;
faceLabels{4} = 10:12;
faceLabels{5} = 13:15;
faceLabels{6} = [16:19 16];
faceLabels{7} = 20;


%% Examine variations from individual PCs
n_pc = 3;
b = sqrt(D(n_pc))*(-3:3);
n_vars = length(b);

% Create some shape variations
P = V(:,n_pc);
shapeVariations = repmat(xBar,1,n_vars) + P*b;
xLim = [floor(min(min(shapeVariations(1:2:end,:)))) ceil(max(max(shapeVariations(1:2:end))))];
yLim = [floor(min(min(shapeVariations(2:2:end,:)))) ceil(max(max(shapeVariations(2:2:end))))];


% Color variation
mew(:,1) = xBar(1:2:end);
mew(:,2) = xBar(2:2:end);
figure, hold on
colors = hsv(n_vars);
for n = 1:n_vars    
    iVar = zeros(size(shapeVariations,1)/2,2);
    iVar(:,1) = shapeVariations(1:2:end,n);
    iVar(:,2) = shapeVariations(2:2:end,n);
    
    % Plot the PC variations
    plot(iVar(:,1),iVar(:,2),'o','color',colors(n,:))
    
    % Connect the dots
    for i = 1:length(faceLabels)
        plot(mew(faceLabels{i},1), mew(faceLabels{i},2), 'k-','linewidth',1)
        plot(iVar(faceLabels{i},1),iVar(faceLabels{i},2), '-','linewidth',1,'color',colors(n,:))
    end
    
end
plot(xBar(1:2:end),xBar(2:2:end),'k.','linewidth',3)
set(gca,'ydir','reverse'), axis square
xlim(xLim), ylim(yLim)
title(sprintf('Variation of PC #%d',n_pc),'fontsize',20)

%% Edge detection using ASMs

% Probably going to need to do some image processing to enhance edges in the images
imDir = './Images/faces_B';
imFile = 'B_49_0.jpg';
im = imread(fullfile(imDir,imFile));
imshow(im), hold on
plot(xBar(1:2:end),xBar(2:2:end),'ro','linewidth',2)

% Loop through each landmark point, calculating the normal vector
n_points = length(xBar);
xy = [xBar(1:2:end) xBar(2:2:end)];
R = [0 -1 1 0]; % Rotate 90 degrees
for n = 2:n_points
    
    
    
    
    plot(xy([n-1 n+1],1),xy([n-1 n+1],2),'o-'), hold on
    a = R*[xy([n-1 n+1],1),xy([n-1 n+1],2)];
    plot(a(:,1),a(:,2),'ro-')

end











%%
% bw = edge(im,'canny',.1);
bw = edge(im,'zerocross');

imshow(bw), hold on
% plot(xBar(1:2:end),xBar(2:2:end),'ro','linewidth',2)

%%
im_thresh = im2bw(im,graythresh(im));
figure, imshow(im_thresh)

[gMag, gDir] = imgradient(im_thresh); % Image gradient
figure, imshow(gMag)








