function plotPrinComp(V,D,xBar,n_pc)
% PLOTPRINCOMP displays the variation along a single principal component from an
% active shapes model
%
%	INPUT
%       V: Principal components (eigen vectors)
%       D: Shape weights (eigenvalues)
%       xBar: Mean shape [2*n_landmarks x 1]
%       n_pc: Which PC do you want to display variations along?
%
%   The shape variations are determined by the equation x = xBar + P*b,
%   where P is a single principal component and b is a vector of weights.
%
%   See also PLACELANDMARKS, PLOTLANDMARKS, ALIGNSHAPES
%
% John W. Miller
% 16-Mar-2017

% Connect dots around the face
faceLabels = cell(7,1);
faceLabels{1} = 1:3;
faceLabels{2} = 4:6;
faceLabels{3} = 7:9;
faceLabels{4} = 10:12;
faceLabels{5} = 13:15;
faceLabels{6} = [16:19 16];
faceLabels{7} = 20;

% Examine variations from individual PCs
b = sqrt(D(n_pc))*(-3:3);
n_vars = length(b);

% Create some shape variations
P = V(:,n_pc);
shapeVariations = repmat(xBar,1,n_vars) + P*b;
xLim = [floor(min(min(shapeVariations(1:2:end,:)))) ceil(max(max(shapeVariations(1:2:end))))];
yLim = [floor(min(min(shapeVariations(2:2:end,:)))) ceil(max(max(shapeVariations(2:2:end))))];

% Change plot color for different weights of the selected PC
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
























end % End of main