function [] = guiPrinComps(V,D,xBar)
% GUIPRINCOMPS
%
%	INPUT
%       V: Principal components (eigenvectors)
%       D: Shape weights (eigenvalues)
%       xBar: Mean shape [2*n_landmarks x 1]
%
%	OUTPUT
%
%
% John W. Miller
% 17-Mar-2017

%https://www.mathworks.com/help/control/ug/build-app-with-interactive-plot-updates.html

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
n_variations = 7;
weights = zeros(3,n_variations);
for n = 1:3
    weights(n,:) = sqrt(D(n))*(-3:3);
end

% Create some shape variations
n_pcs = 2;
P = V(:,n_pcs);
shapeVariations = repmat(xBar,1,n_variations) + P*weights(n_pcs,:);

% https://www.mathworks.com/help/matlab/ref/uicontrol.html

%% Create the GUI
% Create a figure and axes
f = figure();
ax = axes('Parent',f);%,'position',[0.13 0.39  0.77 0.54]);

iVar = zeros(size(shapeVariations,1)/2,2);
iVar(:,1) = shapeVariations(1:2:end,n);
iVar(:,2) = shapeVariations(2:2:end,n);

% Plot the PC variations
mew(:,1) = xBar(1:2:end);
mew(:,2) = xBar(2:2:end);
xLim = [0.8*min(min(shapeVariations(1:2:end,:))) 1.1*max(max(shapeVariations(1:2:end,:)))]';
yLim = [0.8*min(min(shapeVariations(2:2:end,:))) 1.1*max(max(shapeVariations(2:2:end,:)))]';


plot(mew(:,1),mew(:,2),'o','color','k','linewidth',2); hold on
% Connect the dots
for i = 1:length(faceLabels)
    h = plot(iVar(faceLabels{i},1),iVar(faceLabels{i},2), '-','linewidth',2,'color','g');
end

set(gca,'xlim',xLim,'ylim',yLim,'ydir','reverse')

% Create slider
min_val = 1;
max_val = 7;
b1 = uicontrol('Parent',f,'Style','slider','Position',[81,54,419,23],...
    'value',4, 'min',min_val, 'max',max_val, 'callback', @b_callback_pc1,...
    'SliderStep', [1/max_val 1/max_val]);

% b2 = uicontrol('Parent',f,'Style','slider','Position',[81,54,419,23],...
%     'value',4, 'min',min_val, 'max',max_val, 'callback', @b_callback_pc2,...
%     'SliderStep', [1/max_val 1/max_val]);
% 
% b3 = uicontrol('Parent',f,'Style','slider','Position',[81,54,419,23],...
%     'value',4, 'min',min_val, 'max',max_val, 'callback', @b_callback_pc3,...
%     'SliderStep', [1/max_val 1/max_val]);


%% Callback functions
    function b_callback_pc1(source,event)
        val = round(get(source,'value'));
        iVar = zeros(size(shapeVariations,1)/2,2);
        iVar(:,1) = shapeVariations(1:2:end,val);
        iVar(:,2) = shapeVariations(2:2:end,val);       
        
        hold off
        plot(mew(:,1),mew(:,2),'o','color','k','linewidth',2), hold on
        for i = 1:length(faceLabels)
            h = plot(iVar(faceLabels{i},1),iVar(faceLabels{i},2), '-','linewidth',2,'color','g');            
            set(gca,'xlim',xLim,'ylim',yLim,'ydir','reverse')
        end        
    end



end % End of main