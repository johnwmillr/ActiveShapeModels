function guiPrinComps(xBar,V,D)
% GUIPRINCOMPS
%
%	INPUT
%       xBar: Mean shape [2*n_landmarks x 1]
%       V: Principal components (eigenvectors)
%       D: Shape weights (eigenvalues)
%
%	OUTPUT
%
%
% John W. Miller
% 17-Mar-2017

% https://www.mathworks.com/help/control/ug/build-app-with-interactive-plot-updates.html
% https://www.mathworks.com/help/matlab/ref/uicontrol.html

% Connect dots around the face
faceRegions = getFaceRegions();

% Examine variations from individual PCs
n_pcs = 5;
n_variations = 5; %
if mod(n_variations,2) == 0
    n_variations = n_variations + 1;
end
var_step = (n_variations-1)/2;
weights = zeros(var_step,n_variations);
for n = 1:n_pcs
    weights(n,:) = sqrt(D(n))*(-var_step:var_step);
end

% Create some shape variations
n_pc = 2;
P = V(:,n_pc);
shapeVariations = repmat(xBar,1,n_variations) + P*weights(n_pc,:);


%% Generate some shapes
P = V(:,1:n_pcs);
mask_weights = [1]
newShape = xBar + P*weights(:,mask_weights);





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
for i = 1:length(faceRegions)
    h = plot(iVar(faceRegions{i},1),iVar(faceRegions{i},2), '-','linewidth',2,'color','g');
end

set(gca,'xlim',xLim,'ylim',yLim,'ydir','reverse')

% Create slider
min_val = 1;
max_val = n_variations;
b1 = uicontrol('Parent',f,'Style','slider','Position',[81,54,419,23],...
    'value',var_step+1, 'min',min_val, 'max',max_val, 'callback', @b_callback_pc1,...
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
        for i = 1:length(faceRegions)
            h = plot(iVar(faceRegions{i},1),iVar(faceRegions{i},2), '-','linewidth',2,'color','g');            
            set(gca,'xlim',xLim,'ylim',yLim,'ydir','reverse')
        end        
    end



end % End of main