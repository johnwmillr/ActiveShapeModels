function guiPrinComps(xBar,V,D,varargin)
% GUIPRINCOMPS
%
%	INPUT
%       xBar: Mean shape [2*n_landmarks x 1]
%       V: Principal components (eigenvectors)
%       D: Shape weights (eigenvalues)
%
% John W. Miller
% 17-Mar-2017

% https://www.mathworks.com/help/control/ug/build-app-with-interactive-plot-updates.html
% https://www.mathworks.com/help/matlab/ref/uicontrol.html

keys = {'show_image'}; default_values = {0};
[show_image] = parseKeyValuePairs(varargin,keys,default_values);

% Initialization
n_pcs = 3;
slider_values = zeros(n_pcs,1);
faceRegions = getFaceRegions(); % Connecting dots around the face

% Examine variations from individual PCs
n_variations = 11; % Must be an odd number
if mod(n_variations,2) == 0
    n_variations = n_variations + 1;
end
var_step = (n_variations-1)/2;
weights = zeros(var_step,n_variations);
for n = 1:n_pcs
    weights(n,:) = sqrt(D(n))*(-var_step:var_step);
end

% Create some shape variations
P = V(:,1:n_pcs);
shapeVariations = repmat(xBar,1,n_variations) + P*weights(1:n_pcs,:);

%% Initial visualization
f = figure();
if show_image    
    imDir = fullfile(fileparts(which(mfilename)),'/Faces/faces_B');
    imFile = 'B_49_0.jpg';
    im = imread(fullfile(imDir,imFile));
    imshow(im), hold on
end

% Determine mean shape and plot dimensions
mew(:,1) = xBar(1:2:end);
mew(:,2) = xBar(2:2:end);
xLim = [0.8 0; 0 1.1]*[min(min(shapeVariations(1:2:end,:))) max(max(shapeVariations(1:2:end,:)))]';
yLim = [0.8 0; 0 1.3]*[min(min(shapeVariations(2:2:end,:))) max(max(shapeVariations(2:2:end,:)))]';

updatePlot(generateShape()); % Initial plot

%% Create the GUI sliders
[min_val,max_val] = deal(-(n_variations-1)/2,(n_variations-1)/2);
stepsize = (1/range([min_val max_val]))*[1 1];

b1 = uicontrol('Parent',f,'Style','slider','Position',[81,75,419,23],...
    'value',0, 'min',min_val, 'max',max_val, 'callback', @b_callback_pc1,...
    'SliderStep', stepsize);

b2 = uicontrol('Parent',f,'Style','slider','Position',[81,50,419,23],...
    'value',0, 'min',min_val, 'max',max_val, 'callback', @b_callback_pc2,...
    'SliderStep', stepsize);

b3 = uicontrol('Parent',f,'Style','slider','Position',[81,25,419,23],...
    'value',0, 'min',min_val, 'max',max_val, 'callback', @b_callback_pc3,...
    'SliderStep', stepsize);

%% Callback functions
    function b_callback_pc1(source,event) %#ok<*INUSD>
        slider_values(1) = round(get(source,'value'));
        updatePlot(generateShape())
    end

    function b_callback_pc2(source,event)
        slider_values(2) = round(get(source,'value'));
        updatePlot(generateShape())
    end

    function b_callback_pc3(source,event)
        slider_values(3) = round(get(source,'value'));
        updatePlot(generateShape())
    end

%% Generate and plot new shapes

    function newShape = generateShape(newVal)        
        % Calculate the weights for the updated shape
        b = zeros(n_pcs,1);
        for n_pc = 1:n_pcs
            b(n_pc,1) = sqrt(D(n_pc))*(slider_values(n_pc));
        end
        
        % Generate the shape
        x = xBar + P*b;
        newShape(:,1) = x(1:2:end);
        newShape(:,2) = x(2:2:end);
    end

    function updatePlot(newShape)
        hold off
        if show_image imshow(im), hold on, end
        for i = 1:length(faceRegions)
            plot(newShape(faceRegions{i},1),newShape(faceRegions{i},2), '.-','linewidth',2,'color','g'), hold on
        end
        plot(mew(:,1),mew(:,2),'.','color','k','linewidth',2)
        set(gca,'xtick',[],'ytick',[])
        set(gca,'xlim',xLim,'ylim',yLim,'ydir','reverse')
        title('Move the sliders to change the weights on the first 3 PCs','fontsize',FS)
    end

end % End of main