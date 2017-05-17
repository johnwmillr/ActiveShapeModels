function [x_aligned, f] = placeShape(im,x,layout)
% PLACESHAPE 
%
%	INPUT
%
%
%
%	OUTPUT
%
%
% John W. Miller
% 21-Apr-2017

if nargin < 3
    layout = 'standard';
end

% View the image
f = figure('units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
hold on, imshow(im,[],'InitialMagnification','fit')
text(0.07,0.95,'Click on center of nose. Or close by. Test your luck.','fontsize',FS,'units','normalized','color','r')

% Get input from user
[I,J] = ginput(1);

% Center shape on user's point
x_aligned = x;

switch lower(layout)
    case 'standard'
        x_aligned(1:2:end) = x_aligned(1:2:end)-x(27); % Center on middle of nose
        x_aligned(2:2:end) = x_aligned(2:2:end)-x(28);
    case 'muct'
        x_aligned(1:2:end) = x_aligned(1:2:end)-x(135); % Center on middle of nose
        x_aligned(2:2:end) = x_aligned(2:2:end)-x(136);
end
x_aligned(1:2:end) = x_aligned(1:2:end)+I;
x_aligned(2:2:end) = x_aligned(2:2:end)+J;

% Display centered shape
plotLandmarks(x_aligned,'show_lines',1,'hold',1,'layout',layout)

% Varargout
if nargout == 2 % User wants the figure handle
    return
else
    pause(1), close(f)    
end

end % End of main