function [x_aligned, varargout] = placeShape(im,x)
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

% View the image
f = figure; hold on
imshow(im,[])
text(0.1,0.9,'Click on center of nose.','fontsize',FS,'units','normalized','color','r')

% Get input from user
[I,J] = ginput(1);

% Center shape on user's point
x_aligned = x;
x_aligned(1:2:end) = x_aligned(1:2:end)-x(27); % Center on middle of nose
x_aligned(2:2:end) = x_aligned(2:2:end)-x(28);
x_aligned(1:2:end) = x_aligned(1:2:end)+I;
x_aligned(2:2:end) = x_aligned(2:2:end)+J;

% Calculate transform from model space to image space (this may not be useful)
xa = [x_aligned(1:2:end) x_aligned(2:2:end)];
xo = [x(1:2:end) x(2:2:end)];
[~,~,T_model_to_image] = procrustes(xa,xo); 

% Display centered shape
plotLandmarks(x_aligned,'show_lines',1,'hold',1)

% Varargout
if nargout == 2
    varargout{1} = T_model_to_image;
elseif nargout == 3
    varargout{1} = T_model_to_image;
    varargout{2} = f;
    return
end

if nargout ~= 3
    pause(1), close(f)    
end

end % End of main