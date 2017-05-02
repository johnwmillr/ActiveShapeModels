function varargout = xylim(xy)
% XYLIM sets the x and y axis to equal limits
%
% John W. Miller
% 01-May-2017

ax = gca;
if nargin ==0
    
    xLim = xlim(ax); yLim = ylim(ax);
    xy = [min([xLim(1) yLim(1)]) max([xLim(2) yLim(2)])];
end
xlim(xy), ylim(xy)

if nargout == 1
    varargout{1} = xy;
end

end % End of main