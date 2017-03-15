function plotLandmarks(landmarks)
% PLOTLANDMARKS plots all of the aligned landmarks from an active shapes model
%
%	INPUT
%       landmarks: The aligned landmarks from multiple images [2*n_landmarks x n_shapes]
%
% John W. Miller
% 14-Mar-2017

% Mean shape
n_shapes = size(landmarks,2);
if n_shapes > 1        
    meanShape = mean(landmarks,2); % x1, y1, x2, y2, ..., x20, y20        
else
    meanShape = landmarks;
end

% Plot the landmarks for each shape
h = figure; hold on
try colors = parula(n_shapes);
catch
    colors = hsv(n_shapes);
end
for n_shape = 1:n_shapes        
    iShape = [landmarks(1:2:end,n_shape) landmarks(2:2:end,n_shape)];
    plot(iShape(:,1),iShape(:,2),'o','color',colors(n_shape,:),'linewidth',1,'markersize',3)    
end
ax = plot(meanShape(1:2:end),meanShape(2:2:end),'ko','markersize',5,'linewidth',1,'markerfacecolor','k');
legend(ax,{'Mean shape'},'fontsize',FS,'location','northwest')

% Touch up the plot
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]) 
axis off;
set(gca,'YDir','reverse'); 
text(0.05,0.5,sprintf('n=%d',n_shapes),'units','normalized','fontsize',FS,'fontweight','bold')

end % End of main