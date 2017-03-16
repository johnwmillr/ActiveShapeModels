function plotLandmarks(landmarks)
% PLOTLANDMARKS plots all of the aligned landmarks from an active shapes model
%
%	INPUT
%       landmarks: The aligned landmarks from multiple images [2*n_landmarks x n_shapes]
%
% John W. Miller
% 14-Mar-2017

% Face points
faceLabels = cell(7,1);
faceLabels{1} = 1:3;
faceLabels{2} = 4:6;
faceLabels{3} = 7:9;
faceLabels{4} = 10:12;
faceLabels{5} = 13:15;
faceLabels{6} = [16:19 16];
faceLabels{7} = 20;


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
    plot(iShape(:,1),iShape(:,2),'o','color',...
        colors(n_shape,:),'linewidth',2,'markersize',5,...
    'markerfacecolor',colors(n_shape,:))    
    
end
% 
% Connect the dots
mew = [meanShape(1:2:end) meanShape(2:2:end)];
for i = 1:length(faceLabels)
    plot(mew(faceLabels{i},1), mew(faceLabels{i},2), 'k-','linewidth',1)    
end

ax = plot(meanShape(1:2:end),meanShape(2:2:end),'ko','markersize',5,'linewidth',3,'markerfacecolor','k');
legend(ax,{'Mean shape'},'fontsize',FS,'location','northwest')

% Touch up the plot
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]) 
axis off;
set(gca,'YDir','reverse'); 
text(0.05,0.5,sprintf('n=%d',n_shapes),'units','normalized','fontsize',FS,'fontweight','bold')

end % End of main