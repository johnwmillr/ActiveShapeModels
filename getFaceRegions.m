function faceRegions = getFaceRegions()
% GETFACEREGIONS returns a cell array grouping the different regions of a face shape
% (e.g. left eye, right eye, nose, etc.)
%
%	OUTPUT
%       faceRegions: (cell) Arrays corresponding to different regions of the face
%
%
% John W. Miller
% 19-Mar-2017

% Connect dots around the face
faceRegions = cell(7,1);
faceRegions{1} = 1:3;         % Left eye
faceRegions{2} = 4:6;         % Right eye
faceRegions{3} = 7:9;         % Left eyebrow
faceRegions{4} = 10:12;       % Right eyebrow
faceRegions{5} = 13:15;
faceRegions{6} = [16:19 16];
faceRegions{7} = 20;

end % End of main