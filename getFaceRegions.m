function faceRegions = getFaceRegions(layout)
% GETFACEREGIONS returns a cell array grouping the different regions of a face shape
% (e.g. left eye, right eye, nose, etc.)
%
%	OUTPUT
%       faceRegions: (cell) Arrays corresponding to different regions of the face
%
%
% John W. Miller
% 19-Mar-2017

if nargin == 0
    layout = 'standard';    
end

switch lower(layout)
    case 'standard'
        
        % Connect dots around the face
        faceRegions = cell(7,1);
        faceRegions{1} = 1:3;         % Left eye
        faceRegions{2} = 4:6;         % Right eye
        faceRegions{3} = 7:9;         % Left eyebrow
        faceRegions{4} = 10:12;       % Right eyebrow
        faceRegions{5} = 13:15;       % Nose
        faceRegions{6} = [16:19 16];  % Mouth
        faceRegions{7} = 20;          % Chin
    case 'nobrows'
        
        % Connect dots around the face
        faceRegions = cell(9,1);
        faceRegions{1} = 1:3;         % Left eye
        faceRegions{2} = 4:6;         % Right eye
        faceRegions{3} = 7;           % Left side of head
        faceRegions{4} = 8:9;         % Left eyebrow
        faceRegions{5} = 10:11;       % Right eyebrow
        faceRegions{6} = 12;          % Right side of head
        faceRegions{7} = 13:15;
        faceRegions{8} = [16:19 16];
        faceRegions{9} = 20;
    
    case 'ilm' % For the ILM layer in OCT scans
        faceRegions{1} = 1:31;
    case 'both'
        faceRegions{1} = 1:10;
        faceRegions{2} = 11:20;
        faceRegions{3} = 21:51;
    case 'muct'
        faceRegions{1} = 1:76;
        
end


end % End of main