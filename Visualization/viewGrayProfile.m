function viewGrayProfile(grayModel,nr,nl)
% VIEWGRAYPROFILE 
%
%	INPUT
%       grayModel: Comes from BUILDGRAYLEVELMODEL.
%       nr: resolution number
%       nl: landmark number
%
%   See also BUILDGRAYLEVELMODEL
%
% John W. Miller
% 15-May-2017

% Reshape from vector to square
graySquare = reshape(grayModel(nr).meanProfile{nl},grayModel(nr).Info.rc_squaresize-1,grayModel(nr).Info.rc_squaresize-1);

% View the square region
figure, imshow(imresize(graySquare,round(200/length(graySquare))),[])

end % End of main