function fd=dif2(f)
%dif2  Computes two point diference estimate of velocity.
%
%  DIFFERENCE = dif2(SIGNAL)
%
%  Computes the difference at a point as being the average of slopes on
%  either side of the point, returns a vector DIFFERENCE, the same size as
%  in the input SIGNAL
%
%  If SIGNAL is a matrix it takes the difference of the rows in each column
% 
%  To convert the output to dif2 to the proper units multiply by fs
%
% tags: math, calculus, slope
[r c] = size(f);

%May 04, 2010 : JAH, modified to keep dimensions of single vector

if isempty(f)
    fd = f;
    return
end

if r == 1 && c == 1
    error('Dif2 is undefined for a single point')
elseif r >1 && c > 1
    fd=zeros(r,c);
    fd(1,:)=f(2,:)-f(1,:);
    fd(r,:)=f(r,:)-f(r-1,:);
    df=diff(f);
    fd(2:r-1,:)=.5*(df(2:r-1,:)+df(1:r-2,:));  
else %The single vector case
    fd=zeros(r,c);
    fd(1)=f(2)-f(1);
    fd(end)=f(end)-f(end-1);
    df=diff(f);
    fd(2:end-1)=.5*(df(2:end)+df(1:end-1));  
end
end