function  newroi  = cl_seperator( roi )
%function newroi  = cl_seperator( roi )
%   cl_seperators attemps to separate the signal from a second, close-by signal
%   A second run in both directions is nessecary for increased accuracy
%   
%   The function assumes that the true center of the signal is in the area of the pixel in the middle
%
%created by: Manfred Kirchgessner < Manfred.Kirchgessner@web.de>, 
%            Frederik Grüll <Frederik.Gruell@kip.uni-heidelberg.de>%
%%

[DimY DimX]=size(roi);
RadX = (DimX+1)/2;
RadY = (DimY+1)/2;

noise = sqrt(roi);
n_roi = roi + noise;

xframe = [n_roi(:,2:RadX-1),n_roi(:,RadX-1:RadX+1),n_roi(:,RadX+1:DimX-1)];
roi = roi.*(xframe>roi);

yframe = [n_roi(2:RadY-1,:);n_roi(RadY-1:RadY+1,:);n_roi(RadY+1:DimY-1,:)];
roi = roi.*(yframe>roi);

xframe = [n_roi(:,2:RadX-1),n_roi(:,RadX-1:RadX+1),n_roi(:,RadX+1:DimX-1)];
roi = roi.*(xframe>roi);

newroi=roi;

end

