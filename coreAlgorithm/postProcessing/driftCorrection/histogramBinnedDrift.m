function [ templateImage ] = histogramBinnedDrift(positions, pixelSize, xSize, ySize, minX, minY)
%histogramBinnedDrift Binning of subset SMLM images for drift correction
%
% Jan Neumann, 26.09.2017, modified 05.03.18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1 = round(positions(:,3)./pixelSize);
y1 = round(positions(:,2)./pixelSize);
templateImage = sparse(x1-minX+1, y1-minY+1, 1, xSize, ySize);
templateImage = full(templateImage);
end


                