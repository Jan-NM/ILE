function [yPosition, xPosition] = findCorrelationPeak(correlationImage)
%findCorrelationPeak Fits gaussian function to correlation image
%   Detailed explanation goes here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set parameters
paddingFactor = 10; % zero-padding
cropFactor = 50;
% default cropFactor = 50; if you want to adjust the size of the cropped
% region for zero-padding in the fourier domain, adjust this value; minimum crop size is
% 50px
fittingMethod = 'lscurvefit'; % default: lscurvefit; alternative: 'fminsearch'
% initial guess for sigma -> x = -log(y/ymax) with y=ymwx/2 -> 0.69 =
% FWHM/2 -> transfer to sigma: FWHM ~2.35 * sigma
sigma = 0.59;
% size of the cropped region for fitting; default 15 px
fitRegion = 15;
%% zero-padding and initial fit parameter estimation
% cut out center, corelationImage should be a square
[szX, szY] = size(correlationImage);
if ~eq(szX, szY)
    error('Image is not squared!')
end
cropSize = szX / cropFactor;
cropSize = max(50, cropSize);
% crop correlation image
croppedImage = imcrop(correlationImage, [szX/2-cropSize/2 szX/2-cropSize/2 cropSize cropSize]);
% perform zero-padding in fourier domain
interpolatedImage = abs(ifft2((fftshift(fft2(croppedImage))), floor(cropSize * paddingFactor), floor(cropSize * paddingFactor)));
% get peak position and value
[yMax, xMax] = find(interpolatedImage == max(max(interpolatedImage)));
yMax = ceil(yMax / paddingFactor + (szX/2-cropSize/2));
xMax = ceil(xMax / paddingFactor + (szX/2-cropSize/2));
maxVal = correlationImage(yMax, xMax);
% estimate background
bgLevel = mean(mean(correlationImage));

fitVector = [maxVal, xMax, yMax, sigma,  sigma, bgLevel];
%% fitting part
data = imcrop(correlationImage, [xMax-fitRegion yMax-fitRegion 2*fitRegion 2*fitRegion]);

% create meshgrid with same size as data
[x, y] = meshgrid(1:size(data, 2), 1:size(data, 1));

offset = fitVector(2:3) - fitRegion - 1;

fitVector(2) = fitRegion + 1;
fitVector(3) = fitRegion + 1;
 
switch fittingMethod
    case 'lscurvefit'
        xData(:, :, 1) = x;
        xData(:, :, 2) = y;
        [pfit, ~] = lsqcurvefit(@GaussfitLSQ, fitVector, xData, data);
    case 'fminsearch'   
        options = optimset('MaxFunEvals', 500);
        [pfit, ~] = fminsearch(@GaussfitFMinSearch, fitVector, options);     
end

    function rval = GaussfitFMinSearch(pfitGuess)
        mygauss = pfitGuess(1)*exp(-((x-pfitGuess(2)).^2+(y-pfitGuess(3)).^2)/(pfitGuess(4)^2)) + pfitGuess(6);
        rval = sum(sum((data - mygauss).^2));
    end
xPosition = pfit(2) + offset(1, 1);
yPosition = pfit(3) + offset(1, 2);
end

function mygauss = GaussfitLSQ(pfitGuess, data)
    x = data(:, :, 1);
    y = data(:, :, 2);
    mygauss = pfitGuess(1)*exp(-((x-pfitGuess(2)).^2+(y-pfitGuess(3)).^2)/(pfitGuess(4)^2)) + pfitGuess(6);
end