function crossCorr = driftfindshift(img1,img2)
% resShift = driftfindshift(img1,img2)
%
% This function is called by the fastCorrectDriftLocNS and is not intended
% for manual usage.
% driftfindshift determines the respective shift between the two input
% images by a cross-correlation and fitting the peak.
% 
% Input:
%  img1, img2:  The two images that are to be correlated. Both must be
%               equal in size and of class dip_image.
%
% Output:
%  resShift:    Respective shift between the two images in measures of pixels.
%
% Kirti Prakash
% k.prakash@imb-mainz.de
% Martin Hagmann
% martin.hagmann@kip.uni-heidelberg.de
% 07 May 2014, Cremer Group
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%determine dimensionality
dimension = length(size(img1));
if length(size(img2)) ~= dimension; error('dimensionality of input images differ'); end

%convert into matlab to exploit the amazingly fast fftn
switch dimension
    case 2
        crossCorr = abs(fftshift(ifftn(fftn(img1).*conj(fftn(img2))))); 
    case 3
        crossCorr = abs(fftshift(ifftn(fftn(img1).*conj(fftn(img2))))); 
    otherwise
        warning('Unspecified dimension. Implement me!')
        keyboard;
end
% 
% switch dimension
%     case 2
%         sigma = max(size(img1))/5;
%         kernel = permute(dip_image(fspecial('gaussian',size(img1),sigma)),[2 1]);
%     case 3
%         %define 3Dim Gauss kernel
%         sigma = size(img1)/5;
%         kernel = fspecial('gaussian',[size(img1,1) size(img1,2)],max(sigma(1:2)));
%         kernelz = fspecial('gaussian',size(img1,3),sigma(3));
%     	kernelz = squeeze(kernelz(:,floor(size(kernelz,1)/2)));
%         kernel = repmat(kernel,[1 1 size(img1,3)]);
%         kernelz = permute(repmat(kernelz,[1 size(img1,1) size(img1,2)]),[2 3 1]);
%         kernel = kernel.*kernelz;
%         kernel = dip_image(kernel/sum(kernel(:))); %normalize
%         kernel = permute(kernel,[2 1 3]);%dip_image changes x and y
%     otherwise
%         error('not yet implemented for higher dimesions than 3')
% end
%     
% crossCorrW = crossCorr.*kernel; %weighting of the shift, since it is expected to be somehwere in the center
% 
% [~, estimate] = max(crossCorrW);
% cp = size(img1)/2; %centerpixel
% roiRadius = 0.0385*max(size(img1)); %the value as shown to give a reasonable radius
% roiMask = rr(img1)<roiRadius;
% 
% roiMaskShift = circshift(roiMask,round(estimate-cp));
% roi = crossCorr.*roiMaskShift;
% % here we have potential to accelerate the evaluation by taking only the
% % non-zero pixels into account by reducing the roi.
% 
% %center of gravity
% rampX = xx(img1,'corner');
% rampY = yy(img1,'corner');
% cg(1) = sum(roi.*rampX)/sum(roi);
% cg(2) = sum(roi.*rampY)/sum(roi);
% if dimension > 2
%     rampZ = zz(img1,'corner');
%     cg(3) = sum(roi.*rampZ)/sum(roi);
% end
% 
% resShift = cg-cp; %shift is with respect to the 0 shift in the cross correlation which is located at the center pixel

end