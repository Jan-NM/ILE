function [ diff,meanbg,newbgimage ] = addimagetobackground( bgimage, image , bgfactor, oldmeanbg, imgnum )
%function [ diff,meanbg,newbgimage ] = addimagetobackground( bgimage, image )
%   subtract background and add new image to background
%   additionally calculates the mean value of image
%
% to avoid filling the background image with signal peaks a mask is created which
% eliminates big noise signals in the image
%
% Output:
%   diff        -background eliminated image which keeps (in best case) only
%                signal peaks
%   meanbg      -meanvalue over the whole slice to calculate noise strength
%   newbgimage  -returns the new background image for the next slice
%
%
%created by: Manfred Kirchgessner < Manfred.Kirchgessner@web.de>, 
%            Frederik Grüll <Frederik.Gruell@kip.uni-heidelberg.de>%
%%
BGF=bgfactor;

meanbg = round(mean(mean(image)));

cutoff = sqrt(meanbg);

%normally the background lowers over the pictures, if it raises from one
%pic to another over 4 times cutoff there must be sth wrong in the
%pictures, and the backround is adapted to the new level. 
if( (meanbg - oldmeanbg) < ( 4 * cutoff) )

    diff   = image - bgimage;
      
    newbgimage = bgimage + min( diff, cutoff )./BGF;   
    
    diff = max(diff,0);
    
else
    disp(['Warning at img #' num2str(imgnum) ' : big change in background detected. Try to use different start Image: ' num2str(imgnum+1)])
    newbgimage = image;
    diff = 0. * image;  

end

%mesh(image-imagemask);figure(gcf);

%% same in cstyle 
% [sizey sizex] = size(image);
% 
% diff = zeros(sizey,sizex);
% meanbg = 0;
% newbgimage = zeros(sizey,sizex);
% 
% for y = 1:sizey
%     for x = 1: sizex
%         diff(y,x) = max(image(y,x)-bgimage(y,x) , 0);
%         newbgimage(y,x) = (bgimage(y,x) * 15.0 + image(y,x) * 1.0)/16;
%         meanbg = meanbg + image(y,x);
%     end
% end
% 
% meanbg = round(meanbg / (sizex*sizey));

end

