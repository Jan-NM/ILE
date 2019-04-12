function [Orte]= clusterfind(slice , meanbg, slice_no, SPDM_EMCCD_Variance)
%%function [Orte ,Ortef]= clusterfind(slice , meanbg, slice_no)
%
%   clusterfind goes through one slice of a stack and detects signal peaks
%   Signals that have been found will be passed to the gaussian estimator for feature extraction,
%   using a ROI with edge length 2 * ROIRADIUS + 1
%
%   A signal peak is defined to be (Threshold-1) * (noise width) higher than average, where the width of the noise distribution is 
%   assumed to be the square root of the the average intensity meanbg. For increased accuracy twice the width of the noise distribution
%   is subtracted before calculating the center of mass, effectively rounding the ROI.
%
%
%   cutoff,filter,seperator
%  
%   Output:
%       Orte  - Array fo estimations [maxQ,mx,my,dx,dy,sx,sy,Qs,slice_no]
%
%created by: Manfred Kirchgessner < Manfred.Kirchgessner@web.de>, 
%            Frederik Grüll <Frederik.Gruell@kip.uni-heidelberg.de>%
%% 
persistent stencil;

if (isempty(stencil))
    stencil = ones(3,3)./9;
end

global SPDMparamstruct;

if(isempty(SPDMparamstruct.CAMpixelsize))
    SPDMparamstruct.CAMpixelsize = 140;
    disp('No value for camera pixelsize found! Value was now set to 140 nm.');
end
if(isempty(SPDMparamstruct.fastSPDMthreshold))
    SPDMparamstruct.fastSPDMthreshold = 3;
    disp('No value for fastSPDM threshold found! Value was now set to 3.');
end

Pix = SPDMparamstruct.CAMpixelsize;
Orte =[];


ROIRADIUS = 3;

THRESHOLD = (SPDMparamstruct.fastSPDMthreshold-1) * sqrt(meanbg);

cutoff =  2 * sqrt(meanbg);

[sizey sizex] = size(slice); 

firimg = imfilter(slice,stencil);

[ypos,xpos] = find(firimg>THRESHOLD);

if ~isempty(ypos)
    for i=1:size(ypos,1)

        y=ypos(i);
        x=xpos(i);

        if (x > ROIRADIUS && y > ROIRADIUS && x <= sizex - ROIRADIUS && y <= sizey - ROIRADIUS)

           roiyrange = y - ROIRADIUS : y + ROIRADIUS;
           roixrange = x - ROIRADIUS : x + ROIRADIUS;
           roimidrange = ROIRADIUS : ROIRADIUS + 2;
           
           froi = firimg(roiyrange,roixrange); 
           fmid_roi = froi(roimidrange,roimidrange);
           
           if( max(fmid_roi(:))==firimg(y,x) )
           
               roi = slice(roiyrange,roixrange); 
                             
               Qges_old = sum(roi(:));


               %[maxQlx,mxlx,mylx,dxl,dyl,sxlx,sylx,Qlx]= linfit( roi );
               %[maxQl,mxl,myl,dxlx,dylx,sxl,syl,Ql]= levenbergh( roi );
                %[maxQf,mx,my,dx,dy,sx2,sy2,Qf]= levenbergh( roi );
               [maxQf,mxf,myf,dx,dy,sx2,sy2,Qf]= estimator( roi,meanbg, SPDM_EMCCD_Variance );

               roi = cl_seperator(roi);                   

               Q = sum(roi(:));               
               maxQ=max(roi(:));

               roi = max(roi - cutoff,0);

               [maxQf,mx,my,dxf,dyf,sx2f,syf2,Qf]= estimator( roi,meanbg, SPDM_EMCCD_Variance );

               if( Qf>0 && sx2 > 0.005 && sy2 > 0.005 && Q/Qges_old > 0.70)% && dy < 0.7 && dx < 0.7 && Qges_mid/maxQ > 3 )
                    mx = (mx+x-1-ROIRADIUS) * Pix;
                    my = (my+y-1-ROIRADIUS) * Pix;
                    dx = dx * Pix;
                    dy = dy * Pix;
                    sx = sqrt(sx2) * Pix;
                    sy = sqrt(sy2) * Pix;

                    estim = [maxQ,my,mx,dx,dy,sx,sy,Q,slice_no,Q/Qges_old];
                    %disp(sprintf('maxQ %d,my %d,mx %d,dx %d,dy %d,sx %d,sy %d,Q %d,slice_no %d ,Q/Qges_old %d',maxQ,my,mx,dx,dy,sx,sy,Q,slice_no,Q/Qges_old))

                    Orte = [Orte;estim];
               end 
%               if( Ql>0 && sxl > 0.005 && syl > 0.005 )% && dy < 0.7 && dx < 0.7 )%&& Qges_mid/maxQ > 3 )
%                    mxl = (mxl+x-1-ROIRADIUS) * Pix;
%                    myl = (myl+y-1-ROIRADIUS) * Pix;
%                    dxl = dxl * Pix;
%                    dyl = dyl * Pix;
%                    sxl = sqrt(sxl) * Pix;
%                    syl = sqrt(syl) * Pix;
%
%                    fit = [maxQl,myl,mxl,dxl,dyl,sxl,syl,Ql,slice_no,Q/Qges_old];
%
%                    Ortef = [Ortef;fit];
%                end  
               
               
           end
          
        end
    end
end

end