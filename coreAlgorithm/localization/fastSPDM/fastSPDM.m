function fastSPDM(varargin)
%%fastSPDM:
%   Detects signal peaks in the image stack by comparing them to the average pixel
%   value over time. The background is subtracted and the signal is fitted with a 
%   Gaussian esimator.
%
%
%params used from SPDM:
%   global d;             % input:  matrix (3D) with raw microscopy data
%   global startImg;      % input (int):  start evaluation with frame #
%   global Pixelsize;     % input: pixel pitch in the object space (in nm)
%   global Threshold;     % input: see function clusterfind
% 
%   global Orte;          % output: matrix with all localized signals
%
%%created by: Manfred Kirchgessner < Manfred.Kirchgessner@web.de>, 
%            Frederik Gruell <Frederik.Gruell@kip.uni-heidelberg.de>
%%
%profile on;  

    % fprintf('Number of arguments is %d\n', nargin);

    if (nargin == 0)
        filename = 'unknown file';
    else
        filename = varargin{1};
        c = varargin{1,2};
        SPDM_EMCCD_Variance = varargin{3};
    end

    global Orte;
    
    global d;
    global SPDMparamstruct;
        
    if(isempty(SPDMparamstruct.STACKstartFrame)) || (SPDMparamstruct.STACKstartFrame < 1)
        SPDMparamstruct.STACKstartFrame = 1;
    end
        
    if(isempty(d))
        error('what the hell the file d is gone');
    end
    
    mybar = NaN;
        
    wbstr = 'Localization of signals running ...     Please wait ...';
    if(isempty(SPDMparamstruct))
        SPDMparamstruct.outputmode = 'console';
    end
    waitbar = evalin('base','@waitbar');
    switch SPDMparamstruct.outputmode
        case 'console'
            waitbar = @outputDisp;
            fprintf('File: %s - in process.\n',filename);
            fprintf(wbstr);
        case 'silent'
            waitbar = @outputNone;
        otherwise
            mybar = waitbar(0,'Localization of signals running... Please wait.', 'Name', [filename ' - progress']);
    end
    
    % factor for background calculation
    BGF = 8; 
    % old background image
    if mod(SPDMparamstruct.STACKstartFrame, c.simultaneousFrames) == 0
        bgimage = double(d(:,:,1));
    else        
        bgimage = double(d(:,:,mod(SPDMparamstruct.STACKstartFrame, c.simultaneousFrames)));
    end
    % create a background profile over the first 8 images
    try
        bgimage = double(mean(d(:,:,SPDMparamstruct.STACKstartFrame+(0:(BGF-1))),3));
    end
    meanbg = round(mean(bgimage(:)));
    
    tic  
%% Image Handling
    if mod(SPDMparamstruct.STACKstartFrame, c.simultaneousFrames) == 0
        iTmp = 1;
    else
        iTmp = mod(SPDMparamstruct.STACKstartFrame, c.simultaneousFrames);
    end
    
    if SPDMparamstruct.STACKendFrame > c.siz(3)
        staple_size = c.siz(3);
    else
        staple_size = SPDMparamstruct.STACKendFrame;
    end
    % pre-allocate Orte matrix
    nFrames = staple_size - SPDMparamstruct.STACKstartFrame;
    Orte = cell(1, nFrames);
    
    for i = SPDMparamstruct.STACKstartFrame : staple_size
         if (iTmp > c.simultaneousFrames)
            iTmp = 1;
            try
                c.spdm_hFlushFrames(c,i,c.simultaneousFrames);
                % convert data into photons
                d = (d-SPDMparamstruct.CAMoffset).*SPDMparamstruct.CAMconvfactor./SPDMparamstruct.EMgain./SPDMparamstruct.qe;
            catch
                set(SPDMparamstruct.handles.startbut,'Enable','on');
                c.spdm_hCleanUp()
                rethrow(lasterror);
            end
        end
        [diff,meanbg,bgimage] = addimagetobackground( bgimage,double(d(:,:,iTmp)) ,BGF, meanbg,i);  
        % create and subtract backgroundimage from image and calculate mean background
        % mesh(diff);figure(h);
        Orte{i} = clusterfind(diff,meanbg,i, SPDM_EMCCD_Variance);
        % find Signal peaks and estimate meanvalues
        % Orte = [Orte;Orte_slice];
                

        if (mod(i,staple_size/100) == 0)
           % surf(bgimage);figure(h);
            try      
                wbstr = sprintf('Frame # %8d of %8d.', i, staple_size);
                waitbar((i/staple_size), mybar, wbstr);
            catch
                set(SPDMparamstruct.handles.startbut,'Enable','on');
                c.spdm_hCleanUp(c.filehandle)
                rethrow(lasterror);
            end
        end
        iTmp = iTmp + 1;
    end
    Orte = (cell2mat(Orte.'));
    c.spdm_hCleanUp(c.filehandle)
    % profile viewer
%% DIFF IMAGE
%     tic
%     for i=startImg+1:staple_size
% 
%         meanbg = mean(mean(double(d(:,:,i))));
%         diff = max(double(d(:,:,i)-d(:,:,i-1)),0);
% 
%         Orte_slice = clusterfind(diff,meanbg,i);
%         %find Signal peaks and estimate meanvalues
%         Orte = [Orte;Orte_slice]; 
%         
%         try 
%             waitbar((i/staple_size),mybar,['Calculating image nr. ' num2str(i) ' of ' num2str(staple_size) '. \n ' num2str(size(Orte,1)) ' Points found.']);
%         catch e
%             rethrow(e);
%         end       
%        
%     end
    
    
    fprintf('\n\n');
    fprintf('Done evaluating %s\n',filename); 
    fprintf('Execution time: %d\n\n', toc');        % time measurement
    fprintf('Signals found: %g\n',size(Orte,1));
    if(isempty(Orte))
        mean_loc = NaN;
    else
        mean_loc = mean(mean(Orte(:,4:5),2));
    end
    fprintf('Mean localization precision:  %g\n\n', mean_loc);

    close(mybar);
end

function outputNone(fraction,mybar,wbstr)
    ;
end
    
function outputDisp(fraction,mybar,wbstr)
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b')
    fprintf(wbstr)
end