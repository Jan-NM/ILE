function SPDMprocessing()
% SPDMprocessing.m
%   evaluation and visualization script
%   to be called after reading in the data set
%   reads in data stack via SPDMdatahandle
%
%   Input:
%   SPDMparamstruct   global variable from startSPDM
%   files             global variable from startSPDM
%
%   Output:
%   Orte              .mat-file containing the (corrected) list of localizations
%
%   Internal global variables
%   d
%   file_dir
% 
% Cremer Group, Institute of Molecular Biology (IMB), Mainz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% input parameters:
global files
global SPDMparamstruct

% parameters used internally:
global file_dir;
global d

% return parameters:
global Orte
result = {};
                 
%--------------------------------------------------------------------------
%% setup modified filelist

[~, FILESmax] = size(files);

if FILESmax < 10
    jm = ['0' num2str(FILESmax)];
else
    jm = num2str(FILESmax);
end

pn = [SPDMparamstruct.DIR_input,filesep];
sampleID = regexpi(SPDMparamstruct.DIR_input, filesep, 'split');
sampleID=sampleID{end};

%--------------------------------------------------------------------------
%% start looping through files
for j = 1:FILESmax
    if j < 10
        js = ['0' num2str(j)];    
    else
        js = num2str(j);
    end
    SPDMparamstruct.displayStatus(['running Evaluation of File No. ' js ' of ' jm '...']);
    wd = pwd;
    if (~isempty(file_dir))
        fd = 1;
        cd(file_dir);
    end
    fn = char(files(j));
    if isequal(fn,0)||isequal(pn,0)
        SPDMparamstruct.displayStatus('Error - Evaluation not successful');
        SPDMparamstruct.displayError('No File Selected');
        SPDMparamstruct.enableStartButton();
        error('No File Selected');
    end
    %% start accessing file
    file_dir = pn;
    cd(wd);
    switch SPDMparamstruct.orteSel
        case 'image'
            % correct noise model for EMCCD cameras:
            if(SPDMparamstruct.CAMtypeEMCCD)
                SPDM_EMCCD_Variance = sqrt(2);
            else
                SPDM_EMCCD_Variance = 1;
                SPDMparamstruct.EMgain  = 1; % if EMCCD camera is not used, EM gain is set to 1
            end
            % initialize file handle
            c = spdmDataHandle([pn fn]);
            if c.isTIF
                nFiletype = 4;
            elseif c.isH5
                nFiletype = 3;
            else
                disp('Unsupported filetype!');
            end
            fileID = [sampleID '_#' fn(1:end-nFiletype)];
            try
                c = c.getImageDimensions(c);
                if (c.siz(3) < SPDMparamstruct.STACKstartFrame )
                    error('Parameter: "start frame" larger than frame number of image stack. Try a smaller start frame!');   
                end
                DimX = c.siz(1);
                DimY = c.siz(2);
                DimZ = c.siz(3);
                DimE = c.siz(4);
                d = zeros(DimY,DimX,c.simultaneousFrames,'single');
                % reading first image batch
                c.spdm_hFlushFrames(c,SPDMparamstruct.STACKstartFrame,c.simultaneousFrames);
            catch
                set(SPDMparamstruct.handles.startbut,'Enable','on');
                c.spdm_hCleanUp(c.filehandle)
                rethrow(lasterror)
            end
            % convert data into photons
            d = (d-SPDMparamstruct.CAMoffset).*SPDMparamstruct.CAMconvfactor./SPDMparamstruct.EMgain./SPDMparamstruct.qe;
            %% start localization process
            % ##############################################
            % ###        OBJECT FINDING ROUTINE         ####
            % ##############################################
            try
               fastSPDM(fullfile(pn,fn), c, SPDM_EMCCD_Variance)
            catch 
                SPDMparamstruct.enableStartButton();
                warning('User ABORT! or:')
                rethrow(lasterror)
            end
        case 'Orte'
            Orte = load([pn, fn]);
            nameOfField = fieldnames(Orte);
            Orte = Orte.(nameOfField{1});
            fileID = [sampleID '_#' fn(1:end-4)];         
    end
    % initialize correction flag's for result file
    flag_mfp = 'mfp_N';
    flag_filt_PSF = 'filterPSF_N';
    flag_filt_Loc = 'filterLoc_N';
    flag_filt_Phot = 'filterPhot_N';
    flag_drift_corr = 'driftCorr_N';
    %% start correction routine: remove multiframe localizations
    % #################################################################
    % ###         remove multiframe objects from Orte-matrix       ####
    % #################################################################
    meanloc =  mean(mean(Orte(:,4:5),2));
    if SPDMparamstruct.flag_multiframe == 1
        Orte = multiframepoints(Orte, SPDMparamstruct.searchRadius.*meanloc);
        disp('After Multiframe Correction');
        meanloc = mean(mean(Orte(:,4:5),2));
        fprintf('Signals found: %g\n',size(Orte,1));
        fprintf('Mean localization precision:  %g\n\n', meanloc);
        flag_mfp = ['mfp_Y_' num2str(SPDMparamstruct.searchRadius)];
        % fileID = [fileID '_mfp' num2str(SPDMparamstruct.searchRadius)];
    end
    %% start correction routine: filter localizations
    % #################################################################
    % ###         remove overlapping signals                       ####
    % #################################################################
    if SPDMparamstruct.flag_filterSignals == 1
        indexPSFx = Orte(:,6) <= SPDMparamstruct.upperLimitPSF;
        Orte = Orte(indexPSFx, :);
        indexPSFy = Orte(:,7) <= SPDMparamstruct.upperLimitPSF;
        Orte = Orte(indexPSFy, :);
        indexLocx = Orte(:,4) <= SPDMparamstruct.upperLimitLoc;
        Orte = Orte(indexLocx, :);
        indexLocy = Orte(:,5) <= SPDMparamstruct.upperLimitLoc;
        Orte = Orte(indexLocy, :);
        indexPhot = Orte(:,1) >= SPDMparamstruct.limitPhot;
        Orte = Orte(indexPhot, :);
        disp('After Filtering');
        meanloc = mean(mean(Orte(:,4:5),2));
        fprintf('Signals found: %g\n',size(Orte,1));
        fprintf('Mean localization precision:  %g\n\n', meanloc);
        flag_filt_PSF = ['filterPSF_Y_' num2str(SPDMparamstruct.upperLimitPSF)];
        flag_filt_Loc = ['filterLoc_Y_' num2str(SPDMparamstruct.upperLimitLoc)];
        flag_filt_Phot = ['filterPhot_Y_' num2str(SPDMparamstruct.limitPhot)];
        % fileID = [fileID '_filt_' num2str(SPDMparamstruct.upperLimitPSF) '_' num2str(SPDMparamstruct.upperLimitLoc)];
    end
    %% start correction routine: drift correction
    % #################################################################
    % ###         correct Orte matrix for drift                    ####
    % #################################################################
    if SPDMparamstruct.flag_driftCorrection == 1
        Orte = fastCorrectDriftLocNS(Orte, SPDMparamstruct.subsetNo);
        flag_drift_corr = ['driftCorr_Y_SubsetNo: ' num2str(SPDMparamstruct.subsetNo)];
        % fileID = [fileID '_driftCorr_' num2str(SPDMparamstruct.subsetNo)];
    end
    %% save ORTE file (list of localizations)    
    result = [result; {fileID} ,size(Orte,1), meanloc, {flag_mfp},...
    {flag_filt_PSF}, {flag_filt_Loc}, {flag_filt_Phot}, {flag_drift_corr}];
    SPDMparamstruct.OrtePathname = [SPDMparamstruct.DIR_output filesep];
    SPDMparamstruct.OrteFilename = [fileID '_orte.mat'];
    save([SPDMparamstruct.OrtePathname SPDMparamstruct.OrteFilename],'Orte');
    clear('global','d');
    %% create data visualizations / reconstructions
    % #################################################################
    % ###         generating super resolution images               ####
    % #################################################################
    if SPDMparamstruct.OUTPUTscatterplot == 1
       visModule(Orte, 'scatterplot', SPDMparamstruct.OUTPUTpixelsize, SPDMparamstruct, [fileID, '_scatterplot'], SPDMparamstruct.DIR_output);
       % maybe test if images should be shown
    end
    if SPDMparamstruct.OUTPUThistogram == 1
        visModule(Orte, 'histogramBinning', SPDMparamstruct.OUTPUTpixelsize, SPDMparamstruct, [fileID, '_histogramBinning'], SPDMparamstruct.DIR_output);
    end
    if SPDMparamstruct.OUTPUTgaussianBlur == 1
        visModule(Orte, 'gaussianBlur', SPDMparamstruct.OUTPUTpixelsize, SPDMparamstruct, [fileID, '_gaussianBlur'], SPDMparamstruct.DIR_output);
    end
    if SPDMparamstruct.OUTPUTtriangulation == 1
        visModule(Orte, 'triangulation', SPDMparamstruct.OUTPUTpixelsize, SPDMparamstruct, [fileID, '_triangulation'], SPDMparamstruct.DIR_output);        
    end
    if SPDMparamstruct.OUTPUTvoronoi == 1
        visModule(Orte, 'voronoi', SPDMparamstruct.OUTPUTpixelsize, SPDMparamstruct, [fileID, '_voronoi'], SPDMparamstruct.DIR_output);        
    end
    if SPDMparamstruct.OUTPUTdensity == 1
        visModule(Orte, 'densityMap', SPDMparamstruct.OUTPUTpixelsize, SPDMparamstruct, [fileID, '_densityMap'], SPDMparamstruct.DIR_output);        
    end
end
save([SPDMparamstruct.DIR_output filesep 'LocResults_' sampleID '.mat'],'result');
end