function [varargout] = visModule(locTable, visMethod, pixelSize, paramStruct, filename, pathname)
%visModule contains various methods to map localization table into SMLM image
%
% Following visualization methods are supported
% scatterplot           [2D only]
% histogramBinning      [2D only]
% gaussianBlur          [2D only]
% triangulation         [2D only]
% voronoi               [2D only]
% densityMap            [2D only]
%
% input parameters
%   locTable        Matrix that contains point coordinates. Rows should
%                   correspond to detected signals. Columns should be
%                   arranged in the follwing order (column 2 - x coordiante, column 3 - y coordiante, , column 4 - x loc. prec., column 5 - y loc. prec., column 9 - frame number)
%   visMethod       any of the specified visualization methods as string
%                   for example 'histogramBinning'
%   pixelSize       final super-resolution pixel size
%   paramStruct     Struct that contains further parameters:
%       paramStructName.dim             2 or 3 for 2D or 3D visualization
%                                       default: 2
%       paramStructName.zStep           z-Step size for 3D reconstruction
%                                       default: 0
%       paramStructName.saveImage       boolean value if SR images should be
%                                       saved or returned as matrix by the function
%                                       default: 0
%       paramStructName.blurringMethod  specifies the blurring method for
%                                       gaussian blur, possible values:
%                                       'individualBlur' - blurs every
%                                       point with its respective loc.
%                                       prec.
%                                       'globalBlur' - blurs every point
%                                       with the mean loc. prec. of the
%                                       localization table
%                                       'specificBlur' - blurs every point
%                                       with user defined sigma
%                                       default: 'individualBlur'
%       paramStructName.blurringSize    size of sigma (in nm) for gaussianBlur if
%                                       'specificBlur' is specified
%                                       default: 200
%       paramStructName.nPertubation    Number of runs for trinagulation
%                                       and voronoi based visualization to smooth triangles
%                                       default: 0
%       paramStrcutName.acc             Entire localization table is devided into acc parts to accelerate
%                                       the evaluation. Distances are computed on these subsets only.
%                                       If this value is set too high and densities are low,
%                                       horizontal artefacts may occure. acc == 1 always generates
%                                       the correct image, however higher acc values do not necessarily
%                                       generate false images
%                                       default: 10
%       paramStructName.numNN           Number of next neighbours that are taken into account to
%                                       calculate the local density.
%                                       default: 5
%       paramStrcutName.calib           Boolean factor to calibrate the local density to number of
%                                       events per square nm. Result containes more noise. 
%                                       default: false
%   filename    name of image to be saved. if not specified filename is
%               'SRimage' + name of selected visualization method
%   pathname    location where the image will be saved. if not specified
%               image will be saved in current folder
%   
% requires Image Processing Toolbox, Statistics and Machine Learning
% Toolbox
% requires multiWaitbar from Matlab File Exchange
%
%   by Jan Neumann, IMB Mainz, 01.11.2017
%       'triangulation' - adapted from Martin Hagmann and Gerrit Best
%       'voronoi' - adapted from Martin Hagmann, Gerrit Best and Kirti Prakash
%       'densityMap' - adapted from Martin Hagmann, Gerrit Best and Kirti Prakash
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dimError = 0; % for error handling
nSteps = 100; % needed in gaussian blur if individual blurring method is used, specifies binning for loc. prec.

if nargin < 3
    multiWaitbar('CloseAll');
    error('Please specify at least visualization method and pixel size!');
end

if nargin < 4
    dim = 2;
    zStep = 0;
    saveImage = 0; % do not change otherwise recursion for triangulation and voronoi may fail
    % set default blurring type for gaussian blurring
    blurringMethod = 'individualBlur'; % options are 'globalBlur' or 'specificBlur' or 'individualBlur'
    blurringSize = 200; % in nm, only valid for 'specificBlur'
    nPertubation = 0; % do not change this parameter otherwise recursion for triangulation and voronoi may fail
    acc = 10;
    calib = 0;
    numNN = 5;
    if saveImage == 1 && ~exist('filename', 'var') && ~exist('pathname', 'var')
        pathname = pwd;
        filename = ['SRimage_', visMethod];
        disp(['Image will be saved in ', pathname, ' as ', filename, '.']);
    end
else
    try
        dim = paramStruct.dim;
        if dim == 3
            try
                zStep = paramStruct.zStep;
            catch
                dimError = 1;
                msg = ('3D reconstruction requested, but no zStep size specified!');
                error(msg);
            end
        end
        saveImage = paramStruct.saveImage;
        if saveImage == 1 && nargin < 5
            pathname = pwd;
            filename = ['SRimage_', visMethod];
            disp(['Image will be saved in ', pathname, ' as ', filename, '.']);
        end
        switch visMethod
            case 'gaussianBlur'
                % set default blurring type for gaussian blurring
                blurringMethod = paramStruct.blurringMethod; % options are 'globalBlur', 'specificBlur' or 'individualBlur'
                if strcmp(blurringMethod, 'specificBlur')
                    blurringSize = paramStruct.blurringSize; % in nm, only valid for 'specificBlur'
                end
            case {'triangulation', 'voronoi'}
                % parameters for triangulation
                nPertubation = paramStruct.nPertubation;
            case 'densityMap'
                acc = paramStruct.acc;
                calib = paramStruct.calib;
                numNN = paramStruct.numNN;
        end
    catch
        if dimError == 1
            multiWaitbar('CloseAll');
            error(msg)
        else
           multiWaitbar('CloseAll');
           error('You must specifiy a list of parameters containing at minimum the following fields: .dim, .saveImage and the visualization specific parameters. Alternatively do not specify any parameters and run everything on default!');
        end
    end
end
%% start with visualization method
switch visMethod
    case 'scatterplot'
        multiWaitbar('generating scatter plot...', 'Busy');
        SRimage = figure('visible','on'); % if visible is off and figure can not show up when clicking on saved figure
        scatter(locTable(:, 3)./pixelSize, -locTable(:, 2)./pixelSize, 'b.');
        axis([min(locTable(:, 3))./pixelSize max(locTable(:, 3))./pixelSize...
            -max(locTable(:, 2))./pixelSize -min(locTable(:, 2))./pixelSize])
        axis off;
    case 'histogramBinning'
        multiWaitbar('generating histogram binned image...', 'Busy');
        x1 = round(locTable(:, 2)./pixelSize);
        y1 = round(locTable(:, 3)./pixelSize);
        % generate histogram binned image
        SRimage = sparse(x1-min(x1)+1,y1-min(y1)+1,1);
        SRimage = full(SRimage);
    case 'gaussianBlur'
        switch blurringMethod
            case 'globalBlur'
                multiWaitbar('generating gaussian blurred image...', 'Busy');
                % each point is blurred with the average loc. prec.
                globalBlur = mean(mean(locTable(:, 4:5)))/pixelSize;
                x1 = round(locTable(:,2)./pixelSize);
                y1 = round(locTable(:,3)./pixelSize);
                % generate histogram binned image
                SRimage = sparse(x1-min(x1)+1,y1-min(y1)+1,1);
                SRimage = full(SRimage);
                % check version of matlab
                if verLessThan('matlab', '8.5')
                    h = fspecial('gaussian', 2*ceil(2*globalBlur) + 1, globalBlur);
                    SRimage = imfilter(SRimage, h, 'conv');
                    % alternative implement own filter
                    % maskSize = 2*ceil(2*globalBlur) + 1;
                    % ind = -floor(maskSize/2) : floor(maskSize/2);
                    % [X, Y] = meshgrid(ind, ind);
                    % h = exp(-(X.^2 + Y.^2)/(2*globalBlur*globalBlur));
                    % h = h / sum(h(:));
                    % SRimage = conv2(SRimage ,h);
                else
                    SRimage = imgaussfilt(SRimage, globalBlur);
                end
            case 'specificBlur'
                multiWaitbar('generating gaussian blurred image...', 'Busy');
                % each point is blurred with a user defined loc. prec.
                globalBlur = blurringSize/pixelSize;
                x1 = round(locTable(:,2)./pixelSize);
                y1 = round(locTable(:,3)./pixelSize);
                % generate histogram binned image
                SRimage = sparse(x1-min(x1)+1,y1-min(y1)+1,1);
                SRimage = full(SRimage);
                % check version of matlab
                if verLessThan('matlab', '8.5')
                    h = fspecial('gaussian', 2*ceil(2*globalBlur) + 1, globalBlur);
                    SRimage = imfilter(SRimage, h, 'conv');
                    % alternative implement own filter
                    % maskSize = 2*ceil(2*globalBlur) + 1;
                    % ind = -floor(maskSize/2) : floor(maskSize/2);
                    % [X, Y] = meshgrid(ind, ind);
                    % h = exp(-(X.^2 + Y.^2)/(2*globalBlur*globalBlur));
                    % h = h / sum(h(:));
                    % SRimage = conv2(SRimage ,h);
                else
                    SRimage = imgaussfilt(SRimage, globalBlur);
                end
            case 'individualBlur'
                multiWaitbar('generating gaussian blurred image...', 0);
                minPrecision = floor(min(min(locTable(:, 4:5))));
                maxPrecision = ceil(max(max(locTable(:, 4:5))));
                % get dimension of image
                x1 = round(locTable(:,2)./pixelSize);
                y1 = round(locTable(:,3)./pixelSize);
                % generate histogram binned image
                [xSize, ySize] = size(full(sparse(x1-min(x1)+1,y1-min(y1)+1,1)));
                minX = min(x1);
                minY = min(y1);
                SRimage = zeros(xSize, ySize, 'single');
                prevStep = minPrecision;
                steps = linspace(minPrecision, maxPrecision, nSteps);
                steps(1) = [];
                % for waitbar
                prevPercent = 0;
                counter = 1;
                for ii = steps
                    % find points with specific loc prec.
                    idx = (mean(locTable(:, 4:5), 2) <= ii) & (mean(locTable(:, 4:5), 2) > prevStep);
                    tempPoints = locTable(idx, :);
                    globalBlur = ((0.5*(ii - prevStep)) + prevStep) / pixelSize;
                    % generate histogram binned with current selection of points
                    x1 = round(tempPoints(:,2)./pixelSize);
                    y1 = round(tempPoints(:,3)./pixelSize);
                    tempImage = sparse(x1-minX+1, y1-minY+1, 1, xSize, ySize); % check if x and y are correct
                    tempImage = full(tempImage);
                    % images should have fixed size to be abble to be added
                    % check version of matlab
                    if verLessThan('matlab', '8.5')
                        h = fspecial('gaussian', 2*ceil(2*globalBlur) + 1, globalBlur);
                        SRimage = SRimage + imfilter(tempImage, h, 'conv');
                        % alternative implement own filter
                        % maskSize = 2*ceil(2*globalBlur) + 1;
                        % ind = -floor(maskSize/2) : floor(maskSize/2);
                        % [X, Y] = meshgrid(ind, ind);
                        % h = exp(-(X.^2 + Y.^2)/(2*globalBlur*globalBlur));
                        % h = h / sum(h(:));
                        % SRimage = SRimage + conv2(tempImage ,h);
                    else
                        SRimage = SRimage + imgaussfilt(tempImage, globalBlur);
                    end
                    prevStep = ii;
                    currentPercent = fix(100*counter/size(steps, 2));
                    if currentPercent > prevPercent
                        multiWaitbar( 'generating gaussian blurred image...', 'Value', counter/size(steps, 2));
                        prevPercent = currentPercent;
                    end
                    counter = counter + 1; 
                end
        end
    case  'triangulation'% adapted from Orte2TriBild.m
        if nPertubation > 0
             multiWaitbar('generating triangulation image...', 0);
        end 
        if length(nPertubation) > 1 % initialize the slice index if necessary
            sliceInd = 1;
            if ~issorted(nPertubation)
                multiWaitbar('CloseAll');
                error('You specified a vector as 3rd argument that was not sorted. Sort and try again.')
            end
        else
            sliceInd = 0;
        end
        % chop off negative values (might be introduced by jittering)
        locTable(locTable(:, 2) < 0 | locTable(:, 3) < 0,:) = [];
        % get the locations in pixels, not nm (x and y is swapped in the Orte files)
        pos = [locTable(:, 2) locTable(:, 3)] / pixelSize;
        imSize = [ ceil( max(pos(:, 1)) )  ceil( max(pos(:, 2)) )] + 1;
        % generate triangles
        if verLessThan('matlab', '8.1') % check earlier than matlab 2013a
            triangles = DelaunayTri(pos);
            SRimage = zeros(imSize, 'single');
            totnumtriang = length(triangles.Triangulation);
        else
            triangles = delaunayTriangulation(pos);
            SRimage = zeros(imSize, 'single');
            totnumtriang = length(triangles.ConnectivityList);
        end
        trianglenum = 1:totnumtriang;
        % this value should be much higher than the actual number of threads since some runs may take order of magnitudes longer than others depending on ROI size.
        numthreads = 100;
        threadruns = ceil(totnumtriang/numthreads);
        % allocate cells to hold the rois
        rmask = cell(1, numthreads);
        ROIstartx = cell(1, numthreads);
        ROIstarty = cell(1, numthreads);
        ROIendx = cell(1, numthreads);
        ROIendy = cell(1, numthreads);
        multiWaitbar( 'processing triangles...', 0);
        % generate the pixelated cells
        for actrun = 1:threadruns
            starttriang = (actrun-1)*numthreads + 1;
            endtriang = min(actrun*numthreads, totnumtriang);
            acttriangs = trianglenum(starttriang:endtriang);
            actnumthreads = numel(acttriangs);
            % for progress output
            [checkifmod1000, outputPos] = min(mod(acttriangs,1000));
            if checkifmod1000 == 0
                multiWaitbar( 'processing triangles...', 'Value', acttriangs(outputPos) / totnumtriang);
            end
            for threadrun = 1:actnumthreads
                ii = acttriangs(threadrun);
                % extract current triangle indices
                if verLessThan('matlab', '8.1') % check earlier than matlab 2013a
                    triInd = triangles.Triangulation(ii, 1:3);
                else
                    triInd = triangles.ConnectivityList(ii, 1:3);
                end
                % get absolute positions of the cutten triangle in pixels
                apos = pos(triInd, :);
                apos = round(apos);
                ROIstartx{threadrun} = min(apos(:, 1));
                ROIstarty{threadrun} = min(apos(:, 2));
                ROIendx{threadrun} = max(apos(:, 1)) + 1;
                ROIendy{threadrun} = max(apos(:, 2)) + 1;
                ROIsize_x = ROIendx{threadrun}-ROIstartx{threadrun};
                ROIsize_y = ROIendy{threadrun}-ROIstarty{threadrun};
                % current triangles in the respective ROIs
                rpos = zeros(size(apos, 1), 2);% allocate
                rpos(:, 1) = apos(:, 1)-ROIstartx{threadrun} + 1;
                rpos(:, 2) = apos(:, 2)-ROIstarty{threadrun} + 1;
                % map the roi to the image
                rmask{threadrun} = poly2mask(rpos(:, 2),rpos(:, 1), ROIsize_x, ROIsize_y);
                % Area of current triangle
                A = sum(rmask{threadrun}(:));
                if A == 0 % distance was too small to draw a polygon. Print pixel at the center instead
                    x = round(mean(rpos(:,1)));
                    y = round(mean(rpos(:,2)));
                    rmask{threadrun}(x,y) = 1;
                    A = 1;
                end
                % brightness of triangle is 1/Area
                rmask{threadrun} = rmask{threadrun}/A;
            end
            % add rois from threads to bild array
            for ii = 1:length(acttriangs)
                SRimage(ROIstartx{ii}+1:ROIendx{ii},...
                    ROIstarty{ii}+1:ROIendy{ii})=SRimage(ROIstartx{ii}+1:ROIendx{ii},...
                    ROIstarty{ii}+1:ROIendy{ii})+rmask{ii};
            end
        end
        % call the function recursively
        if any(nPertubation) > 0
            if any(nPertubation < 3)
                warning('You specified pertubation = %i. Such low amount of perturbation should be used with caution because of the low and thus bad statistics.\n', nPertubation(1))
            end
            locTableP = zeros([size(locTable, 1) 3 nPertubation(end)]);
            for ii = 1:nPertubation(end)
                locTableP(:, 2, ii) = normrnd(locTable(:,2), locTable(:, 4));
                locTableP(:, 3, ii) = normrnd(locTable(:,3), locTable(:, 5));
                % Fluorophroes at the edge might slip beyond and become negative which
                % is why we omit them in the triangulation (where the locTable is 2
                % dimensional - at the beginning of this function)
                % the size may change on the different sets, which is why we calculate
                % it again
                pos = [locTableP(:, 2, ii) locTableP(:, 3, ii)] / pixelSize; % get the locations in pixels, not nm (x and y is swapped in the locTable files)
                imSizeTemp = [ ceil( max(pos(:,1)) )  ceil( max(pos(:,2)) )] + 1;
                imSize = max(imSize, imSizeTemp); % and we take whichever values are bigger
            end
            % expanding already generated SRimage to match the new size
            SRimage = cat(2,SRimage,zeros(size(SRimage,1),imSize(2)-size(SRimage,2)));
            SRimage = cat(1,SRimage,zeros(imSize(1)-size(SRimage,1),imSize(2)));
            % extend imSize and the image by a third dimension if nPertubation is a vector
            if sliceInd > 0
                imSize = [imSize length(nPertubation)];
                SRimage = cat(3,SRimage,zeros(imSize(1),imSize(2),imSize(3)-1));
            end
            for ii = 1 : nPertubation(end)
                multiWaitbar( 'generating triangulation image...', 'Value', (ii-1) / nPertubation(end));
                % now of course without perturbation, we wouldn't finish ever otherwise
                % with that we can call Orte2TriBild pert times with the new Orte sets
                bildNew = visModule(locTableP(:,:,ii), 'triangulation',  pixelSize);
                % (1,1), so top left should be absolute in every image and we can
                % concatenate zeros at the bottom right (matlab convention)
                bildNew = cat(2,bildNew,zeros(size(bildNew,1),imSize(2)-size(bildNew,2)));
                bildNew = cat(1,bildNew,zeros(imSize(1)-size(bildNew,1),imSize(2)));
                
                % summing up the result of the recursions
                if sliceInd > 0 %pert is a vector
                    % due to matlab's stupid indexing we have to subtract 1 for the indexing, aarr!
                    if any(ii-1 == nPertubation(1:end-1)) %except the last, because then we're done already
                        SRimage(:,:,sliceInd+1) = SRimage(:,:,sliceInd);%copy last generated image to new slice
                        sliceInd = sliceInd+1;%increment slice index
                    end
                    SRimage(:,:,sliceInd) = SRimage(:,:,sliceInd) + bildNew;
                else % pert is a scalar
                    % no need to care about storing the images and we just add the
                    % results
                    SRimage = SRimage + bildNew;
                end
            end
        end
    case 'voronoi' % adapted from Orte2CellBild
        if nPertubation > 0
             multiWaitbar('generating voronoi image...', 0);
        end 
        if length(nPertubation) > 1
            sliceInd = 1;
            if ~issorted(nPertubation)
                multiWaitbar('CloseAll');
                error('You specified a vector as 3rd argument that was not sorted. Sort and try again.')
            end
        else
            sliceInd = 0;
        end % initialize the slice index if necessary
        % get the locations in pixels, not nm (x and y is swapped in the Orte files)
        pos = [locTable(:, 2) locTable(:, 3)] / pixelSize;
        imSize = [ ceil( max(pos(:,1)) )  ceil( max(pos(:,2)) )] + 1;
        pos = unique(pos, 'rows'); %remove duplicates
        [vertices, regions] = voronoin([pos(:,1) pos(:,2)]);
        % run the cell generation
        SRimage = zeros( imSize , 'single');
        totnumcells = length(regions);
        cellnum = 1:totnumcells;
        % this value should be much higher than the actual number of threads since some runs may take order of magnitudes loger than others depending on ROI size.
        numthreads = 100;
        threadruns = ceil(totnumcells/numthreads);
        % allocate cells to hold the rois
        rmask = cell(1,numthreads);
        ROIstartx = cell(1,numthreads);
        ROIstarty = cell(1,numthreads);
        ROIendx = cell(1,numthreads);
        ROIendy = cell(1,numthreads);
        multiWaitbar( 'processing voronoi cells...', 0);
        % generate the pixelated cells
        for actrun = 1:threadruns
            starttriang = (actrun-1)*numthreads + 1;
            endtriang = min(actrun*numthreads, totnumcells) ;
            acttriangs = cellnum(starttriang:endtriang);
            actnumthreads = numel(acttriangs);
            % for progress output
            [checkifmod1000, outputPos] = min(mod(acttriangs, 1000));
            if checkifmod1000 == 0
                multiWaitbar( 'processing voronoi cells...', 'Value', acttriangs(outputPos) / totnumcells);
            end
            for threadrun = 1:actnumthreads
                ii=acttriangs(threadrun);
                % extract vertices of the current cell ii
                apos = vertices(regions{ii},:);
                if length(apos(:, 1)) <= 2 % if there are less than 2 vertices present for the current cell
                    continue; % with the next loop iteration, since we cannot define an area
                end
                % check if one point is outside image
                below_x = apos(:, 1) < 1;
                below_y = apos(:, 2) < 1;
                beyond_x = apos(:, 1) > imSize(1);
                beyond_y = apos(:, 2) > imSize(2);
                remtup = below_x | below_y | beyond_x | beyond_y;
                if any(remtup(:)) > 0
                    continue; % with the next loop iteration, since we cannot map cells intercepting the border
                end
                apos = round(apos);
                % called beyond x and beyond y why to check for crossing image border??
                ROIstartx{threadrun} = max(1          ,   min(apos(:, 1)));
                ROIstarty{threadrun} = max(1          ,   min(apos(:, 2)));
                ROIendx{threadrun} = min(imSize(1),   max(apos(:, 1)));
                ROIendy{threadrun} = min(imSize(2),   max(apos(:, 2)));
                ROIsize_x = ROIendx{threadrun} - ROIstartx{threadrun} + 1;
                ROIsize_y = ROIendy{threadrun} - ROIstarty{threadrun} + 1;
                % current cells in the respecitve ROIs
                rpos = zeros(size(apos, 1), 2); % allocate
                rpos(:, 1) = apos(:,1)-ROIstartx{threadrun} + 1;
                rpos(:, 2) = apos(:,2)-ROIstarty{threadrun} + 1;
                if numel(rpos(:, 2)) <= 2 %only 2 vertices found for current cell
                    continue; % with the next cell
                end
                
                % map the roi to the image
                rmask{threadrun} = poly2mask(rpos(:, 2), rpos(:, 1), ROIsize_x, ROIsize_y);
                A = sum(rmask{threadrun}(:)); % area of current ws-cell
                if A == 0 % distance was too small to draw a polygon. Print pixel at the center instead
                    x = round(mean(rpos(:, 1)));
                    y = round(mean(rpos(:, 2)));
                    rmask{threadrun}(x,y) = 1;
                    A = 1;
                end
                % brightness of ws-cell is 1/Area
                rmask{threadrun} = rmask{threadrun}/A;
            end
            % add rois from threads to bild array
            for ii = 1:length(acttriangs)
                if ~isempty(rmask{ii})
                    SRimage(ROIstartx{ii}:ROIendx{ii},ROIstarty{ii}:ROIendy{ii})=SRimage(ROIstartx{ii}:ROIendx{ii},ROIstarty{ii}:ROIendy{ii})+rmask{ii};
                end
            end
        end
        % call the function recursively
        if any(nPertubation) > 0
            if any(nPertubation < 3)
                warning('You specified pertubation = %i. Such low amount of perturbation should be used with caution because of the low and thus bad statistics.\n', nPertubation(1))
            end
            locTableP = zeros([size(locTable, 1) 3 nPertubation(end)]);
            for ii = 1:nPertubation(end)
                locTableP(:, 2, ii) = normrnd(locTable(:, 2), locTable(:, 4));
                locTableP(:, 3, ii) = normrnd(locTable(:, 3), locTable(:, 5));
                % the size may change on the different sets, which is why we calculate
                % it again
                pos = [locTableP(:, 2, ii) locTableP(:, 3, ii)] / pixelSize; % get the locations in pixels, not nm (x and y is swapped in the Orte files)
                imSizeTemp = [ ceil( max(pos(:, 1)) )  ceil( max(pos(:, 2)) )] + 1;
                imSize = max(imSize, imSizeTemp); % and we take whichever values are bigger
            end
            % expanding already generated cellBild to match the new size
            SRimage = cat(2, SRimage, zeros(size(SRimage,1), imSize(2)-size(SRimage, 2)));
            SRimage = cat(1, SRimage, zeros(imSize(1)-size(SRimage, 1), imSize(2)));
            %extend imSize and the image by a third dimension if pert is a vector
            if sliceInd > 0
                imSize = [imSize length(nPertubation)];
                SRimage = cat(3, SRimage, zeros(imSize(1), imSize(2), imSize(3)-1));
            end
            for ii = 1 : nPertubation(end)
                multiWaitbar( 'generating voronoi image...', 'Value', (ii-1) / nPertubation(end));
                % now of course without perturbation, we wouldn't finish ever otherwise
                % with that we can call Orte2CellBild pert times with the new Orte sets
                bildNew = visModule(locTableP(:, :, ii), 'voronoi', pixelSize);
                % (1,1), so top left should be absolute in every image and we can
                % concatenate zeros at the bottom right (matlab convention)
                bildNew = cat(2, bildNew, zeros(size(bildNew,1), imSize(2)-size(bildNew, 2)));
                bildNew = cat(1, bildNew, zeros(imSize(1)-size(bildNew, 1), imSize(2)));
                % summing up the result of the recursions
                if sliceInd > 0 % nPertubation is a vector
                    % due to matlab's stupid indexing we have to subtract 1 for the indexing, aarr!
                    if any(ii-1 == nPertubation(1:end-1)) %except the last, because then we're done already
                        SRimage(:,:,sliceInd + 1) = SRimage(:, :, sliceInd); % copy last generated image to new slice
                        sliceInd = sliceInd + 1; % increment slice index
                    end
                    SRimage(:, :, sliceInd) = SRimage(:, :, sliceInd) + bildNew;
                else % nPertubation is a scalar
                    % no need to care about storing the images and we just add the
                    % results
                    SRimage = SRimage + bildNew;
                end
            end
        end
    case 'densityMap' % from Orte2DensMap
        multiWaitbar('Starting super resolution image reconstruction', 'Relabel', 'generating density map...');
        multiWaitbar('generating density map...', 0);
        % extract the localizations, working with the entries directly might save memory.
        % Using it anyways, though, for the sake of readability
        locVals = locTable(:, 2:3) / pixelSize;
        mapSize = ceil(max(locVals));
        SRimage = nan(mapSize);
        % for waitbar
        prevPercent = 0;
        counter = 1;
        for line = 1:mapSize(1)
            % compute boundaries of current fluorophores (this takes out a stripe
            % of the locations to accelerate processing
            lowBound = max(line-ceil(mapSize(1)/acc), 1);% take only positive values
            uppBound = min(mapSize(1), line+ceil(mapSize(1)/acc));% prevent exceeding image boundaries
            currLocVals = locVals(locVals(:, 1) > lowBound & locVals(:, 1) < uppBound,:);% use a subset only (if acc>1)
            pixMat = [repmat(line,mapSize(2),1) linspace(1,mapSize(2),mapSize(2))']-.5;% coordinates of the pixels
            if calib
                largestDistances = max(pdist2(currLocVals, pixMat, 'euclidean', 'Smallest', numNN));% compute numNN closest distances of one line and return the maximum, type help PDIST2 for weighting
                if isempty(largestDistances); largestDistances = nan(1, size(pixMat, 1)); end % workaround for the case if no values are found (not necessary for mean).
                SRimage(line, :) = largestDistances.^2;
            else
                SRimage(line, :) = mean(pdist2(currLocVals, pixMat, 'euclidean', 'Smallest', numNN));% compute numNN closest distances of one line and build the mean, type help PDIST2 for weighting
            end
            currentPercent = fix(100*counter/mapSize(1));
            if currentPercent > prevPercent
                multiWaitbar( 'generating density map...', 'Value', counter/mapSize(1));
                prevPercent = currentPercent;
            end
            counter = counter + 1;
        end
        if calib
            SRimage(SRimage == 0) = NaN; % prevent division by 0
            SRimage = numNN./(SRimage.*pi*pixelSize^2);% density is numNN/area where area is pi*(largestDistance*pixelsize)^2
        else
            SRimage = 1./(SRimage.*pixelSize);% every pixel value is 1/(the mean distance in nm) to its next numNN neighbours
        end
        % The function might be faster if PDIST2 is not called with a line but with a square ROI
    otherwise
        multiWaitbar('CloseAll');
        error([visMethod,' not available!']); 
end
%% save SRimage routine
% save image as Tiff file, create Tiff Object
% except scatterplot - will be saved as .fig
if saveImage == 1
    switch visMethod % first differ if its figure object or image
        case 'scatterplot'
            saveas(SRimage, [pathname, '\', filename, '.fig'])
            close(SRimage)     
        otherwise
            imageTiff = Tiff([pathname, '\', filename, '.tif'], 'w');
            % set standard tag values
            setTag(imageTiff, 'ImageLength', size(SRimage, 1))
            setTag(imageTiff, 'ImageWidth', size(SRimage, 2))
            setTag(imageTiff, 'Photometric', Tiff.Photometric.MinIsBlack)
            switch visMethod
                case 'histogramBinning'
                    setTag(imageTiff, 'BitsPerSample', 16) % save as 16 bit
                case {'gaussianBlur', 'triangulation', 'voronoi', 'densityMap'}
                    setTag(imageTiff, 'BitsPerSample', 32) % save as 32 bit
            end
            setTag(imageTiff, 'SamplesPerPixel', 1)
            switch visMethod
                case 'histogramBinning'
                    setTag(imageTiff, 'SampleFormat', 1)
                case {'gaussianBlur', 'triangulation', 'voronoi', 'densityMap'}
                    setTag(imageTiff, 'SampleFormat', 3)
            end
            setTag(imageTiff, 'Compression', Tiff.Compression.LZW)
            setTag(imageTiff, 'PlanarConfiguration', Tiff.PlanarConfiguration.Chunky)
            setTag(imageTiff, 'ResolutionUnit', Tiff.ResolutionUnit.Centimeter);
            setTag(imageTiff, 'XResolution', round(1/(pixelSize*10^(-7))));
            setTag(imageTiff, 'YResolution', round(1/(pixelSize*10^(-7))));
            setTag(imageTiff, 'Software', 'MATLAB')
            switch visMethod
                case 'histogramBinning'
                    write(imageTiff, uint16(SRimage)); % transform image data to 16 bit
                case {'gaussianBlur', 'triangulation', 'voronoi', 'densityMap'}
                    write(imageTiff, single(SRimage)); % transform image data to single
            end
            close(imageTiff);
    end
else
    varargout{1} = SRimage;
end
multiWaitbar('CloseAll');
end