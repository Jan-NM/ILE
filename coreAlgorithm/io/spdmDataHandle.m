function c = spdmDataHandle(filename,varargin)
% c = spdmDataHandle(filename,varargin)
%
% Input:
%    [filename]  <1xN char> full or relative path including extension to a 
%                    supported image file format
%                    (presently .tif and .h5).
% 
% Output:
%     [c]       struct containing the following fields
%                 .filename: (char) filename used to open the image file
%                 .isTIF:    (bool)
%                 .simultaneousFrames: (int/double) number of slices 
%                            simultaneously kept in memory 
%                 .nextFrame (int/double)  framenumber of next image slice to read
% 
%    future changes: spdmDataHandle - data stack access class constructor
%
% Author Udo Birk, IMB Mainz gGmbH.
%
% See also FASTSPDMFLUSH

if nargin==1
    c.filename = filename;
    c.filetype = '';
    c.isH5 = false;
    c.isTIF = false;
    c.nextFrame = 0;
    c.filehandle = [];
    c.H5data_address = [];
    [pn,fn,ext] = fileparts(filename);
    c.fn = fn;
    
    switch(lower(ext))
        case {'.tif','.tiff'}
            c.filetype = 'tif';
            c.isTIF = true;
        case '.h5'
            c.filetype = 'h5';
            c.isH5 = true;
        otherwise
            error('unknown file format');
    end
    
    % reduce number of simultaneousFrames to evaluate on older PCs
    c.simultaneousFrames = 500;
    
    c.getImageDimensions = @getImageDimensions;
end

    
function TifFlushFrames(c,startFrame,numFrames)
global d
    if(startFrame + numFrames > c.siz(3))
        numFrames = c.siz(3) - startFrame + 1;
    end
    warning('off', 'all')
    for i=1:numFrames
       c.filehandle.setDirectory(startFrame+i-1);
       d(:,:,i) = c.filehandle.read();
    end
    warning('on', 'all')


function TifCleanUp(fileHandle)
    fileHandle.close();

function H5FlushFrames(c,startFrame,numFrames)
global d
    if(startFrame + numFrames > c.siz(3))
        numFrames = c.siz(3) - startFrame + 1;
    end
    for i=1:numFrames
        d(:,:,i)=(h5read(c.filehandle,['/' c.H5data_address]',[1 1 startFrame+i-1],[c.siz(1) c.siz(2) 1]))';
    end

function H5CleanUp(fileHandle)
    ;
        
function TifStdReadFrames(c,startFrame,numFrames)
global d
    if(startFrame + numFrames > c.siz(3))
        numFrames = c.siz(3) - startFrame + 1;
    end
    for i=1:numFrames
        d(:,:,i)=imread(c.filehandle,'Index',startFrame+i-1, 'Info', c.info);
    end    
    
function TifStdCleanUp(fileHandle)
    ;

    
%%
function c = getImageDimensions(c)

    if (c.isTIF)

         if(verLessThan('matlab', '7.12'))
            c.filehandle  = c.filename;
            c.spdm_hFlushFrames = @TifStdReadFrames;
            c.spdm_hCleanUp     = @TifStdCleanUp;
         else
             c.filehandle  = Tiff(c.filename,'r');
             c.spdm_hFlushFrames = @TifFlushFrames;
             c.spdm_hCleanUp     = @TifCleanUp;
         end
        
        mybar = waitbar(0,'Reading tif header', 'Name', [c.filename ' - progress']);
        c.info = imfinfo(c.filename);
        close(mybar);

        DimZ = length(c.info); % size(info,1);
        DimX = c.info(1).Width;
        DimY = c.info(1).Height;
        DimE = 1;

        siz = [DimX,DimY,DimZ,DimE];
        
    elseif (c.isH5)
        
        c.filehandle = c.filename;
        c.spdm_hFlushFrames = @H5FlushFrames;
        c.spdm_hCleanUp = @H5CleanUp;
        
        info = h5info(c.filename);
        c.H5data_address = info.Datasets(2).Name;
        DimZ = info.Datasets(2).Dataspace.Size(3);
        DimY = info.Datasets(2).Dataspace.Size(2);
        DimX = info.Datasets(2).Dataspace.Size(1);
        DimE = 1;
        
        siz = [DimX,DimY,DimZ,DimE];
        
    else
        error('only tif or h5 files supported');            
    end
    c.siz = siz;
