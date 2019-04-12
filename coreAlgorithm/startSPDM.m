function varargout = startSPDM(varargin)
% STARTSPDM GUI script for fastSPDM evaluation.
%   
%   
% Cremer Group, Institute of Molecular Biology (IMB), Mainz

% Last Modified by GUIDE v2.5 16-Feb-2018 21:58:20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @startSPDM_OpeningFcn, ...
                   'gui_OutputFcn',  @startSPDM_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


%% Function to load content into the listbox 'filebox'
function load_filebox(dir_path, handles)

filter_str = get(handles.fn,'String');   
dir_struct= dir([dir_path '/' filter_str]); 
[sorted_names,sorted_index] = sortrows({dir_struct.name}');

handles.file_names = sorted_names;
handles.is_dir = [dir_struct.isdir];
handles.sorted_index = sorted_index;
guidata(handles.figure1,handles)
set(handles.filebox,'String',handles.file_names,...
 'Value',1)


%% Function to load content into the listbox 'selbox' from the 'filebox'
function load_selbox(fnames, handles)

[~, x] = size(fnames);
[sorted_names,sorted_index] = sortrows(fnames);
handles.file_names = sorted_names;
handles.sorted_index = sorted_index;
guidata(handles.figure1,handles)
set(handles.filenum,'String',[num2str(x) ' files selected'])
set(handles.selbox,'String',handles.file_names,...
 'Value',1)


%% Callback: openingFcn
% --- Executes just before startSPDM is made visible.
function startSPDM_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to startSPDM (see VARARGIN)

% Choose default command line output for startSPDM
global SPDM_GUI
global files;

DIR_input = pwd;
DIR_output = [pwd filesep 'out'];

files ={};
set(handles.displayDIR_input,'String',DIR_input);
set(handles.displayDIR_output,'String',DIR_output);
load_filebox(DIR_input,handles);

handles.output = hObject;

loadState(handles)

guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = startSPDM_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in selDIR_input.
function selDIR_input_Callback(hObject, eventdata, handles)
DIR_input = get(handles.displayDIR_input,'String');
if(isempty(DIR_input))
    startdir = pwd;
else
    startdir = DIR_input;
end
newDIR_input = uigetdir(startdir,'Select Input Directory');
if(newDIR_input == 0) % user pressed cancel
    return;
else
    DIR_input = newDIR_input;
end
DIR_output = [DIR_input filesep 'out'];
set(handles.filebox,'String',DIR_input);
guidata(hObject, handles);
set(handles.displayDIR_input,'String',DIR_input);
set(handles.displayDIR_output,'String',DIR_output);
load_filebox(DIR_input,handles);


% --- Executes on selection change in filebox.
function filebox_Callback(hObject, eventdata, handles)
global files;
contents = cellstr(get(hObject,'String'));
if(isempty(contents))
    return;
end
newfilename = contents{get(hObject,'Value')};
myfiles = get(handles.selbox,'String');
x = strmatch(newfilename, myfiles);  %% check if file is already selected
if(isempty(x))                  
    files(end+1) ={newfilename};
    load_selbox(files, handles)
end


% --- Executes during object creation, after setting all properties.
function filebox_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in selbox.
function selbox_Callback(hObject, eventdata, handles)
global files;
[a b ]= size(files);

contents = cellstr(get(hObject,'String'));
if(isempty(contents))
    return;
end
filenames = contents{get(hObject,'Value')};
newfiles = {};
j=1;
for i = 1:b    
    if (strcmp(files(i),filenames)==0)
        newfiles(j) = files(i);     
        j=j+1;
    end
end
files = newfiles;
load_selbox(files, handles)


% --- Executes during object creation, after setting all properties.
function selbox_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function fn_Callback(hObject, eventdata, handles)
DIR_input = get(handles.displayDIR_input,'String');
load_filebox(DIR_input,handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function fn_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in refresh.
function refresh_Callback(hObject, eventdata, handles)
DIR_input = get(handles.displayDIR_input,'String');
load_filebox(DIR_input,handles);


% --- Executes on button press in loadall.
function loadall_Callback(hObject, eventdata, handles)
global files;
files = cellstr(get(handles.filebox,'String'))';
load_selbox(files, handles);


% --- Executes on button press in selDIR_output.
function selDIR_output_Callback(hObject, eventdata, handles)
% global DIR_output
DIR_output = get(handles.displayDIR_output,'String');
newDIR_output = uigetdir('','Select Output Directory');
if(newDIR_output == 0)
    return;
else
    DIR_output = newDIR_output;
end
set(handles.displayDIR_output,'String',DIR_output);
% hObject    handle to selDIR_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%% Main routine: START BUTTON PRESSED
% #######################################################
% ###         start main evaluation routine          ####
% #######################################################
% --- Executes on button press in startbut.
function startbut_Callback(hObject, eventdata, handles)

global files
global SPDMparamstruct

set(handles.startbut,'Enable','off');
drawnow;
set(handles.errors,'String','no Problems so far...');
%--------------------------------------------------------------------------
% initializing spdm paramters
SPDMparamstruct = struct('CAMpixelsize', 102, ...
                     'CAMconvfactor',2, ...
                     'OUTPUTpixelsize',20);

% initialize callbacks
SPDMparamstruct.displayStatus = @(text) displayStatus(handles.displaystatus,text);
SPDMparamstruct.displayError  = @(text) displayStatus(handles.errors,text);
SPDMparamstruct.enableStartButton = @() enableStartButton(handles.startbut);

% data subset   
SPDMparamstruct.STACKstartFrame = str2num(get(handles.startFrameNo,'String')); 
SPDMparamstruct.STACKendFrame   = str2num(get(handles.endFrameNo,'String'));
if(isempty(SPDMparamstruct.STACKendFrame))
    SPDMparamstruct.STACKendFrame = Inf;
end

% data properties
SPDMparamstruct.CAMpixelsize   = str2num(get(handles.dataCAMpixelsize,'String'));
SPDMparamstruct.CAMconvfactor  = str2num(get(handles.dataPhotonsPerCount,'String'));
SPDMparamstruct.qe  = str2num(get(handles.qe, 'String'));
SPDMparamstruct.EMgain  = str2num(get(handles.EMgain, 'String'));
SPDMparamstruct.CAMtypeEMCCD   = get(handles.dataEMCCDcheckbox,'Value');
SPDMparamstruct.CAMoffset = str2num(get(handles.dataCountsOffset, 'String'));

% localization
SPDMparamstruct.fastSPDMthreshold = str2num(get(handles.thfactor,'String'));

% corrections
SPDMparamstruct.flag_multiframe = (get(handles.multiframe,'value'));
SPDMparamstruct.searchRadius = str2num(get(handles.searchRadius,'String'));
SPDMparamstruct.flag_filterSignals = (get(handles.filterSignals,'value'));
SPDMparamstruct.upperLimitPSF = str2num(get(handles.filterMeanPSF,'String'));
SPDMparamstruct.upperLimitLoc = str2num(get(handles.filterMeanLoc,'String'));
SPDMparamstruct.limitPhot = str2num(get(handles.filterPhot,'String'));
SPDMparamstruct.flag_driftCorrection = get(handles.driftCorrection,'value');
SPDMparamstruct.subsetNo = str2num(get(handles.subsetNo,'String'));

% visualization
SPDMparamstruct.OUTPUTscatterplot = (get(handles.visScatterplot, 'value'));
SPDMparamstruct.OUTPUThistogram = (get(handles.visHistBinning, 'value'));
SPDMparamstruct.OUTPUTgaussianBlur = (get(handles.visGaussianBlur, 'value'));
switch get(get(handles.blurringMethod, 'SelectedObject'), 'Tag')
    case 'specificBlur'
        SPDMparamstruct.blurringMethod = 'specificBlur';
    case 'globalBlur'
        SPDMparamstruct.blurringMethod = 'globalBlur';
    case 'individualBlur'
        SPDMparamstruct.blurringMethod = 'individualBlur';
end
SPDMparamstruct.blurringSize = str2num(get(handles.blurringSize, 'String'));
SPDMparamstruct.OUTPUTtriangulation = (get(handles.visTriangulation, 'value'));
SPDMparamstruct.nPertubation = str2num(get(handles.numPertubation, 'String'));
SPDMparamstruct.OUTPUTvoronoi = (get(handles.visVoronoi, 'value'));
SPDMparamstruct.OUTPUTdensity = (get(handles.visDensity, 'value'));
SPDMparamstruct.calib = (get(handles.visCalibration, 'value'));
SPDMparamstruct.acc = str2num(get(handles.visPartitions, 'String'));
SPDMparamstruct.numNN = str2num(get(handles.visNN, 'String'));

SPDMparamstruct.dim = 2;
SPDMparamstruct.saveImage = 1;
SPDMparamstruct.zStep = 0;
SPDMparamstruct.OUTPUTpixelsize = str2num(get(handles.bildpixel,'String'));

SPDMparamstruct.driftCorrPlot = (get(handles.driftCorrPlot,'value'));

set(handles.displaystatus,'String','Evaluation running...');
%--------------------------------------------------------------------------
%% - create outputfolder if not existing
SPDMparamstruct.DIR_input = get(handles.displayDIR_input,'String');
SPDMparamstruct.DIR_output = get(handles.displayDIR_output,'String');
files = cellstr(get(handles.selbox,'String'))';
if(strcmp(get(handles.filenum,'String'),'no files selected'))
    files = {};
end
if(isempty(files))
    set(handles.displaystatus,'String','Error - Evaluation not successful');
    set(handles.errors,'String','No File Selected');
    set(handles.startbut,'Enable','on');
    error('No files selected');
end

if(~isdir(SPDMparamstruct.DIR_output))  
    mkdir(SPDMparamstruct.DIR_output); 
    disp('output folder created')
end
%% - start evaluation - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SPDMparamstruct.handles = handles;
SPDMparamstruct.outputmode = 'gui';
if get(handles.ortesel,'value') == 1
    SPDMparamstruct.orteSel = 'Orte';
else
    SPDMparamstruct.orteSel = 'image';
end

SPDMprocessing
saveState(handles, SPDMparamstruct.DIR_output);
set(handles.displaystatus,'String','Evaluation finished!');
set(handles.startbut,'Enable','on');


%% parameter settings - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dataCAMpixelsize_Callback(hObject, eventdata, handles)
input = str2num(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function dataCAMpixelsize_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function dataPhotonsPerCount_Callback(hObject, eventdata, handles)
input = str2num(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function dataPhotonsPerCount_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function bildpixel_Callback(hObject, eventdata, handles)
input = str2num(get(hObject,'String'));
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function bildpixel_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function darstellung_CreateFcn(hObject, eventdata, handles)
% hObject    handle to darstellung (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in visGaussianBlur.
function visGaussianBlur_Callback(hObject, eventdata, handles)
% hObject    handle to visGaussianBlur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(get(hObject,'Value'))
    set(handles.specificBlur, 'Enable', 'on')
    set(handles.globalBlur, 'Enable', 'on')
    set(handles.individualBlur, 'Enable', 'on')
    set(handles.blurringSize, 'Enable', 'on')
    set(handles.text59,'Enable','on')
else
    set(handles.specificBlur, 'Enable', 'off')
    set(handles.globalBlur, 'Enable', 'off')
    set(handles.individualBlur, 'Enable', 'off')
    set(handles.blurringSize, 'Enable', 'off')
    set(handles.text59,'Enable','off')
end


% --- Executes on button press in multiframe.
function multiframe_Callback(hObject, eventdata, handles)
% hObject    handle to multiframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of multiframe
if(get(hObject,'Value'))
    set(handles.searchRadius,'enable','on');
    set(handles.text58,'Enable','on')
else
    set(handles.searchRadius,'enable','off');
    set(handles.text58,'Enable','off')
end


% --- Executes when selected object is changed in filetypegroup.
function filetypegroup_SelectionChangeFcn(hObject, eventdata, handles)
switch get(eventdata.NewValue,'Tag') 
    case 'tifsel'
        set(handles.fn,'String','*.tif');
        set(handles.text30,'Enable','on')
        set(handles.text45,'Enable','on')
        set(handles.startFrameNo,'Enable','on')
        set(handles.endFrameNo,'Enable','on')
        set(handles.text19,'Enable','on')
        set(handles.text44,'Enable','on')
        set(handles.text8,'Enable','on')
        set(handles.dataCAMpixelsize,'Enable','on')
        set(handles.dataCountsOffset,'Enable','on')
        set(handles.dataPhotonsPerCount,'Enable','on')
        set(handles.dataEMCCDcheckbox,'Enable','on')
        set(handles.text28,'Enable','on')
        set(handles.thfactor,'Enable','on')
        set(handles.qe,'Enable','on')
        set(handles.text62,'Enable','on')
        if get(handles.dataEMCCDcheckbox,'Value')
            set(handles.EMgain,'Enable','on')
            set(handles.text63,'Enable','on')
        end
    case 'h5sel'
        set(handles.fn,'String','*.h5');
        set(handles.text30,'Enable','on')
        set(handles.text45,'Enable','on')
        set(handles.startFrameNo,'Enable','on')
        set(handles.endFrameNo,'Enable','on')
        set(handles.text19,'Enable','on')
        set(handles.text44,'Enable','on')
        set(handles.text8,'Enable','on')
        set(handles.dataCAMpixelsize,'Enable','on')
        set(handles.dataCountsOffset,'Enable','on')
        set(handles.dataPhotonsPerCount,'Enable','on')
        set(handles.dataEMCCDcheckbox,'Enable','on')
        set(handles.text28,'Enable','on')
        set(handles.thfactor,'Enable','on')
        set(handles.qe,'Enable','on')
        set(handles.text62,'Enable','on')
        if get(handles.dataEMCCDcheckbox,'Value')
            set(handles.EMgain,'Enable','on')
            set(handles.text63,'Enable','on')
        end
    case 'ortesel'
        set(handles.fn,'String','*.mat');
        set(handles.text30,'Enable','off')
        set(handles.text45,'Enable','off')
        set(handles.startFrameNo,'Enable','off')
        set(handles.endFrameNo,'Enable','off')
        set(handles.text44,'Enable','off')
        set(handles.text8,'Enable','off')
        set(handles.text19,'Enable','off')
        set(handles.dataCAMpixelsize,'Enable','off')
        set(handles.dataCountsOffset,'Enable','off')
        set(handles.dataPhotonsPerCount,'Enable','off')
        set(handles.dataEMCCDcheckbox,'Enable','off')
        set(handles.text28,'Enable','off')
        set(handles.thfactor,'Enable','off')
        set(handles.qe,'Enable','off')
        set(handles.EMgain,'Enable','off')
        set(handles.text62,'Enable','off')
        set(handles.text63,'Enable','off')
end
DIR_input = get(handles.displayDIR_input,'String');
load_filebox(DIR_input,handles);
% hObject    handle to the selected object in filetypegroup 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function filetypegroup_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to filetypegroup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function filetypegroup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filetypegroup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


function thfactor_Callback(hObject, eventdata, handles)
% hObject    handle to thfactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thfactor as text
%        str2double(get(hObject,'String')) returns contents of thfactor as a double


% --- Executes during object creation, after setting all properties.
function thfactor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thfactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function startFrameNo_Callback(hObject, eventdata, handles)
% hObject    handle to startFrameNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of startFrameNo as text
%        str2double(get(hObject,'String')) returns contents of startFrameNo as a double


% --- Executes during object creation, after setting all properties.
function startFrameNo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to startFrameNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function uimenuEnableStartButton_Callback(hObject, eventdata, handles)
% hObject    handle to uimenuEnableStartButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.startbut,'Enable','on');


% --- Executes on button press in debugbutton.
function debugbutton_Callback(hObject, eventdata, handles)
% hObject    handle to debugbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
keyboard

 
function displayStatus(hObject,text)
%    call this function with hObject = handles.displaystatus
    set(hObject,'String',text);


function enableStartButton(hObject)
%    call this function with hObject = handles.startbut
    set(hObject,'Enable','on');


% --- Executes on button press in visNLBild.
function visNLBild_Callback(hObject, eventdata, handles)
% hObject    handle to visNLBild (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of visNLBild


function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit22 as text
%        str2double(get(hObject,'String')) returns contents of edit22 as a double


% --- Executes during object creation, after setting all properties.
function edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit23_Callback(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit23 as text
%        str2double(get(hObject,'String')) returns contents of edit23 as a double


% --- Executes during object creation, after setting all properties.
function edit23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox10.
function listbox10_Callback(hObject, eventdata, handles)
% hObject    handle to listbox10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox10 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox10


% --- Executes during object creation, after setting all properties.
function listbox10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function order_Callback(hObject, eventdata, handles)
% hObject    handle to order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of order as text
%        str2double(get(hObject,'String')) returns contents of order as a double


% --- Executes during object creation, after setting all properties.
function order_CreateFcn(hObject, eventdata, handles)
% hObject    handle to order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function dataCountsOffset_Callback(hObject, eventdata, handles)
% hObject    handle to dataCountsOffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dataCountsOffset as text
%        str2double(get(hObject,'String')) returns contents of dataCountsOffset as a double


% --- Executes during object creation, after setting all properties.
function dataCountsOffset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dataCountsOffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in dataEMCCDcheckbox.
function dataEMCCDcheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to dataEMCCDcheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dataEMCCDcheckbox
if(get(hObject,'Value'))
    set(handles.EMgain, 'Enable','on')
    set(handles.text63, 'Enable','on')
else
    set(handles.EMgain, 'Enable','off')
    set(handles.text63, 'Enable','off')
end

function endFrameNo_Callback(hObject, eventdata, handles)
% hObject    handle to endFrameNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of endFrameNo as text
%        str2double(get(hObject,'String')) returns contents of endFrameNo as a double


% --- Executes during object creation, after setting all properties.
function endFrameNo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to endFrameNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function searchRadius_Callback(hObject, eventdata, handles)
% hObject    handle to searchRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of searchRadius as text
%        str2double(get(hObject,'String')) returns contents of searchRadius as a double


% --- Executes during object creation, after setting all properties.
function searchRadius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to searchRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in filterSignals.
function filterSignals_Callback(hObject, eventdata, handles)
% hObject    handle to filterSignals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of filterSignals
if(get(hObject,'Value'))
    set(handles.filterMeanPSF,'Enable','on')
    set(handles.text51,'Enable','on')
    set(handles.filterMeanLoc,'Enable','on')
    set(handles.text52,'Enable','on')
    set(handles.filterPhot,'Enable','on')
    set(handles.text57,'Enable','on')
else
    set(handles.filterMeanPSF,'Enable','off')
    set(handles.text51,'Enable','off')
    set(handles.filterMeanLoc,'Enable','off')
    set(handles.text52,'Enable','off')
    set(handles.filterPhot,'Enable','off')
    set(handles.text57,'Enable','off')
end


function filterMeanPSF_Callback(hObject, eventdata, handles)
% hObject    handle to filterMeanPSF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filterMeanPSF as text
%        str2double(get(hObject,'String')) returns contents of filterMeanPSF as a double


% --- Executes during object creation, after setting all properties.
function filterMeanPSF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filterMeanPSF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function filterMeanLoc_Callback(hObject, eventdata, handles)
% hObject    handle to filterMeanLoc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filterMeanLoc as text
%        str2double(get(hObject,'String')) returns contents of filterMeanLoc as a double


% --- Executes during object creation, after setting all properties.
function filterMeanLoc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filterMeanLoc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function subsetNo_Callback(hObject, eventdata, handles)
% hObject    handle to subsetNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of subsetNo as text
%        str2double(get(hObject,'String')) returns contents of subsetNo as a double


% --- Executes during object creation, after setting all properties.
function subsetNo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to subsetNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in visHistBinning.
function visHistBinning_Callback(hObject, eventdata, handles)
% hObject    handle to visHistBinning (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of visHistBinning


% --- Executes on button press in driftCorrection.
function driftCorrection_Callback(hObject, eventdata, handles)
% hObject    handle to driftCorrection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of driftCorrection
if(get(hObject,'Value'))
    set(handles.subsetNo,'Enable','on')
    set(handles.text53,'Enable','on')
    set(handles.driftCorrPlot,'Enable','on')
else
    set(handles.subsetNo,'Enable','off')
    set(handles.text53,'Enable','off')
    set(handles.driftCorrPlot,'Enable','off')
end


% --- Executes on button press in cropOrte.
function cropOrte_Callback(hObject, eventdata, handles)
% hObject    handle to cropOrte (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cropImage;


% --- Executes on button press in helpbutton.
function helpbutton_Callback(hObject, eventdata, handles)
% hObject    handle to helpbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
winopen('instructions.pdf');


% --- Executes on button press in driftCorrPlot.
function driftCorrPlot_Callback(hObject, eventdata, handles)
% hObject    handle to driftCorrPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of driftCorrPlot


function numPertubation_Callback(hObject, eventdata, handles)
% hObject    handle to numPertubation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numPertubation as text
%        str2double(get(hObject,'String')) returns contents of numPertubation as a double


% --- Executes during object creation, after setting all properties.
function numPertubation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numPertubation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function filterPhot_Callback(hObject, eventdata, handles)
% hObject    handle to filterPhot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filterPhot as text
%        str2double(get(hObject,'String')) returns contents of filterPhot as a double


% --- Executes during object creation, after setting all properties.
function filterPhot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filterPhot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in visTriangulation.
function visTriangulation_Callback(hObject, eventdata, handles)
% hObject    handle to visTriangulation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of visTriangulation
if(get(hObject,'Value'))
    set(handles.numPertubation,'Enable','on')
    set(handles.text56,'Enable','on')
else
    if ~get(handles.visVoronoi, 'Value')
        set(handles.numPertubation,'Enable','off')
        set(handles.text56,'Enable','off')
    end
end

% --- Executes on button press in visScatterplot.
function visScatterplot_Callback(hObject, eventdata, handles)
% hObject    handle to visScatterplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of visScatterplot



function blurringSize_Callback(hObject, eventdata, handles)
% hObject    handle to blurringSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of blurringSize as text
%        str2double(get(hObject,'String')) returns contents of blurringSize as a double


% --- Executes during object creation, after setting all properties.
function blurringSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to blurringSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in visVoronoi.
function visVoronoi_Callback(hObject, eventdata, handles)
% hObject    handle to visVoronoi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of visVoronoi
if(get(hObject,'Value'))
    set(handles.numPertubation,'Enable','on')
    set(handles.text56,'Enable','on')
else
    if ~get(handles.visTriangulation, 'Value')
        set(handles.numPertubation,'Enable','off')
        set(handles.text56,'Enable','off')
    end
end

% --- Executes on button press in visDensity.
function visDensity_Callback(hObject, eventdata, handles)
% hObject    handle to visDensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of visDensity
if(get(hObject,'Value'))
    set(handles.visPartitions,'Enable','on')
    set(handles.visNN,'Enable','on')
    set(handles.visCalibration,'Enable','on')
    set(handles.text61,'Enable','on')
    set(handles.text60,'Enable','on')
else
    set(handles.visPartitions,'Enable','off')
    set(handles.visNN,'Enable','off')
    set(handles.visCalibration,'Enable','off')
    set(handles.text61,'Enable','off')
    set(handles.text60,'Enable','off')
end

% --- Executes on button press in visCalibration.
function visCalibration_Callback(hObject, eventdata, handles)
% hObject    handle to visCalibration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of visCalibration



function visNN_Callback(hObject, eventdata, handles)
% hObject    handle to visNN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of visNN as text
%        str2double(get(hObject,'String')) returns contents of visNN as a double


% --- Executes during object creation, after setting all properties.
function visNN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to visNN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function visPartitions_Callback(hObject, eventdata, handles)
% hObject    handle to visPartitions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of visPartitions as text
%        str2double(get(hObject,'String')) returns contents of visPartitions as a double


% --- Executes during object creation, after setting all properties.
function visPartitions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to visPartitions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in specificBlur.
function specificBlur_Callback(hObject, eventdata, handles)
% hObject    handle to specificBlur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of specificBlur


% --- Executes on button press in globalBlur.
function globalBlur_Callback(hObject, eventdata, handles)
% hObject    handle to globalBlur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of globalBlur


% --- Executes on button press in individualBlur.
function individualBlur_Callback(hObject, eventdata, handles)
% hObject    handle to individualBlur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of individualBlur



function qe_Callback(hObject, eventdata, handles)
% hObject    handle to qe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of qe as text
%        str2double(get(hObject,'String')) returns contents of qe as a double


% --- Executes during object creation, after setting all properties.
function qe_CreateFcn(hObject, eventdata, handles)
% hObject    handle to qe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EMgain_Callback(hObject, eventdata, handles)
% hObject    handle to EMgain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EMgain as text
%        str2double(get(hObject,'String')) returns contents of EMgain as a double


% --- Executes during object creation, after setting all properties.
function EMgain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EMgain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
saveState(handles) % saves current state of figure if figure is closed
delete(hObject);


function saveState(handles, location)
    % data subset   
    state.STACKstartFrame = get(handles.startFrameNo,'String'); 
    state.STACKendFrame = get(handles.endFrameNo,'String');
    
    % data properties
    state.CAMpixelsize = get(handles.dataCAMpixelsize,'String');
    state.CAMconvfactor = get(handles.dataPhotonsPerCount,'String');
    state.qe = get(handles.qe, 'String');
    state.EMgain = get(handles.EMgain, 'String');
    state.CAMtypeEMCCD = get(handles.dataEMCCDcheckbox, 'Value');
    state.CAMoffset = (get(handles.dataCountsOffset, 'String'));
    % localization
    state.fastSPDMthreshold = get(handles.thfactor,'String');

    % corrections
    state.flag_multiframe = get(handles.multiframe,'value');
    state.searchRadius = get(handles.searchRadius,'String');
    state.flag_filterSignals = get(handles.filterSignals,'value');
    state.upperLimitPSF = get(handles.filterMeanPSF,'String');
    state.upperLimitLoc = get(handles.filterMeanLoc,'String');
    state.limitPhot = get(handles.filterPhot,'String');
    state.flag_driftCorrection = get(handles.driftCorrection,'value');
    state.subsetNo = get(handles.subsetNo,'String');

    % visualization
    state.OUTPUTscatterplot = get(handles.visScatterplot, 'value');
    state.OUTPUThistogram = get(handles.visHistBinning, 'value');
    state.OUTPUTgaussianBlur = get(handles.visGaussianBlur, 'value');
    state.blurringMethod = get(get(handles.blurringMethod, 'SelectedObject'), 'Tag');
    state.blurringSize = get(handles.blurringSize, 'String');
    state.OUTPUTtriangulation = get(handles.visTriangulation, 'value');
    state.nPertubation = get(handles.numPertubation, 'String');
    state.OUTPUTvoronoi = get(handles.visVoronoi, 'value');
    state.OUTPUTdensity = get(handles.visDensity, 'value');
    state.calib = get(handles.visCalibration, 'value');
    state.acc = get(handles.visPartitions, 'String');
    state.numNN = get(handles.visNN, 'String');

    state.dim = 2;
    state.saveImage = 1;
    state.zStep = 0;
    state.OUTPUTpixelsize = str2num(get(handles.bildpixel,'String'));

    state.driftCorrPlot = (get(handles.driftCorrPlot,'value'));
    if nargin > 1
        save([location, '/', 'out_state.mat'], 'state')
    else
        save([mfilename('fullpath') '_state.mat'], 'state')
    end
        

function loadState(handles)
    fileName = [mfilename('fullpath') '_state.mat'];
    if exist(fileName)
        load(fileName)
        
        % data subset   
        set(handles.startFrameNo,'String', state.STACKstartFrame); 
        set(handles.endFrameNo,'String', state.STACKendFrame);

        % data properties
        set(handles.dataCAMpixelsize,'String', state.CAMpixelsize);
        set(handles.dataPhotonsPerCount,'String', state.CAMconvfactor);
        set(handles.qe, 'String', state.qe);
        set(handles.EMgain, 'String', state.EMgain);
        set(handles.dataEMCCDcheckbox, 'Value', state.CAMtypeEMCCD);
        set(handles.dataCountsOffset, 'String', state.CAMoffset);
        
        % localization
        set(handles.thfactor,'String', state.fastSPDMthreshold);

        % corrections
        set(handles.multiframe,'value', state.flag_multiframe);
        set(handles.searchRadius,'String', state.searchRadius);
        set(handles.filterSignals,'value', state.flag_filterSignals);
        set(handles.filterMeanPSF,'String', state.upperLimitPSF);
        set(handles.filterMeanLoc,'String', state.upperLimitLoc);
        set(handles.filterPhot,'String', state.limitPhot);
        set(handles.driftCorrection,'value', state.flag_driftCorrection);
        set(handles.subsetNo,'String', state.subsetNo);

        % visualization
        set(handles.visScatterplot, 'value', state.OUTPUTscatterplot);
        set(handles.visHistBinning, 'value', state.OUTPUThistogram);
        set(handles.visGaussianBlur, 'value', state.OUTPUTgaussianBlur);
        set(handles.blurringSize, 'String', state.blurringSize);
        set(handles.visTriangulation, 'value', state.OUTPUTtriangulation);
        set(handles.numPertubation, 'String', state.nPertubation);
        set(handles.visVoronoi, 'value', state.OUTPUTvoronoi);
        set(handles.visDensity, 'value', state.OUTPUTdensity);
        set(handles.visCalibration, 'value', state.calib);
        set(handles.visPartitions, 'String', state.acc);
        set(handles.visNN, 'String', state.numNN);

        state.dim = 2;
        state.saveImage = 1;
        state.zStep = 0;
        set(handles.bildpixel,'String', state.OUTPUTpixelsize);

        set(handles.driftCorrPlot,'value', state.driftCorrPlot);
        
        delete(fileName) 
    end
