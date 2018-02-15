function varargout = LSRanalyze(varargin)
%LSRANALYZE M-file for LSRanalyze.fig
%      LSRANALYZE, by itself, creates a new LSRANALYZE or raises the existing
%      singleton*.
%
%      H = LSRANALYZE returns the handle to a new LSRANALYZE or the handle to
%      the existing singleton*.
%
%      LSRANALYZE('Property','Value',...) creates a new LSRANALYZE using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to LSRanalyze_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      LSRANALYZE('CALLBACK') and LSRANALYZE('CALLBACK',hObject,...) call the
%      local function named CALLBACK in LSRANALYZE.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LSRanalyze

% Last Modified by GUIDE v2.5 23-Dec-2013 17:05:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LSRanalyze_OpeningFcn, ...
                   'gui_OutputFcn',  @LSRanalyze_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before LSRanalyze is made visible.
function LSRanalyze_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for LSRanalyze
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes LSRanalyze wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = LSRanalyze_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function FilesPath_Callback(hObject, eventdata, handles)
% hObject    handle to FilesPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FilesPath as text
%        str2double(get(hObject,'String')) returns contents of FilesPath as a double


% --- Executes during object creation, after setting all properties.
function FilesPath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FilesPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fishSet_Callback(hObject, eventdata, handles)
% hObject    handle to fishSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fishSet as text
%        str2double(get(hObject,'String')) returns contents of fishSet as a double
if (strcmp(get(hObject,'String'),''))
    set(hObject,'String','[:]');
end

% --- Executes during object creation, after setting all properties.
function fishSet_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fishSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in FilesBut.
function FilesBut_Callback(hObject, eventdata, handles)
% hObject    handle to FilesBut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.FilesPath,'String','No files selected...');
[FileName,PathName] = uigetfile( ...
{  '*.mat','mat files (*.mat)'; ...
   '*.*',  'All Files (*.*)'}, ...
   'Pick the movies you want to analyze', ...
   'MultiSelect', 'off');
fileList = fullfile(PathName, FileName);
set(handles.FilesPath,'String',fileList);
set(handles.FileName,'String',FileName);


% --- Executes on button press in goBut.
function goBut_Callback(hObject, eventdata, handles)
% hObject    handle to goBut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject,'Enable','off');
set(handles.FilesBut,'Enable','off');
set(handles.FilesPath,'Enable','off');
set(handles.OutputBut,'Enable','off');
set(handles.OutputPath,'Enable','off');

runAnalysis3(handles);

set(hObject,'Enable','on');
set(handles.FilesBut,'Enable','on');
set(handles.FilesPath,'Enable','on');
set(handles.OutputBut,'Enable','on');
set(handles.OutputPath,'Enable','on');

% --- Executes during object creation, after setting all properties.
function fileList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fileList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if (str2double(get(hObject,'String')) < 0 || str2double(get(hObject,'String')) >1)
    set(hObject,'String',0.05);
end


function noiseThresh_Callback(hObject, eventdata, handles)
% hObject    handle to noiseThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of noiseThresh as text
%        str2double(get(hObject,'String')) returns contents of noiseThresh as a double


% --- Executes during object creation, after setting all properties.
function noiseThresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to noiseThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ppoTime_Callback(hObject, eventdata, handles)
% hObject    handle to ppoTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ppoTime as text
%        str2double(get(hObject,'String')) returns contents of ppoTime as a double
if (strcmp(get(hObject,'String'),''))
    set(hObject,'String','[1:100]');
end

% --- Executes during object creation, after setting all properties.
function ppoTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ppoTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plotNoise.
function plotNoise_Callback(hObject, eventdata, handles)
% hObject    handle to plotNoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotNoise


% --- Executes on button press in plotGMVOT.
function plotGMVOT_Callback(hObject, eventdata, handles)
% hObject    handle to plotGMVOT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotGMVOT


% --- Executes on button press in plotGPTOT.
function plotGPTOT_Callback(hObject, eventdata, handles)
% hObject    handle to plotGPTOT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotGPTOT


% --- Executes on button press in plotHeatmap.
function plotHeatmap_Callback(hObject, eventdata, handles)
% hObject    handle to plotHeatmap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotHeatmap


% --- Executes on button press in plotIntensities.
function plotIntensities_Callback(hObject, eventdata, handles)
% hObject    handle to plotIntensities (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotIntensities


% --- Executes on button press in plotHistograms.
function plotHistograms_Callback(hObject, eventdata, handles)
% hObject    handle to plotHistograms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotHistograms


% --- Executes on button press in plotPPO.
function plotPPO_Callback(hObject, eventdata, handles)
% hObject    handle to plotPPO (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotPPO



function OutputPath_Callback(hObject, eventdata, handles)
% hObject    handle to OutputPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of OutputPath as text
%        str2double(get(hObject,'String')) returns contents of OutputPath as a double


% --- Executes during object creation, after setting all properties.
function OutputPath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OutputPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in OutputBut.
function OutputBut_Callback(hObject, eventdata, handles)
% hObject    handle to OutputBut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
outputPath = uigetdir(pwd);
if (outputPath == 0)
    outputPath = 'No directory selected...';
end
set(handles.OutputPath,'String',outputPath);



function FileName_Callback(hObject, eventdata, handles)
% hObject    handle to FileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FileName as text
%        str2double(get(hObject,'String')) returns contents of FileName as a double


% --- Executes during object creation, after setting all properties.
function FileName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plotColormap.
function plotColormap_Callback(hObject, eventdata, handles)
% hObject    handle to plotColormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotColormap


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over goBut.
function goBut_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to goBut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in plotGQPTOT.
function plotGQPTOT_Callback(hObject, eventdata, handles)
% hObject    handle to plotGQPTOT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotGQPTOT


% --- Executes on button press in plotGMQOT.
function plotGMQOT_Callback(hObject, eventdata, handles)
% hObject    handle to plotGMQOT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotGMQOT
