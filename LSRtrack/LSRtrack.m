function varargout = LSRtrack(varargin)
%%%%%%%%%%%%%% VERSION 1.2

% LSRTRACK M-file for LSRtrack.fig
%      LSRTRACK, by itself, creates a new LSRTRACK or raises the existing
%      singleton*.
%
%      H = LSRTRACK returns the handle to a new LSRTRACK or the handle to
%      the existing singleton*.
%
%      LSRTRACK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LSRTRACK.M with the given input arguments.
%
%      LSRTRACK('Property','Value',...) creates a new LSRTRACK or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LSRtrack_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LSRtrack_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LSRtrack

% Last Modified by GUIDE v2.5 10-Oct-2013 11:58:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LSRtrack_OpeningFcn, ...
                   'gui_OutputFcn',  @LSRtrack_OutputFcn, ...
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


% --- Executes just before LSRtrack is made visible.
function LSRtrack_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LSRtrack (see VARARGIN)

% Choose default command line output for LSRtrack
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes LSRtrack wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = LSRtrack_OutputFcn(hObject, eventdata, handles) 
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


% --- Executes on button press in FilesBut.
function FilesBut_Callback(hObject, eventdata, handles)
% hObject    handle to FilesBut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.FilesPath,'String','No files selected...');
[files, filesPath] = uigetfile( ...
{  '*.mp4;*.avi','MP4-,AVI-Files (*.mp4, *.avi)'; ...
   '*.*',  'All Files (*.*)'}, ...
   'Pick the movies you want to analyze', ...
   'MultiSelect', 'on');
set(handles.directory,'String',filesPath);
set(handles.fileList,'String',files);
if iscell(files)
    set(handles.FilesPath,'String',strcat(int2str(size(files,2)),' files selected'));
elseif ischar(files)
    set(handles.FilesPath,'String',files);
end

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




function wellThresh_Callback(hObject, eventdata, handles)
% hObject    handle to wellThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wellThresh as text
%        str2double(get(hObject,'String')) returns contents of wellThresh as a double
if (str2double(get(hObject,'String')) < 0 || str2double(get(hObject,'String')) >200000)
    set(hObject,'String',500);
end
set(handles.reAlign,'String','True');


% --- Executes during object creation, after setting all properties.
function wellThresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wellThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fishThresh_Callback(hObject, eventdata, handles)
% hObject    handle to fishThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fishThresh as text
%        str2double(get(hObject,'String')) returns contents of fishThresh as a double
if (str2double(get(hObject,'String')) < 0 || str2double(get(hObject,'String')) >10)
    set(hObject,'String',5);
end

% --- Executes during object creation, after setting all properties.
function fishThresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fishThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trackingThresh_Callback(hObject, eventdata, handles)
% hObject    handle to trackingThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trackingThresh as text
%        str2double(get(hObject,'String')) returns contents of trackingThresh as a double
if (str2double(get(hObject,'String')) < 0 || str2double(get(hObject,'String')) >1)
    set(hObject,'String',0.75);
end

% --- Executes during object creation, after setting all properties.
function trackingThresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trackingThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function watchWell_Callback(hObject, eventdata, handles)
% hObject    handle to watchWell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of watchWell as text
%        str2double(get(hObject,'String')) returns contents of watchWell as a double
if (str2double(get(hObject,'String')) < 1 || str2double(get(hObject,'String')) >96 || isnan(str2double(get(hObject,'String'))))
    set(hObject,'String','1');
end
axes(handles.WatchWellFig);
cla();
text(40,60,sprintf('%s %s','Will watch well',get(handles.watchWell,'String')));

% --- Executes during object creation, after setting all properties.
function watchWell_CreateFcn(hObject, eventdata, handles)
% hObject    handle to watchWell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in displayChk.
function displayChk_Callback(hObject, eventdata, handles)
% hObject    handle to displayChk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of displayChk
if (get(hObject,'Value') == 1.0)
    set(handles.watchWell,'Enable','on');
    axes(handles.WatchWellFig);
    cla();
    text(40,60,sprintf('%s %s','Will watch well',get(handles.watchWell,'String')));
else
    set(handles.watchWell,'Enable','off');
    %figure(handles.WatchWellFig);
    axes(handles.WatchWellFig);
    cla();
    text(40,60,'No well to watch');
end    

% --- Executes on button press in goBut.
function goBut_Callback(hObject, eventdata, handles)
% hObject    handle to goBut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject,'Enable','off');
set(handles.FilesBut,'Enable','off');
set(handles.OutputBut,'Enable','off');
set(handles.FilesPath,'Enable','off');
set(handles.OutputPath,'Enable','off');

runTracking3(handles);

set(hObject,'Enable','on');
set(handles.FilesBut,'Enable','on');
set(handles.OutputBut,'Enable','on');
set(handles.FilesPath,'Enable','on');
set(handles.OutputPath,'Enable','on');


function scaleFactor_Callback(hObject, eventdata, handles)
% hObject    handle to scaleFactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scaleFactor as text
%        str2double(get(hObject,'String')) returns contents of scaleFactor as a double
if (str2double(get(hObject,'String')) < 0 || str2double(get(hObject,'String')) >1)
    set(hObject,'String',0.9);
end
set(handles.reAlign,'String','True');

% --- Executes during object creation, after setting all properties.
function scaleFactor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scaleFactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minimumMovement_Callback(hObject, eventdata, handles)
% hObject    handle to minimumMovement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minimumMovement as text
%        str2double(get(hObject,'String')) returns contents of minimumMovement as a double
if (str2double(get(hObject,'String')) < 0 || str2double(get(hObject,'String')) >7)
    set(hObject,'String',1);
end

% --- Executes during object creation, after setting all properties.
function minimumMovement_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minimumMovement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function fileList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fileList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in listbox3.
function listbox3_Callback(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox3


% --- Executes during object creation, after setting all properties.
function listbox3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in TextInfo.
function TextInfo_Callback(hObject, eventdata, handles)
% hObject    handle to TextInfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns TextInfo contents as cell array
%        contents{get(hObject,'Value')} returns selected item from TextInfo


% --- Executes during object creation, after setting all properties.
function TextInfo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TextInfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function status_Callback(hObject, eventdata, handles)
% hObject    handle to status (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of status as text
%        str2double(get(hObject,'String')) returns contents of status as a double


% --- Executes during object creation, after setting all properties.
function status_CreateFcn(hObject, eventdata, handles)
% hObject    handle to status (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function alignFreq_Callback(hObject, eventdata, handles)
% hObject    handle to alignFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alignFreq as text
%        str2double(get(hObject,'String')) returns contents of alignFreq as a double
if (str2double(get(hObject,'String')) < 0 || str2double(get(hObject,'String')) >100)
    set(hObject,'String',10);
end

% --- Executes during object creation, after setting all properties.
function alignFreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alignFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function WatchWellFig_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WatchWellFig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate WatchWellFig
%hObject.Ylim = [0 70];
%hObject.Color = [179 179 179];
text(40,60,'No well to watch');


% --- Executes on button press in wellUp.
function wellUp_Callback(hObject, eventdata, handles)
% hObject    handle to wellUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.watchWell,'String',int2str(uint8(str2double(get(handles.watchWell,'String'))+1)));

% --- Executes on button press in wellDown.
function wellDown_Callback(hObject, eventdata, handles)
% hObject    handle to wellDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.watchWell,'String',int2str(uint8(str2double(get(handles.watchWell,'String'))-1)));


% --- Executes on button press in wTrack.
function wTrack_Callback(hObject, eventdata, handles)
% hObject    handle to wTrack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.wQuant,'Value',0);
% Hint: get(hObject,'Value') returns toggle state of wTrack


% --- Executes on button press in wQuant.
function wQuant_Callback(hObject, eventdata, handles)
% hObject    handle to wQuant (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.wTrack,'Value',0);



function quantingThresh_Callback(hObject, eventdata, handles)
% hObject    handle to quantingThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of quantingThresh as text
%        str2double(get(hObject,'String')) returns contents of quantingThresh as a double
if (str2double(get(hObject,'String')) < 0)
    set(hObject,'String',30);
end


% --- Executes during object creation, after setting all properties.
function quantingThresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to quantingThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
