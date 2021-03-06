function varargout = AdvancedSettingsGUI(varargin)
% ADVANCEDSETTINGSGUI MATLAB code for AdvancedSettingsGUI.fig
%      ADVANCEDSETTINGSGUI, by itself, creates a new ADVANCEDSETTINGSGUI or raises the existing
%      singleton*.
%
%      H = ADVANCEDSETTINGSGUI returns the handle to a new ADVANCEDSETTINGSGUI or the handle to
%      the existing singleton*.
%
%      ADVANCEDSETTINGSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ADVANCEDSETTINGSGUI.M with the given input arguments.
%
%      ADVANCEDSETTINGSGUI('Property','Value',...) creates a new ADVANCEDSETTINGSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AdvancedSettingsGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AdvancedSettingsGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AdvancedSettingsGUI

% Last Modified by GUIDE v2.5 08-Sep-2015 16:21:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AdvancedSettingsGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @AdvancedSettingsGUI_OutputFcn, ...
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


% --- Executes just before AdvancedSettingsGUI is made visible.
function AdvancedSettingsGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AdvancedSettingsGUI (see VARARGIN)

% Choose default command line output for AdvancedSettingsGUI
handles.output = hObject;

%Accept Settings from MainGUI or Load Settings.mat
if isempty(varargin)
    stemp=load('Settings.mat');
    Settings = stemp.Settings;
    clear stemp
else
    Settings = varargin{1};
end
handles.PrevSettings = Settings;

%HROIM Method
HROIMMethodList = {'Simulated-Kinematic','Simulated-Dynamic','Real-Grain Ref','Real-Single Ref'};
set(handles.HROIMMethod, 'String', HROIMMethodList);
if strcmp(Settings.HROIMMethod,'Simulated')
    SetPopupValue(handles.HROIMMethod,'Simulated-Kinematic');
elseif strcmp(Settings.HROIMMethod,'Dynamic Simulated')
    SetPopupValue(handles.HROIMMethod,'Simulated-Dynamic');
else
    if Settings.RefImageInd == 0
        SetPopupValue(handles.HROIMMethod,'Real-Grain Ref');
    else
        SetPopupValue(handles.HROIMMethod,'Real-Single Ref');
    end
end

%Standard Deviation
set(handles.StandardDeviation,'String',num2str(Settings.StandardDeviation));
%Misorientation Tolerance
set(handles.MisoTol,'String',num2str(Settings.MisoTol));
%Grain Ref Type
GrainRefImageTypeList = {'Min Kernel Avg Miso','IQ > Fit > CI'};
set(handles.GrainRefType, 'String', GrainRefImageTypeList);
SetPopupValue(handles.GrainRefType,Settings.GrainRefImageType);
%Calculate Dislocation Density
set(handles.DoDD,'Value', Settings.CalcDerivatives);
%Do Split DD
set(handles.DoSplitDD,'Value',Settings.DoDDS);
%Number of Skip Points
set(handles.SkipPoints,'String',Settings.NumSkipPts);
%IQ Cutoff
set(handles.IQCutoff,'String',num2str(Settings.IQCutoff));
%SplitDD Method
DDSList = {'Nye-Kroner', 'Nye-Kroner (Pantleon)','Distortion Matching'};
set(handles.SplitDDMethod,'String',DDSList);
SetPopupValue(handles.SplitDDMethod,Settings.DDSMethod);
%Kernel Avg Miso
if iscell(Settings.KernelAvgMisoPath)
    Settings.KernelAvgMisoPath = Settings.KernelAvgMisoPath{1};
end
if exist(Settings.KernelAvgMisoPath,'file')
    [path,name,ext] = fileparts(Settings.KernelAvgMisoPath);
    set(handles.KAMname,'String',[name ext]);
    set(handles.KAMpath,'String',path);
else
    set(handles.KAMname,'String','No File Selected');
    set(handles.KAMpath,'String','No File Selected');
end
%Calculation Options
set(handles.EnableProfiler,'Value',Settings.EnableProfiler);

%Set Position
if length(varargin) > 1
    MainSize = varargin{2};
    set(hObject,'Units','pixels');
    GUIsize = get(hObject,'Position');
    set(hObject,'Position',[MainSize(1)-GUIsize(3)-20 MainSize(2) GUIsize(3) GUIsize(4)]);
    movegui(hObject,'onscreen');
end

% Update handles structure
handles.Settings = Settings;

%Update Components
HROIMMethod_Callback(handles.HROIMMethod,eventdata,handles);
handles = guidata(hObject);
DoDD_Callback(handles.DoDD, eventdata, handles);
handles = guidata(hObject);

guidata(hObject, handles);

% UIWAIT makes AdvancedSettingsGUI wait for user response (see UIRESUME)
uiwait(handles.AdvancedSettingsGUI);


% --- Outputs from this function are returned to the command line.
function varargout = AdvancedSettingsGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.Settings;
delete(handles.AdvancedSettingsGUI);

% --- Executes when user attempts to close AdvancedSettingsGUI.
function AdvancedSettingsGUI_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to AdvancedSettingsGUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if strcmp(get(hObject,'waitstatus'),'waiting')
    uiresume(hObject);
else
    delete(hObject);
end

% --- Executes on button press in SaveButton.
function SaveButton_Callback(hObject, eventdata, handles)
% hObject    handle to SaveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
AdvancedSettingsGUI_CloseRequestFcn(handles.AdvancedSettingsGUI, eventdata, handles);

% --- Executes on button press in CancelButton.
function CancelButton_Callback(hObject, eventdata, handles)
% hObject    handle to CancelButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Settings = handles.PrevSettings;
guidata(hObject,handles);
AdvancedSettingsGUI_CloseRequestFcn(handles.AdvancedSettingsGUI, eventdata, handles);


% --- Executes on selection change in HROIMMethod.
function HROIMMethod_Callback(hObject, eventdata, handles)
% hObject    handle to HROIMMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns HROIMMethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from HROIMMethod
contents = cellstr(get(hObject,'String'));
HROIMMethod = contents{get(hObject,'Value')};
switch HROIMMethod
    case 'Simulated-Kinematic'
        set(handles.HROIMlabel,'String','Iteration Limit');
        set(handles.HROIMedit,'String',num2str(handles.Settings.IterationLimit));
        set(handles.HROIMedit,'Enable','on');
        set(handles.GrainRefType,'Enable','off');
        set(handles.SelectKAM,'Enable','off');
        handles.Settings.HROIMMethod = 'Simulated';
    case 'Simulated-Dynamic'
        %Check for EMsoft
        EMsoftPath = GetEMsoftPath;
        if isempty(EMsoftPath)
            HROIMMethod = 'Simulated';
            SetPopupValue(hObject,'Simulated-Kinematic');
            handles.Settings.HROIMMethod = HROIMMethod;
        end
        
        %Validate EMsoft setup
        if ~isempty(EMsoftPath)  && strcmp(HROIMMethod,'Simulated-Dynamic')
            %Check if Monte-Carlo Simulation data has been generated
            valid = 1;
            if isfield(handles.Settings,'Phase')
                mats = unique(handles.Settings.Phase);
            elseif strcmp(handles.Settings.Material,'Auto-detect')
                valid = 0;
                warndlgpause({'Material must be specified before selecting Simulation-Dynamic','Resetting to kinematic simulation'},'Select Material');
                SetPopupValue(hObject,'Simulated-Kinematic');
                handles.Settings.HROIMMethod = 'Simulated';
                mats = {};
            else
                mats = handles.Settings.Material;
            end
            EMdataPath = fullfile(fileparts(EMsoftPath),'EMdata');
            EMsoftMats = dir(EMdataPath);
            EMsoftMats = {EMsoftMats(~cellfun(@isempty,strfind({EMsoftMats.name},'EBSDmaster'))).name}';
            EMsoftMats = cellfun(@(x) x(1:strfind(x,'_EBSDmaster')-1),EMsoftMats,'UniformOutput',false);
            inlist = ismember(mats,EMsoftMats);
            if ~all(inlist)
                valid = 0;
                msg = {['No master EBSD files for: ' strjoin(mats(~inlist),', ')], ['Search path: ' EMdataPath],'Resetting to kinematic simulation'};
                warndlgpause(msg,'No EMsoft data found');
                SetPopupValue(hObject,'Simulated-Kinematic');
                handles.Settings.HROIMMethod = 'Simulated';
            end
            
            %Set method to simualated dynamic
            if valid
                set(handles.HROIMlabel,'String','Iteration Limit');
                set(handles.HROIMedit,'String',num2str(handles.Settings.IterationLimit));
                set(handles.HROIMedit,'Enable','off');
                set(handles.GrainRefType,'Enable','off');
                set(handles.SelectKAM,'Enable','off');
                handles.Settings.HROIMMethod = 'Dynamic Simulated';
            end
        end
    case 'Real-Grain Ref'
        set(handles.HROIMlabel,'String','Ref Image Index');
        handles.Settings.RefImageInd = 0;
        set(handles.HROIMedit,'String',num2str(handles.Settings.RefImageInd));
        set(handles.HROIMedit,'Enable','off');
        set(handles.GrainRefType,'Enable','on');
        GrainRefType_Callback(handles.GrainRefType, eventdata, handles);
        handles.Settings.HROIMMethod = 'Real';
    case 'Real-Single Ref'
        set(handles.HROIMlabel,'String','Ref Image Index');
        if handles.Settings.RefImageInd == 0
            handles.Settings.RefImageInd = 1;
        end
        set(handles.HROIMedit,'String',num2str(handles.Settings.RefImageInd));
        set(handles.HROIMedit,'Enable','on');
        set(handles.GrainRefType,'Enable','on');
        GrainRefType_Callback(handles.GrainRefType, eventdata, handles);
        handles.Settings.HROIMMethod = 'Real';
end
guidata(hObject,handles);

function warndlgpause(msg,title)
h = warndlg(msg,title);
uiwait(h,7);
if isvalid(h); close(h); end;


% --- Executes during object creation, after setting all properties.
function HROIMMethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HROIMMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function HROIMedit_Callback(hObject, eventdata, handles)
% hObject    handle to HROIMedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of HROIMedit as text
%        str2double(get(hObject,'String')) returns contents of HROIMedit as a double
contents = cellstr(get(handles.HROIMMethod,'String'));
HROIMMethod = contents{get(handles.HROIMMethod,'Value')};
switch HROIMMethod
    case 'Simulated-Kinematic'
        handles.Settings.IterationLimit = str2double(get(hObject,'String'));
    case 'Simulated-Dynamic'
        handles.Settings.IterationLimit = str2double(get(hObject,'String'));
    case 'Real-Grain Ref'
        handles.Settings.RefImageInd = 0;
    case 'Real-Single Ref'
        input = str2double(get(hObject,'String'));
        if isfield(handles.Settings,'ScanLength') && ...
                input > 0 && input <= handles.Settings.ScanLength
            handles.Settings.RefImageInd = input;
        else
            msgbox(['Invalid input. Must be between 1 and ' num2str(handles.Settings.ScanLength) '.'],'Invalid Image Index');
            set(hObject,'String',num2str(handles.Settings.RefImageInd));
        end
end
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function HROIMedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HROIMedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function StandardDeviation_Callback(hObject, eventdata, handles)
% hObject    handle to StandardDeviation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of StandardDeviation as text
%        str2double(get(hObject,'String')) returns contents of StandardDeviation as a double
handles.Settings.StandardDeviation = str2double(get(hObject,'String'));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function StandardDeviation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StandardDeviation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MisoTol_Callback(hObject, eventdata, handles)
% hObject    handle to MisoTol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MisoTol as text
%        str2double(get(hObject,'String')) returns contents of MisoTol as a double
handles.Settings.MisoTol = str2double(get(hObject,'String'));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function MisoTol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MisoTol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in GrainRefType.
function GrainRefType_Callback(hObject, eventdata, handles)
% hObject    handle to GrainRefType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns GrainRefType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from GrainRefType
contents = cellstr(get(hObject,'String'));
GrainRefType = contents{get(hObject,'Value')};
switch GrainRefType
    case 'Min Kernel Avg Miso'
        set(handles.SelectKAM,'Enable','on');
    case 'IQ > Fit > CI'
        set(handles.SelectKAM,'Enable','off');
end
handles.Settings.GrainRefImageType = GrainRefType;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function GrainRefType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GrainRefType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in DoDD.
function DoDD_Callback(hObject, eventdata, handles)
% hObject    handle to DoDD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DoDD
handles.Settings.CalcDerivatives = get(hObject,'Value');
if get(hObject,'Value')
    set(handles.DoSplitDD,'Enable','on');
    set(handles.SkipPoints,'Enable','on');
    set(handles.IQCutoff,'Enable','on');
    DoSplitDD_Callback(handles.DoSplitDD, eventdata, handles);
    handles = guidata(hObject);
else
    set(handles.DoSplitDD,'Enable','off');
    set(handles.SkipPoints,'Enable','off');
    set(handles.IQCutoff,'Enable','off');
    set(handles.SplitDDMethod,'Enable','off');
end
guidata(hObject,handles);


% --- Executes on button press in DoSplitDD.
function DoSplitDD_Callback(hObject, eventdata, handles)
% hObject    handle to DoSplitDD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DoSplitDD
Settings = handles.Settings;
valid = 0;
j = 1;
allMaterials = unique(Settings.Phase);
for i = 1:length(allMaterials)
    M = ReadMaterial(allMaterials{i});
    if isfield(M,'SplitDD')
        valid = 1;
    else
        valid = 0;
        invalidInd(j) = i;
        j = j + 1;
    end
end
if get(hObject,'Value')
    if valid
        set(handles.SplitDDMethod,'Enable','on');
    else
        warndlg(['Split Dislocation data not available for ' allMaterials{invalidInd(1)}],'OpenXY');
        set(hObject,'Value',0);
        set(handles.SplitDDMethod, 'Enable', 'off');
    end
else
    set(hObject,'Value',0);
    set(handles.SplitDDMethod,'Enable','off');
end
handles.Settings.DoDDS = get(hObject,'Value');
guidata(hObject,handles);


function SkipPoints_Callback(hObject, eventdata, handles)
% hObject    handle to SkipPoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SkipPoints as text
%        str2double(get(hObject,'String')) returns contents of SkipPoints as a double
%Input validation
Settings = handles.Settings;
UserInput = get(hObject,'String');
if isempty(str2num(UserInput))
    if ~strcmp(UserInput, 'a') && ~strcmp(UserInput, 't')
        set(hObject, 'String', Settings.NumSkipPts);
        warndlg('Input must be numerical, "a" or "t"');
    end
elseif str2double(UserInput) < 0
    set(hObject, 'String', Settings.NumSkipPts);
    warndlg('Input must be positive');
else
    set(hObject, 'String', round(str2double(UserInput)));
end
handles.Settings.NumSkipPts = get(hObject,'String');

%Updates handles object
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function SkipPoints_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SkipPoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function IQCutoff_Callback(hObject, eventdata, handles)
% hObject    handle to IQCutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of IQCutoff as text
%        str2double(get(hObject,'String')) returns contents of IQCutoff as a double
handles.Settings.IQCutoff = str2double(get(hObject,'String'));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function IQCutoff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IQCutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in SplitDDMethod.
function SplitDDMethod_Callback(hObject, eventdata, handles)
% hObject    handle to SplitDDMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns SplitDDMethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SplitDDMethod
contents = cellstr(get(hObject,'String'));
handles.Settings.DDSMethod = contents{get(hObject,'Value')};
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function SplitDDMethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SplitDDMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SelectKAM.
function SelectKAM_Callback(hObject, eventdata, handles)
% hObject    handle to SelectKAM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
w = pwd;
if exist(handles.Settings.ScanFilePath,'file')
    path = fileparts(handles.Settings.ScanFilePath);
else
    path = pwd;
end
cd(path);
[name, path] = uigetfile('*.txt','OIM Map Data');
set(handles.KAMname,'String',name);
set(handles.KAMname,'TooltipString',name);
set(handles.KAMpath,'String',path);
set(handles.KAMpath,'TooltipString',path);
handles.Settings.KernelAvgMisoPath = fullfile(path,name);
guidata(hObject,handles);
cd(w);


function SetPopupValue(Popup,String)
String = num2str(String);    
List = get(Popup,'String');
IndList = 1:length(List);
Value = IndList(strcmp(List,String));
if isempty(Value); Value =1; end;
set(Popup, 'Value', Value);



% --- Executes on button press in EnableProfiler.
function EnableProfiler_Callback(hObject, eventdata, handles)
% hObject    handle to EnableProfiler (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of EnableProfiler
handles.Settings.EnableProfiler = get(hObject,'Value');


