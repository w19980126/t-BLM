function varargout = ViewAmp(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ViewAmp_OpeningFcn, ...
    'gui_OutputFcn',  @ViewAmp_OutputFcn, ...
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


% --- Executes just before ViewAmp is made visible.
function ViewAmp_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = ViewAmp_OutputFcn(hObject, eventdata, handles)

varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function axes3_CreateFcn(hObject, eventdata, handles)

guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function axes4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes4
guidata(hObject, handles);


% --- Executes on slider movement.
function ImageSlider_Callback(hObject, eventdata, handles)
% hObject    handle to ImageSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% Figure out which image to show
index = ceil(get(hObject, 'Value'));

% Update existing image object in the GUI using this image data
set(handles.image, 'CData', handles.intensity{index});


% --- Executes during object creation, after setting all properties.
function ImageSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImageSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function TimeInterval_Callback(hObject, eventdata, handles)
% hObject    handle to TimeInterval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TimeInterval as text
%        str2double(get(hObject,'String')) returns contents of TimeInterval as a double

% get 'Time interval (s)'
TimeInterval = str2double(get(hObject, 'String'));
handles.TimeInterval = TimeInterval;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function TimeInterval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TimeInterval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SampleRate_Callback(hObject, eventdata, handles)
% hObject    handle to SampleRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SampleRate as text
%        str2double(get(hObject,'String')) returns contents of SampleRate as a double

% get 'Sample rate (fps)'
Fs = str2double(get(hObject, 'String'));
handles.Fs = Fs;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function SampleRate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SampleRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in EquipNiUSB.
function EquipNiUSB_Callback(hObject, eventdata, handles)
% hObject    handle to EquipNiUSB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%---------------------------Ni 采集设置-----------------------------

s = daq.createSession('ni');
[ch1, idx1] = addAnalogInputChannel(s, 'Dev1', 0, 'Voltage');
[ch2, idx2] = addAnalogInputChannel(s, 'Dev1', 1, 'Voltage');
[ch3, idx3] = addAnalogInputChannel(s, 'Dev1', 7, 'Voltage');
handles.s = s;
guidata(hObject, handles);

warndlg('Ni-Daq is ready!')

% --- Executes on button press in StartTimer.
function StartTimer_Callback(hObject, eventdata, handles)
% hObject    handle to StartTimer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Start timer
s = handles.s;
s.Rate = 1000;
s.DurationInSeconds = handles.NiDaqTime;
% lh = addlistener(s,'DataAvailable', @(src,event) plot(event.TimeStamps, event.Data));
data = startForeground(s);
TimeStamps = (1:size(data, 1))/s.Rate;
% 1秒1000个数据，除以Rate后曲线横轴单位便是秒
guidata(hObject, handles);
% delete(lh);
figure('color', 'w');
plot(TimeStamps, data);

[folder_structure, current_folder] = fileparts(handles.DirectoryName);
if length(data) > 0
    mkdir([folder_structure '\MAT']);
    mkdir([folder_structure '\MAT\' current_folder]);
end

SavePath = [folder_structure '\MAT\' handles.expName '_timer.mat'];
save(SavePath, 'data', 'TimeStamps');



function AcqTime_Callback(hObject, eventdata, handles)
% hObject    handle to AcqTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AcqTime as text
%        str2double(get(hObject,'String')) returns contents of AcqTime as a double

% get 'Ni Acquisition time (s)'
NiDaqTime = str2double(get(hObject, 'String'));
handles.NiDaqTime = NiDaqTime;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function AcqTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AcqTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PicsPath.
function PicsPath_Callback(hObject, eventdata, handles)
% hObject    handle to PicsPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

DirectoryName = uigetdir();
handles.DirectoryName = DirectoryName;
handles.roiNumber = 0;daq.getVendors
guidata(hObject,handles);


function PreviewNum_Callback(hObject, eventdata, handles)
% hObject    handle to PreviewNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PreviewNum as text
%        str2double(get(hObject,'String')) returns contents of PreviewNum as a double

% Preview Number
PreviewNum = str2double(get(hObject, 'String'));
handles.PreviewNum = PreviewNum;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function PreviewNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PreviewNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Load.
function Load_Callback(hObject, ~, handles)
% hObject    handle to Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Load Images

BeginPoint = 1;
if length(dir([handles.DirectoryName '\*.raw'])) == 0
    % *.tiff images
    PreviewList = dir([handles.DirectoryName '\*.tiff']);
    if length(PreviewList) == 0
        PreviewList = dir([handles.DirectoryName '\*.tif']);
    end
    
    % The first image 'intensity0'
    intensity0 = GetTiffIntensity(handles.DirectoryName, PreviewList, 1);
    
    % Check if amount of images is enough for processing
    FileList = dir([handles.DirectoryName '\*.tiff']);
    if length(FileList) == 0
        FileList = dir([handles.DirectoryName '\*.tif']);
    end
    if length(FileList) <= handles.PreviewNum
        pause(handles.TimeInterval);
    end
    
    % All images' intensity
    intensity = cell(handles.PreviewNum + 1, 1);
    for file = BeginPoint:(handles.PreviewNum + BeginPoint)
        temp = GetTiffIntensity(handles.DirectoryName, PreviewList, file);
        intensity{file-BeginPoint+1, 1} = temp - intensity0;
    end
    
else
    % *.raw images
    PreviewList = dir([handles.DirectoryName '\*.raw']);
    
    % The first image 'intensity0'
    intensity0 = GetRawIntensity(handles.DirectoryName, PreviewList, 1);
    
    % Check if amount of images is enough for processing
    FileList = dir([handles.DirectoryName '\*.raw']);
    if length(FileList) <= handles.PreviewNum
        pause(handles.TimeInterval);
    end
    
    % All images' intensity
    intensity = cell(handles.PreviewNum + 1, 1);
    for file = BeginPoint:(handles.PreviewNum + BeginPoint)
        temp = GetRawIntensity(handles.DirectoryName, PreviewList, file);
        intensity{file-BeginPoint+1, 1} = temp - intensity0;
    end
    
end
handles.intensity = intensity;
handles.intensity0 = intensity0;

% Display the first one and store the graphics handle to the imshow object
axes(handles.axes1);
handles.image = imshow(handles.intensity{1}, 'Parent', handles.axes1);
handles.hZoom = [];

set(handles.ImageSlider, 'Min', 1, 'Max', (handles.PreviewNum + 1), ...
    'SliderStep', [1 1]/(handles.PreviewNum+1 - 1), 'Value', 1)

handles.BeginPoint = BeginPoint+handles.PreviewNum;
guidata(hObject, handles);

% --- Get Raw images Intensity
function RawIntensity = GetRawIntensity(DirectoryName, DirectoryFileList, num)
eight_bit = 0; % select 8 bit  or 16 bit raw files; default is 16 bit

fid = fopen([DirectoryName '\' DirectoryFileList(num).name]);

A = fread(fid, 'uint8=>uint8');
fclose(fid);
% allign bits
E = double(A(1:2:end));
F = double(A(2:2:end));
G = 64*E+F/4;
if eight_bit == 1
    % Low speed [640 480] High speed [320 240]
    RawIntensity = reshape(E, [640 480]); % ======================
elseif eight_bit == 0
    RawIntensity = reshape(G, [640 480])';% ======================
end

% --- Get Tiff images Intensity
function TiffIntensity = GetTiffIntensity(DirectoryName, DirectoryFileList, num)
TiffIntensity = double(imread([DirectoryName '\' DirectoryFileList(num).name]));


% --- Executes on button press in LoadMore.
function LoadMore_Callback(hObject, eventdata, handles)
% hObject    handle to LoadMore (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

BeginPoint = handles.BeginPoint;
intensity0 = handles.intensity0;

if length(dir([handles.DirectoryName '\*.raw'])) == 0
    % *.tiff images
    PreviewList = dir([handles.DirectoryName '\*.tiff']);
    if length(PreviewList) == 0
        PreviewList = dir([handles.DirectoryName '\*.tif']);
    end
    
    % All images' intensity
    intensity = cell(handles.PreviewNum + 1, 1);
    for file = BeginPoint:(handles.PreviewNum + BeginPoint)
        temp = GetTiffIntensity(handles.DirectoryName, PreviewList, file);
        intensity{file-BeginPoint+1, 1} = temp - intensity0;
    end
    
else
    % *.raw images
    PreviewList = dir([handles.DirectoryName '\*.raw']);
    
    % All images' intensity
    intensity = cell(handles.PreviewNum + 1, 1);
    for file = BeginPoint:(handles.PreviewNum + BeginPoint)
        temp = GetRawIntensity(handles.DirectoryName, PreviewList, file);
        intensity{file-BeginPoint+1, 1} = temp - intensity0;
    end
    
end
handles.intensity = intensity;


% Display the first one and store the graphics handle to the imshow object
handles.image = imshow(handles.intensity{1}, 'Parent', handles.axes1);

% Update the slider to accomodate all of the images
set(handles.ImageSlider, 'Min', 1, 'Max', (handles.PreviewNum + 1), ...
    'SliderStep', [1 1]/(handles.PreviewNum+1 - 1), 'Value', 1)

handles.BeginPoint = BeginPoint+handles.PreviewNum;
guidata(hObject, handles);


function ACfrequency_Callback(hObject, eventdata, handles)
% hObject    handle to ACfrequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ACfrequency as text
%        str2double(get(hObject,'String')) returns contents of ACfrequency as a double

% get 'AC frequency (Hz)'
ACfrequency = str2double(get(hObject, 'String'));
handles.ACfrequency = ACfrequency;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function ACfrequency_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ACfrequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function expName_Callback(hObject, eventdata, handles)
% hObject    handle to expName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of expName as text
%        str2double(get(hObject,'String')) returns contents of expName as a double

% get 'Experiment No.'
expName = get(hObject, 'String');
handles.expName = expName;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function expName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to expName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ACtime_Callback(hObject, eventdata, handles)
% hObject    handle to ACtime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ACtime as text
%        str2double(get(hObject,'String')) returns contents of ACtime as a double


% get 'AC time (s)'
ACtime = str2double(get(hObject, 'String'));
handles.ACtime = ACtime;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function ACtime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ACtime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CreateROI.
function CreateROI_Callback(hObject, eventdata, handles)
% hObject    handle to CreateROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Create ROI
axes(handles.axes1);
roi = imrect;
mask = createMask(roi);
c = length(find(mask(:)~=0));
handles.mask = mask;
guidata(hObject, handles);

TimeInterval = handles.TimeInterval;
intensity0 = handles.intensity0;
Fs = handles.Fs;


BeginPoint = 1;
IntensityOfROI = zeros(BeginPoint+Fs*TimeInterval, 1);
if length(dir([handles.DirectoryName '\*.raw'])) == 0
    % *.tiff images
    FileList = dir([handles.DirectoryName '\*.tiff']);
    if length(FileList) == 0
        FileList = dir([handles.DirectoryName '\*.tif']);
    end
    Amp = zeros(fix(length(FileList)/Fs), 1);
    
    if length(FileList) <= Fs*TimeInterval
        pause(TimeInterval);
    end
    
    
    while BeginPoint < length(FileList)
        EndPoint = BeginPoint+Fs*TimeInterval;
        
        if EndPoint > length(FileList)
            EndPoint = length(FileList);
        end
        % All images' intensity
        for file = BeginPoint:EndPoint
            intensity = GetTiffIntensity(handles.DirectoryName, FileList, file);
            
            % Intensity of ROI (axes2)
            temp = (intensity-intensity0).*mask;
            IntensityOfROI(file, 1) = sum(temp(:))/c;
        end
        
        % Intensity of ROI (axes2)
        axes(handles.axes2);
        Lm = length(IntensityOfROI);
        handles.roiIntensity_Plot = plot((1:Lm)', IntensityOfROI);
        xlim([0 Lm])
        xlabel('Frames in interval')
        ylabel('Intensity (a.u.)')
        
        % Amplitude spectrum (axes3)
        Y = fft(IntensityOfROI(BeginPoint:EndPoint));
        L = EndPoint - BeginPoint + 1;
        P2 = abs(Y/L);
        % Amp(num, 1) = max(2*P2(2:end-1));
        P1 = P2(1:ceil(L/2)+1);
        P1(2:end-1) = 2*P1(2:end-1);
        f = Fs*(0:ceil(L/2))/L;
        axes(handles.axes3);
        plot(f, P1)
        title('Single-Sided Amplitude Spectrum of X(t)')
        xlim([1 fix(Fs/2)])
        xlabel('f (Hz)')
        ylabel('|P1(f)|')
        
        % Average Intensity vs time(axes4)
        axes(handles.axes4);
        num = fix(EndPoint/Fs);
        % Amp(num, 1) = 2*std(IntensityOfROI(BeginPoint:EndPoint)); % a wrong data
        Amp(num, 1) = 1.414*max(P1(3:end-1));
        TestTime = (1:length(Amp))';
        plot(TestTime, Amp, '.');
        xlim([0 handles.ACtime])
        xlabel('t (s)')
        ylabel('Oscillation intensity (a.u.)')
        
        BeginPoint = BeginPoint+Fs*TimeInterval;
        FileList = dir([handles.DirectoryName '\*.tiff']);
        if length(FileList) == 0
            FileList = dir([handles.DirectoryName '\*.tif']);
        end
        
        pause(0.05)
    end
else
    % *.raw images
    FileList = dir([handles.DirectoryName '\*.raw']);
    Amp = zeros(fix(length(FileList)/Fs), 1);
    
    if length(FileList) <= Fs*TimeInterval
        pause(TimeInterval);
    end
    
    IntensityOfROI = zeros(BeginPoint+Fs*TimeInterval, 1);
    while BeginPoint < length(FileList)
        EndPoint = BeginPoint+Fs*TimeInterval;
        
        if EndPoint > length(FileList)
            EndPoint = length(FileList);
        end
        
        for file = BeginPoint:EndPoint
            intensity = GetRawIntensity(handles.DirectoryName, FileList, file);
            
            % Intensity of ROI (axes2)
            temp = (intensity-intensity0).*mask;
            IntensityOfROI(file, 1) = sum(temp(:))/c;
            
        end
        
        % Intensity of ROI (axes2)
        axes(handles.axes2);
        Lm = length(IntensityOfROI);
        handles.roiIntensity_Plot = plot((1:Lm)', IntensityOfROI);
        xlim([0 Lm])
        xlabel('Frames in interval')
        ylabel('Intensity (a.u.)')
        
        % Amplitude spectrum (axes3)
        Y = fft(IntensityOfROI(BeginPoint:EndPoint));
        L = EndPoint - BeginPoint + 1;
        P2 = abs(Y/L);
        % Amp(num, 1) = max(2*P2(2:end-1));
        P1 = P2(1:ceil(L/2)+1);
        P1(2:end-1) = 2*P1(2:end-1);
        f = Fs*(0:ceil(L/2))/L;
        axes(handles.axes3);
        plot(f, P1)
        title('Single-Sided Amplitude Spectrum of X(t)')
        xlim([1 fix(Fs/2)])
        xlabel('f (Hz)')
        ylabel('|P1(f)|')
        
        % Average Intensity vs time(axes4)
        axes(handles.axes4);
        num = fix(EndPoint/Fs);
        % Amp(num, 1) = 2*std(IntensityOfROI(BeginPoint:EndPoint)); % a wrong data
        Amp(num, 1) = 1.414*max(P1(3:end-1));
        TestTime = (1:length(Amp))';
        plot(TestTime, Amp, '.');
        xlim([0 handles.ACtime])
        xlabel('t (s)')
        ylabel('Oscillation intensity (a.u.)')
        
        BeginPoint = BeginPoint+Fs*TimeInterval;
        FileList = dir([handles.DirectoryName '\*.raw']);
        
        pause(0.02)
    end
end

handles.IntensityOfROI = IntensityOfROI;
handles.Amp = Amp;
handles.TestTime = TestTime;
handles.roiNumber = handles.roiNumber + 1;
guidata(hObject, handles);



% --- Executes on button press in Delete_Save.
function Delete_Save_Callback(hObject, eventdata, handles)
% hObject    handle to Delete_Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% delete ROI and save


IntensityOfROI = handles.IntensityOfROI;
roiNumber = handles.roiNumber;
mask = handles.mask;
Amp = handles.Amp;
TestTime = handles.TestTime;
ACtime = handles.ACtime;
% fitresult = handles.fitresult;
% gof = handles.gof;
expName = handles.expName;

[folder_structure, current_folder] = fileparts(handles.DirectoryName);
if length(TestTime) > 0
    mkdir([folder_structure '\MAT']);
    mkdir([folder_structure '\MAT\' current_folder  ]);
end

h = getframe(gcf);
GUIsaved = [folder_structure '\MAT\' current_folder '\' expName '_roi' num2str(roiNumber) '_GUIsaved.tif'];
imwrite(h.cdata, GUIsaved);

GUIsaved = [folder_structure '\MAT\' current_folder '\' expName '_roi' num2str(roiNumber) '_Mask.tif'];
imwrite(mask, GUIsaved);


SavePath = [folder_structure '\MAT\' current_folder '\' expName '_roi' num2str(roiNumber) '.mat'];
save(SavePath, 'IntensityOfROI', 'Amp', 'TestTime', 'ACtime');
% % save(SavePath, 'IntensityOfROI', 'Amp', 'TestTime', 'ACtime', 'fitresult', 'gof');


% --- Executes on button press in ZoomOn.
function ZoomOn_Callback(hObject, eventdata, handles)
% hObject    handle to ZoomOn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.hZoom)
    handles.hZoom = zoom;
    guidata(hObject, handles);
end
if ~strcmp(get(handles.hZoom,'Enable'), 'on')
    set(handles.hZoom, 'Enable', 'on');
end


% --- Executes on button press in ZoomOff.
function ZoomOff_Callback(hObject, eventdata, handles)
% hObject    handle to ZoomOff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isempty(handles.hZoom)
    set(handles.hZoom, 'Enable', 'off');
end


% --- Executes on button press in Adjust.
function Adjust_Callback(hObject, eventdata, handles)
% hObject    handle to Adjust (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% I = getimage(handles.axes1);
% J = imadjust(I);
% axes(handles.axes1);
% imshow(J, 'DisplayRange',[], 'InitialMagnification', 'fit')

for ii = 1:(handles.PreviewNum+1)
   I = handles.intensity{ii};
   temp = I - ones(size(I))*min(I(:));
   handles.intensity{ii} = uint8(temp/max(temp(:))*255);
end
handles.image = imshow(handles.intensity{1}, 'Parent', handles.axes1);
set(handles.ImageSlider, 'Min', 1, 'Max', (handles.PreviewNum + 1), ...
    'SliderStep', [1 1]/(handles.PreviewNum+1 - 1), 'Value', 1)
guidata(hObject, handles);


% --- Executes on button press in Debug.
function Debug_Callback(hObject, eventdata, handles)
% hObject    handle to Debug (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
keyboard


% --- Executes on button press in Fit_Curve.
function Fit_Curve_Callback(hObject, eventdata, handles)
% hObject    handle to Fit_Curve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Amp = handles.Amp;
TestTime = handles.TestTime;
ACtime = handles.ACtime;

axes(handles.axes4);
[x, ~] = ginput(2);
intensity = Amp(round(x(1)): round(x(2))); % get intensity in the selected region
% qt = 1./intensity; % from intensity to qt
t = TestTime(round(x(1)): round(x(2))); % x is t.

absorptionEqn = '(1+a*x)/(a*b*x)'; % a = k*qe, b = qe
% absorptionEqn = 'ae^x/(1+b*e^x)';
[fitobject, gof] = fit(t,intensity, absorptionEqn);
coeffvals = coeffvalues(fitobject); % get a & b;
qe1 = 1/Amp(round(x(2))); % selected fit
k1 = coeffvals(1)./qe1;
qe2 = coeffvals(2); % direct cfit
k2 = coeffvals(1)./coeffvals(2);
plot(TestTime, Amp, 'b.');
hold on
t1 = TestTime(round(x(1))-10: round(x(2))+10);
intensity1 = (1+coeffvals(1).*t1)./(coeffvals(2).*coeffvals(1).*t1);
plot(t1, intensity1, 'r-')
xlim([0 ACtime])
xlabel('t (s)')
ylabel('Oscillation intensity (a.u.)')
hold off

output = cell(7, 1);
output{1, 1} = ['SSE: ' num2str(gof.sse)];
output{2, 1} = ['R-square: ' num2str(gof.rsquare)];
output{3, 1} = ['DFE: ' num2str(gof.dfe)];
output{4, 1} = ['adjR-square: ' num2str(gof.adjrsquare)];
output{5, 1} = ['RMSE: ' num2str(gof.rmse)];
output{6, 1} = ['if qe = 1/y2, then qe = ' num2str(qe1) ', k = ' num2str(k1) '.'];
output{7, 1} = ['if qe follows the cfit, then qe = ' num2str(qe2) ', k = ' num2str(k2) '.'];
msgbox(output, 'Directly fit curve: Goodness of fit & ''k''');

handles.fitresult = fitobject;
handles.gof = gof;
% handles.qe = qe;
handles.k1 = k1;
handles.k2 = k2;
guidata(hObject, handles);


% --- Executes on button press in Unfit.
function Unfit_Callback(hObject, eventdata, handles)
% hObject    handle to Unfit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes4);
plot(handles.TestTime, handles.Amp, '.');
xlim([0 handles.ACtime])
xlabel('t (s)')
ylabel('Oscillation intensity (a.u.)')


% --- Executes on button press in smooth_curve.
function smooth_curve_Callback(hObject, eventdata, handles)
% hObject    handle to smooth_curve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Amp = handles.Amp;
TestTime = handles.TestTime;
ACtime = handles.ACtime;

axes(handles.axes4);
Amp_smooth = smooth(TestTime, Amp, 0.1,'rloess');
plot(TestTime, Amp, 'b.',TestTime, Amp_smooth,'r-');
xlim([0 ACtime])
xlabel('t (s)')
ylabel('Oscillation intensity (a.u.)')

handles.Amp_smooth = Amp_smooth;
guidata(hObject, handles);


% --- Executes on button press in Smooth_fit.
function Smooth_fit_Callback(hObject, eventdata, handles)
% hObject    handle to Smooth_fit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Amp = handles.Amp;
TestTime = handles.TestTime;
ACtime = handles.ACtime;

axes(handles.axes4);
Amp_smooth = smooth(TestTime, Amp, 0.1,'rloess');
plot(TestTime, Amp, 'b.',TestTime, Amp_smooth,'r-');
xlim([0 ACtime])
xlabel('t (s)')
ylabel('Oscillation intensity (a.u.)')

[x, ~] = ginput(2);
intensity = Amp_smooth(round(x(1)): round(x(2))); % get intensity in the selected region
% qt = 1./intensity; % from intensity to qt
t = TestTime(round(x(1)): round(x(2))); % x is t.

absorptionEqn = '(1+a*x)/(a*b*x)'; % a = k*qe, b = qe
% absorptionEqn = 'ae^x/(1+b*e^x)';
[fitobject, gof] = fit(t,intensity, absorptionEqn);
coeffvals = coeffvalues(fitobject); % get a & b;
qe1 = 1/Amp(round(x(2))); % selected fit
k1 = coeffvals(1)./qe1;
qe2 = coeffvals(2); % direct cfit
k2 = coeffvals(1)./coeffvals(2);
plot(TestTime, Amp, 'b.');
t1 = TestTime(round(x(1))-10: round(x(2))+10);
intensity1 = (1+coeffvals(1).*t1)./(coeffvals(2).*coeffvals(1).*t1);
hold on
plot(t1, intensity1, 'r-')
xlim([0 ACtime])
xlabel('t (s)')
ylabel('Oscillation intensity (a.u.)')
hold off

output = cell(7, 1);
output{1, 1} = ['SSE: ' num2str(gof.sse)];
output{2, 1} = ['R-square: ' num2str(gof.rsquare)];
output{3, 1} = ['DFE: ' num2str(gof.dfe)];
output{4, 1} = ['adjR-square: ' num2str(gof.adjrsquare)];
output{5, 1} = ['RMSE: ' num2str(gof.rmse)];
output{6, 1} = ['if qe = 1/y2, then qe = ' num2str(qe1) ', k = ' num2str(k1) '.'];
output{7, 1} = ['if qe follows the cfit, then qe = ' num2str(qe2) ', k = ' num2str(k2) '.'];
msgbox(output, 'Smoothed curve fit: Goodness of fit & ''k''');

handles.fitresult = fitobject;
handles.gof = gof;
% handles.qe = qe;
handles.k1 = k1;
handles.k2 = k2;
guidata(hObject, handles);
