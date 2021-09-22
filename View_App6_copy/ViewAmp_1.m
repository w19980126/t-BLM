function varargout = ViewAmp_1(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ViewAmp_1_OpeningFcn, ...
                   'gui_OutputFcn',  @ViewAmp_1_OutputFcn, ...
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

% --- Executes just before ViewAmp_1 is made visible.
function ViewAmp_1_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;

guidata(hObject, handles);

function varargout = ViewAmp_1_OutputFcn(hObject, eventdata, handles)

varargout{1} = handles.output;

%------------------------------------------------------------------------------


% -------------------------閺夆晜鐟�?Σ鍛婄▔閿熸枻鎷锋繛鍫ユ涧閼荤喐绋夊鍥╁弨闂侇剚鎸诲﹢渚�?1�?7宕妷褎鏆忛柣銊ュ閸ら亶寮敓鏂ゆ�?------------------
function text6_CreateFcn(hObject, eventdata, handles)

function text6_DeleteFcn(hObject, eventdata, handles)

function text7_CreateFcn(hObject, eventdata, handles)

function axes1_text_CreateFcn(hObject, eventdata, handles)

function axes2_text_CreateFcn(hObject, eventdata, handles)

% -------------------------閻犲鍟弳锝夊炊閸撗冾暬閻庝絻顫夐惁顔芥償閿熸枻鎷�?1�?7--------------------------------
function Adjust_Callback(hObject, eventdata, handles)

hf = figure;
subplot(121)
histogram(get(handles.hf1,'CData'));
title('The histogram of axes1','fontsize',15,'fontweight','bold');
subplot(122)
histogram(get(handles.hf2,'CData'));
title('The histogram of axes2','fontsize',15,'fontweight','bold');

prompt = {'axes1:min','axes1:max','axes2:min','axes2:max'};
dlgtitle = 'caxis';
dims = [1 35];
definput = {'-1000' '1000' '0' '10000'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
handles.color_axis = [str2num(answer{1}) str2num(answer{2}) str2num(answer{3}) str2num(answer{4})];
axes(handles.axes1)
color_axis = handles.color_axis;
caxis([color_axis(1) color_axis(2)]);
axes(handles.axes2)
caxis([color_axis(3) color_axis(4)]);
close(hf)

guidata(hObject,handles);

%---------------------------闁猴拷锟介幆褋浜ｇ紓鍌楁櫅閻拷锟�1�?7-----------------------------------

function Zoom_Off_Callback(hObject, eventdata, handles)

if ~isempty(handles.hZoom)
    set(handles.hZoom, 'Enable', 'off');
end

function Zoom_On_Callback(hObject, eventdata, handles)

if isempty(handles.hZoom)
    handles.hZoom = zoom;
    guidata(hObject, handles);
end
if ~strcmp(get(handles.hZoom,'Enable'), 'on')
    set(handles.hZoom, 'Enable', 'on');
end

%-------------------------------婵犲﹥鍨垫慨鈺呭级閿熸枻鎷�1�?7--------------------------------

function axes1_slider_Callback(hObject, eventdata, handles)

temp1 = floor(get(hObject,'Value'));
set(handles.axes1_text,'String',num2str(temp1));
axes(handles.axes1);
Intensity = handles.Intensity;
Intensity0 = handles.Intensity0;

if get(handles.Sub,'Value') == 0
    temp2 = Intensity(:,:,temp1) + Intensity0;
else
    temp2 = Intensity(:,:,temp1);
end

handles.hf1 = imagesc(temp2);
colormap(handles.map)
axis off

names = fields(handles);
temp_num = zeros(length(names),1);
for ii = 1:length(names)
    if strcmp(names{ii},'color_axis')
        temp_num(ii) = 1;
    else
        temp_num(ii) = 0;
    end
end
if sum(temp_num) == 1
    color_axis = handles.color_axis;
    caxis([color_axis(1) color_axis(2)]);
end

guidata(hObject,handles);

function axes1_slider_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function axes2_slider_Callback(hObject, eventdata, handles)

temp1 = floor(get(hObject,'Value'));
set(handles.axes2_text,'String',num2str(temp1));
axes(handles.axes2);
Intensity_FFT = handles.Intensity_FFT;
temp2 = Intensity_FFT(:,:,temp1);
handles.hf2 = imagesc(temp2);
colormap(handles.map)
axis off

names = fields(handles);
temp_num = zeros(length(names),1);
for ii = 1:length(names)
    if strcmp(names{ii},'color_axis')
        temp_num(ii) = 1;
    else
        temp_num(ii) = 0;
    end
end
if sum(temp_num) == 1
    color_axis = handles.color_axis;
    caxis([color_axis(3) color_axis(4)]);
end

guidata(hObject,handles);

function axes2_slider_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

%-------------------鐎殿噯鎷烽敓鑺ユ叏鐎ｎ喖娅氶梻鍡楁閺嗙喖骞戦鑹板珯閻忓繐鎽糰ta闁告粌鐣篿g濞ｅ洦绻傞悺銊╁礆閻х剾mer濞戞搫鎷烽敓锟斤�?---------------

function Start_Acq_Callback(hObject, eventdata, handles)

s = handles.s;
s.Rate = 10000;     
s.DurationInSeconds = handles.DAQ_Duration;
data = startForeground(s);

t = (1:size(data, 1))/s.Rate;
figure('color', 'w');
hf = plot(t, data);
xlabel('t(s)','fontsize',15,'fontweight','bold');
ylabel('Intensity(a.u.)','fontsize',15,'fontweight','bold');
set(gca,'fontsize',15,'fontweight','bold');
axis tight
Fs = handles.SampleRate;
handles.begin = TriggerTime(data,t,Fs);

[folder_structure, current_folder] = fileparts(handles.Directory_Name);
SavePath = [folder_structure '\Timer'];
if exist(SavePath)==0
    mkdir([folder_structure '\Timer']);
end
save([SavePath '\' handles.Expname '.mat'], 'data', 't','-v7.3');
saveas(gcf,[SavePath '\' handles.Expname '.fig']);

guidata(hObject, handles);

%------------------------TriggerTime----------------------------

function begin = TriggerTime(data,~,Fs)
[~, begin.Pike] = max(diff(data(:,1)));
begin.Cs = FindEdge(data(:,2));
if begin.Pike >= begin.Cs
    begin.Frame = 1;
else
    begin.Frame = floor((begin.Cs - begin.Pike)/10000*Fs);
end

function loc = FindEdge(data)

temp1 = mean(abs(data(1:100)));
loc = 0;
for ii = 1:length(data)-1
    temp2 = abs(data(ii + 1) - data(ii));
    if temp2 > 25*temp1
        loc = ii+1;
        return
    end
end

%------------------------Ni 閻犱礁澧介悿锟斤�?----------------------------

function Ni_Ready_Callback(hObject, eventdata, handles)

s = daq.createSession('ni');
[ch1, idx1] = addAnalogInputChannel(s, 'Dev1', 0, 'Voltage');
[ch2, idx2] = addAnalogInputChannel(s, 'Dev1', 1, 'Voltage');
[ch3, idx3] = addAnalogInputChannel(s, 'Dev1', 7, 'Voltage');
handles.s = s;
guidata(hObject, handles);

warndlg('Ni-Daq is ready!')

%----------------------------闁搞儱澧芥晶鏍ㄦ償韫囨挸鐏欓悹渚灠缁讹拷锟�?1�?7----------------------------

function Raw_Path_Callback(hObject, eventdata, handles)

persistent startDir

if ~ischar(startDir)
    startDir = [];
end

Directory_Name = uigetdir(startDir);
startDir = Directory_Name;
if Directory_Name ~= 0
    handles.Directory_Name = Directory_Name;
end

guidata(hObject,handles);

%-------------------------------AC_Time-------------------------------

function AC_Time_Callback(hObject, eventdata, handles)

AC_Time = str2double(get(hObject,'String'));
handles.AC_Time = AC_Time;
guidata(hObject,handles);

function AC_Time_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%-------------------------------AC_Freq-------------------------------

function AC_Freq_Callback(hObject, eventdata, handles)

AC_Freq = str2double(get(hObject,'String'));
handles.AC_Freq = AC_Freq;

names = fields(handles);
temp_num = zeros(length(names),1);
for ii = 1:length(names)
    if strcmp(names{ii},'Intensity_FFT')
        temp_num(ii) = 1;
    else
        temp_num(ii) = 0;
    end
end
if sum(temp_num) == 1
    temp1 = handles.Time_Res;
    Intensity = handles.Intensity;
    Intensity0 = handles.Intensity0;
    n = temp1;
    Interval = floor(n*(1/handles.AC_Freq*handles.SampleRate));     % Translate 'Interval' Frames to one fig by FFT
    FFT_Frames = floor(handles.PreviewNum/Interval);   % after trans,we have FFT_length Frames new figs,or Amp figs
    handles.FFT_Frames = FFT_Frames;
    handles.Interval = Interval;
    delta_f = handles.SampleRate/Interval;
    Intensity_FFT = zeros(size(Intensity0,1),size(Intensity0,2),FFT_Frames);
    for ii = 1:size(Intensity0,1)
        for jj = 1:size(Intensity0,2)
            for kk = 1:FFT_Frames
                temp1 = Intensity(ii,jj,1 + (kk - 1)*Interval:...
                    kk*Interval);
                temp2 = fft(temp1);
                Intensity_FFT(ii,jj,kk) = 2*abs(temp2(1 + floor(handles.AC_Freq/delta_f)))/Interval;
                Intensity_FFT_complex(ii,jj,kk) = 2*temp2(1 + floor(handles.AC_Freq/delta_f))/Interval;
            end
        end
    end
    handles.Intensity_FFT = Intensity_FFT;
    handles.Intensity_FFT_complex = Intensity_FFT_complex;
    axes(handles.axes2);
    handles.hf2 = imagesc(Intensity_FFT(:,:,1));
    colormap(gray(256));
    axis off
    set(handles.axes2_slider,'Max',FFT_Frames);
    set(handles.axes2_slider,'Value',1);
end

guidata(hObject,handles);

function AC_Freq_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%-----------------------------DAQ_Duration-------------------------------

function DAQ_Duration_Callback(hObject, eventdata, handles)

DAQ_Duration = str2double(get(hObject,'String'));
handles.DAQ_Duration = DAQ_Duration;
guidata(hObject,handles);

function DAQ_Duration_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%------------------------- SampleRate-------------------------------

function Sample_Rate_Callback(hObject, eventdata, handles)

SampleRate = str2num(get(hObject, 'String'));
handles.SampleRate = SampleRate;
guidata(hObject, handles);

function Sample_Rate_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------ROI闁烩晝顭堥崣锟斤�?---------------------------------
function Save_ROI_Callback(hObject, eventdata, handles)

handles.pin = handles.pin + 1;
Directory_Name = handles.Directory_Name;
names = split_(Directory_Name,'\');
temp1 = strcmp(names,'TIFF');
temp2 = join_(names(1:find(temp1 == 1) - 1),'\');

if isempty(temp2{1}) == 1
    temp2 = join_(names(1:end - 1),'\');
end

temp3 = [temp2{1} '\Result\Temp Processing\' handles.Expname];
if exist(temp3) == 0
    mkdir(temp3)
end
ROI_data.ROI_Mean_y = get(findobj(handles.axes3,'type','line'),'YData');
ROI_data.ROI_Mean_x = get(findobj(handles.axes3,'type','line'),'XData');
ROI_data.FFT_y = get(findobj(handles.axes4,'type','line'),'YData');
ROI_data.FFT_x = get(findobj(handles.axes4,'type','line'),'XData');
ROI_data.Amp_y = get(findobj(handles.axes5,'type','line'),'YData');
ROI_data.Amp_x = get(findobj(handles.axes5,'type','line'),'XData');
ROI_data.Mask = handles.Mask;
saveroute = [temp3 '\ROI_' num2str(handles.pin) '.mat'];
save(saveroute,'ROI_data','-v7.3');

guidata(hObject,handles);

function Create_ROI_Callback(hObject, eventdata, handles)

axes(handles.axes1)
ROI = imrect;
Mask = createMask(ROI);
handles.Mask = Mask;
delete(ROI)
Intensity = handles.Intensity;
ROI_Mean = zeros(size(Intensity,3),1);

for ii = 1:length(ROI_Mean)
    ROI_Mean(ii) = sum(sum(Mask.*Intensity(:,:,ii)))/sum(Mask(:));
end
% 濞戞柨顑嗘晶宥嗙閵夈儱甯ラ柣銏㈤ケxes4闁挎稑鏈�?Σ鍛婄▔鏉炴壆鍟婇梺顒�?1�?7鐏濋崢顥窸I_Mean濞戞搩鍘惧▓鎴︽儎鐎涙銈﹂柛鎺戞閸ｇ瘞FT闁告艾姘︾涵锟斤拷闁绘凹鍠栭妵濠冨緞濮橆偆绉靛Λ鐗堝灥閸ㄥ酣鏌岃箛搴♀枏鐎电増顨婇悵顔斤紣閿熸枻鎷�?1�?7
% 濞ｅ洠鍓濇导鍛偖椤愩垹绔剧紓鍌楁櫇濠�1�?7鍛▔�?�ュ棛顏�?1�?7
axes(handles.axes4)
Fs = handles.SampleRate;
temp1 = abs(fft(ROI_Mean))/length(ROI_Mean);
temp2 = temp1(1:floor(length(ROI_Mean)/2) - 1);
temp2(2:end) = 2*temp2(2:end);
f = Fs/length(ROI_Mean)*(0:length(temp2)-1);
handles.hf4 = plot(f,temp2,'k','linewidth',1);
title('ROI\_Mean FFT','fontsize',10,'fontweight','bold');
xlabel('Freq (Hz)')
ylabel('Amplitude (a.u.)');
set(gca,'fontsize',10,'fontweight','bold');
xlim([0 Fs/2])

if get(handles.Sub,'Value') == 0
    ROI_Mean = ROI_Mean + sum(sum(Mask.*handles.Intensity0))/sum(Mask(:));
end
axes(handles.axes3)
handles.hf3 = plot(ROI_Mean,'k','linewidth',1);
title('ROI\_Mean','fontsize',10,'fontweight','bold');
xlabel('Frame')
ylabel('Intensity (a.u.)');
set(gca,'fontsize',10,'fontweight','bold');
xlim([0 length(ROI_Mean)]);

axes(handles.axes5)
Intensity_FFT_complex = handles.Intensity_FFT_complex;
temp3 = zeros(1,size(Intensity_FFT_complex,3));
for ii = 1:length(temp3)
    temp3(ii) = abs(sum(sum(Mask.*Intensity_FFT_complex(:,:,ii))))/sum(Mask(:));
end
Interval = handles.Interval;
t = (1:handles.FFT_Frames)*Interval/Fs;
if length(temp3) == 1
    handles.hf5 = plot(t,temp3,'ko');
else
    handles.hf5 = plot(t,temp3,'k','linewidth',1);
end

title('Amp-Time','fontsize',10,'fontweight','bold');
xlabel('t (s)')
ylabel('Amplitude (a.u.)');
set(gca,'fontsize',10,'fontweight','bold');
xlim([0 max(t)+0.5]);

guidata(hObject,handles);

%-----------------------鐎殿噯鎷烽敓鑺ユ叏鐎ｎ亜顫ｉ弶鐐舵濞存﹢鎮ч敓鏂ゆ嫹------------------------------

function Load_Callback(hObject, eventdata, handles)

Frame = handles.begin.Frame;
handles.hZoom = [];
handles.pin = 0;
handles.map = gray(256);

% 閻犲洩顕ч崣鍡涘炊閸撗冾暬
All_Raw_List = dir([handles.Directory_Name '\*.raw']);
if size(All_Raw_List,1) ~= 0
    
    if length(All_Raw_List) < Frame + handles.PreviewNum
        handles.PreviewNum = length(All_Raw_List) - Frame;
    end
        
    handles.Raw_Num = size(All_Raw_List,1);
    File_Path = handles.Directory_Name;
    Intensity0 = Get_Raw_Intensity(File_Path,All_Raw_List(Frame).name);
    Intensity = zeros(size(Intensity0,1),size(Intensity0,2),handles.PreviewNum);
    % 閻犲洩顕ч崣鍡涘储閻斿娼楅柛銉ュ⒔婢ф牕顕ｉ崫鍕唺闁轰胶澧楀畵锟斤拷
    for ii = 1:handles.PreviewNum
        temp = Get_Raw_Intensity(File_Path,All_Raw_List(Frame + ii - 1).name);
        Intensity(:,:,ii) = temp - Intensity0;
    end
    handles.Intensity0 = Intensity0;
    handles.Intensity = Intensity;
    
    % 閻庝絻顫夐惁鈩冪▔椤忓嫬鍓肩紒杈╁Х閸嬶絾娼诲☉婊庢斀FFT
    n = handles.Time_Res;
    Interval = floor(n*(1/handles.AC_Freq*handles.SampleRate));     % Translate 'Interval' Frames to one fig by FFT
    FFT_Frames = floor(handles.PreviewNum/Interval);   % after trans,we have FFT_length Frames new figs,or Amp figs
    handles.FFT_Frames = FFT_Frames;
    handles.Interval = Interval;
    delta_f = handles.SampleRate/Interval;
    Intensity_FFT = zeros(size(Intensity0,1),size(Intensity0,2),FFT_Frames);
    for ii = 1:size(Intensity0,1)
        for jj = 1:size(Intensity0,2)
            for kk = 1:FFT_Frames
                temp1 = Intensity(ii,jj,1 + (kk - 1)*Interval:...
                    kk*Interval);
                temp2 = fft(temp1);
                Intensity_FFT(ii,jj,kk) = 2*abs(temp2(1 + floor(handles.AC_Freq/delta_f)))/Interval;
                Intensity_FFT_complex(ii,jj,kk) = 2*temp2(1 + floor(handles.AC_Freq/delta_f))/Interval;
            end
        end
    end
    handles.Intensity_FFT = Intensity_FFT;
    handles.Intensity_FFT_complex = Intensity_FFT_complex;
    axes(handles.axes1);
    if get(handles.Sub,'Value') == 1
        handles.hf1 = imagesc(Intensity(:,:,2));
        colormap(gray)
    else 
        handles.hf1 = imagesc(Intensity(:,:,2) + Intensity0);
        colormap(gray)
    end
    axis off
    set(handles.axes1_slider,'Value',2);
    set(handles.axes1_text,'String',2);
    axes(handles.axes2);
    handles.hf2 = imagesc(Intensity_FFT(:,:,1));
    colormap(gray)
    axis off
    set(handles.axes2_slider,'Value',1);
    set(handles.axes1_text,'String',1);

    % 缂備胶鍠嶇粩瀛樻交濞戞粠�?界憸鐗堝笂缁�?挳宕犻弽鐢电閻㈩垰鏈﹢婊堟嚄閽樺妾褏鍋涙慨鐐哄炊閸撗冾暬闁汇劌瀚顔夹掗弬鍨唺
    % 闁搞儱澧芥晶鏍焾閽樺韬�1�?71濞戞柨顑夊Λ鍧�?及閸撗佷粵闁挎稑鐭傛导鈺呭礂閿熸枻鎷峰ù锝呯Т濞存﹢宕撹箛搴ｇ憹闁煎疇濮ゅΟ澶�?矆閾忚鐣遍梻鍌ゅ櫍椤ｏ拷锟�?1�?7
   
else
    
    if length(dir([handles.Directory_Name '\*.tif'])) == 0
        All_Tiff_List = dir([handles.Directory_Name '\*.tiff']);
        handles.Tiff_Num = size(All_Tiff_List,1);
    else
        All_Tiff_List = dir([handles.Directory_Name '\*.tif']);
        handles.Tiff_Num = size(All_Tiff_List,1);
    end
    
    if length(All_Tiff_List) < Frame + handles.PreviewNum
        handles.PreviewNum = length(All_Tiff_List) - Frame;
    end
    
    File_Path = handles.Directory_Name;
    % 閻犲洩顕ч崣鍡涘储閻斿娼楅柛銉ュ⒔婢ф牕顕ｉ崫鍕唺闁轰胶澧楀畵锟斤拷
    Intensity0 = Get_Tiff_Intensity(File_Path,All_Tiff_List(Frame).name);
    Intensity = zeros(size(Intensity0,1),size(Intensity0,2),handles.PreviewNum);
    for ii = 1:size(Intensity,3)
        temp = Get_Tiff_Intensity(File_Path,All_Tiff_List(Frame + ii - 1).name);
        Intensity(:,:,ii) = temp - Intensity0;
    end
    handles.Intensity0 = Intensity0;
    handles.Intensity = Intensity;

    % 閻庝絻顫夐惁鈩冪▔椤忓嫬鍓肩紒杈╁Х閸嬶絾娼诲☉婊庢斀FFT
    n = handles.Time_Res;
    Interval = floor(n*(1/handles.AC_Freq*handles.SampleRate));    % Translate 'Interval' Frames to one fig by FFT
    FFT_Frames = floor(handles.PreviewNum/Interval);   % after trans,we have FFT_length Frames new figs,or Amp figs
    handles.FFT_Frames = FFT_Frames;
    handles.Interval = Interval;
    delta_f = handles.SampleRate/Interval;
    Intensity_FFT = zeros(size(Intensity0,1),size(Intensity0,2),FFT_Frames);
    for ii = 1:size(Intensity0,1)
        for jj = 1:size(Intensity0,2)
            for kk = 1:FFT_Frames
                temp1 = Intensity(ii,jj,1 + (kk - 1)*Interval:...
                    kk*Interval);
                temp2 = fft(temp1);
                Intensity_FFT(ii,jj,kk) = 2*abs(temp2(1 + floor(handles.AC_Freq/delta_f)))/Interval;
                Intensity_FFT_complex(ii,jj,kk) = 2*temp2(1 + floor(handles.AC_Freq/delta_f))/Interval;
            end
        end
    end
    handles.Intensity_FFT = Intensity_FFT;
    handles.Intensity_FFT_complex = Intensity_FFT_complex;
    axes(handles.axes1);
    if get(handles.Sub,'Value') == 1
        handles.hf1 = imagesc(Intensity(:,:,2));
        colormap(gray)
    else 
        handles.hf1 = imagesc(Intensity(:,:,2) + Intensity0);
        colormap(gray)
    end
    axis off
    set(handles.axes1_slider,'Value',2);
    set(handles.axes1_text,'String',2);
    axes(handles.axes2);
    handles.hf2 = imagesc(Intensity_FFT(:,:,1));
    colormap(gray)
    axis off
    set(handles.axes2_slider,'Value',1);
    set(handles.axes1_text,'String',1);

end

guidata(hObject, handles);
    
%-----------------------閻犲洩顕цぐ鍢w閹兼潙绻愰崹顏堝椽鐎涚FF閹兼潙绻愰崹顏勵嚕閸濆嫬顔�1�?7-------------------------------

function Raw_Intensity = Get_Raw_Intensity(File_Path,File_Name)
fid = fopen([File_Path '\' File_Name]);
A = fread(fid,'uint8=>uint8');
fclose(fid);
E = double(A(1:2:end));
F = double(A(2:2:end));
G = 64*E+F/4;
Raw_Intensity = reshape(G,[640 480])';

function TIFF_Intensity = Get_Tiff_Intensity(File_Path,File_Name)
TIFF_Intensity = double(imread([File_Path '\' File_Name]));

%------------------------------Time_Resolution---------------------------------

function Time_Resolution_Callback(hObject, eventdata, handles)

temp1 = get(hObject,'String');
n = str2num(temp1);
handles.Time_Res = n;

names = fields(handles);
temp_num = zeros(length(names),1);
for ii = 1:length(names)
    if strcmp(names{ii},'PreviewNum')
        temp_num(ii) = 1;
    else
        temp_num(ii) = 0;
    end
end

if sum(temp_num) == 1
    Res_max = floor(handles.PreviewNum*handles.AC_Freq/handles.SampleRate);
else 
    msgbox('PreviewNum is required');
    return
end

if n > Res_max
    msgbox(['Time Resolution must small than ' num2str(Res_max)]);
    n = Res_max;
    set(hObject,'String',num2str(n));
end

clear temp_num
names = fields(handles);
temp_num = zeros(length(names),1);
for ii = 1:length(names)
    if strcmp(names{ii},'Intensity_FFT')
        temp_num(ii) = 1;
    else
        temp_num(ii) = 0;
    end
end

if sum(temp_num) == 1
    Intensity = handles.Intensity;
    Intensity0 = handles.Intensity0;
    Interval = floor(n*(1/handles.AC_Freq*handles.SampleRate));     % Translate 'Interval' Frames to one fig by FFT
    FFT_Frames = floor(handles.PreviewNum/Interval);   % after trans,we have FFT_length Frames new figs,or Amp figs
    handles.FFT_Frames = FFT_Frames;
    handles.Interval = Interval;
    delta_f = handles.SampleRate/Interval;
    Intensity_FFT = zeros(size(Intensity0,1),size(Intensity0,2),FFT_Frames);
    for ii = 1:size(Intensity0,1)
        for jj = 1:size(Intensity0,2)
            for kk = 1:FFT_Frames
                temp1 = Intensity(ii,jj,1 + (kk - 1)*Interval:...
                    kk*Interval);
                temp2 = fft(temp1);
                Intensity_FFT(ii,jj,kk) = 2*abs(temp2(1 + floor(handles.AC_Freq/delta_f)))/Interval;
                Intensity_FFT_complex(ii,jj,kk) = 2*temp2(1 + floor(handles.AC_Freq/delta_f))/Interval;
            end
        end
    end
    handles.Intensity_FFT = Intensity_FFT;
    handles.Intensity_FFT_complex = Intensity_FFT_complex;
    axes(handles.axes2);
    handles.hf2 = imagesc(Intensity_FFT(:,:,1));
    colormap(gray(256));
    axis off
    set(handles.axes2_slider,'Value',1);
    set(handles.axes2_slider,'Max',FFT_Frames);
    set(handles.axes2_text,'Value',1);
end

guidata(hObject, handles);

function Time_Resolution_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------PreviewNum----------------------------------

function Preview_Num_Callback(hObject, eventdata, handles)

PreviewNum = str2num((get(hObject,'String')));
handles.PreviewNum = PreviewNum;

set(handles.axes1_slider,'Min',1,'Max',handles.PreviewNum,'SliderStep',...
    [1/handles.PreviewNum 20/handles.PreviewNum]);
set(handles.axes2_slider,'Min',1,'Max',...
    floor(handles.PreviewNum/handles.SampleRate),...
    'SliderStep',...
    [1/floor(handles.PreviewNum/handles.SampleRate) ...
    2/floor(handles.PreviewNum/handles.SampleRate)]);

guidata(hObject,handles);

function Preview_Num_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------exp_name----------------------------------

function Exp_Name_Callback(hObject, eventdata, handles)

Expname = get(hObject,'String');
handles.Expname = Expname;

guidata(hObject,handles);

function Exp_Name_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------闁告劕绉存慨鐐存姜閽樺娈ょ�?1�?7殿喚濮村ù姗�1�?7鎮ч敓鏂ゆ�?----------------------------------
function uipanel14_CreateFcn(hObject, eventdata, handles)

% --- Executes on button press in Sample_Time.
function Sample_Time_Callback(hObject, eventdata, handles)

function Connected_Time_Callback(hObject, eventdata, handles)

function Sub_Callback(hObject, eventdata, handles)

temp1 = get(hObject,'Value');
if temp1 == 0
    axes(handles.axes1)
    temp2 = ceil(get(handles.axes1_slider,'Value'));
    Intensity = handles.Intensity;
    Intensity0 = handles.Intensity0;
    imagesc(Intensity(temp2) + Intensity0);
    axis off
    colormap(handles.map)
end
guidata(hObject,handles)

% ------------------------------缁绢収鍠栭悾楣冨礉閻樺灚鏆╅柡鍐ㄧ埣濡拷锟�?1�?7------------------------------
function Timer_Callback(hObject, eventdata, handles)

names = split_(handles.Directory_Name,'\');
temp1 = strcmp(names,'TIFF');
temp2 = join_(names(1:find(temp1 == 1) -1),'\');

if isempty(temp2{1}) == 1
    temp2 = join_(names(1:end - 1),'\');
end

if exist([temp2{1} '\Timer\' handles.Expname '.mat']) ~= 0
    temp3 = load([temp2{1} '\Timer\' handles.Expname '.mat']);
    figure
    plot((temp3.t)',temp3.data);
    xlabel('t (s)','fontsize',15,'fontweight','bold');
    ylabel('Intensity(a.u.)','fontsize',15,'fontweight','bold');
    set(gca,'fontsize',15,'fontweight','bold');
    axis tight
    Fs = handles.SampleRate;
    handles.begin = TriggerTime(temp3.data,temp3.t,Fs);
    set(handles.Pike,'String',num2str(handles.begin.Pike));
    set(handles.Cs,'String',num2str(handles.begin.Cs));
    set(handles.Frame,'String',num2str(handles.begin.Frame));
    guidata(hObject,handles);
else
    msgbox('The Timer_data you want to analyze doesn''t exist');
    handles.begin.Pike = 0;
    handles.begin.Cs = 0;
    handles.begin.Frame = 1;
    set(handles.Pike,'String',num2str(0));
    set(handles.Cs,'String',num2str(0));
    set(handles.Frame,'String',num2str(1));
    guidata(hObject,handles);
end

guidata(hObject,handles);

% ---------------------------濞戞挸锕ｉ崥澶庛亹閳哄嫬顥�1�?7------------------------------------
function PseudoColor_Callback(hObject, eventdata, handles)

    prompt = {'colormap'};
    dlgtitle = 'Colormap';
    dims = [1 35];
    definput = {'parula(256)'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    if size(answer,1) == 0
        return
    end
    map = eval(answer{1});
    handles.map = map;
    
    axes(handles.axes1)
    colormap(map)
    axes(handles.axes2)
    colormap(map)

guidata(hObject,handles);

%-------------------------闁规惌浜滃ù姗�?1�?7宕戝鎰佹綊濡府鎷烽敓锟斤�?-----------------------------------
function Video_Callback(hObject, eventdata, handles)

%--------------闁告瑥鍊归弳鐔烘媼閸撗呮�?-------------------
prompt = {'From','Frames','FrameRate','Name','colormap'};
dlgtitle = '閻熸瑥妫濋�?�鍫曞矗閸屾稒娈堕悹浣稿⒔閻わ拷锟�1�?7';
dims = [1 35];
definput = {'1','1000','50','A1','parula(256)'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
if size(answer,1) == 0
    return
end
From = str2num(answer{1});
Frames = str2num(answer{2});
FR = str2num(answer{3});
name = answer{4};
map = answer{5};
if (From + Frames - 1) > handles.PreviewNum
    CreateStruct.Interpreter = 'tex';
    CreateStruct.WindowStyle = 'modal';
    msgbox('\fontsize{20}Frames閻℃帒鎳撶换鍐炊閸撗冾暬鐎殿喚濮甸弳锟斤拷',CreateStruct);
    Frames = handles.PreviewNum - From + 1
end
%------------------闁告帗甯掕ぐ鍥�?枎闊厾绠介柣锝嗙懃鐏忣垶宕洪敓鏂ゆ�?------------------
axes(handles.axes1)
roi = imrect;
p = floor(getPosition(roi));
delete(roi);
%----------------閻犱礁澧介悿鍡樼┍濠靛棛鎽犻悹渚灠缁讹拷锟�?1�?7----------------------
Directory_Name = handles.Directory_Name;
names = split_(Directory_Name,'\');
temp1 = strcmp(names,'TIFF');
temp2 = join_(names(1:find(temp1 == 1) - 1),'\');

if isempty(temp2{1}) == 1
    temp2 = join_(names(1:end - 1),'\');
end

temp3 = [temp2{1} '\Result\Temp Processing\' handles.Expname];
if exist(temp3) == 0
    mkdir(temp3)
end
saveroute = [temp3 '\' name '.avi'];
%-------------------闂佹彃娲▔锕傚极閻�?牆绁�?1�?7--------------------
Intensity = handles.Intensity;
Intensity0 = handles.Intensity0;
temp = imcrop(Intensity0,p);
crop_fig = zeros(size(temp,1),size(temp,2),Frames);
for ii = 1:Frames
    if get(handles.Sub,'Value') == 0
        crop_fig(:,:,ii) = imcrop((Intensity(:,:,From + ii - 1) + Intensity0),p);
    else
        crop_fig(:,:,ii) = imcrop(Intensity(:,:,From + ii - 1),p);
    end
end
Max_ = max(crop_fig(:));
Min_ = min(crop_fig(:));
crop_fig = im2uint8((crop_fig - Min_)/(Max_ - Min_));
%-------------------闂佹彃娲▔锔炬喆閸℃侗鏆�1�?7--------------------
v = VideoWriter(saveroute,'Indexed AVI');
v.FrameRate = FR;
v.Colormap = eval(map);
open(v)
hwait = waitbar(0,'0');
for ii = 1:Frames

    writeVideo(v,crop_fig(:,:,ii));
    waitbar(ii/Frames,hwait,[num2str(floor(ii/Frames*100)) '%']);

end
delete(hwait)

guidata(hObject,handles);

%---------------------Timer相关---------------------------------

function Pike_Callback(hObject, eventdata, handles)

Pike = str2num(get(hObject,'String'));
handles.begin.Pike = Pike;

guidata(hObject,handles)

function Cs_Callback(hObject, eventdata, handles)

Cs = str2num(get(hObject,'String'));
handles.begin.Cs = Cs;
guidata(hObject,handles)

function Frame_Callback(hObject, eventdata, handles)

Frame = str2num(get(hObject,'String'));
handles.begin.Frame = Frame;
guidata(hObject,handles)

function Pike_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Cs_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Cs_DeleteFcn(hObject, eventdata, handles)

function Frame_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Cs__Callback(hObject, eventdata, handles)

persistent CS_Dir

if ~ischar(CS_Dir)
    CS_Dir = [];
    [CS_name,CS_path] = uigetfile;
    CS_Dir = CS_path;
else
    [CS_name,CS_path] = uigetfile(fullfile(CS_Dir,'*.cor'));
end

if CS_path == 0
    return
end

CS_data = importdata(fullfile(CS_path,CS_name));
V = CS_data.data(:,1);
I = CS_data.data(:,2);
t = CS_data.data(:,3);

figure
subplot(121)
plot(V,I);
title('Potential-Curernt');
xlabel('Potential(V vs. Ag/AgCl)');
ylabel('Current(A)');
subplot(122)
hAx = plotyy(t,V,t,I);
ha = findobj(gcf,'type','axes');
xlabel('Time (s)');
ylabel(hAx(1),'Potential (V vs. Ag/AgCl)');
ylabel(hAx(2),'Current (A)');
title('Time-Potential-Curernt');
set(ha,'fontsize',15,'fontweight','bold','titlefontweight','bold');

    
    
