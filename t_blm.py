# t-BLM模块用来处理纳米孔道实验相关的实验数据

######################################################################################
def int_fun():
    '''
    初始化函数，用来生成filepath字符串和files列表，以及mask，这里的mask是通过二值化生成的，所以很局限

    '''
    
    import os
    import tkinter as tk
    from tkinter import filedialog
    import numpy as np
    import cv2
    import matplotlib.pyplot as plt
    
    global filepath,files,mask
    
    window = tk.Tk()
    window.withdraw()
    filepath = filedialog.askdirectory()
    files = os.listdir(filepath)   
    
    # 生成mask 
    filename0 = filepath + "/" + files[0]
    tiff0_uint8 = cv2.imread(filename0,cv2.IMREAD_GRAYSCALE)
    hist,_ = np.histogram(tiff0_uint8.ravel(),bins = 256,range = (0,256))    
    plt.figure()
    plt.plot(hist)
    plt.show()
    plt.pause(5)
    threshold = int(input('请输入分割阈值：\n'))
    ret,thre = cv2.threshold(tiff0_uint8,threshold,255,cv2.THRESH_BINARY)
    mask = np.empty((480,640),np.uint8)
    mask[:] = (thre>1)
    kernel = np.ones(9)
    mask = cv2.erode(cv2.dilate(mask,kernel),kernel)
    
    return filepath,files,mask
    
######################################################################################    
def data_sum(Raw_Data,Period):    
    '''
    
    Parameters
    ----------
    Raw_Data : 原始周期数据
    Period : 每个周期包含的数据点个数

    Returns
    -------
    Data_Out : 通过周期电平对系统进行扰动，然后将每个周期的信号叠加到一起，以达到压制噪声提高信噪比的目的

    '''
    
    import numpy as np
    Cycle = round(np.size(Raw_Data)/Period)
    temp = np.empty((Period,Cycle))
    for ii in range(Cycle):
        temp[:,ii] = Raw_Data[ii*Period:(ii+1)*Period]
    Data_Out = np.sum(temp,axis = 1)/Cycle
    
    return Data_Out

######################################################################################
def get_roi_intensity(rawtiffs1,rawtiffs2,filepath,files,start = True,end = True):
    '''

    Parameters
    ----------
    rawtiffs1 / rawtiffs2 : 两个tiff矩阵，因为把所有的图片都存在同一个矩阵里面太占内存，
    所有拆分成两个矩阵进行输入
    filepath/files : 用于读取图片选取roi.
    start : 开始帧数，默认为0
    end : 结束帧数，默认为np.size(rawtiffs1,2)

    Returns
    -------
    intensity : 所选取的roi中强度平均值时间序列

    '''
    
    import numpy as np
    import cv2
    from tqdm import tqdm
    import matplotlib.pyplot as plt
    import time
    
    time_start = time.time()
    img = cv2.imread(filepath+'/'+files[0],cv2.IMREAD_GRAYSCALE)
    img = np.around(((img-np.min(img))/(np.max(img)-np.min(img))*255)).astype(np.uint8)
    if start == True:
        start = 0
    if end == True:
        end = np.size(rawtiffs1,axis = 2)
    
            
# -------------定义OnMouseAction函数-------------------   
    def OnMouseAction(event, x, y, flags, param):
        
        nonlocal img
        global position1,position2  
        image = img.copy()   
        
        if event == cv2.EVENT_LBUTTONDOWN:                                          #按下左键
            position1 = (x,y)                                                     #获取鼠标的坐标(起始位置)     
        elif event == cv2.EVENT_MOUSEMOVE and flags == cv2.EVENT_FLAG_LBUTTON:      #按住左键拖曳不放开
            cv2.rectangle(image, position1, (x,y), (0,255,0), 3)                    #画出矩形选定框
            cv2.imshow('ROI_Slection', image)    
        elif event == cv2.EVENT_LBUTTONUP:                                          #放开左键
            position2 = (x,y)                                                       #获取鼠标的最终位置
            cv2.rectangle(image, position1, position2, (0,0,255), 3)                #画出最终的矩形 
            cv2.imshow('ROI_Slection', image)
# --------------通过SetMouseCallback调用函数------------
    cv2.imshow('ROI_Slection',img)
    cv2.setMouseCallback('ROI_Slection',OnMouseAction)
    cv2.waitKey(0)
    cv2.destroyWindow('ROI_Slection')
# --------------通过返回的position读取roi强度---------------
    row,col = [],[]
    col.append(min(position1[0],position2[0]))
    col.append(max(position1[0],position2[0]))
    row.append(min(position1[1],position2[1]))
    row.append(max(position1[1],position2[1]))
    NormalizationFctor = (col[1]-col[0])*(row[1]-row[0])
    
    length = end - start
    intensity = np.empty(length,np.float64) 
    if row[0]<240 and row[1]<240:
        for ii in tqdm(range(length)):
            jj = ii + start
            intensity[ii] = np.sum(rawtiffs1[row[0]:row[1],col[0]:col[1],jj])/NormalizationFctor
    elif row[0]>=240 and row[1]>=240:
        for ii in tqdm(range(length)):
            jj = ii + start
            intensity[ii] = np.sum(rawtiffs2[row[0]-240:row[1]-240,col[0]:col[1],jj])/NormalizationFctor
    else:
         for ii in tqdm(range(length)):
            jj = ii + start
            intensity[ii] = (np.sum(rawtiffs1[row[0]:,col[0]:col[1],jj]) + \
                             np.sum(rawtiffs1[0:row[1]-240,col[0]:col[1],jj]))/NormalizationFctor
        
# -----------通过返回intensity----------------------        
    plt.figure()
    plt.plot(start+np.arange(length),intensity) 
    plt.xlabel('Frames')
    plt.ylabel('Intensity')
    time_end = time.time()
    print('\n耗时' + str(time_end-time_start) + 's' )
    
    return intensity

######################################################################################
def tiff_imread(filepath,files):
    '''

    Parameters
    ----------
    filepath : tiff文件存储路径
    files : tiff图片名组成的列表，用来索引各个图片
    输出两个3维的矩阵，每个矩阵大小为240*640*len(files)，因为凑成一个大
    矩阵太占内存，以至于无法导入，所以分成两个小矩阵输出

    Returns
    -------
    rawtiffs1 / rawtiffs2 : 输出两个3维的矩阵，每个矩阵大小为240*640*len(files)，因为凑成一个大
    矩阵太占内存，以至于无法导入，所以分成两个小矩阵输出

    '''
    
    import time
    from tqdm import tqdm
    import numpy as np
    import cv2
    
    frames = len(files)
    rawtiffs1 = np.empty((240,640,frames),np.uint16)
    rawtiffs2 = np.empty((240,640,frames),np.uint16)
    start = time.time()
    for kk in tqdm(range(frames)):
        filename = filepath + '/' + files[kk]
        rawtiffs1[:,:,kk] = cv2.imread(filename,cv2.IMREAD_UNCHANGED)[:240,:]
        rawtiffs2[:,:,kk] = cv2.imread(filename,cv2.IMREAD_UNCHANGED)[240:,:]
    end = time.time()
    print('\n耗时' + str(end - start) + 's') 
    
    return rawtiffs1,rawtiffs2

######################################################################################
def tiff_merge(rawtiffs1,rawtiffs2,period,start = True,end = True):
    '''
    将多个周期的数据合并成一个周期的数据以达到压制噪声的目的
    
    Parameters
    ----------
    rawtiffs1 : 上半tiff图片矩阵
    rawtiffs2 : 上半tiff图片矩阵
    period : 一个周期包含的帧数
    start : 起始帧，默认为0       
    end : 终止帧，默认为最后一帧
    
    Returns
    -------
    tiff_out : 返回一个480*640*period大小的矩阵，用于后续的stft
    '''
    import numpy as np
    from t_blm import data_sum
    from tqdm import tqdm
# ---------------参数初始化----------------------    
    tiff_out = np.empty((480,640,period))
    if start == True:
        start = 0
    if end == True:
        end = np.size(rawtiffs1,axis = 2)
# -----------------数据合并-----------------------    
    for ii in tqdm(range(480)):
        for jj in range(640):
            Raw_Data = (ii<240)*rawtiffs1[ii-ii*(ii>=240),jj,start:end] +\
                (ii>=240)*rawtiffs2[ii-240,jj,start:end]
            tiff_out[ii,jj,:] = data_sum(Raw_Data, period)
# -----------------返回-----------------------
    return tiff_out

######################################################################################
def stft(Raw_Data,window,noverlap):
    '''
    Parameters
    ----------
    Raw_Data : 要求是已经进行周期合并后的时间序列
    window : FFT窗口大小，即多少帧进行一次fft
    noverlap : fft窗口之间的重叠帧数

    Returns
    -------
    Out_Data : 该矩阵是一个频谱图，行是频率，列是时间，强度已经进行过归一化
    '''    
    import numpy as np
    
    interval = int(window-noverlap) # 间隔interval帧进行一次fft
    col_len = int(np.floor((len(Raw_Data)-window)/interval))    # 输出矩阵的列数，即压缩后的时间轴数据点个数       
    temp_matrix = np.empty((window,int(col_len)))   # 临时存放输出数据的矩阵
    for ii in np.arange(col_len):
        temp_matrix[:,ii] = Raw_Data[interval*ii:interval*ii+window]
    Out_Data = np.fft.fft(temp_matrix,axis = 0)       
    Out_Data = abs(Out_Data[0:int(window/2)])*2/window  # 将频域的振幅转换成振幅强度
    Out_Data[0,:] = Out_Data[0,:]/2
    return Out_Data

######################################################################################
def extr_sig_in_BF(tiff_merged,Fs,BF,window,noverlap):
    '''
    
    函数名解释：在基频(BF)处提取信号：extract signal in BaseFrequency
    
    Parameters
    ----------
    tiff_merged : 要求是已经进行周期合并后的图像序列
    Fs : 采样频率，即相机采样速度
    BF : BaseFrequency，基频率，即给系统的扰动频率
    window : 进行stft的窗口大小
    noverlap : 进行stft的窗口之间的重叠数
    
    Returns
    -------
    out_data : 对每个像素点的时间序列挑去基频处的振幅并存在输出矩阵中
    out_data ：(480,640,时间)
    '''

    import numpy as np
    from tqdm import tqdm
    from t_blm import stft
    
    period = np.size(tiff_merged,axis = 2)
    interval = int(window-noverlap) # 间隔interval帧进行一次fft
    t_points = int(np.floor((period-window)/interval))  # 处理后的时间轴上的数据点个数
    
    f_res = Fs/window     # 频率分辨率
    BFloc = int(BF/f_res)   # 基频位置
    
    out_data = np.empty((480,640,t_points))
    for ii in tqdm(range(480)):
        for jj in range(640):
            temp = stft(tiff_merged[ii,jj,:],window,noverlap)
            out_data[ii,jj,:] = temp[BFloc,:]
    
    return out_data 
