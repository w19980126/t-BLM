# t-BLM模块用来处理纳米孔道实验相关的实验数据

######################################################################################
def Data_Sum(Raw_Data,Period):
    
    """ 
    通过周期电平对系统进行扰动，然后将每个周期的信号
    叠加到一起，以达到压制噪声提高信噪比的目的
    Raw_Data：原始周期数据
    Period：每个周期包含的数据点个数
    """
    
    import numpy as np
    Cycle = round(np.size(Raw_Data)/Period)
    print(Cycle)
    Data_Out = np.zeros(Period)
    
    for ii in np.arange(Period):
        for jj in np.arange(Cycle):
            Data_Out[ii] += Raw_Data[jj*Period + ii]
            
    Data_Out /= Cycle
    return Data_Out

######################################################################################
def STFT(Raw_Data,window):
    
    """
    STFT，所谓短时傅里叶变换，所不同的是，这里每隔一个信号就在窗内进行一个FFT，以提高时间分辨率
    Raw_Data：原始数据
    window：窗口大小
    """
    
    import numpy as np
    from scipy.fftpack import fft
    import copy
    
    Pre_Length = int(window/2)
    End_Length = window - 1 - Pre_Length
    x = np.zeros(len(Raw_Data) + window)
    x[np.arange(Pre_Length,Pre_Length+len(Raw_Data))] = copy.deepcopy(Raw_Data)
    x[0:Pre_Length] = 0
    x[-End_Length:] = 0
    
    Out_Data = np.empty([int(window/2),len(Raw_Data)])
    for ii in np.arange(len(Raw_Data)):
        temp = fft(x[ii:ii+window])
        Out_Data[:,ii] = 2/window*abs(temp[np.arange(window/2,dtype = int)])
        
    return Out_Data 
