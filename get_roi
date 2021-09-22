
    
    #--------------------------------------
    
def get_roi(figs):
    
    import numpy as np
    import cv2
    from tqdm import tqdm
    import matplotlib.pyplot as plt
    import time
    
    time_start = time.time()
    img = (np.round(A2[:,:,0])).astype(np.uint8)
            
# -------------定义OnMouseAction函数-------------------   
    def OnMouseAction(event, x, y, flags, param):
        
        # nonlocal img
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
    cv2.namedWindow('ROI_Slection')
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
    
    length = np.size(A2,2)
    intensity = np.empty(length,np.float64) 
    for ii in tqdm(range(length)):
        intensity[ii] = np.sum(A2[row[0]:row[1],col[0]:col[1],ii])/NormalizationFctor

# -----------通过返回intensity----------------------        
    plt.figure()    
    plt.plot(np.arange(length),intensity) 
    plt.xlabel('Frames')
    plt.ylabel('Amplitude')
    time_end = time.time()
    print('\n耗时' + str(time_end-time_start) + 's' )
    
    return intensity
