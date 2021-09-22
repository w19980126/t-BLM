# 工作流程
import sys
sys.path.append(r'D:\Data_Processing\pycode')
import t_blm
import imageio
import os
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np

# 初始化
expname = 'B3'
filepath,files,mask = t_blm.int_fun()

# 读入图片
rawtiffs1,rawtiffs2 = t_blm.tiff_imread(filepath, files)
saveroute = 'D:/work/t_LBM/Optical Experiment/20210903_t-BLM/Result/'+expname

# 对时间序列做分析
plt.imshow(plt.imread(filepath+'/'+files[10]))
I = t_blm.get_roi_intensity(rawtiffs1, rawtiffs2, filepath,files)
plt.plot(I)
plt.title('The light intensity curve on t-BLM')
os.makedirs(saveroute)
plt.savefig(saveroute+'/'+'膜上信号曲线')
# 各周期数据堆在一起看
temp = I[678:40678]
plt.figure()
for ii in range(10):
    plt.plot(temp[4000*ii:4000*(ii+1)],label = str(ii))
plt.legend(loc = 1)
plt.xlabel('Frames')
plt.ylabel('Intensity')
plt.title('The light intensity curve on t-BLM')
plt.savefig(saveroute + '/' + '膜上不同周期信号曲线叠加')

# 单独一个周期的信号分析
I_sum = I[16620:18620]
plt.figure()
ax1 = plt.subplot2grid((6,4),(0,3),rowspan=2)
plt.plot(I_sum)
plt.xlabel('Frames',size=8)
plt.ylabel('Intensity',size=8)
plt.xticks(size=8)
plt.yticks(size=8)
plt.title('The 9th cycle light intensity \n curve on t-BLM',size=8)
ax2 = plt.subplot2grid((6,4),(1,0),rowspan=4,colspan=2)
temp = t_blm.stft(I_sum, 200, 175)
# temp[0,:] = 0
plt.imshow(temp,aspect='auto')
plt.xlabel('Time(*0.25s)',size=8)
plt.ylabel('Frequency(*0.5Hz)',size=8)
plt.xticks(size=8)
plt.yticks(size=8)
plt.title('The spectrum of the 9th cycle time \n domain signal on t-BLM',size=8)
ax3 = plt.subplot2grid((6,4),(4,3),rowspan=2)
plt.plot(temp[0,:])
plt.xlabel('Time(*0.25s)',size=8)
plt.ylabel('Frequency(*0.5Hz)',size=8)
plt.xticks(size=8)
plt.yticks(size=8)
plt.title('The amplitude of the 9th cycle time \ndomain signal on t-BLM at 0Hz',size=8)
plt.savefig(saveroute + '/' + '膜上第9圈信号分析')

I_sum = t_blm.data_sum(I[678:40678], 2000)
plt.figure()
ax1 = plt.subplot2grid((6,4),(0,3),rowspan=2)
plt.plot(I_sum)
plt.xlabel('Frames',size=8)
plt.ylabel('Intensity',size=8)
plt.xticks(size=8)
plt.yticks(size=8)
plt.title('The sum_averaged light intensity \n curve on t-BLM',size=8)
ax2 = plt.subplot2grid((6,4),(1,0),rowspan=4,colspan=2)
temp = t_blm.stft(I_sum, 200, 175)
# temp[0:2,:] = 0
plt.imshow(temp,aspect='auto')
plt.xlabel('Time(*0.25s)',size=8)
plt.ylabel('Frequency(*0.5Hz)',size=8)
plt.xticks(size=8)
plt.yticks(size=8)
plt.title('The spectrum of the sum_averaged time \n domain signal on t-BLM',size=8)
ax3 = plt.subplot2grid((6,4),(4,3),rowspan=2)
plt.plot(temp[0,:])
plt.xlabel('Time(*0.25s)',size=8)
plt.ylabel('Frequency(*0.5Hz)',size=8)
plt.xticks(size=8)
plt.yticks(size=8)
plt.title('The amplitude of the sum_averaged time \ndomain signal on the t-BLM at 0Hz',size=8)
plt.savefig(saveroute + '/' + '处理后的膜上时域信号的频谱图')

# 伪彩色并转为gif
tiff_merged = t_blm.tiff_merge(rawtiffs1, rawtiffs2, 2000,start = 678,end = 40678)

B3 = t_blm.extr_sig_in_BF(tiff_merged, 100, 0, 200, 175)
B = []
for ii in range(72):
    B.append(B3[:,:,ii]-B3[:,:,0])
a1,a2,a3 = plt.hist(B3[:,:,40].ravel(),bins = 256)
vmin = -100
vmax = 150
savepath = saveroute+'/'+'ForGif'
os.makedirs(savepath)
plt.figure()
for ii in tqdm(range(72)):
    I = plt.imshow(B[ii],vmin = vmin,vmax = vmax)
    plt.title('%.2f s' %(ii*0.25))
    plt.xticks([])
    plt.yticks([])
    plt.colorbar()
    plt.savefig(savepath+'/'+str(ii))
    plt.close()

gif = []
for ii in tqdm(range(np.size(B3,2))):
    gif.append(imageio.imread(savepath + '/' + str(ii) + '.png'))
gifpath = saveroute + '/B3.gif'
imageio.mimsave(gifpath,gif)
    
# 保存数据
import pickle
A1 = {}
A1['B3_stfted_pics'] = B3
A1['filepath'] = filepath
A1['files'] = files
A1['mask'] = mask
pickle.dump(A1,open(saveroute+'/'+'B3','wb'))

