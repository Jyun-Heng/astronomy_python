import astropy
import numpy as np
import glob,sys,os
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.visualization import ZScaleInterval

imglist = glob.glob(os.path.join('/Users/linjyunheng/Documents/ad_obs_astrophysics/M6/Flat','*.fit'))
# dark_path = glob.glob(os.path.join('/Users/linjyunheng/Documents/ad_obs_astrophysics/M6/bias_dark/dark','*.fit'))
bias_path = glob.glob(os.path.join('/Users/linjyunheng/Documents/ad_obs_astrophysics/M6/bias_dark/bias','*.fit'))
# print(imglist)
flatlist = []
darklist = []
bias_list = []
exptime = []
imglist = sorted(imglist)
# dark_path = sorted(dark_path)
bias_path = sorted(bias_path)

# for i in dark_path:
#     imgdata,imgheader = fits.getdata(i,header=True)
#     exptime = imgheader['EXPTIME']
#     darklist.append(imgdata/exptime)
# dark = np.median(np.array(darklist),axis=0)

for i in bias_path:
    imgdata = fits.getdata(i)
    bias_list.append(imgdata)
bias = np.mean(np.array(bias_list),axis=0)
print(bias)

for i in imglist:
    imgdata,imghdr = fits.getdata(i,header=True)
    if imghdr['FILTER'] == 'V':
        imgdata = imgdata - bias
        flatlist.append(imgdata/np.nanmean(imgdata,axis=0))
    else:
        continue
flat = np.median(np.array(flatlist),axis=0)
print(flat)

zsc = ZScaleInterval()
Zmin, Zmax = zsc.get_limits(flat)
plt.imshow(flat, vmin=Zmin, vmax=Zmax, cmap='gray',origin='lower')
plt.colorbar()
plt.show()
plt.close()
fits.writeto('/Users/linjyunheng/Documents/ad_obs_astrophysics/M6/Flat/master_flat_V.fits',flat,overwrite=True)
# fits.writeto('/Users/linjyunheng/Documents/LOT_data/SR12_C/20230524/flat/flat.new',flat_master,overwrite=True)

