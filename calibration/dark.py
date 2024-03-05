from astropy.io import fits
import numpy as np
import glob,os,sys
import matplotlib.pyplot as plt
from astropy.visualization import ZScaleInterval

dark_path = glob.glob(os.path.join('/Users/linjyunheng/Documents/ad_obs_astrophysics/M22/bias_dark/dark/n_10deg','*8min*'))
# imglist = glob.glob(os.path.join('/Users/linjyunheng/Documents/ad_obs_astrophysics/M6/bias_dark/dark','*.fit'))
# imglist = sorted(imglist)
darklist = []
for i in dark_path:
    imgdata,imghdr = fits.getdata(i,header=True)
    exptime = imghdr['EXPTIME']
    darklist.append(imgdata)    
print(darklist)
dark = np.median(np.array(darklist),axis=0)
print(np.shape(dark))
fits.writeto('/Users/linjyunheng/Documents/ad_obs_astrophysics/M22/bias_dark/master_dark_8min_n10deg.fits',dark,overwrite=True)
zsc = ZScaleInterval()
Zmin, Zmax = zsc.get_limits(dark)
plt.imshow(dark, vmin=Zmin, vmax=Zmax, cmap='gray')
plt.colorbar()
plt.show()
plt.close()
     