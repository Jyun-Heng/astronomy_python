import numpy as np
import matplotlib.pyplot as plt
import glob,os,sys,time
import astropy.io.fits as fits
from astropy.visualization import ZScaleInterval

imgpath = glob.glob(os.path.join('/Users/linjyunheng/Documents/ad_obs_astrophysics/M6/Science','*.fit'))
imgpath = sorted(imgpath)

def dark_master(deg,exptime):
    if deg == -15:
        if exptime ==300:
            dark_master = fits.getdata('/Users/linjyunheng/Documents/ad_obs_astrophysics/M6/bias_dark/master_dark_5min.fits')
            return dark_master
        elif exptime == 480:
            dark_master = fits.getdata('/Users/linjyunheng/Documents/ad_obs_astrophysics/M6/bias_dark/master_dark_8min.fits')
            return dark_master
        elif exptime == 40:
            dark_master = fits.getdata('/Users/linjyunheng/Documents/ad_obs_astrophysics/M6/bias_dark/master_dark_40s.fits')
            return dark_master
        elif exptime == 30:
            dark_master = fits.getdata('/Users/linjyunheng/Documents/ad_obs_astrophysics/M6/bias_dark/master_dark_30sec.fits')
            return dark_master
        elif exptime == 5:
            dark_master = fits.getdata('/Users/linjyunheng/Documents/ad_obs_astrophysics/M6/bias_dark/master_dark_5s.fits')
            return dark_master
        elif exptime == 4:
            dark_master = fits.getdata('/Users/linjyunheng/Documents/ad_obs_astrophysics/M6/bias_dark/master_dark_4s.fits')
            return dark_master
        elif exptime == 3:
            dark_master = fits.getdata('//Users/linjyunheng/Documents/ad_obs_astrophysics/M6/bias_dark/master_dark_3s.fits')
            return dark_master
    if deg == -10:
        if exptime  == 300:
            dark_master = fits.getdata('/Users/linjyunheng/Documents/ad_obs_astrophysics/M6/Flat/master_flat_V.fits')
            return dark_master
        elif exptime == 480:
            dark_master = fits.getdata('/Users/linjyunheng/Documents/dark_150s.fits')
            return dark_master
    if deg == -5:
        if exptime == 300:
            dark_master = fits.getdata('/Users/linjyunheng/Documents/ad_obs_astrophysics/M6/Flat/master_flat_V.fits')
            return dark_master


        
count = 1
for i in imgpath:
    img_data,imghdr = fits.getdata(i,header=True)
    if imghdr['FILTER'] == 'V':
        flatmaster = fits.getdata('/Users/linjyunheng/Documents/ad_obs_astrophysics/M6/Flat/master_flat_V.fits')
        darkmaster = dark_master(-15,imghdr['EXPTIME'])
        img_data = (img_data-darkmaster)/flatmaster
    elif imghdr['FILTER'] == 'B':
        flatmaster = fits.getdata('/Users/linjyunheng/Documents/ad_obs_astrophysics/M6/Flat/master_flat_B.fits')
        darkmaster = dark_master(-15,imghdr['EXPTIME'])
        img_data = (img_data-darkmaster)/flatmaster
    filter = imghdr['FILTER']
    exptime = int(imghdr['EXPTIME'])
    fits.writeto(f'/Users/linjyunheng/Documents/ad_obs_astrophysics/M6/calibrated/M6_{filter}_{exptime}s_{count:003}.fits',img_data,imghdr,overwrite=True)
    count = count + 1


# zsc = ZScaleInterval()
# Zmin, Zmax = zsc.get_limits(img_data)
# plt.imshow(img_data, vmin=Zmin, vmax=Zmax, cmap='gray')
# plt.colorbar()
# # plt.show()
# plt.close()