from ccdproc import cosmicray_lacosmic
from astropy.io import fits
import glob,sys,os
import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import ZScaleInterval

imglist = glob.glob(os.path.join("/Users/linjyunheng/Documents/LOT_data/VV_ser/VV_ser/VV_ser_new", '*.fits'))
imgpath = sorted(imglist)
count = 1
for i in imgpath:
    imgarr = fits.getdata(i)
    header = fits.getheader(i)
    newimg,mask = np.array(cosmicray_lacosmic(imgarr,readnoise=8.5, gain=2., sigclip=5, cleantype='medmask', verbose=True))
    fits.writeto(f'/Users/linjyunheng/Documents/LOT_data/VV_ser/VV_ser/VV_ser_new/20230524_VV_ser-{count:003}ha.fits',newimg,header,overwrite=True)
    count =count+1
# zsc = ZScaleInterval() 
# Zmin, Zmax = zsc.get_limits(newimg)
# plt.imshow(newimg, vmin=Zmin, vmax=Zmax, cmap='gray')
# plt.colorbar()
# plt.show()
# plt.close()

