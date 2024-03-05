from astropy.io import fits
import numpy as np
import glob,os,sys
import matplotlib.pyplot as plt
from astropy.visualization import ZScaleInterval
from ccdproc import cosmicray_lacosmic

def master_dark(darkpath,exptime,savefile = False,newfilepath = None):
    imglist = glob.glob(os.path.join(darkpath,'*'))
    # imglist = sorted(imglist)
    darklist = []
    for i in imglist:
        imgdata,imghdr = fits.getdata(i,header=True)
        if np.shape(imgdata) != (2048,2048):
            continue
        if imghdr['IMAGETYP'] != 'DARK':
            print(imghdr['IMAGETYP'])
            continue
        if exptime == imghdr['EXPTIME']:
            print(imghdr['EXPTIME'])
            darklist.append(imgdata)
    dark = np.median(np.array(darklist),axis=0)
    # print(dark)
    if savefile:
        fits.writeto(newfilepath,dark,overwrite=True)
    return dark


def master_dark_for_flat(darkpath,exptime,savefile = False,newfilepath = None):
    imglist = glob.glob(os.path.join(darkpath,'*'))
    # print(imglist)
    # imglist = sorted(imglist)
    darklist = []
    for i in imglist:
        imgdata,imghdr = fits.getdata(i,header=True)
        dark_exptime = imghdr['EXPTIME']
        if imghdr['IMAGETYP'] != 'DARK':
            # print(imghdr['IMAGETYP'])
            continue
        if np.shape(imgdata) != (2048,2048) or dark_exptime != exptime:
            continue     
        darklist.append(imgdata)
    dark = np.median(np.array(darklist),axis=0)
    # print(dark)
    if savefile:
        fits.writeto(newfilepath,dark,overwrite=True)
    return dark

def master_bias(biaspath,savefile = False,newfilepath = None):
    imglist = glob.glob(os.path.join(biaspath,'*'))
    biaslist = []
    for i in imglist:
        imgdata,imghdr = fits.getdata(i,header=True)
        if np.shape(imgdata) != (2048,2048):
            continue
        if imghdr['IMAGETYP'] != 'BIAS':
            continue
        biaslist.append(imgdata)
    bias = np.median(np.array(biaslist),axis=0)
    # print(bias)
    if savefile:
        fits.writeto(newfilepath,bias,overwrite=True)
    return bias

def master_flat(flatpath,darkpath,biaspath,filter,savefile = False,newfliepath = None):
    imglist = glob.glob(os.path.join(flatpath,'*'))
    # print(imglist)
    flatlist = []
    count = 0
    for i in imglist:
        imgdata,header = fits.getdata(i,header=True)
        if np.shape(imgdata) != (2048,2048) :
            print(np.shape(imgdata),'image shape wrong')
            continue 
        if header['IMAGETYP'] != 'FLAT':
            print(header['IMAGETYP'],'image type wrong')
            continue
        if header['FILTER'] != filter:
            print(header['FILTER'],'filter wrong')
            continue
        if np.nanmean(imgdata) > 49000:
            print(np.nanmean(imgdata),'flat is out of linear range')
            continue
        # print(np.max(imgdata))
        exptime=header['exptime']
        if count == 0:
            dark = master_dark_for_flat(darkpath,exptime)
            print('******************dark') 
            bias = master_bias(biaspath)      
            print('******************bias')
        elif exptime != header['exptime']:
            exptime = header['exptime']
            dark = master_dark_for_flat(darkpath,exptime)
            print('******************dark')
        imgdata = imgdata - dark - bias
        flatlist.append(imgdata/np.nanmedian(imgdata))
        print(flatlist)
        # print(imgdata/np.median(imgdata))
        print('wrong' if (imgdata/np.median(imgdata,axis=0)).any() < 0 else 'right')
        count = count + 1
    flat = np.median(np.array(flatlist),axis=0)
    # print(flat)

    if savefile:
        fits.writeto(newfliepath,flat,overwrite=True)
    return flat
 
def cosmicray(imgarr):
    newimg,mask = np.array(cosmicray_lacosmic(imgarr,readnoise=8.5, gain=2., sigclip=4.5, cleantype='medmask', verbose=True)) 
    return newimg,mask

def calibration(darkpath,flatpath,biaspath,imgpath,filter,calibrated_path,file_name,f_dark = False,sigmaclip  = False,flat_other_path = None,flat_other_dark_path = None):
    imglist = glob.glob(os.path.join(imgpath,'*'))
    # print(imglist)
    imglist = sorted(imglist)
    # print(exptime)
    if flat_other_path==None:
        if f_dark:
            flat = master_flat(flatpath,flatpath,biaspath,filter=filter)
        else:
            flat = master_flat(flatpath,darkpath,biaspath,filter=filter)
    else:
        if f_dark:
            flat = master_flat(flat_other_path,flat_other_path,biaspath,filter=filter)
        else:
            flat = master_flat(flat_other_path,flat_other_dark_path,biaspath,filter=filter)
    count = 1
    for i in imglist:
        rawimg,imghdr = fits.getdata(i,header=True)
        exptime = imghdr['EXPTIME'] if count == 1 else exptime
        dark = master_dark(darkpath,exptime) if count == 1 else dark
        if exptime != imghdr['EXPTIME']:
            exptime = imghdr['EXPTIME']
            dark = master_dark(darkpath,exptime)
        # print(rawimg)
        imgdata = (rawimg - dark)/flat
        cleanimg = imgdata
        if sigmaclip:
            cleanimg,mask = cosmicray(imgdata)
        else:
            cleanimg = imgdata
        print('wrong') if (imgdata/np.median(imgdata,axis=0)).any() < 0 else print('right')
        fits.writeto(calibrated_path + f'{file_name}-{count:003}ha{int(exptime)}s.fits',cleanimg,imghdr,overwrite=True)
        count = count + 1
    return None

file_name = '/LOT20240217_GJ3470'
darkpath = '/Users/linjyunheng/Documents/LOT_data/calibration_data/LOT20240217/bias-dark'
flatpath = '/Users/linjyunheng/Documents/LOT_data/calibration_data/LOT20240217/flat'
biaspath = '/Users/linjyunheng/Documents/LOT_data/calibration_data/LOT20240217/bias-dark'
imgpath = '/Users/linjyunheng/Documents/LOT_data/GJ3470/LOT20240217/ntnuphy'
calibrated_path = '/Users/linjyunheng/Documents/LOT_data/GJ3470/LOT20240217/calibrated'
flat_other_path = '/Users/linjyunheng/Documents/LOT_data/calibration_data/LOT20240202/flat'
flat_other_dark_path = '/Users/linjyunheng/Documents/LOT_data/calibration_data/LOT20240202/bias-dark'

'''
darkpath: the path of dark folder
flatpath: the path of flat folder
imgpath: the path of raw image folder
filter: the filter of the image we use. (see header of the image)
calibrated_path: the path of calibrated image folder
flat_other_path: the path of flat folder on other day
flat_other_dark_path: the path of dark folder on other day

if the dark of the flat is save in the same folder with the flat, f_dark = True; otherwise, f_dark = False. Defult is False.
if you want to use sigmaclip to remove the cosmic ray, sigmaclip = True; otherwise, sigmaclip = False. Defult is False.
'''
calibration(darkpath=darkpath,flatpath=flatpath,biaspath=biaspath,imgpath=imgpath,filter='Ha_6563_30',file_name = file_name,
            calibrated_path=calibrated_path,f_dark = False,sigmaclip  = True,flat_other_path = None,flat_other_dark_path = None)

# imglist = glob.glob(os.path.join(imgpath,'*.fts'))
# for i in imglist:
#     data = fits.getdata(i)
#     print(max(data))


