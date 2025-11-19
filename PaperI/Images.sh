#!/bin/bash
#Program to make images for paper I
source /home/silvia/miriad/MIRRC.sh

export MIRIADPATH=/home/silvia/miriad/

#Cutout of GPM data.
if [[ ! -e SNR_G16.0+0.7_MEERKAT_2.fits ]]
then
	python Cutout.py /home/silvia/all_epochs_wavg.fits SNR_G16.0+0.7_GPM.fits
	python Cutout_meerkat.py G017.5+0.0IFx_Mosaic.fits SNR_G16.0+0.7_MEERKAT.fits
	python /home/silvia/Desktop/Papers/Paper_I/Axes_meerkat.py SNR_G16.0+0.7_MEERKAT.fits SNR_G16.0+0.7_MEERKAT_2.fits
	python ../New_fk5_template.py  --ra=274.2906556 --dec=-14.6538405 --pixscale=4.16666700000E-04 --fwhm=2.22219992429E-03 --nx=1200 --ny=1200 --output=Template.fits --overwrite
fi

#Download Wise image of 2 degrees in Sin projection from skyview (size = number of pixels times the internal pixel size of an image).
if [[ ! -e SNR_G16.0+0.7_WISE34.fits ]]
then
	perl /home/silvia/Desktop/Papers/Paper_I/skvbatch_wget.pl file=SNR_G16.0+0.7_WISE12.fits position='18 17 09, -14 39 13' Survey='WISE12' Coordinates='J2000' Projection=Sin Size=0.5 Pixels=1309
	perl /home/silvia/Desktop/Papers/Paper_I/skvbatch_wget.pl file=SNR_G16.0+0.7_WISE22.fits position='18 17 09, -14 39 13' Survey='WISE22' Coordinates='J2000' Projection=Sin Size=0.5 Pixels=1309
	perl /home/silvia/Desktop/Papers/Paper_I/skvbatch_wget.pl file=SNR_G16.0+0.7_WISE34.fits position='18 17 09, -14 39 13' Survey='WISE3.4' Coordinates='J2000' Projection=Sin Size=0.5 Pixels=1309
fi 

#Download RACS image of 2 degrees.
if [[ ! -e SNR_G16.0+0.7_RACS_2.fits ]]
then
	python Download_RACS.py 18:17:09.7573 -14:39:13.826
	python /home/silvia/Desktop/Papers/Paper_I/Axes_racs.py SNR_G16.0+0.7_RACS.fits SNR_G16.0+0.7_RACS_2.fits
fi

#Cutting images:
#if [[ ! -e SNR_G16.0+0.7_RACS.fits ]]
#then
#	python3 Cutout_wise.py /home/silvia/DS9_file/Summary_Plot/Wise/19:01:56.0932:+6:04:49.977_Band3.fits SNR_G16.0+0.7_WISE.fits
#	python3 Cutout_racs.py /home/silvia/DS9_file/Summary_Plot/RAC/19:01:56.0932:+6:04:49.977_RAC.fits SNR_G16.0+0.7_RACS.fits 
#fi

#Regridding to match same projection (I want to use Sin projection):
if [[ ! -e SNR_G16.0+0.7_MEERKAT_rg.fits ]]
then
	fits op=xyin in=SNR_G16.0+0.7_MEERKAT_2.fits out=SNR_G16.0+0.7_MEERKAT.im
	#fits op=xyin in=SNR_G16.0+0.7_RACS_2.fits out=SNR_G16.0+0.7_RACS.im
	fits op=xyin in=SNR_G16.0+0.7_GPM.fits out=SNR_G16.0+0.7_GPM.im
	fits op=xyin in=SNR_G16.0+0.7_WISE12.fits out=SNR_G16.0+0.7_WISE.im	#Using this to regrid others
	fits op=xyin in=Template.fits out=Template.im
	
	regrid tin=Template.im in=SNR_G16.0+0.7_MEERKAT.im out=MEERKAT.im
	#regrid tin=SNR_G16.0+0.7_WISE.im in=SNR_G16.0+0.7_RACS.im out=RACS.im
	regrid tin=SNR_G16.0+0.7_WISE.im in=SNR_G16.0+0.7_GPM.im out=GPM.im
	
	fits op=xyout in=MEERKAT.im out=SNR_G16.0+0.7_MEERKAT_rg.fits
	#fits op=xyout in=RACS.im out=SNR_G16.0+0.7_RACS_rg.fits
	fits op=xyout in=GPM.im out=SNR_G16.0+0.7_GPM_rg.fits
fi

#Plotting the figure:
python Plot-GPM-Wise-Radio.py SNR_G16.0+0.7_GPM_rg.fits SNR_G16.0+0.7_WISE12.fits SNR_G16.0+0.7_WISE22.fits SNR_G16.0+0.7_WISE34.fits SNR_G16.0+0.7_MEERKAT_rg.fits