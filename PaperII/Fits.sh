#!/bin/bash
source /home/silvia/miriad/MIRRC.sh

export MIRIADPATH=/home/silvia/miriad/

#Creo un cutout di MWA, CHIPASS e SPASS della stessa dimensione (3.5 deg x 3.5 deg)
if [[ ! -e SNR_G321.3-3.9_MWA_rg.fits ]]
then
	python3 Cutout.py MWA.fits SNR_G321.3-3.9_MWA_rg.fits
fi

if [[ ! -e SNR_G321.3-3.9_SUMSS.fits ]]
then
	python3 Cutout.py SUMSS.fits SNR_G321.3-3.9_SUMSS.fits
fi

if [[ ! -e SPASS_model_subtracted_nb.fits ]]
then
	python3 Cutout.py SPASS.fits SPASS_model_subtracted_nb.fits
fi

if [[ ! -e SPASS_q_rg_mod.fits ]]
then
	python3 Cutout.py SPASS_q.fits SPASS_q_rg_mod.fits
fi

if [[ ! -e SPASS_u_rg_mod.fits ]]
then
	python3 Cutout.py SPASS_u.fits SPASS_u_rg_mod.fits
fi

if [[ ! -e CHIPASS_rg_subtracted_nb.fits ]]
then
	python3 Cutout.py CHIPASS.fits CHIPASS_rg_subtracted_nb.fits
fi

if [[ ! -e SPASS_q_rg_mod_nb.fits ]]
then
	python3 Blank-remover.py --infits SPASS_q_rg_mod.fits
	python3 Blank-remover.py --infits SPASS_u_rg_mod.fits
fi

if [[ ! -e Halpha_rg_nb.fits ]]
then
	python3 Blank-remover.py --infits Halpha_rg.fits
fi

if [[ ! -e SNR_G321.3-3.9_EMU_rg.fits ]]
then
	#python3 Cutout.py SNR_G321.3-3.9_EMU.fits SNR_G321.3-3.9_EMU_cut.fits
	fits op=xyin in=SNR_G321.3-3.9_EMU.fits out=SNR_G321.3-3.9_EMU.im
	fits op=xyin in=SNR_G321.3-3.9_SUMSS.fits out=SNR_G321.3-3.9_SUMSS.im
	regrid tin=SNR_G321.3-3.9_SUMSS.im axes=1,2 in=SNR_G321.3-3.9_EMU.im out=EMU_rg.im
	fits op=xyout in=EMU_rg.im out=SNR_G321.3-3.9_EMU_cut_rg.fits
	python3 /home/silvia/DS9_file/Summary_Plot/Candidates/G321.3-3.9/SpectralIndex/Blank-remover.py --infits SNR_G321.3-3.9_EMU_cut_rg.fits
fi