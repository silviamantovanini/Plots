#!/home/silvia/PhD/bin/python3
#Paper 1: same as Statistic.py but using HII region catalog instead of candidates from other radio works.

import matplotlib.pyplot as plt   
import numpy as np
import sys, os
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
import seaborn as sns
import pandas as pd


plt.rcParams.update({
    "font.size": 15,
    "font.sans-serif": ["Helvetica"]})

cm = 1/2.54  # centimeters in inches

def mwrap(l):
	if l>180:
		l-=360
	return l
	
def posv(b):
	if b<0:
		b=-np.log10(-b)
	else:
		b=np.log10(b)
	return b

#Faccio un istogramma per vedere le dimensioni delle mie sorgenti.
dataS=np.loadtxt('Silvia_candidates.txt', dtype=str, delimiter=',', skiprows=1)
r_sec=np.float64(dataS[:,2])
r_min=r_sec/60

n=0
r_min_small=[]
r_min_big=[]
for r in r_min:
	if r_min[n] > 40: 
		r_min_big.append(r_min[n])
	elif r_min[n] < 40: 
		r_min_small.append(r_min[n])
	n=n+1
	
dataK=np.loadtxt('Green.txt', dtype=str, delimiter=',', skiprows=1)
r_sec_snrs=np.float64(dataK[:,4])
r_min_snrs=r_sec_snrs/60

n=0
r_min_snrs_small=[]
r_min_snrs_big=[]
for r in r_min_snrs:
	if r_min_snrs[n] > 40: 
		r_min_snrs_big.append(r_min_snrs[n])
	elif r_min_snrs[n] < 40: 
		r_min_snrs_small.append(r_min_snrs[n])
	n=n+1

dataW=np.loadtxt('wise_hii_V2.txt', dtype=str, delimiter=',', skiprows=1, usecols=(0,1,2,3,4))
r_sec_hii=np.float64(dataW[:,4])
r_min_hii=r_sec_hii/60

n=0
r_min_hii_small=[]
r_min_hii_big=[]
for r in r_min_hii:
	if r_min_hii[n] > 40: 
		r_min_hii_big.append(r_min_hii[n])
	elif r_min_hii[n] < 40: 
		r_min_hii_small.append(r_min_hii[n])
	n=n+1

fig=plt.figure(figsize=(45*cm, 25*cm))
ax1=fig.add_axes([0.1,0.1,0.45,0.45])
#ax4=fig.add_axes([0.55,0.1,0.4,0.4])

sns.histplot(data=r_min_hii, fill=False, color='maroon', label='HII regions', ax=ax1, log_scale=True, stat="count", element='step')
sns.histplot(data=r_sec_snrs, fill=False, color='navy', label='Green SNRs', ax=ax1, log_scale=True, stat="count", element='step')
sns.histplot(data=r_min, fill=True, color='cyan', label='This work', ax=ax1, log_scale=True, stat="count", element='step')

ax1.legend()
ax1.set(xlabel="Angular size / arcmin", ylabel="Normalised count")

#ax4.legend()
#ax4.set(xlabel="Angular size / arcmin", ylabel="Count")

fig.savefig('Size_hii.png', bbox_inches='tight')
fig.savefig('Size_hii.pdf', bbox_inches='tight')

##########################################################################################################
#Faccio un confronto fra i miei candidati e i candidati di GLEAM (vedere se trovo altri paper con candidati MWA) e Green. Plotto la flux density (tratteggoata per i miei perche lower limit) e le sizes
dataG=np.loadtxt('natasha.txt', dtype=str, delimiter=',', usecols=(0,1,2,3,4))
r_min_gleam=np.float64(dataG[:,2])

FluxS=np.float64(dataS[:,3])
Flux1GHz=np.float64(dataK[:,5])
index=np.float64(dataK[:,6])
FluxG=np.float64(dataG[:,3])

FluxK=[]
for n in range(0,224):
	FluxK.append(Flux1GHz[n]*pow(200.0/1000.0, -index[n]))

figurea=plt.figure(figsize=(45*cm, 30*cm))
figureb=plt.figure(figsize=(45*cm, 30*cm))
axa=figurea.add_axes([0.1,0.1,0.5,0.5])
axb=figureb.add_axes([0.1,0.1,0.5,0.5])

bins=np.linspace(-1, 4.2, num=20, endpoint=True)
sns.histplot(data=FluxK, fill=False, color='navy', label='Green SNRs', ax=axa, log_scale=True, stat="probability", element='step', bins=bins)
sns.histplot(data=FluxG, fill=False, color='orange', label='GLEAM', ax=axa, log_scale=True, stat="probability", element='step', bins=bins)
sns.histplot(data=FluxS, fill=True, color='cyan', label='This work', ax=axa, log_scale=True, stat="probability", element='step', bins=bins, linestyle='--', edgecolor='darkcyan')

bins=np.linspace(0, 2.7, num=20, endpoint=True)
sns.histplot(data=r_sec_snrs, fill=False, color='navy', label='Green SNRs', ax=axb, log_scale=True, stat="probability", element='step', bins=bins)
sns.histplot(data=r_min_gleam, fill=False, color='orange', label='GLEAM', ax=axb, log_scale=True, stat="probability", element='step', bins=bins)
sns.histplot(data=r_min*2, fill=True, color='cyan', label='This work', ax=axb, log_scale=True, stat="probability", element='step', bins=bins)

axa.legend()
axb.legend()
axa.set(xlabel="Flux density (200 MHz) / Jy", ylabel="Normalised count")
axb.set(xlabel="Angular size / arcmin", ylabel="Normalised count")

figurea.savefig('Flux.png', bbox_inches='tight')
figurea.savefig('Flux.pdf', bbox_inches='tight')
figureb.savefig('Dimension.png', bbox_inches='tight')
figureb.savefig('Dimension.pdf', bbox_inches='tight')

##########################################################################################################
#Rappresento il piano galattico per vedere come sono distribuite. 
Ra=dataS[:,0]
Dec=dataS[:,1]

n=0
l=[]
b=[]
l1=[]
for line in Ra:
	coord=SkyCoord(Ra[n], Dec[n], unit=(u.hourangle, u.deg), frame='fk5')
	l.append(coord.galactic.l.value)
	b.append(coord.galactic.b.value)
	l1.append(mwrap(np.float64(l[n])))
	n+=1

Ra1=dataK[:,2]
Dec1=dataK[:,3]

n=0
l4=[]
b4=[]
l5=[]
for line in Ra1:
	coord=SkyCoord(Ra1[n], Dec1[n], unit=(u.hourangle, u.deg), frame='fk5')
	l4.append(coord.galactic.l.value)
	b4.append(coord.galactic.b.value)
	l5.append(mwrap(np.float64(l4[n])))
	n+=1
	
Ra2=dataW[:,2]
Dec2=dataW[:,3]

n=0
l2=[]
b2=[]
l3=[]
for line in Ra2:
	coord=SkyCoord(Ra2[n], Dec2[n], unit=(u.deg, u.deg), frame='galactic')
	l2.append(coord.galactic.l.value)
	b2.append(coord.galactic.b.value)
	l3.append(mwrap(np.float64(l2[n])))
	n+=1

fig=plt.figure(figsize=(45*cm, 25*cm))
ax2=fig.add_axes([0.1, 0.1, 0.45, 0.45])

bins=np.linspace(-250, 250, num=20, endpoint=True)
sns.histplot(data=l3, fill=False, color='maroon', label='HII regions', ax=ax2, stat="probability", common_norm=False, element='step', bins=bins)
sns.histplot(data=l5, fill=False, color='navy', label='Green SNRs', ax=ax2, stat="probability", common_norm=False, element='step', bins=bins)
sns.histplot(data=l1, fill=True, color='cyan', label='This work', ax=ax2, stat="probability", common_norm=False, element='step', bins=bins)

ax2.set(xlabel="Galactic longitude / deg", ylabel="Normalised count")
ax2.legend()
ax2.set_xlim(-220, 220)
ax2.invert_xaxis()

fig.savefig('Distribution_hii.png', bbox_inches='tight')
fig.savefig('Distribution_hii.pdf', bbox_inches='tight')

fig=plt.figure(figsize=(45*cm, 25*cm))
ax5=fig.add_axes([0.1, 0.1, 0.45, 0.45])

bins=np.linspace(-10, +10, num=40, endpoint=True)
sns.histplot(data=b2, fill=False, color='maroon', label='HII regions', ax=ax5, common_norm=False, stat="probability", element='step', bins=bins)
sns.histplot(data=b4, fill=False, color='navy', label='Green SNRs', ax=ax5, common_norm=False, stat="probability", element='step', bins=bins)
sns.histplot(data=b, fill=True, color='cyan', label='This work', ax=ax5, common_norm=False, stat="probability", element='step', bins=bins)

ax5.set(xlabel="Galactic latitude / deg", ylabel="Normalised count")
ax5.legend()
ax5.set_xlim(-6.0, 6.0)

above=0
below=0
n=0
for l in b:
	if b[n]>0:
		above=above+1
	else:
		below=below+1
	n+=1
	
print(above, below)

fig.savefig('Latitude_hii.png', bbox_inches='tight')
fig.savefig('Latitude_hii.pdf', bbox_inches='tight')


###########################################################################################################
#Rappresento il piano galattico inserendo sia i miei candidati che quelli gia' noti in modo da poter effettuare un confronto.
array=np.loadtxt('Silvia_candidates.txt', dtype=str, delimiter=',', skiprows=1)
RaS=array[:,0]
DecS=array[:,1]

n=0
l=[]
b=[]
l1=[]
b1=[]
for line in RaS:
	coord=SkyCoord(Ra[n], Dec[n], unit=(u.hourangle, u.deg), frame='fk5')
	l.append(coord.galactic.l.value)
	b.append(coord.galactic.b.value)
	l1.append(mwrap(np.float64(l[n])))
	b1.append(posv(np.float64(b[n])))
	n+=1
	
n=0
l_snrs=[]
b_snrs=[]
l1_snrs=[]
b1_snrs=[]
for line in Ra1:
	coord=SkyCoord(Ra1[n], Dec1[n], unit=(u.hourangle, u.deg), frame='fk5')
	l_snrs.append(coord.galactic.l.value)
	b_snrs.append(coord.galactic.b.value)
	l1_snrs.append(mwrap(np.float64(l_snrs[n])))
	b1_snrs.append(posv(np.float64(b_snrs[n])))
	n+=1

n=0
l_hii=[]
b_hii=[]
l1_hii=[]
b1_hii=[]
for line in Ra2:
	coord=SkyCoord(Ra2[n], Dec2[n], unit=(u.deg, u.deg), frame='galactic')
	l_hii.append(coord.galactic.l.value)
	b_hii.append(coord.galactic.b.value)
	l1_hii.append(mwrap(np.float64(l_hii[n])))
	b1_hii.append(posv(np.float64(b_hii[n])))
	n+=1
		
fig=plt.figure(figsize=(40*cm, 12*cm))		
ax3=fig.add_axes([0.1, 0.1, 0.45, 0.45])

ax3.grid(True)
ax3.scatter(l1_hii, b1_hii, color='maroon', alpha=0.6, label='HII regions', marker='*')
ax3.scatter(l1_snrs, b_snrs, color='navy', alpha=0.8, label='Green SNRs', marker='+')
ax3.scatter(l1, b, color='cyan', alpha=0.6, label='This work')
ax3.set_xlabel('Galactic longitude \n l / deg')
ax3.set_ylabel('Galactic latitude \n b / deg')
ax3.invert_xaxis()
ax3.legend()
ax3.set_ylim(-5.0, 5.0)

fig.savefig('Gal_hii.png', bbox_inches='tight')
fig.savefig('Gal_hii.pdf', bbox_inches='tight')