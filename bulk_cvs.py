import tessreduce as tr
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import lightkurve as lk
from astropy.coordinates import SkyCoord
from astropy import units as u

import os
dirname = os.path.dirname(__file__)

file = os.path.join(dirname,'/data/cataclysmic_variables.csv')
cvs = pd.read_csv(file)

# don't want to deal with the crowded tuc, Pav, or Sgr zones for now
ind = (cvs['GCVS'].values == 'Tuc') | (cvs['GCVS'].values == 'Pav') | (cvs['GCVS'].values == 'Sgr')
cvs = cvs.iloc[~ind]


for i in range(len(cvs)):
	cv = cvs.iloc[i]

	ra = cv['RAJ2000']
	dec = cv['DEJ2000']

	c = SkyCoord(ra=float(ra)*u.degree, dec=float(dec) *u.degree, frame='icrs')

	tess = lk.search_tesscut(c,sector=None)

	if len(tess) > 0:
		lcs = []
		zps = []
		err = []
		sectors = []
		trends1 = []
		trends2 = []

		if len(tess) > 1:
		    tpfs = []
		    for t in tess:
		        tpf = t.download(cutout_size=90)
		        aper_b18 = np.zeros(tpf.shape[1:], dtype=bool)
		        aper_b18[44:48, 44:47] = True
		        res = tr.Quick_reduce(tpf,aper=aper_b18)
		        lcs += [res['lc']]
		        err += [res['err']]
		        zps += [res['zp']]
		        sectors += [tpf.sector]
		        trends1 += [tr.Remove_stellar_variability(lcs[i],err[i],variable=True)]
		        trends2 +=  [tr.Remove_stellar_variability(lcs[i],err[i],variable=False)]

	        name = cv['Names']
	        plt.figure(figsize=(6.5,8))
			plt.subplot(311)
			plt.title(name)
			for i in range(len(lcs)):
			    plt.plot(lcs[i][0],lcs[i][1],label='S ' + sectors[i])
			plt.legend()

			plt.subplot(312)
			plt.title('trend method 1')
			for i in range(len(lcs)):
			    plt.fill_between(lcs[i][0],lcs[i][1]-trends1[i]-err[i],lcs[i][1]-trends1[i]+err[i],alpha=.5)
			    plt.plot(lcs[i][0],lcs[i][1]-trends1[i])
		    plt.subplot(313)
		    plt.title('trend method 2')
			for i in range(len(lcs)):
			    plt.fill_between(lcs[i][0],lcs[i][1]-trends2[i]-err[i],lcs[i][1]-trends2[i]+err[i],alpha=.5)
			    plt.plot(lcs[i][0],lcs[i][1]-trends2[i])
		    savename = name.replace('/',' ').replace(' ','_')
		    plt.savefig('./figs/{}.pdf'.format(savename))

		    # save to cvs
	    	mjd = lcs[0][0].copy()
			flux = lcs[0][1].copy()
			e = err[0].copy()
			t1 = trends1[0].copy()
			t2 = trends2[0].copy()
			z = np.ones(len(lcs[0][0])) * zps[0]
			s = np.ones(len(lcs[0][0])) * sectors[0]
			for i in range(len(lcs)-1):
			    i += 1
			    mjd = np.append(mjd,lcs[i][0])
			    flux = np.append(flux,lcs[i][1])
			    e = np.append(e,err[i])
			    t1 = np.append(t1,trends1[i])
			    t2 = np.append(t2,trends2[i])
			    
			    zz = np.ones(len(lcs[i][0])) * zps[i]
			    ss = np.ones(len(lcs[i][0])) * sectors[i]
			    z = np.append(z,zz)
			    s = np.append(s,ss)
			df['mjd'] = mjd
			df['flux'] = flux
			df['err'] = e
			df['trend1'] = t1
			df['trend2'] = t2
			df['zp'] = z
			df['sector'] = s

			df.to_csv('./lcs/{}.csv'.format(savename),index=False)

			print('finished {}'.format(name))




		        