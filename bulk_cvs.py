import matplotlib
matplotlib.use('Agg')
import tessreduce as tr
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import lightkurve as lk
from astropy.coordinates import SkyCoord
from astropy import units as u

import os
dirname = os.path.dirname(__file__)

#where we're going we dont need warnings!!
import warnings
warnings.filterwarnings("ignore")

file = os.path.join(dirname,'/data/cataclysmic_variables.csv')
cvs = pd.read_csv('./data/cataclysmic_variables.csv')

# don't want to deal with the crowded Tuc, Pav, or Sgr zones for now
ind = (cvs['GCVS'].values == 'Tuc      ') | (cvs['GCVS'].values == 'Pav      ') | (cvs['GCVS'].values == 'Sgr      ')
cvs = cvs.iloc[~ind]


for j in range(len(cvs)):
	cv = cvs.iloc[j]
	print('NAME: ',cv['Names'])
	ra = cv['RAJ2000']
	dec = cv['DEJ2000']

	c = SkyCoord(ra=float(ra)*u.degree, dec=float(dec) *u.degree, frame='icrs')

	tess = lk.search_tesscut(c,sector=None)
	try:
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
					tpf = t.download(cutout_size=50)
					#aper_b18 = np.zeros(tpf.shape[1:], dtype=bool)
					#aper_b18[44:48, 44:47] = True
					res = tr.Quick_reduce(tpf,calibrate=False)#,aper=aper_b18)
					lcs += [res['lc']]
					err += [res['err']]
					zps += [res['zp']]
					sectors += [tpf.sector]
					try:
						trends1 += [tr.Remove_stellar_variability(lcs[-1],err[-1],variable=True)]
						trends2 +=  [tr.Remove_stellar_variability(lcs[-1],err[-1],variable=False)]
					except:
						print('trend error in {} sector {}'.format(cv['Names'],tpf.sector))
						filler = np.nan * np.ones(len(err[-1]))
						trends1 += [filler]
						trends2 += [filler]
					

				name = cv['Names']
				plt.figure(figsize=(6.5,8))
				plt.subplot(311)
				plt.title(name)
				for i in range(len(lcs)):
					plt.plot(lcs[i][0],lcs[i][1],label='S ' + str(sectors[i]))
				plt.legend()
				plt.ylabel('Counts')

				plt.subplot(312)
				plt.title('trend method 1')
				for i in range(len(lcs)):
					#plt.fill_between(lcs[i][0],lcs[i][1]-trends1[i]-err[i],lcs[i][1]-trends1[i]+err[i],alpha=.5)
					plt.plot(lcs[i][0],lcs[i][1]-trends1[i])
				plt.ylabel('Counts')
				
				plt.subplot(313)
				plt.title('trend method 2')
				for i in range(len(lcs)):
					#plt.fill_between(lcs[i][0],lcs[i][1]-trends2[i]-err[i],lcs[i][1]-trends2[i]+err[i],alpha=.5)
					plt.plot(lcs[i][0],lcs[i][1]-trends2[i])
				plt.ylabel('Counts')
				plt.xlabel('MJD')
				plt.tight_layout()

				savename = name.replace('/',' ').replace(' ','_')
				plt.savefig('./figs/{}.pdf'.format(savename))

				# save to cvs
				mjd = lcs[0][0].copy()
				flux = lcs[0][1].copy()
				e = err[0].copy()
				t1 = trends1[0].copy()
				t2 = trends2[0].copy()
				z = np.ones(len(lcs[0][0])) * zps[0][0]
				s = np.ones(len(lcs[0][0])) * sectors[0]
				for i in range(len(lcs)-1):
					i += 1
					mjd = np.append(mjd,lcs[i][0])
					flux = np.append(flux,lcs[i][1])
					e = np.append(e,err[i])
					t1 = np.append(t1,trends1[i])
					t2 = np.append(t2,trends2[i])

					zz = np.ones(len(lcs[i][0])) * zps[i][0]
					ss = np.ones(len(lcs[i][0])) * sectors[i]
					z = np.append(z,zz)
					s = np.append(s,ss)
				df = pd.DataFrame(columns=['mjd','flux','err','trend1','trend2','zp','sector'])
				df['mjd'] = mjd
				df['flux'] = flux
				df['err'] = e
				df['trend1'] = t1
				df['trend2'] = t2
				df['zp'] = z
				df['sector'] = s
				try:
					t1events = tr.Event_isolation(flux-t1,err=e,sig=3)
					for i in range(len(t1events)):
						df['t1event' + str(i+1)] = t1events[i]
				except:
					pass
				try:
					t2events = tr.Event_isolation(flux-t2,err=e,sig=3)
					for i in range(len(t2events)):
						df['t2event' + str(i+1)] = t2events[i]
				except:
					pass
				
				df.to_csv('./lcs/{}.csv'.format(savename),index=False)

				print('finished {}'.format(name))
	except:
		pass