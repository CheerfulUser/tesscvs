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
	j += 4
    cv = cvs.iloc[j]
    print('NAME: ',cv['Names'])
    ra = cv['RAJ2000']
    dec = cv['DEJ2000']

    obs = tr.spacetime_lookup(ra=ra, dec=dec)
    lcs = []
    t1 = []
    t2 = []
    sectors = []
    try:
	    for ob in obs:
	        try:
	            tt = tr.tessreduce(obs_list=ob,reduce=True)
	            tt.to_flux()
	            lcs += [tt.lc]
	            sectors += [tt.tpf.sector]
	        except:
	            print('Failed for ', ob)


	    name = cv['Names']
	    print('MAKE FIGURE')
	    plt.figure(figsize=(6.5,8))
	    plt.title(name)
	    for i in range(len(lcs)):
	        plt.plot(lcs[i][0],lcs[i][1],label='S ' + str(sectors[i]))
	    plt.legend()
	    plt.ylabel('mJy')
	    plt.xlabel('MJD')
	    plt.tight_layout()

	    savename = name.replace('/',' ').replace(' ','_')
	    plt.savefig('./figs/{}.pdf'.format(savename))

	    # save to cvs
	    print('SAVE TO CSV')
	    mjd = lcs[0][0].copy()
	    flux = lcs[0][1].copy()
	    e = lcs[0][2].copy()
	    s = np.ones(len(lcs[0][0])) * sectors[0]
	    for i in range(len(lcs)-1):
	        i += 1
	        mjd = np.append(mjd,lcs[i][0])
	        flux = np.append(flux,lcs[i][1])
	        e = np.append(e,lcs[i][2])

	        ss = np.ones(len(lcs[i][0])) * sectors[i]
	        s = np.append(s,ss)
	    df = pd.DataFrame(columns=['mjd','flux','err','trend1','trend2','sector'])
	    df['mjd'] = mjd
	    df['flux'] = flux
	    df['err'] = e
	    df['sector'] = s
	    
	    df.to_csv('./lcs/{}.csv'.format(savename),index=False)

	    print('finished {}'.format(name))
    except:
    	print('eh')