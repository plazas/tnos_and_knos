from astroplan.plots import plot_airmass 
from astroplan import Observer 

from astroplan import FixedTarget
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time 

import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

import matplotlib.pyplot as plt

# Las Campanas Observatory
#lco = Observer.at_site("LCO", timezone="America/Santiago") 
apo = Observer.at_site("APO", timezone="US/Mountain") 

observatory = apo


ra_max = np.array([28.77793, 15.54044, 28.31347, 9.53626, 19.42361, 329.83334,
49.50065, 59.44437, 59.77828, 54.05925, 50.04861, 43.06966])

ra_min = np.array([26.42629, 13.1782, 25.33274, 5.91576, 16.81335, 326.57642, 45.77923, 56.83081, 
56.93357, 52.33194, 48.64039, 41.36328 ]) 

dec_max = np.array([5.35481, 2.44644, 7.32161, -40.52014, 6.84466, -44.20068, -11.04898, -16.94137,
-22.91617, -17.72259, -9.22384, 16.47803])

dec_min = np.array([4.42221, 1.56007, 6.17567, -42.04645, 5.69949, -45.27301, -12.8212, -18.13655, 
-24.6034, -18.65246, -9.91562, -17.39446])

target_names=[r'2007 TB$_{418}$', r'2010 RF$_{188}$', r'2013 RJ$_{124}$', r'2013 RV$_{124}$', r'2013 SA$_{100}$', r'2013 SR$_{102}$', r'2014 RV$_{86}$', r'2014 SQ$_{403}$', r'2014 ST$_{373}$', r'2014 US$_{277}$', r'2014 UZ$_{224}$', r'2014 YC$_{92}$']

print (target_names)
#stop

ra_mid = 0.5*(ra_max + ra_min)
dec_mid = 0.5*(dec_max + dec_min)



pp = PdfPages ("./airmassTargets-APO.pdf")

times = []
for i in range(1,13):
    if i < 10:
        times += [f'2022-0{i}-01']
    else:
        times += [f'2022-{i}-01']

print ("Times: ", times)

months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
print ("names: ", target_names)


for i, time in enumerate(times):
    fig = plt.figure()
    #ax = fig.add_subplot(111)
    if i in [3,4,5,6]:
        loc = 'upper right'
    else:
        loc = 'upper left'
    for j, (ra_midpoint, dec_midpoint) in enumerate(zip (ra_mid, dec_mid)):
        ax = fig.add_subplot(111)
        coord = SkyCoord(ra=ra_midpoint*u.deg, dec=dec_midpoint*u.deg, frame='icrs') 
        target = FixedTarget(name=target_names[j], coord=coord)
        time = Time(time)
        print (i, j, target_names[j])
        plot_airmass(target, observatory, time, use_local_tz=True, brightness_shading=True, altitude_yaxis=True)
        ax.legend(shadow=True, loc=loc)
    plt.tight_layout()
    fig.suptitle(months[i])
    pp.savefig()
pp.close()
