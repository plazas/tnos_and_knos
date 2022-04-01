import math
import astropy.io.fits as pf
import numpy as np
import matplotlib.pyplot as plt

# from source code of http://www.lco.cl/?epkb_post_type_1=direct-ccd-exposure-time-calculator-1-0
#  <SCRIPT language="Javascript" src="http://www.lco.cl/wp-content/uploads/data.js"> </SCRIPT>
#  <SCRIPT language="Javascript" src="http://www.lco.cl/wp-content/uploads/functions.js"> </SCRIPT>


class Instrument: 
    def __init__(self,name,scale,rdnoise,dark,pixsize):
        self.name = name		#not used! (only for clearness)
        self.scale = scale
        self.rdnoise = rdnoise	# e-
        self.dark = dark		# e-/s
        self.pixsize = pixsize	# microns


class Filter:
    def __init__(self,name,extinction,mag,sky0,sky1,sky2):
        self.name = name		#not used! (only for clearness)
        self.extinction = extinction	#(mag/airmass)
        self.mag = mag		#star mag zero point (for 1 e-/sec @ 1 airmass)
        self.sky0 = sky0		#sky mag (zero order term)
        self.sky1 = sky1		#sky mag (1st order term)
        self.sky2 = sky2		#sky mag (2nd order term)


filter_dict={"g": Filter("g'",0.18,27.65,21.6,-0.0626,-0.00902),
             "r": Filter("r'",0.10,27.65,20.7,-0.0202,-0.00440),
             "i": Filter("i'",0.04,27.45,20.0,-0.0079,-0.00392),
             "z": Filter("z'",0.02,26.90,18.4,-0.0003,-0.00364)}

instrument = Instrument ("IMACS f/4",7.41,5.0,0.0,15.0)

def etc_imacs_f4 (magnitude = 24.5, phase=7, airmass=1,
                  binning = 1, snr = 10, seeing = 0.6, filter="r"):

    phase2 = phase*phase

    filter = filter_dict[filter]

    sky     = filter.sky0 + filter.sky1 * phase + filter.sky2 * phase2
    scount  = math.pow(10,(0.4*(1-airmass)*filter.extinction)) #star count rate
    pixel   = (instrument.pixsize/1000) * instrument.scale * binning
    pixel2  = pixel*pixel
    aux1    = 0.4*(filter.mag-sky)
    aux2    = math.pow(10,aux1)     
    b       = scount*pixel2*aux2	 	#sky rate per pixel
    aux3    = (seeing/pixel)
    d       = math.ceil(1.4*aux3*aux3)   	#npix

    if (d<9):
        d = 9

    a = scount*math.pow(10,(0.4*(filter.mag-magnitude)))
    c = instrument.dark
    e = instrument.rdnoise*instrument.rdnoise
    f = snr 
    g =  (a*a) / (f*f) 

    aux4 = math.sqrt ( (math.pow((-a-b*d-c*d),2))+ (4*g*d*e))
    aux5 = a + (b*d) + (c*d)

    aux6 = ((f*f)/(2*(a*a)))*(aux5+aux4)
    time = round(aux6)
    tscounts = a*aux6
    tskypixel = b*aux6
    starnoise = math.sqrt (tscounts)
    skynoise = math.sqrt (tskypixel*d)
    ccdnoise = math.sqrt (d*e)

    return time, tscounts, tskypixel, starnoise, skynoise, ccdnoise



# Des Y6 TNO's
file_name = "/Users/plazas/Documents/TNO_KPO/y6_res.fits"

data = pf.open(file_name)[1].data

mask_detached = (data['CLASS'] == 'Detached')
mask_scattered = (data['CLASS'] == 'Scattering')
mask_classical = (data['CLASS'] == 'Classical')
mask_resonant = (data['CLASS'] == 'Resonant')

detached = data[mask_detached]
scattered = data[mask_scattered]
classical = data[mask_classical]
resonant = data[mask_resonant]

#detached


test = etc_imacs_f4 (magnitude = 24, phase=0, airmass=1.2, snr = 20, seeing = 0.6, filter="z")
print (test)

stop


def total_time_per_class (data_class):
    mag_g, mag_r, mag_z = data_class["m_g"], data_class["m_r"], data_class["m_z"]
    mag_i = data_class["m_i"]
    phase_vec = [0,7,14]
    snr_vec = [20, 50]
    seeing_vec  = [0.6, 0.8]
    airmass_vec = [1.0, 1.4, 2.0]
    for airmass in airmass_vec:
        for phase in phase_vec:
            for snr in snr_vec:
                for seeing in seeing_vec:
                    total_time=0
                    for g in mag_g:
                        time, _, _, _, _, _ = etc_imacs_f4 (magnitude = g, phase=phase, airmass=airmass, snr = snr, seeing = seeing, filter="g")
                        total_time+=time
                    for r in mag_r:
                        time, _, _, _, _, _ = etc_imacs_f4 (magnitude = r, phase=phase, airmass=airmass, snr = snr, seeing = seeing, filter="r")
                        total_time+=time
                    for i in mag_i:
                        time, _, _, _, _, _ = etc_imacs_f4 (magnitude = i, phase=phase, airmass=airmass, snr = snr, seeing = seeing, filter="i")
                        total_time+=time
                    for z in mag_z:
                        time, _, _, _, _, _ = etc_imacs_f4 (magnitude = z, phase=phase, airmass=airmass, snr = snr, seeing = seeing, filter="z")
                        total_time+=time
                    total_time/=3600
                    print (f"Phase: {phase}, Airmass: {airmass}, SNR: {snr},  Seeing: {seeing}, Total time (g+r+i+z): {total_time:.2f} hours")


print ("Detached", len(detached))
total_time_per_class (detached)
print (" ")
print ("Scattered", len(scattered))
total_time_per_class (scattered)

