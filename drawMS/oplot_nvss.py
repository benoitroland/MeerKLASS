import numpy as np
from astropy.io import fits

def give_nvss(rarad,decrad,radius,smin):
    #rarad,decrad,radius in radians

    cat=fits.open('/home/tasse/drawMS/CATALOG41.FIT')[1]
    cat.data.dtype=[('RA', '>f8'), ('DEC', '>f8'), ('PEAKINT', '>f4'), ('MAJORAX', '>f4'), ('MINORAX', '>f4'), ('POSANGLE', '>f4'), ('QCENTER', '>f4'), ('UCENTER', '>f4'), ('PFLUX', '>f4'), ('IRMS', '>f4'), ('POLRMS', '>f4'), ('RESRMS', '>f4'), ('RESPEAK', '>f4'), ('RESFLUX', '>f4'), ('CENTERX', '>f4'), ('CENTERY', '>f4'), ('FIELD', '|S8'), ('JDPROCESSED', '>i4')]

    ind=np.where(cat.data.PEAKINT > smin)[0]
    cat=cat.data[ind]

    ra=cat.RA*np.pi/180.
    dec=cat.DEC*np.pi/180.

    rac=rarad
    decc=decrad

    ra0,ra1=rac-radius,rac+radius
    dec0,dec1=decc-radius,decc+radius
    pi2=np.pi/2.
    cosd=np.cos(pi2-dec)*np.cos(pi2-decc)+np.sin(pi2-dec)*np.sin(pi2-decc)*np.cos(ra-rac)
    d=np.arccos(cosd)

    ind=np.where((d<radius))[0]
    cat=cat[ind]

    

    return cat.RA*np.pi/180.,cat.DEC*np.pi/180.,cat.PEAKINT
