from pygsm import GlobalSkyModel,GlobalSkyModel2016
import numpy as np
from astropy.coordinates import SkyCoord
import healpy
import pint.models as models
from astropy import units as u
import astropy
from optparse import OptionParser
import os,sys

def SkyCoord_Tsky(source, freq=350*u.MHz, model='2008'):
    """
    T=SkyCoord_Tsky(source, freq=350*u.MHz)
    takes an astropy SkyCoord object and computes the sky temperature for that position
    using pyGSM
    scaled to the given frequency

    If frequency has no units, MHz is assumed

    model='2008' or '2016' for GSM2008 (Oliveira-Costa et al., 2008) or GSM2016 (Zheng et al., 2016)

    returns sky temperature in K
    """
    assert str(model) in ['2008','2016']

    if not isinstance(freq, astropy.units.quantity.Quantity):
        # assume MHz
        freq=freq*u.MHz

    if str(model)=='2008':
        gsm = GlobalSkyModel()
    elif str(model)=='2016':
        gsm=GlobalSkyModel2016()
    map=gsm.generate(freq.to(u.MHz).value)
    T=healpy.pixelfunc.get_interp_val(map,
                                      source.galactic.l.value,
                                      source.galactic.b.value,
                                      lonlat=True)
    return T*u.K



def pulsar_Tsky(parfile, freq=350*u.MHz, model='2008'):
    """
    T=pulsar_Tsky(parfile, freq=350*u.MHz)
    takes a parfile and computes the sky temperature for that pulsar
    using pyGSM
    scaled to the given frequency

    If frequency has no units, MHz is assumed

    model='2008' or '2016' for GSM2008 (Oliveira-Costa et al., 2008) or GSM2016 (Zheng et al., 2016)

    returns sky temperature in K
    """

    m=models.get_model(parfile)
    try:
        psr=SkyCoord(m.RAJ.value,m.DECJ.value,unit='deg')
    except:
        try:
            psr=SkyCoord(m.ELONG.value*u.deg,m.ELAT.value*u.deg,frame='pulsarecliptic')
        except:
            raise KeyError,'Cannot find RAJ,DECJ or ELONG,ELAT in:\n%s' % m.as_parfile()
    T=source_Tsky(psr, freq=freq, model=model)
    return T

def main():
    
    usage="Usage: %prog -f <freq> -m <model> <parfile> [<parfile2>]\n"
    usage+='\tEvaluates Global Sky Model temperature for the pulsar specified in <parfile>'
    parser = OptionParser(usage=usage)
    parser.add_option('-f','--freq',dest='freq',default=350,type='float',
                      help='Frequency to compute sky temperature (MHz) [default=%default]')
    parser.add_option('-m','--model',dest='model',default='2008',choices=['2008','2016'],
                      help='GSM model, 2008 or 2016 [default=%default]')
    
    (options, args) = parser.parse_args()
    for parfile in args:
        if not os.path.exists(parfile):
            raise IOError, 'No such file %s' % parfile
        T=pulsar_Tsky(parfile, freq=options.freq, model=options.model)
        m=models.get_model(parfile)
        print 'PSR %s: %.1f K (at %d MHz)' % (m.PSR.value,
                                              T.value,
                                              options.freq)
        
######################################################################

if __name__=="__main__":
    main()

