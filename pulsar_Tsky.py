from pygsm import GlobalSkyModel,GlobalSkyModel2016
import numpy as np
from astropy.coordinates import SkyCoord
import healpy
import pint.models as models
from astropy import units as u
import astropy
from optparse import OptionParser
import os,sys

class GSModel:
    """
    class interface to GSM model
    so that structures/data stay in memory
    gsm=GSModel(freq=freq, model=model)

    to use on astropy SkyCoord:
    gsm.SkyCoord_Tsky(source)

    to use on pulsar parfile:
    gsm.pulsar_Tsky(parfile)
    
    """

    def __init__(self, freq=350*u.MHz, model='2008'):
        assert str(model) in ['2008','2016']

        if not isinstance(freq, astropy.units.quantity.Quantity):
            # assume MHz
            freq=freq*u.MHz

        if str(model)=='2008':
            self.gsm = GlobalSkyModel()
        elif str(model)=='2016':
            self.gsm=GlobalSkyModel2016()
        self.map=self.gsm.generate(freq.to(u.MHz).value)
        self.freq=freq
        self.model=model

    def SkyCoord_Tsky(self, source):
        T=healpy.pixelfunc.get_interp_val(self.map,
                                          source.galactic.l.value,
                                          source.galactic.b.value,
                                          lonlat=True)
        return T*u.K

    def pulsar_Tsky(self, parfile):
        m=models.get_model(parfile)
        try:
            psr=SkyCoord(m.RAJ.value,m.DECJ.value,unit='deg')
        except:
            try:
                psr=SkyCoord(m.ELONG.value*u.deg,m.ELAT.value*u.deg,frame='pulsarecliptic')
            except:
                raise KeyError,'Cannot find RAJ,DECJ or ELONG,ELAT in:\n%s' % m.as_parfile()
        return self.SkyCoord_Tsky(psr)


##############################
# convenience functions
# provide a single interface
##############################
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

    gsm=GSModel(freq=freq, model=model)
    return gsm.SkyCoord_Tsky(source)


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
    gsm=GSModel(freq=freq, model=model)
    return gsm.pulsar_Tsky(parfile)

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

