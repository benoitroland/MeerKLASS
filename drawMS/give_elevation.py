
import sys
import numpy
import pyrap.measures as pm
import pyrap.tables as pt
import pyrap.quanta as qa
import pylab
import ephem

class WritableObject:
    def __init__(self):
        self.content = []
    def write(self, string):
        self.content.append(string)
 


def give_elevation(timevals,ra,dec,xarray,yarray,zarray):

    me = pm.measures()
    position = me.position('wgs84',xarray,yarray,zarray)
    me.doframe(position)
    ra_qa  = qa.quantity( ra, 'rad' )
    dec_qa  = qa.quantity( dec, 'rad' )
    pointing = me.direction('j2000',ra_qa,dec_qa)

    el = numpy.zeros(len(timevals))
    az = numpy.zeros(len(timevals))

    #stdout_backup = sys.stderr
    #sys.stderr = open("/dev/null", "w")

    #original_stderr = sys.stderr
    #f = open("/tmp/stderr.txt", "w")
    #sys.stderr = f

    for i in range(len(timevals)):
        tt = qa.quantity(timevals[i],'s')
        tt1 = me.epoch('utc',tt)
        me.doframe(tt1)
        a=me.measure(pointing,'azel')
        el[i]=a['m1']['value']
        az[i]=a['m0']['value']

    #sys.stderr = original_stderr
    #f.close() 
    #sys.stderr = stdout_backup

    return az,el
