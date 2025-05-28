import numpy as np

def give_str_time(t):
    #str('%03d'%20)
    sign=""
    if t<0.:
        sign="-"
        t=np.abs(t)
    t=float(t)
    h=int(t)
    f=1000.
    m=int(round((f*60*(t-float(h))))/f)
    s=abs(round(f*3600.*(t-h-m/60.))/f)
    sout=sign+str('%02i'%h)+':'+str('%02i'%m)+':'+str('%05.2f'%s)
    return sout
