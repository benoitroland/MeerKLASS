
import sys
import pylab
import numpy as np
import scipy.fftpack
import time as time_mod
import oplot_nvss
import os
try:
    import ephem
    from pyrap.tables import table
except ImportError:
    print('you have to type "use LofIm" and "use Pythonlibs"')
    exit()
    
# # ##############################
# # Catch numpy warning
# np.seterr(all='raise')
# import warnings
# warnings.catch_warnings()
# warnings.filterwarnings('error')
# # ##############################


def init(MSname, DoReadMS, ColName='DATA'):
    global na,nbl,reffreq, wavelength, wavelength_chan, rarad,decrad, radeg, decdeg, vis, uvw_all, time_all, timewindow, MaxNfoundPeaks, ntheta_solve, SlimNVSS, bl_selection, fill_factor
    global blcut_min, blcut_max, wmax, oldmax, ntimes, mrot, vl, vm, vn, vl2, vm2, vn2, R, snrcut, mask_freq, freq_cut, Nchan, MultiPlotTime, RadPlotNVSS, ranvss,decnvss,snvss, A0,A1
    global table_all, work_colname, antenna_list,progress_str, time_slots_all, delta_time_bins, time_window_bins, xarray, yarray, zarray, time_range_min, time_range_max
    progress_str="|"

    ta=table(MSname+'/ANTENNA',ack=False)
    
    na=ta.getcol('POSITION').shape[0]
    nbl=int((na*(na-1))/2+na)
    ta.close()
    ta=table(MSname+'/SPECTRAL_WINDOW/',ack=False)
    reffreq=ta.getcol('REF_FREQUENCY')[0]
    wavelength=3.e8/reffreq
    wavelength_chan=3.e8/ta.getcol('CHAN_FREQ')[0]
    Nchan=wavelength_chan.shape[0]
    ta.close()
    ta=table(MSname+'/FIELD/',ack=False)
    rarad,decrad=ta.getcol('DELAY_DIR')[0][0]
    rarad,decrad=ta.getcol('PHASE_DIR')[0][0]
    ta.close()
    radeg=rarad*180./np.pi
    decdeg=decrad*180./np.pi
    #print("Pointing center: ",deg2hmsdms.give_str_time(radeg/15.),deg2hmsdms.give_str_time(decdeg))
    
    # extract visibilities, coordinates ...,color=colors[col]
    work_colname=ColName
    if DoReadMS:
        print("  Reading MS ...")
        table_all=table(MSname,ack=False)#,readonly=False)
        tb = table(table_all.getkeyword('ANTENNA'),ack=False)
        pos = tb.getcol('POSITION')

        import pyrap.quanta as qa
        xarray = qa.quantity(pos[0,0],'m')
        yarray = qa.quantity(pos[0,1],'m')
        zarray = qa.quantity(pos[0,2],'m')

        vis=table_all.getcol(work_colname)
        uvw_all=table_all.getcol('UVW')
        uvw_all/=wavelength_chan[-1]
        A0=table_all.getcol('ANTENNA1')
        A1=table_all.getcol('ANTENNA2')
        time_all=table_all.getcol("TIME")
        time_slots_all=time_all[0::nbl]
        flag=table_all.getcol("FLAG")
        #table_all.close()
        ind=np.where(flag==True)[0]
        vis[ind]=0.
        ntimes=np.int(np.floor(uvw_all.shape[0]/float(nbl)))
        print("  ... Done!")
        table_all.close()

    blcut_min=100
    blcut_max=100000
    snrcut=5.
    wmax=1000000.
    oldmax=0.
    delta_time_bins=500
    time_window_bins=500
    mask_freq=5
    freq_cut=-1#1
    MaxNfoundPeaks=7
    MultiPlotTime=False
    ntheta_solve=20
    antenna_list=range(na)
    time_range_min=0.
    time_range_max=1.
    RadPlotNVSS=0.#7.*wavelength/70.
    SlimNVSS=0.6
    fill_factor=1.
    bl_selection=range(nbl)

    if RadPlotNVSS!=0.:
        ranvss,decnvss,snvss=oplot_nvss.give_nvss(rarad,decrad,RadPlotNVSS,SlimNVSS)
    # Declare rotation Matrix
    cos=np.cos
    sin=np.sin
    mrot=np.matrix([[cos(rarad)*cos(decrad), sin(rarad)*cos(decrad),sin(decrad)],[-sin(rarad),cos(rarad),0.],[-cos(rarad)*sin(decrad),-sin(rarad)*sin(decrad),cos(decrad)]]).T
    vm=np.matrix([[0.,0.,1.]]).T
    vl=np.matrix([[0.,1., 0.]]).T
    vn=np.matrix([[1., 0, 0.]]).T
    vl2=mrot*vl
    vm2=mrot*vm
    vn2=mrot*vn
    R=np.matrix([[cos(decrad)*cos(rarad),cos(decrad)*sin(rarad),sin(decrad)]]).T

    if DoReadMS==True:
        print()
        print(" MS PROPERTIES: ")
        print("   - Frequency = ",reffreq/1e6," MHz")
        print("   - Wavelength = ",wavelength," meters")
        print("   - Time bin = ",time_all[nbl+1]-time_all[0]," seconds")
        print("   - Total Integration time = ",(time_all[-1]-time_all[0])/3600.," hours")
        print("   - Number of antenna  = ",na)
        print("   - Number of baseline = ",nbl)
        print("   - Number of channels = ",Nchan)

    # if Nchan>10:
    #     print("   I think you have too many channels, and I don't manage chunks yet.")
    #     sys.exit()




def find_source(DoReadMS,noplot=False):
    #init("/data/scratch/tasse/L36383_SAP000_SB060_uv.MS.dppp",DoReadMS)
    #init("/data/scratch/tasse/L36383_SAP001_SB152_uv.MS.dppp",DoReadMS)
    #init("/media/6B5E-87D0/MS/L34239_SAP001_SB121_uv.MS.dppp.dppp",DoReadMS)
    #init("/media/6B5E-87D0/MS/data/scratch/tasse/L29689cyril.ms",DoReadMS)
    #init("/media/6B5E-87D0/MS/L36383_SAP000_SB102_uv.MS.dppp.dppp",DoReadMS)

    global catalog_lines
    init_catalog_lines()
    ss=list(range(0,ntimes,delta_time_bins))
    off=0.

    if ss[-1]!=(ntimes-1): ss.append(ntimes-1)
    
    if len(ss)<2:
        print("Time bin too long, there was no selected time windows")
        exit()

    if noplot==False:
        if MultiPlotTime:
            for i in range(len(ss)-1):
                fig=pylab.figure(i+1)
                pylab.clf()
        else:
            fig=pylab.figure(1,figsize=(8, 4))
            pylab.clf()
        
        pylab.figure(1)
        pylab.clf()
        pylab.figure(3,figsize=(8,8))
        pylab.clf()


    coutnbl=0

    print("  Calculating the lines ... ")
    toolbar_width = int(na)
    sys.stdout.write("[%s]" % (" " * (toolbar_width-1)))
    sys.stdout.flush()
    sys.stdout.write("\b" * (toolbar_width)) # return to start of line, after '['
    incr=int(nbl/float(toolbar_width))
    tot=0

    index_cat=0
    for bl in bl_selection:#range(0,nbl):
        #print("Process baseline: ", bl,"/",nbl)
        if tot==int(incr):
            tot=0
            sys.stdout.write(progress_str)
            sys.stdout.flush()
        tot+=1

        for i in range(len(ss)-1):
            if np.random.rand(1)[0]>fill_factor: continue

            #delta_time_bins=100
            #time_window_bins=400
            time_bin_start=int(np.max([ss[i]-time_window_bins/2.,0.]))
            time_bin_end=int(np.min([ss[i]+time_window_bins/2.,ntimes-1]))


            x,y,freq_list,Rout,theta_rad,snr=test_find(bl,ss[i],ss[i+1],off)

            if x.shape!=(1,1):
                #print("will plot...")
                coutnbl+=1
                #number of peaks detected in the fft
                for ii in range(len(x)):
                    #print("will plot 2...")
                    #number of radial solutions per peak
                    #print("freq=",freq_list[ii])
                    for iii in range(2):
                        line=0
                        cont=True
                        while cont==True:
                            RR=np.sqrt(x[ii][iii]**2+y[ii][iii]**2)
                            mask,cont=give_mask(np.abs(RR-0.555)>0.0001,line)

                            line+=1
                            #if freq_list[ii]==0.:
                            #    mask[:]=True
                            if cont==True:
                                if noplot==False: 
                                    pylab.figure(3)
                                    pylab.plot(x[ii][iii][mask],y[ii][iii][mask],lw=0.3,ls=':',color='grey',marker=".",ms=.8,mew=.8)#,alpha=np.min(np.array([1.,(snr[ii]-snrcut)/10.])))

                                catalog_lines.line.x[index_cat][mask] = np.float32(x[ii][iii][mask])
                                catalog_lines.line.y[index_cat][mask] = np.float32(y[ii][iii][mask])
                                catalog_lines.line.theta  [index_cat][mask] = np.float32(theta_rad[mask])
                                catalog_lines.line.r      [index_cat][mask] = np.float32(Rout[ii][iii][mask])
                                catalog_lines.freq        [index_cat]           = freq_list[ii]
                                catalog_lines.bln         [index_cat]           = bl
                                catalog_lines.timeslot    [index_cat]           = i
                                index_cat+=1

                                #print("will plot 4...")
                                if noplot==False:
                                    ra,dec=lm2radec(x[ii][iii][mask],y[ii][iii][mask])
                                    if MultiPlotTime:
                                        pylab.figure(i+1)
                                    else:
                                        pylab.figure(1)
                                    plot_polar(ra,dec,col=ii)
    print("\n")
    sys.stdout.flush()


    ind=np.where(catalog_lines.timeslot!=-1.)[0]
    cat=catalog_lines[ind]
    np.save("catalog_lines",cat)
    del(cat,catalog_lines)

    if noplot==False:
        init_plot_lm()
        if MultiPlotTime:
            for i in range(len(ss)-1):
                pylab.figure(i+1)
                init_plot_polar()
        else:
            pylab.figure(1)
            init_plot_polar()

    #pylab.show()
    


    
def test_find(bl,s0,s1,off):
    
    # number of antenna, baselines...
    uvw=uvw_all[bl::nbl]
    NCols=uvw.shape[0]
    time=time_all[bl::nbl]
    vv=vis[bl::nbl,:,0]+vis[bl::nbl,:,3]

    wmean=np.abs(np.mean(uvw[:,2]))

    # Baselines
    B=np.sqrt(uvw[:,1]**2+uvw[:,0]**2)
    Bmean=np.mean(B)
    ant0=A0[bl]
    ant1=A1[bl]
    tbin=(s0+s1)/2.
    condtime=(tbin<time_range_min*NCols)|(tbin>time_range_max*NCols)
   # print
   # print(ant0,ant1,antenna_list,type(antenna_list))
    #print(np.int16([int(ant0),int(ant1)]))
    #print(np.int16(antenna_list[:]))

    if (wmean>wmax)|(uvw[0,0]==0.)|(Bmean<blcut_min)|(Bmean>blcut_max)|condtime:#|(np.intersect1d_nu(np.int16([ant0,ant1]),np.int16(antenna_list)).shape[0]==0):
        #print("exit")
        return np.array([[0.]]),np.array([[0]]),[0],[0],[0],[0]


    #s0=0
    #s1=100
    uvw=uvw[s0:s1]
    time=time[s0:s1]
    vv=vv[s0:s1]

    #select baseline, compute fft
    vsel=np.array(vv,dtype=np.complex)#[ind]
    vsel[np.isnan(vsel)]=0.
    

    import scipy.interpolate
    fft=scipy.fftpack.fftshift(scipy.fftpack.fft(vsel[:,-1]))#/np.abs(vsel)))
    ff0=wavelength_chan[-1]*scipy.fftpack.fftshift(scipy.fftpack.fftfreq(vsel.shape[0],d=time[1]-time[0]))
    #for ch in range(0,Nchan-1):


    ich0=np.max([0,Nchan-1-100])
    ich1=Nchan-1
    for ch in range(ich0,ich1):
        fft_ch=scipy.fftpack.fftshift(scipy.fftpack.fft(vsel[:,ch]))
        ff_ch=wavelength_chan[ch]*scipy.fftpack.fftshift(scipy.fftpack.fftfreq(vsel.shape[0],d=time[1]-time[0]))
        inter=scipy.interpolate.interp1d(ff_ch,fft_ch)
        fft+=inter(ff0)

    absfft=np.abs(fft)
    #std=find_rms(absfft)#np.std(absfft)
    import scipy.stats
    std=scipy.stats.median_absolute_deviation(fft)
    mean=np.median(absfft)
    ff=ff0/wavelength_chan[-1]
    #stop

    # pylab.figure(0)
    # pylab.clf()
    # fft=np.abs(fft)
    # pylab.plot(ff,absfft,marker=".")
    # #stop

    u=uvw[:,0]
    v=uvw[:,1]
    w=uvw[:,2]
    du=u[0:-1]-u[1:]
    dv=v[0:-1]-v[1:]
    dw=w[0:-1]-w[1:]
    dt=time[0:-1]-time[1:]
    du_dt=np.median(du/dt)
    dv_dt=np.median(dv/dt)
    dw_dt=np.median(dw/dt)
    
    theta_rad=np.array(range(ntheta_solve))/(ntheta_solve-1.)*2.*np.pi
    #theta_rad=np.random.rand(ntheta_solve)*4.*np.pi
    ind=theta_rad.argsort()
    theta_rad=theta_rad[ind]


    maxabsfft=np.max(absfft)
    if (maxabsfft<snrcut*std)|(maxabsfft==0.)|(std<1.e-6):
        return np.array([[0.]]),np.array([[0]]),[0],[0],[0],[0]

    iall=0
    acoef=-0.5*dw_dt
    bcoef=(du_dt*np.cos(theta_rad)+dv_dt*np.sin(theta_rad))
    xout=[]
    yout=[]
    NfoundPeaks=0
    freq_list=[]
    Rout=[]
    snr=[]
    while True:
        theta_rad+=+np.random.rand(1)[0]*2.*np.pi
        NfoundPeaks+=1
        #RFI FT level high
        if NfoundPeaks>MaxNfoundPeaks: break
        maxabsfft=np.max(absfft)
        if maxabsfft<mean+snrcut*std: break
        
        # angles, angle increment ...
        ind=np.where(absfft==maxabsfft)[0]
        freq=ff[ind][0]
        #pylab.plot(ff[ind-mask_freq:ind+mask_freq],absfft[ind-mask_freq:ind+mask_freq],marker=".",color='black')
        #pylab.draw()

        absfft[ind[0]-mask_freq:ind[0]+mask_freq]=0.
        #pylab.show()
        #time_mod.sleep(0.1)

        if np.abs(freq)<freq_cut: continue

        ccoef=-freq

        Rs_sol=np.zeros((2,theta_rad.shape[0]),dtype=np.float)-0.555
        if acoef!=0.:
            D=bcoef**2-4.*acoef*ccoef
            ind0=np.where(D>0.)[0]
            if ind0.shape[0]>0:
                Rs_sol[0,ind0]=-(bcoef[ind0]-np.sqrt(D[ind0]))/(2.*acoef)
                Rs_sol[1,ind0]=-(bcoef[ind0]+np.sqrt(D[ind0]))/(2.*acoef)
        else:
            Rs_sol[0,:]=-ccoef/bcoef

        

        #print("====")
        #print(Rs_sol)
        #pylab.figure(0)
        #pylab.clf()
        #cc=["red","blue"]
        for sol in range(2):
            #ind=np.arange(Rs_sol[sol,:].shape[0])
            ind=np.where((np.abs(Rs_sol[sol,:])>-0.3)&(np.abs(Rs_sol[sol,:]+0.555)>0.00001))[0]
            #ind0=np.where((np.abs(Rs_sol[sol,:])<-0.3))[0]
        #!!!!! un peu bancal...
            replace=0
            #x=Rs_sol[sol,ind0]*np.cos(theta_rad[ind0])
            #y=Rs_sol[sol,ind0]*np.sin(theta_rad[ind0])
            #pylab.plot(x,y,color="blue",marker="*",ls='')
            if len(ind)>0:
                Rs_sol[sol,ind]=-0.555
                Rs=Solve_Radius_brent(du_dt,dv_dt,dw_dt,freq,theta_rad[ind])
                #Rs=Solve_Radius_brent(du_dt,dv_dt,dw_dt,0.,theta_rad[ind])
                ind0=np.where(Rs[0,:]!=-0.555)[0]
                if len(ind0)>0:
                    Rs_sol[0,ind]=Rs[0,:]
                    #Rs_sol[0,ind0]=Rs[0,ind0]
                    replace=1
                ind0=np.where(Rs[1,:]!=-0.555)[0]
                if len(ind0)>0:
                    Rs_sol[1,ind]=Rs[1,:]
                    #Rs_sol[0,ind0]=Rs[1,ind0]
                    replace=2
                #x=Rs_sol[sol,ind]*np.cos(theta_rad[ind])
                #y=Rs_sol[sol,ind]*np.sin(theta_rad[ind])
                #pylab.plot(x,y,linestyle='',color=cc[sol],marker='+')
        #pylab.plot([0],[0],linestyle='',color="green",marker='*',)
        #pylab.xlim([-1.,1.])
        #pylab.ylim([-1.,1.])
        #pylab.draw()
        #pylab.show()
        #time_mod.sleep(0.5)
            #print(replace)
            #print(Rs_sol)
            #time_mod.slRs_soleep(0.1)
            
        ind=np.where(Rs_sol>1.)[0]
        if len(ind)>0: stop
        #ind=np.isnan(Rs_sol)
        #Rs_sol[ind]=-0.555

        #ind=np.where(np.abs(Rs_sol[0]+0.555)>0.001)[0]
        #print("max0= ",np.max(np.abs(Rs_sol[0,ind])))
        #if np.max(np.abs(Rs_sol[0,ind]))<0.3:
        #    print(Rs_sol)
        #    stop
        #ind=np.where(np.abs(Rs_sol[1]+0.555)>0.001)[0]
        #print("max1= ",np.max(np.abs(Rs_sol[1,ind])))
        #if np.max(np.abs(Rs_sol[1,ind]))<0.3:
        #    print(Rs_sol)
        #    stop
        x=Rs_sol*np.cos(theta_rad)
        y=Rs_sol*np.sin(theta_rad)
        #pylab.figure(0)
        #pylab.plot(x,y,marker='.',ls='')
        #pylab.draw()
        xout.append(x)
        yout.append(y)
        Rout.append(Rs_sol)
        freq_list.append(freq)
        snr.append(maxabsfft/std)
#    print(xout,yout)
                                 
    if xout==0: return np.array([[0.]]),np.array([[0]]),[0],[0],[0],[0]

    return np.array(xout),np.array(yout),freq_list,Rout,theta_rad,snr

def order_par_spline(x,y):
    import scipy.interpolate
    x0=np.mean(x)
    y0=np.mean(y)
    z=x-x0+1j*(y-y0)
    ang=np.angle(z)
    ind=ang.argsort()
    return ind
    # xn=x[ind]
    # yn=y[ind]
    # theta_rad=theta_rad[ind]

def lm2radec(l_list,m_list):

    ra_list=np.zeros(l_list.shape,dtype=np.float)
    dec_list=np.zeros(l_list.shape,dtype=np.float)
    
    for i in range(l_list.shape[0]):
        l=l_list[i]
        m=m_list[i]
        if (l_list[i]==0.)&(m_list[i]==0.): continue
        Rp=R+vl2*l+vm2*m-(1.-np.sqrt(1.-l**2-m**2))*vn2
        dec_list[i]=np.arcsin(Rp[2])
        #ra_list[i]=np.arccos(Rp[0]/np.cos(dec_list[i]))
        ra_list[i]=np.arctan(Rp[1]/Rp[0])
        if Rp[0]<0.: ra_list[i]+=np.pi
    
    return ra_list,dec_list

def plot_polar(ra,dec,col=0):

    colors=['red', 'blue', 'green', 'grey']
    colors=['grey']#, 'blue', 'green', 'grey']
    alpha=np.arange(0.2,1.1,0.25)[::-1]
    alpha[:]=1.#0.3
    if col>len(colors)-1: col=len(colors)-1
    col=0
    #racut,deccut=cut_array(ra,dec)

    ind=np.where(dec>0.)[0]
    if ind.shape[0]>0:
        pylab.subplot(1,2,1)
        decsel=dec[ind]
        theta=ra[ind]
        radius=(np.pi/2.-decsel)*180./np.pi
        x=radius*np.cos(theta)
        y=radius*np.sin(theta)
        pylab.plot(-x,y,lw=0.3,ls=':',color=colors[col],alpha=alpha[col])
    ind=np.where(dec<0.)[0]
    if ind.shape[0]>0:
        pylab.subplot(1,2,2)
        decsel=dec[ind]
        theta=ra[ind]
        radius=(np.pi/2.+decsel)*180./np.pi
        x=radius*np.cos(theta)
        y=radius*np.sin(theta)
        pylab.plot(-x,y,lw=0.3,ls=':',color=colors[col],alpha=alpha[col])#,color=colors[col])





def load_Ateam():
    global Ateam
    Ateam = [ {'name' : 'CasA', 'ra' : 6.123487680622104,  'dec' : 1.0265153995604648, 'tag': 0},
              {'name' : 'CygA', 'ra' : 5.233686575770755,  'dec' : 0.7109409582180791, 'tag': 0},
              {'name' : 'TauA', 'ra' : 1.4596748493730913, 'dec' : 0.38422502335921294, 'tag': 0},
              {'name' : 'HerA', 'ra' : 4.4119087330382163, 'dec' : 0.087135562905816893, 'tag': 0},
              {'name' : 'VirA', 'ra' : 3.276086511413598,  'dec' : 0.21626589533567378, 'tag': 0}]
    
    import pyrap.quanta as qa
    import pyrap.measures as pm
    time_start = qa.quantity(time_all[0], 's')
    me = pm.measures()
    dict_time_start_MDJ = me.epoch('utc', time_start)
    time_start_MDJ=dict_time_start_MDJ['m0']['value']
    JD=time_start_MDJ+2400000.5-2415020
    d=ephem.Date(JD)
    sun, jup = ephem.Sun(), ephem.Jupiter()
    sun.compute(d)
    jup.compute(d)
    Ateam.append({'name': 'Sun', 'ra': sun.ra, 'dec': sun.dec, 'tag': 1})
    Ateam.append({'name': 'Jupiter', 'ra': jup.ra, 'dec': jup.dec, 'tag': 1})
    


def Solve_Radius_brent(du_dt,dv_dt,dw_dt,freq,theta):
    
    import scipy.interpolate
    import scipy.optimize
    RadS=np.linspace(0.,0.999,10)#np.arange(0.,1.1,0.1)
    Rout=np.zeros((2,theta.shape[0]),dtype=np.float)-0.555

    #pylab.figure(0)
    #pylab.clf()
    for i in range(theta.shape[0]):
        yy=Freq_lmn(RadS,du_dt,dv_dt,dw_dt,freq,theta[i])
        #if (np.min(yy)>0.)|(np.max(yy)<0.): continue
        yyp=Freq_lmn_prim(RadS,du_dt,dv_dt,dw_dt,freq,theta[i])

        #pylab.plot(RadS,yy)
        #pylab.plot([0,1],[0,0],ls='--')

        cond=[yyp<0.,yyp>0.]
        for ii in range(2):
            ind,end=give_mask(cond[ii],0)
            #print(ii,ind,end)
            if end==False: continue
            RadSsel=RadS[ind]
            yysel=yy[ind]
            #pylab.plot(RadSsel,yysel)
            if (np.min(yysel)<0.)&(np.max(yysel)>0.):
                Rout[ii,i]=scipy.optimize.brentq(Freq_lmn, np.min(RadSsel), np.max(RadSsel), args=(du_dt,dv_dt,dw_dt,freq,theta[i]), xtol=1e-3, rtol=1e-3, maxiter=100, full_output=False, disp=True)
                #pylab.plot(Rout[ii,i],[0],marker='*')
        #pylab.show()
        #time_mod.sleep(0.1)
        #

    #pylab.draw()
    ind=np.where(Rout==-0.555)
    if freq==0.:
        Rout[ind]=0.
    # if np.max(np.abs(Rout))<0.6:
    #     pylab.savefig("lala.png")
    #     print(Rout)
    #     print()
    #     stop
    return Rout
        

def Freq_lmn(rad,du_dt,dv_dt,dw_dt,freq,theta): 
    return rad*(du_dt*np.cos(theta)+dv_dt*np.sin(theta))+dw_dt*(np.sqrt(1.-rad**2)-1.)-freq

def Freq_lmn_prim(rad,du_dt,dv_dt,dw_dt,freq,theta):
    return du_dt*np.cos(theta)+dv_dt*np.sin(theta)-dw_dt*rad/np.sqrt(1.-rad**2)






      
def cut_array0(A):
    Aout=[]
    AoutTmp=[]
    for i in range(A.shape[0]):
        if A[i]!=0.:
            AoutTmp.append(A[i])
            
        if (A[i]==0.)|(i==A.shape[0]-1):
            if len(AoutTmp)>0:
                Aout.append(np.array(AoutTmp))
                AoutTmp=[]

 
    return Aout
                
def cut_array(ra,dec):
    ind=np.where(dec>0.)[0]
    RAout=[]
    DECout=[]
    if ind.shape[0]>0.:
        RAoutTmp=[ra[ind[0]]]
        DECoutTmp=[dec[ind[0]]]
        indold=ind[0]
        for i in range(1,ind.shape[0]):
            if ind[i]==indold+1:
                RAoutTmp.append(ra[ind[i]])
                DECoutTmp.append(dec[ind[i]])
            if (ind[i]!=indold+1)|(i==ind.shape[0]-1):
                RAout.append(np.array(RAoutTmp))
                DECout.append(np.array(DECoutTmp))
                RAoutTmp=[ra[ind[i]]]
                DECoutTmp=[dec[ind[i]]]
            indold=ind[i]

    ind=np.where(dec<0.)[0]
    if ind.shape[0]>0.:

        RAoutTmp=[ra[ind[0]]]
        DECoutTmp=[dec[ind[0]]]
        indold=ind[0]
        for i in range(1,ind.shape[0]):
            if ind[i]==indold+1:
                RAoutTmp.append(ra[ind[i]])
                DECoutTmp.append(dec[ind[i]])
            if (ind[i]!=indold+1)|(i==ind.shape[0]-1):
                RAout.append(np.array(RAoutTmp))
                DECout.append(np.array(DECoutTmp))
                RAoutTmp=[ra[ind[i]]]
                DECoutTmp=[dec[ind[i]]]
            indold=ind[i]

    return RAout,DECout
                

      
def give_mask(Mask_bool,n):
    
    #Mask_bool=A>0.
    indn=-1
    ind=0
    Aold=False
    mask_out=np.zeros(Mask_bool.shape,dtype=np.bool)
    while (Mask_bool[ind]!=True)|(indn!=n):
        if (Mask_bool[ind]==True)&(Aold==False):
            indn+=1
        if ind==Mask_bool.shape[0]-1:
            return mask_out,False

        Aold=Mask_bool[ind]
        ind+=1


    ind-=1
    while Mask_bool[ind]==True:
        mask_out[ind]=True
        if ind==Mask_bool.shape[0]-1: return mask_out,True
        ind+=1
        
    return mask_out,True
    

    
#========================================

def init_plot_lm():

    pylab.figure(3)
    ax = pylab.gca()
    ax.format_coord = lambda x,y : "RA=%2.2ih %2.2im %2.2is DEC=%s%2.2id %2.2im %2.2is" % lm2radec_scalar(x,y)
    load_Ateam()
    colAteam=['red','green']
    markerAteam=['+','*']
    SmarkerAteam=[1.,1.5]
    for i in range(len(Ateam)):
        pi2=np.pi/2.
        ra=Ateam[i]['ra']
        dec=Ateam[i]['dec']
        cosd=np.cos(pi2-dec)*np.cos(pi2-decrad)+np.sin(pi2-dec)*np.sin(pi2-decrad)*np.cos(ra-rarad)
        d=np.arccos(cosd)
        if d<np.pi/2.:
            l,m=radec2lm_scalar(ra,dec)
            name=Ateam[i]['name']
            pylab.plot([l],[m],marker=markerAteam[Ateam[i]['tag']],color=colAteam[Ateam[i]['tag']],
                       mew=2*SmarkerAteam[Ateam[i]['tag']],ms=5*SmarkerAteam[Ateam[i]['tag']], mec=colAteam[Ateam[i]['tag']])
            pylab.text(l,m,name,color=colAteam[Ateam[i]['tag']])
        
    theta=np.arange(0.,2.*np.pi,0.01)
    Rbeam=wavelength/70.
    l,m=Rbeam*np.cos(theta),Rbeam*np.sin(theta)
    pylab.plot(Rbeam*np.cos(theta),Rbeam*np.sin(theta),color='blue')
    pylab.plot(np.cos(theta),np.sin(theta),color='black',ls="-")
    pylab.title("Image plane projection")
    pylab.xlim([1.,-1.])
    pylab.ylim([-1.,1.])
    if RadPlotNVSS!=0.:
        for i in range(decnvss.shape[0]):
            l,m=radec2lm_scalar(ranvss[i],decnvss[i])
            pylab.plot([l],[m],marker=".",color='blue', mew=1,ms=1, mec="blue")
            ss=("%6.3f") % snvss[i]
            pylab.text(l,m,ss,fontsize=7)
    

def lm2radec_scalar(l,m):
    import deg2hmsdms
    Rp=R+vl2*l+vm2*m-(1.-np.sqrt(1.-l**2-m**2))*vn2
    dec=np.arcsin(Rp[2])
    ra=np.angle(Rp[0]+1j*Rp[1])
    #if Rp[0]<0.: ra+=np.pi
    if ra<0.: ra+=np.pi*2.
    strra=deg2hmsdms.give_str_time((24./360)*(ra)*180./(np.pi)).split(":")
    strdec=deg2hmsdms.give_str_time(np.abs(dec)*180./(np.pi)).split(":")
    sg="-"
    if np.sign(dec)>0.: sg="+"
    return (int(strra[0]), int(strra[1]), float(strra[2]), sg, int(strdec[0]), int(strdec[1]), float(strdec[2]))


def radec2lm_scalar(ra,dec):
    l = np.cos(dec) * np.sin(ra - rarad)
    m = np.sin(dec) * np.cos(decrad) - np.cos(dec) * np.sin(decrad) * np.cos(ra - rarad)
    return l,m
    

def init_plot_polar():

    theta=np.arange(0.,2.*np.pi,0.01)
    x0=(np.pi/2.-decrad)*180./np.pi*np.cos(rarad)
    y0=(np.pi/2.-decrad)*180./np.pi*np.sin(rarad)

    # pylab.subplot(1,2,1)
    # rad=15.*(8.+03/60.+56./3600)*np.pi/180.
    # decd=49+14/60.+07./3600
    # r=90.-decd
    # print(rad,decd)
    # pylab.plot([r*np.cos(rad)],[r*np.sin(rad)],marker='.',color='blue')
    

    # rad=15.*(8.+13/60.+34./3600)
    # decd=48+13/60.+57./3600
    # print(rad,decd)
    # pylab.plot([rad],[decd],marker='*',color='blue')

    r=np.arange(0,91.,10)
    h=np.arange(0,24.,1)
    dl=-2
    lwgrid=1.
    for i in range(h.shape[0]):
        ang=h[i]*2.*np.pi/24.
        pylab.subplot(1,2,1)
        ax = pylab.gca()
        #ax.format_coord = lambda x,y : "x=%g y=%g" % (np.sqrt(x**2+y**2), np.arctan(y/np.sqrt(x**2+y**2)))
        #ax.format_coord = lambda x,y : "RA=%5.2f DEC=%5.2f" % ((np.angle(x-1j*y)+np.pi)*12./np.pi,90.-np.abs(x+1j*y))
        ax.format_coord = lambda x,y : "RA=%2.2ih %2.2im %2.2is DEC=%2.2id %2.2im %2.2is" % (np.floor((np.angle(x-1j*y)+np.pi)*12./np.pi),
                                                                                             np.floor(60.*((np.angle(x-1j*y)+np.pi)*12./np.pi-np.floor((np.angle(x-1j*y)+np.pi)*12./np.pi))),
                                                                                             60.*(60.*((np.angle(x-1j*y)+np.pi)*12./np.pi-np.floor((np.angle(x-1j*y)+np.pi)*12./np.pi))-
                                                                                                  np.floor(60.*((np.angle(x-1j*y)+np.pi)*12./np.pi
                                                                                                                -np.floor((np.angle(x-1j*y)+np.pi)*12./np.pi)))),
                                                                                             np.floor(90.-np.abs(x+1j*y)),
                                                                                             np.floor(60.*(90.-np.abs(x+1j*y)-np.floor(90.-np.abs(x+1j*y)))),
                                                                                             60.*(60.*(90.-np.abs(x+1j*y)-np.floor(90.-np.abs(x+1j*y)))-
                                                                                                  np.floor(60.*(90.-np.abs(x+1j*y)
                                                                                                                -np.floor(90.-np.abs(x+1j*y))))))
        pylab.title("Northern Celestial Sphere")
        pylab.plot([0,-92.*np.cos(ang)],[0,92.*np.sin(ang)],ls='--',lw=lwgrid,color='black')
        if i>0:
            pylab.text(-92.*np.cos(ang)+dl,92.*np.sin(ang)+dl,str(int(h[i])),fontsize=8)
            pylab.text(-42.*np.cos(ang)+dl,42.*np.sin(ang)+dl,str(int(h[i])),fontsize=8)
        pylab.subplot(1,2,2)
        ax = pylab.gca()
        #ax.format_coord = lambda x,y : "x=%g y=%g" % (np.sqrt(x**2+y**2), np.arctan(y/np.sqrt(x**2+y**2)))
        #ax.format_coord = lambda x,y : "RA=%5.2f DEC=%5.2f" % ((np.angle(x-1j*y)+np.pi)*12./np.pi,-90.+np.abs(x+1j*y))
        ax.format_coord = lambda x,y : "RA=%2.2ih %2.2im %2.2is DEC=-%2.2id %2.2im %2.2is" % (np.floor((np.angle(x-1j*y)+np.pi)*12./np.pi),
                                                                                             np.floor(60.*((np.angle(x-1j*y)+np.pi)*12./np.pi-np.floor((np.angle(x-1j*y)+np.pi)*12./np.pi))),
                                                                                             60.*(60.*((np.angle(x-1j*y)+np.pi)*12./np.pi-np.floor((np.angle(x-1j*y)+np.pi)*12./np.pi))-
                                                                                                  np.floor(60.*((np.angle(x-1j*y)+np.pi)*12./np.pi
                                                                                                                -np.floor((np.angle(x-1j*y)+np.pi)*12./np.pi)))),
                                                                                             np.floor(90.-np.abs(x+1j*y)),
                                                                                             np.floor(60.*(90.-np.abs(x+1j*y)-np.floor(90.-np.abs(x+1j*y)))),
                                                                                             60.*(60.*(90.-np.abs(x+1j*y)-np.floor(90.-np.abs(x+1j*y)))-
                                                                                                  np.floor(60.*(90.-np.abs(x+1j*y)
                                                                                                                -np.floor(90.-np.abs(x+1j*y))))))
        pylab.title("Sourtern Celestial Sphere")
        pylab.plot([0,-92.*np.cos(ang)],[0,92.*np.sin(ang)],ls='--',lw=lwgrid,color='black')
        if i>0:
            pylab.text(-92.*np.cos(ang)+dl,92.*np.sin(ang)+dl,str(int(h[i])),fontsize=8)
            pylab.text(-42.*np.cos(ang)+dl,42.*np.sin(ang)+dl,str(int(h[i])),fontsize=8)

    for i in range(r.shape[0]):
        pylab.subplot(1,2,1)
        x=r[i]*np.cos(theta)
        y=r[i]*np.sin(theta)
        pylab.plot(-x,y,lw=lwgrid,ls='--',color='black')
        if decrad>0.: pylab.plot([-x0],[y0],marker='*',ms=10,color='red')
        if i>0: pylab.text(-(r[i]+1),0.,str(int(90.-r[i])),fontsize=8)
        pylab.plot([0],[0],marker='+',ms=10,color='black')
        pylab.subplot(1,2,2)
        pylab.plot(-x,y,lw=lwgrid,ls='--',color='black')
        if decrad<0.: pylab.plot([-x0],[y0],marker='*',ms=10,color='red')
        if i>0: pylab.text(-(r[i]+1),0,str(int(r[i]-90.)),fontsize=8)


    if RadPlotNVSS!=0.:
        for i in range(decnvss.shape[0]):
            pylab.subplot(1,2,2)
            if decnvss[i]>0.: pylab.subplot(1,2,1)
            x=(np.pi/2.-decnvss[i])*180./np.pi*np.cos(ranvss[i])
            y=(np.pi/2.-decnvss[i])*180./np.pi*np.sin(ranvss[i])
            pylab.plot([-x],[y],marker=".",color='blue', mew=1,ms=1, mec="blue")
            ss=("%6.1f") % snvss[i]
            pylab.text(-x,y,ss,fontsize=7)

    load_Ateam()
    colAteam=['red','green']
    markerAteam=['+','*']
    SmarkerAteam=[1.,1.5]
    for i in range(len(Ateam)):
        pylab.subplot(1,2,1)
        ra=Ateam[i]['ra']
        dec=Ateam[i]['dec']
        if dec<0.:
            dec=-dec
            pylab.subplot(1,2,2)
        x=(np.pi/2.-dec)*180./np.pi*np.cos(ra)
        y=(np.pi/2.-dec)*180./np.pi*np.sin(ra)
        name=Ateam[i]['name']
        pylab.plot([-x],[y],marker=markerAteam[Ateam[i]['tag']],color=colAteam[Ateam[i]['tag']],
                   mew=2*SmarkerAteam[Ateam[i]['tag']],ms=5*SmarkerAteam[Ateam[i]['tag']], mec=colAteam[Ateam[i]['tag']])
        pylab.text(-x,y,name,color=colAteam[Ateam[i]['tag']])
        
    pylab.subplot(1,2,1)
    Rbeam=wavelength/70.
    l,m=Rbeam*np.cos(theta),Rbeam*np.sin(theta)
    ra,dec=lm2radec(l,m)
    r=90.-dec*180./np.pi
    pylab.plot(-r*np.cos(ra),r*np.sin(ra),color='blue')

    


def init_catalog_lines():
    
    dtype_line=[('x',np.float32,(ntheta_solve,)), ('y',np.float32,(ntheta_solve,)), ('r',np.float32,(ntheta_solve,)), ('theta',np.float32,(ntheta_solve,))]
    global catalog_lines
    catalog_lines=np.zeros((int(3.*nbl*(ntimes/delta_time_bins+1.)*MaxNfoundPeaks),),dtype=[('bln',np.int32),('timeslot',np.float),('line',dtype_line),('freq',np.float32),('flag',np.bool)])
    catalog_lines=catalog_lines.view(np.recarray)
    catalog_lines.timeslot=-1.
    catalog_lines.line.x=-1.
    catalog_lines.line.y=-1.
    catalog_lines.flag=False
    

#==============================

def plot_gridded():

    cat=np.load("catalog_lines.npy")
    cat=cat.view(np.recarray)

    import scipy.ndimage as ndi
    
    pylab.figure(5,figsize=(10,10))
    pylab.clf()
    ax = pylab.gca()
    ax.format_coord = lambda x,y : "RA=%2.2ih %2.2im %2.2is DEC=%s%2.2id %2.2im %2.2is" % lm2radec_scalar(-x,y)

    npix=4000
    img = np.ones((npix,npix))
    N_interp_line=100
    dens_interp_line=1e-4
    consize=2

    print("  Grid the lines ... ")
    toolbar_width = 50
    sys.stdout.write("[%s]" % (" " * (toolbar_width)))
    sys.stdout.flush()
    sys.stdout.write("\b" * (toolbar_width+1)) # return to start of line, after '['
    incr=int(cat.shape[0]/float(toolbar_width))
    tot=0

    for c in range(cat.shape[0]):
        if tot==int(incr):
            tot=0
            sys.stdout.write(progress_str)
            sys.stdout.flush()
        tot+=1

        ind=np.where((cat.line.x[c]!=-1.)&(cat.line.x[c]!=np.nan)&(cat.line.y[c]!=np.nan)&(cat.line.x[c]!=1.)&(cat.line.y[c]!=-1.)&(cat.line.y[c]!=1.))
        if ind[0].shape[0]>0:
            for i in range(ind[0].shape[0]-1):
                xm=cat.line.x[c][ind[0]][i:i+2]
                ym=cat.line.y[c][ind[0]][i:i+2]
                dist=np.sqrt((xm[0]-xm[1])**2+(ym[0]-ym[1])**2)
                indu=xm.argsort()
                xm=xm[indu]
                ym=ym[indu]
                if dist==0: continue
                N_interp_line=dist/dens_interp_line
                xs=np.arange(np.min(xm),np.max(xm),(np.max(xm)-np.min(xm))/(N_interp_line))
                ys=np.interp(xs, xm, ym)
                x=npix*(xs+1.)/2.
                y=npix*(ys+1.)/2.
                indnan=np.where((x!=np.nan)&(y!=np.nan))[0]
                img[np.int16(x[indnan]), np.int16(y[indnan])] += 1

    sys.stdout.flush()
    print("\n")

    img = np.log10(img)
    #img = ndi.gaussian_filter(img, (consize,consize))

    pylab.imshow(-img.T[::-1,::-1],extent=(-1.,1.,-1.,1.),cmap=pylab.cm.Greys)#,interpolation="nearest")
    pylab.xlabel("L")
    pylab.ylabel("M")
    pylab.title("Image plane")
    theta=np.arange(0.,2.*np.pi,0.01)
    pylab.plot(np.cos(theta),np.sin(theta),color='black')
    


    load_Ateam()
    colAteam=['red','green']
    markerAteam=['+','*']
    SmarkerAteam=[1.,1.5]
    for i in range(len(Ateam)):
        pi2=np.pi/2.
        ra=Ateam[i]['ra']
        dec=Ateam[i]['dec']
        cosd=np.cos(pi2-dec)*np.cos(pi2-decrad)+np.sin(pi2-dec)*np.sin(pi2-decrad)*np.cos(ra-rarad)
        d=np.arccos(cosd)
        if d<np.pi/2.:
            l,m=radec2lm_scalar(ra,dec)
            name=Ateam[i]['name']
            pylab.plot([-l],[m],marker=markerAteam[Ateam[i]['tag']],color=colAteam[Ateam[i]['tag']],
                       mew=2*SmarkerAteam[Ateam[i]['tag']],ms=5*SmarkerAteam[Ateam[i]['tag']], mec=colAteam[Ateam[i]['tag']])
            pylab.text(-l,m,name,color=colAteam[Ateam[i]['tag']])

    if RadPlotNVSS!=0.:
        for i in range(decnvss.shape[0]):
            l,m=radec2lm_scalar(ranvss[i],decnvss[i])
            pylab.plot([-l],[m],marker=".",color='blue', mew=1,ms=1, mec="blue")
            ss=("%6.3f") % snvss[i]
            pylab.text(-l,m,ss,fontsize=7,color="blue")

    #pylab.colorbar()
    theta=np.arange(0.,2.*np.pi,0.01)
    Rbeam=wavelength/70.
    l,m=Rbeam*np.cos(theta),Rbeam*np.sin(theta)
    pylab.plot(Rbeam*np.cos(theta),Rbeam*np.sin(theta),color='red')#,ls="--")#,lw=.5)

    pylab.draw()

    # import pyfits
    # hdu = pyfits.PrimaryHDU(img)
    # hdulist = pyfits.HDUList([hdu])
    # import os
    # os.system('rm -Rf new.fits')
    # hdulist.writeto('new.fits')

    # from pyrap.images import image
    # im=image('new.fits')
    # os.system('rm -Rf new.imms')
    # im.saveas('new.imms')

def movie_maker(namein):
    from matplotlib.backends.backend_pdf import PdfPages
    import pyrap.quanta as qa
    import scipy.ndimage as ndi
    import pyrap.measures as pm
    import give_elevation
    import gc
    gc.enable()

    cata=np.load("catalog_lines.npy")
    cata=cata.view(np.recarray)
    pdfname=namein+".pdf"
    #pdf = PdfPages(pdfname)
    dirname="dMSprods."+namein+"/"
    dirnamepng=dirname+"png/"
    os.system("rm -Rf "+dirname)
    os.system("mkdir "+dirname)
    os.system("mkdir "+dirnamepng)

    npix=1500
    N_interp_line=100
    dens_interp_line=1e-3
    consize=1

    ss=range(0,ntimes,delta_time_bins)

    load_Ateam()
    az_list=[]
    el_list=[]
    name_list=[]
    for i in range(len(Ateam)):
        if Ateam[i]['name']!='Sun': continue
        ra=Ateam[i]['ra']
        dec=Ateam[i]['dec']
        az,el=give_elevation.give_elevation(time_slots_all[ss],ra,dec,xarray,yarray,zarray)
        azf,elf=give_elevation.give_elevation(time_slots_all[ss],rarad,decrad,xarray,yarray,zarray)
        az_list.append(az)
        el_list.append(el)
        name_list.append(Ateam[i]['name'])


    toolbar_width = 50
    print("  Grid the lines for the movie ... ")
    sys.stdout.write("[%s]" % (" " * (toolbar_width)))
    sys.stdout.flush()
    sys.stdout.write("\b" * (toolbar_width+1)) # return to start of line, after '['
    incr=int(cata.shape[0]/float(toolbar_width))
    tot=0

    img_stack = np.zeros((npix,npix))

    for ii in range(len(ss)-1):
        time_slots_all
        img = np.ones((npix,npix))
        indt=np.where(cata.timeslot==ii)[0]
        cat=cata[indt]
        
        time_center = qa.quantity(time_slots_all[ss[ii]], 's')
        me = pm.measures()
        dict_time_center_MDJ = me.epoch('utc', time_center)
        time_center_MDJ=dict_time_center_MDJ['m0']['value']
        JD=time_center_MDJ+2400000.5-2415020
        d=ephem.Date(JD)
        eldeg=el[ii]*180./np.pi
        streldeg="%5.1f" % eldeg
        # global strdate
        # strdate=str(d)+"  Sun elevation = "+streldeg+" deg"

        for c in range(cat.shape[0]):
            if tot==int(incr):
                tot=0
                sys.stdout.write(progress_str)
                sys.stdout.flush()
            tot+=1
    
            ind=np.where((cat.line.x[c]!=-1.)&(cat.line.x[c]!=np.nan)&(cat.line.y[c]!=np.nan)&(cat.line.x[c]!=1.)&(cat.line.y[c]!=-1.)&(cat.line.y[c]!=1.))
            if ind[0].shape[0]>0:
                # x=cat.line.x[c][ind[0]]
                # y=cat.line.y[c][ind[0]]
                # ind0=order_par_spline(x,y)
                # cat.line.x[c][ind[0]]=x[ind0]
                # cat.line.y[c][ind[0]]=y[ind0]
                for i in range(ind[0].shape[0]-1):
                    xm=cat.line.x[c][ind[0]][i:i+2]
                    ym=cat.line.y[c][ind[0]][i:i+2]
                    dist=np.sqrt((xm[0]-xm[1])**2+(ym[0]-ym[1])**2)
                    dist=np.sqrt((xm[0]-xm[1])**2+(ym[0]-ym[1])**2)
                    N_interp_line=dist/dens_interp_line
                    N_interp_line
                    indu=xm.argsort()
                    xm=xm[indu]
                    ym=ym[indu]
                    xs=np.arange(np.min(xm),np.max(xm),(np.max(xm)-np.min(xm))/N_interp_line)
                    ys=np.interp(xs, xm, ym)
                    x=npix*(xs+1.)/2.
                    y=npix*(ys+1.)/2.
                    indnan=np.where((x!=np.nan)&(y!=np.nan))[0]
                    img[np.int16(x[indnan]), np.int16(y[indnan])] += 1
    
        sys.stdout.flush()
    
        img = np.log10(img)
        img_stack+=img
        
        fig=pylab.figure(5,figsize=(16,8))
        plot2panels(fig,img,noNVSS=True)

        #pdf.savefig()
        fname = namein+'%03d.png'%ii
        fig.savefig(dirnamepng+fname)
        #files.append(fname)

        pylab.close(fig)
        del(cat)
        gc.collect()
    print("\n")
    gc.collect()
    #print("  Saving movie ... "+pdfname)
    #pdf.close()
    #time_mod.sleep(2.)
    #print("  Converting png to animated gif ..."+gifname)
    #os.system("convert _tmp_png/"+namein+"*.png "+gifname)
    print("  Images put in "+dirnamepng)
    print("  Converting png to mpg ... "+dirname+namein+".mpg")
    opt="vbitrate=6160000:mbd=2:keyint=132:v4mv:vqmin=3:lumi_mask=0.07:dark_mask=0.2:mpeg_quant:scplx_mask=0.1:tcplx_mask=0.1:naq"
    #os.system("mencoder 'mf://_tmp_png/"+namein+"*.png' -mf type=png:fps=10 -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o "+namein+".mpg > lala 2>&1 ") 
    os.system("mencoder -ovc lavc -lavcopts vcodec=mpeg4:vpass=1:"+opt+" -mf type=png:fps=10 -nosound -o "+dirname+namein+".mpg mf://"+dirnamepng+"\*.png > lala 2>&1" )

    print("  Stacked image ... "+dirname+namein+".stack.png")
    fig=pylab.figure(5,figsize=(32,16))
    pylab.clf()
    plot2panels(fig,img_stack)
    fig.savefig(dirname+namein+".stack.png")



def plot2panels(fig,img,noNVSS=False,noElevation=False):

    pylab.clf()
    pylab.subplot(1,2,1)
    #pannel=fig.add_subplot(1,2,1)
    #pannel.set_position([0.1,0.1,0.8,0.8])
    pylab.imshow(-img.T[::-1,::-1],extent=(-1.,1.,-1.,1.),cmap=pylab.cm.Greys)#,interpolation="nearest")
    pylab.xlabel("L")
    pylab.ylabel("M")
    pylab.title("Image plane")
    theta=np.arange(0.,2.*np.pi,0.01)
    pylab.plot(np.cos(theta),np.sin(theta),color='black')
    #pylab.annotate(strdate,(0.3,0.9),xycoords='figure fraction', fontsize=20)

    colAteam=['red','green']
    markerAteam=['+','*']
    SmarkerAteam=[1.,1.5]
    for i in range(len(Ateam)):
        pi2=np.pi/2.
        ra=Ateam[i]['ra']
        dec=Ateam[i]['dec']
        cosd=np.cos(pi2-dec)*np.cos(pi2-decrad)+np.sin(pi2-dec)*np.sin(pi2-decrad)*np.cos(ra-rarad)
        d=np.arccos(cosd)
        if d<np.pi/2.:
            l,m=radec2lm_scalar(ra,dec)
            name=Ateam[i]['name']
            pylab.plot([-l],[m],marker=markerAteam[Ateam[i]['tag']],color=colAteam[Ateam[i]['tag']],
                       mew=2*SmarkerAteam[Ateam[i]['tag']],ms=5*SmarkerAteam[Ateam[i]['tag']], mec=colAteam[Ateam[i]['tag']])
            pylab.text(-l,m,name,color=colAteam[Ateam[i]['tag']])


    theta=np.arange(0.,2.*np.pi,0.01)
    Rbeam=wavelength/70.
    l,m=Rbeam*np.cos(theta),Rbeam*np.sin(theta)
    pylab.plot(Rbeam*np.cos(theta),Rbeam*np.sin(theta),color='red')#,ls="--")#,lw=.5)

    pylab.subplot(1,2,2)
    pylab.imshow(-img.T[::-1,::-1],extent=(-1.,1.,-1.,1.),cmap=pylab.cm.Greys)#,interpolation="nearest")
    pylab.xlabel("L")
    pylab.ylabel("M")
    pylab.title("Image plane")
    del(img)
    theta=np.arange(0.,2.*np.pi,0.01)
    pylab.plot(np.cos(theta),np.sin(theta),color='black')
    theta=np.arange(0.,2.*np.pi,0.01)
    Rbeam=wavelength/70.
    l,m=Rbeam*np.cos(theta),Rbeam*np.sin(theta)
    pylab.plot(Rbeam*np.cos(theta),Rbeam*np.sin(theta),color='red')#,ls="--")#,lw=.5)
    xlim=[-0.3,0.3]
    ylim=[-0.3,0.3]
    if (RadPlotNVSS!=0.)&(noNVSS==False):
        mkmin=1
        mkmax=10
        lvec,mvec=radec2lm_scalar(ranvss,decnvss)
        ind=np.where((lvec<xlim[1])&(lvec>xlim[0])&(mvec<ylim[1])&(mvec>ylim[0]))[0]
        
        snvssc=snvss[ind]#np.log10(snvss)
        lvec=lvec[ind]
        mvec=mvec[ind]
        smin=np.min(snvssc)
        smax=10.*smin#np.median(snvssc)+1.*np.std(snvssc)
        
        for i in range(snvssc.shape[0]):
            l,m=lvec[i],mvec[i]#radec2lm_scalar(ranvss[i],decnvss[i])
            mksize=np.min([mkmin+snvssc[i]*(mkmax-mkmin)/(smax-smin),mkmax])
            pylab.plot([-l],[m],marker=".",color='blue', mew=1,ms=mksize, mec="blue")
            sss=("%6.3f") % snvss[i]
            #pylab.text(-l,m,sss,fontsize=7,color="blue")
    pylab.xlim(xlim)
    pylab.ylim(ylim)


def set_options(options):
    print
    print(" DRAWING OPTIONS: ")
    print("   - Column Name: ",options.ColName)
    global blcut_min, blcut_max, snrcut, wmax, mask_freq, delta_time_bins, MaxNfoundPeaks, ntheta_solve, RadPlotNVSS, time_range_min, time_range_max
    blcut_min=float(options.uvrange.split(",")[0])
    blcut_max=float(options.uvrange.split(",")[1])
    print("   - Baseline selection uv-range  = ",blcut_min,blcut_max)
    snrcut=float(options.snrcut)
    print("   - SNR Cut for fringe searching = ",snrcut)
    wmax=float(options.wmax)
    print("   - WMax = ",wmax)
    time_range_min=float(options.timerange.split(",")[0])
    time_range_max=float(options.timerange.split(",")[1])
    print("   - Timerange (fraction of observation time)  = ",time_range_min,time_range_max)
    mask_freq=int(options.maskfreq)
    print("   - Mask Frequencies  = ",mask_freq)
    delta_time_bins=int(options.timestep)
    print("   - Time step   = ",delta_time_bins," bins and ",delta_time_bins*(time_all[nbl+1]-time_all[0])/60.," minutes")
    time_window_bins=int(options.timewindow)
    if time_window_bins==0:
        time_window_bins=delta_time_bins
    print("   - Time window = ",time_window_bins," bins and ",time_window_bins*(time_all[nbl+1]-time_all[0])/60.," minutes")
    RadPlotNVSS=float(options.RadNVSS)
    print("   - Plot NVSS = ",RadPlotNVSS>0)
    SlimNVSS=float(options.SlimNVSS)
    if RadPlotNVSS>0:
        print("     - Radius Plot NVSS (in beams diameters) = ",RadPlotNVSS)
        print("     - Limiting NVSS flux = ",SlimNVSS)
        
        global ranvss,decnvss,snvss
        ranvss,decnvss,snvss=oplot_nvss.give_nvss(rarad,decrad,RadPlotNVSS*wavelength/70.,SlimNVSS)
    MaxNfoundPeaks=int(options.MaxNPeaks)
    print("   - Maximum number of peaks = ",MaxNfoundPeaks)
    ntheta_solve=int(options.NTheta)
    print("   - Number of angles to solve for = ",ntheta_solve)
    fill_f=float(options.FillFactor)
    if fill_f<1.: 
        global fill_factor
        fill_factor=float(options.FillFactor)
        print("   - Filling factor = ",fill_factor)
    if len(options.AntList)>0:
        global antenna_list
        antenna_list=np.int16(options.AntList.split(","))
        print("   - List of antennas = ",antenna_list)
    print

def display_time(stepname,init=False):
    global tstart
    if init==True:
        tstart=time_mod.time()
        return
    print("  * "+stepname+" computation time: " , time_mod.time()-tstart," seconds")
    print()
    tstart=time_mod.time()


def print_logo():

    os.system('clear')
    print
    print("""               =========== drawMS ===========""")
    print

def find_rms(m):
    
    rmsold=np.std(m)
    diff=1e-1
    cut=3.
    bins=np.arange(np.min(m),np.max(m),(np.max(m)-np.min(m))/30.)
    #print(rmsold)
    med=np.median(m)
    #pylab.hist(m,bins=bins,cumulative=True)
    for i in range(10):
        ind=np.where(np.abs(m-med)<rmsold*cut)[0]
        rms=np.std(m[ind])
        #print("rms  =",rms)
        #print("diff =",(rms-rmsold)/rmsold)
        if np.abs((rms-rmsold)/rmsold)<diff: break
        rmsold=rms

#    bins=np.arange(np.min(m[ind]),np.max(m[ind]),(np.max(m[ind])-np.min(m[ind]))/30.)
    #pylab.hist(m[ind],bins=bins,alpha=0.8,cumulative=True)

    return rms
