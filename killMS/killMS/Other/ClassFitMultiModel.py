from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as np
from DDFacet.Other import logger
log=logger.getLogger("ClassFit_1D")
import killMS.Array.ModLinAlg
K=8.4479745e9
import scipy.sparse
from DDFacet.Other import ClassTimeIt
from scipy.optimize import least_squares
import scipy.ndimage
from . import ClassFitTEC
from . import ClassFitAmp
import scipy.signal
import scipy

logger.setSilent(["ClassFit_1D"])

def TECToPhase(TEC,freq):
    phase=K*TEC*(1./freq)
    return phase

def TECToZ(TEC,ConstPhase,freq):
    if ConstPhase is None: ConstPhase=0
    return np.exp(1j*(TECToPhase(TEC,freq)+ConstPhase))

def Dot(*args):
    P=1.
    for M in args:
        #P=np.dot(np.complex128(P),np.complex128(M))
        P=np.dot(P,M)
    return P

# it=208; iDir=14; S=np.load("L229509_merged.npz"); G=S["Sols"]["G"][it,:,:,iDir,0,0]; f=S["FreqDomains"].mean(axis=1)

def Norm(G,iRef=0):
    nf,na=G.shape
    for iFreq in range(nf):
        g0=G[iFreq,iRef]
        G[iFreq]*=g0.conj()/np.abs(g0)


def doFit_2D(GIn2D,nu,Mode=["Amp:Unit","Phase:TEC+CPhase"]):#,InParsGuess={"Phase:TEC+CPhase":(0,0)}):
        if isinstance(Mode,str):
            Mode=[Mode]
        # NSols,nChan
        GIn2D0=GIn2D.copy()
        G=GIn2D[-1,:]
        Flag=np.zeros(GIn2D0.shape,bool)
        dnu_Hz=None
        if nu.size>1:
            dnu=nu.flat[1:]-nu.flat[:-1]
            dnu_Hz=np.median(np.abs(dnu))
                
                
                
            
        for ThisMode in Mode:
            #print(ThisMode,G)
            if ":" in ThisMode:
                Op,Filter=ThisMode.split(":")
            else:
                Op=ThisMode
                
            if Op.startswith("MedianClip"):
                nx,ny,Th=Op.split("_")[1:]
                nx=int(nx)
                ny=int(ny)
                Th=float(Th)
                A=np.abs(GIn2D)
                minsig=0.01
                med=np.median(A)
                sig=1.4*np.median(np.abs(A-med))
                sig=np.max([minsig,sig])
                minsig=sig
                if not (nx==0 and ny==0):
                    med=scipy.ndimage.median_filter(A, (nx,ny))
                    sig=1.4*scipy.ndimage.median_filter(np.abs(A-med), (nx,ny))
                    sig[sig<minsig]=minsig
                indx,indy=np.where(np.abs(A-med)>Th*sig)
                Flag[indx,indy]=1
                G[Flag[-1,:]]=1
            elif Op.startswith("Clip"):
                #print("DoClip",Op)
                Th0,Th1=Op.split("_")[1:]
                Th0=float(Th0)
                Th1=float(Th1)
                A=np.abs(GIn2D)
                indx,indy=np.where((np.abs(A)<Th0)|(np.abs(A)>Th1))
                Flag[indx,indy]=1
                G[Flag[-1,:]]=1
            elif Op=="Amp" or Op=="Re" or Op=="Im":
                if Op=="Amp":
                    Proj=np.abs
                    backProj=lambda g,gf: gf*g/np.abs(g)
                elif Op=="Re":
                    Proj=lambda g: np.real(g)
                    backProj=lambda g,gf: np.real(gf) + 1j * g.imag
                elif Op=="Im":
                    Proj=lambda g: np.imag(g)
                    backProj=lambda g,gf: g.real + np.real(gf) * 1j
                else:
                    stop
                    
                if Filter=="Unit":
                    def FitFunc(gp):
                        return np.ones_like(gp)
                    G[:]=backProj(G[:],FitFunc(Proj(G[:])))
                elif Filter.startswith("PolyFreq"):
                    Order=int(Filter.split("_")[1])
                    # GIn[-1]=10
                    def FitFunc(gp):
                        FitMachine=ClassFitAmp.ClassFitAmp_1d(gp,
                                                              nu,
                                                              Order=Order,
                                                              LogMode=0,
                                                              RemoveMedianAmp=False)
                        gf=FitMachine.doSmooth()
                        return gf
                    G[:]=backProj(G[:],FitFunc(Proj(G[:])))
                    
                elif Filter.startswith("GaussFilter"):
                    Sx,Sy=Filter.split("_")[1:]
                    Sx=float(Sx)
                    if "MHz" in Sy and dnu is not None:
                        Sy=Sy.replace("MHz","")
                        Sy=float(Sy)*1e6/dnu_Hz
                        
                    Sy=float(Sy)
                    #print("SSS",Sy,dnu_Hz)
                    if Sx==0: Sx=1e-2
                    if Sy==0: Sy=1e-2
                    nx=np.max([3,int(3*Sx)])
                    ny=np.max([3,int(3*Sy)])
                    #gf=scipy.ndimage.gaussian_filter(np.abs(GIn2D), (int(nx),int(ny)))
                    x,y = np.mgrid[-3*nx:3*nx+1,-3*ny:3*ny+1]
                    gaussian = np.exp(-x**2/(2*Sx**2)-y**2/(2*Sy**2))
                    def FitFunc(GIn2Dp):
                        A=GIn2Dp
                        A[Flag]=0
                        gf=scipy.signal.fftconvolve(A, gaussian, mode='same')
                        Ones=np.ones_like(A)
                        Ones[Flag]=0
                        gfn=scipy.signal.fftconvolve(Ones, gaussian, mode='same')
                        gfn[gfn==0]=1
                        gf/=gfn
                        gf=gf[-1,:]
                        return gf
                    G[:]=backProj(G[:],FitFunc(Proj(GIn2D[:,:])))
                else:
                    stop
                    
                # import pylab
                # pylab.clf()
                # op0=np.abs
                # op1=np.angle
                # op2=np.real
                # op3=np.imag

                # pylab.figure(0)
                # pylab.subplot(2,2,1)
                # pylab.imshow(op0(GIn2D0),interpolation="nearest")
                # pylab.subplot(2,2,2)
                # pylab.imshow(op1(GIn2D0),interpolation="nearest")
                # pylab.subplot(2,2,3)
                # pylab.imshow(op2(GIn2D0),interpolation="nearest")
                # pylab.subplot(2,2,4)
                # pylab.imshow(op3(GIn2D0),interpolation="nearest")
                
                # pylab.figure(1)
                # pylab.subplot(2,2,1)
                # pylab.imshow(op0(GIn2D),interpolation="nearest")
                # pylab.subplot(2,2,2)
                # pylab.imshow(op1(GIn2D),interpolation="nearest")
                # pylab.subplot(2,2,3)
                # pylab.imshow(op2(GIn2D),interpolation="nearest")
                # pylab.subplot(2,2,4)
                # pylab.imshow(op3(GIn2D),interpolation="nearest")

                # pylab.figure(3)
                # pylab.subplot(2,2,1)
                # pylab.scatter(nu/1e6,op0(GIn2D0[-1,:]))
                # pylab.plot(nu/1e6,op0(GIn2D[-1,:]))
                # pylab.title(op0.__name__)
                # pylab.subplot(2,2,2)
                # pylab.scatter(nu/1e6,op1(GIn2D0[-1,:]))
                # pylab.plot(nu/1e6,op1(GIn2D[-1,:]))
                # pylab.title(op1.__name__)
                # pylab.subplot(2,2,3)
                # pylab.scatter(nu/1e6,op2(GIn2D0[-1,:]))
                # pylab.plot(nu/1e6,op2(GIn2D[-1,:]))
                # pylab.title(op2.__name__)
                # pylab.subplot(2,2,4)
                # pylab.scatter(nu/1e6,op3(GIn2D0[-1,:]))
                # pylab.plot(nu/1e6,op3(GIn2D[-1,:]))
                # pylab.title(op3.__name__)
                
                # pylab.draw()
                # pylab.show(block=False)
                # pylab.pause(1)
                
            elif Op=="Phase":
                if Filter=="Zero":
                    G[:]=np.abs(G[:])
                else:
                    #print(G.shape)
                    FitMachine=ClassFitTEC.ClassFitTEC_1D(G,nu,Mode=Filter.split("+"))
                    _,t0,c0=ClassFitTEC.EstimateThisTECTime(G.flatten(),nu)
                    # X0=InParsGuess.get(ThisMode,None)
                    
                    TEC0CPhase0=np.zeros((len(Filter.split("+")),1),np.float32)
                    TEC0CPhase0[0,0]=t0#+np.random.randn(1)[0]*0.001
                    if "CPhase" in Mode:
                        TEC0CPhase0[1,0]=c0
                    FitMachine.setX0(TEC0CPhase0)
                    X=FitMachine.doFit()
                    G[:]=FitMachine.X2G(X)
        
        # if G.max()>2:
        #     stop
        return G[:]
                



        
def doFit_1D(GIn,nu,Mode=["Amp:Unit","Phase:TEC+CPhase"]):#,InParsGuess={"Phase:TEC+CPhase":(0,0)}):
        for ThisMode in Mode:
            Op,Filter=ThisMode.split(":")
            if Op=="Amp":
                if Filter=="Unit":
                    GOut=GIn/np.abs(GIn)
                elif Filter.startswith("Poly"):
                    Order=int(Filter.split("_")[1])
                    #GIn[-1]=10
                    g=GIn
                    FitMachine=ClassFitAmp.ClassFitAmp_1d(np.abs(GIn.copy()),
                                                          nu,
                                                          Order=Order,
                                                          LogMode=0,
                                                          RemoveMedianAmp=False)
                    gf=FitMachine.doSmooth()
                    # print "Done %i"%iDir
                    gf=gf*g/np.abs(g)
                    GOut=gf
                else:
                    stop
                    
                # import pylab
                # pylab.clf()
                # print(GIn)
                # pylab.plot(np.abs((GIn)))
                # pylab.plot(np.abs(gf),color="black")
                # pylab.plot([0,3],[0,0],color="black",ls=":")
                # #pylab.plot(np.abs((GOut)))
                # pylab.draw()
                # pylab.show(block=False)
                # pylab.pause(0.1)
                
            elif Op=="Phase":
                if Filter=="Zero":
                    GOut=GIn/np.abs(GIn)
                else:
                    FitMachine=ClassFitTEC.ClassFitTEC_1D(GIn,nu,Mode=Filter.split("+"))
                    _,t0,c0=ClassFitTEC.EstimateThisTECTime(GIn.flatten(),nu)
                    # X0=InParsGuess.get(ThisMode,None)
                    
                    TEC0CPhase0=np.zeros((len(Mode),1),np.float32)
                    TEC0CPhase0[0,0]=t0#+np.random.randn(1)[0]*0.001
                    if "CPhase" in Mode:
                        TEC0CPhase0[1,0]=c0
                    FitMachine.setX0(TEC0CPhase0)
                    X=FitMachine.doFit()
                    GOut=FitMachine.X2G(X)
            GIn=GOut
        return GIn
                

