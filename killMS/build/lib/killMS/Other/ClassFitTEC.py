from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as np
from DDFacet.Other import logger
log=logger.getLogger("ClassFitTEC")
import killMS.Array.ModLinAlg
K=8.4479745e9
import scipy.sparse
from DDFacet.Other import ClassTimeIt
from scipy.optimize import least_squares

logger.setSilent(["ClassFitTEC"])

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
    

def test(G,f):
    # nf,na=G.shape
    # #na=3
    # t=np.random.randn(na)*0.01
    # c=np.random.randn(na)*np.pi/10
    # G=TECToZ(t.reshape((1,-1)),c.reshape((1,-1)),f.reshape((-1,1)))

    TECMachine=ClassFitTEC(G,f)
    #TECMachine.DoFit()
    TECMachine.findX0()
    TECMachine.doFit()


NTEC=101
NConstPhase=51
TECGridAmp=0.1
TECGrid,CPhase=np.mgrid[-TECGridAmp:TECGridAmp:NTEC*1j,-np.pi:np.pi:NConstPhase*1j]
Z=None
CentralFreqs=None

def EstimateThisTECTime(g,nu):
    g0=g/np.abs(g)

    W=np.ones(g0.shape,np.float32)
    W[g==1.]=0

    global Z,CentralFreqs
    if Z is None or np.allclose(CentralFreqs,nu.ravel()):
        CentralFreqs=nu.ravel()
        Z=TECToZ(TECGrid.reshape((-1,1)),CPhase.reshape((-1,1)),CentralFreqs.reshape((1,-1)))
        
    for iTry in range(5):
        R=(g0.reshape((1,-1))-Z)*W.reshape((1,-1))
        Chi2=np.sum(np.abs(R)**2,axis=1)
        iTec=np.argmin(Chi2)
        rBest=R[iTec]
        if np.max(np.abs(rBest))==0: break
        Sig=np.sum(np.abs(rBest*W))/np.sum(W)
        ind=np.where(np.abs(rBest*W)>5.*Sig)[0]
        if ind.size==0: break
        W[ind]=0

    t0=TECGrid.ravel()[iTec]
    c0=CPhase.ravel()[iTec]
    
    gz=np.abs(g)*TECToZ(t0,c0,nu)
    return gz,t0,c0


class ClassFitTEC_1D():
    def __init__(self,GIn,nu,
                 Tol=5e-2,Incr=1,
                 #Mode=["TEC","CPhase"],
                 Mode=["TEC"]):

        self.CentralFreqs=nu.flatten()
        na=1
        nch=nu.size
        if GIn.size!=nch: stop
        
        self.GIn=GIn.copy().reshape((nch,1))
        self.Mode=Mode
        self.CrossMode=0
        self.TEC0CPhase0=np.zeros((len(Mode),),np.float32)

    def setX0(self,TEC0CPhase0):
        self.TEC0CPhase0=TEC0CPhase0.ravel()

    def doFit(self):
        G=self.GIn.T.copy()
        na,nch=G.shape
        
        if self.CrossMode:
            A0,A1=np.mgrid[0:na,0:na]
            gg_meas=G[A0.ravel(),:]*G[A1.ravel(),:].conj()
            gg_meas_reim=np.array([gg_meas.real,gg_meas.imag]).ravel()[::self.incrCross]
        else:
            self.incrCross=1
            A0,A1=np.mgrid[0:na],None
            gg_meas=G[A0.ravel(),:]
            gg_meas_reim=np.array([gg_meas.real,gg_meas.imag]).ravel()[::self.incrCross]

        iIter=np.array([0])
        tIter=np.array([0],np.float64)
        def _f_resid(TecConst,A0,A1,ggmeas,iIter,tIter):
            T2=ClassTimeIt.ClassTimeIt("resid")
            T2.disable()
            TEC,CPhase=TecConst.reshape((2,na))
            GThis=TECToZ(TEC.reshape((-1,1)),CPhase.reshape((-1,1)),self.CentralFreqs.reshape((1,-1)))
            #T2.timeit("1")
            if self.CrossMode:
                gg_pred=GThis[A0.ravel(),:]*GThis[A1.ravel(),:].conj()
            else:
                gg_pred=GThis[A0.ravel(),:]

            #T2.timeit("2")
            gg_pred_reim=np.array([gg_pred.real,gg_pred.imag]).ravel()[::self.incrCross]
            #T2.timeit("3")
            r=(ggmeas-gg_pred_reim).ravel()
            #print r.shape
            #T2.timeit("4")
            #return np.angle((ggmeas-gg_pred).ravel())
            #print np.mean(np.abs(r))
            iIter+=1
            #tIter+=T2.timeit("all")
            #print iIter[0]
            return r

        #print _f_resid(TEC0CPhase0,A0,A1,ggmeas)
        
        Sol=least_squares(_f_resid,
                          self.TEC0CPhase0.ravel(),
                          #method="trf",
                          method="lm",
                          args=(A0,A1,gg_meas_reim,iIter,tIter),
                          ftol=1e-2,gtol=1e-2,xtol=1e-2)#,ftol=1,gtol=1,xtol=1,max_nfev=1)
        #Sol=leastsq(_f_resid, TEC0CPhase0.ravel(), args=(A0,A1,gg_meas_reim,iIter),ftol=1e-2,gtol=1e-2,xtol=1e-2)
        #T.timeit("Done %3i %3i %5i"%(it,iDir,iIter[0]))
        #print "total time f=%f"%tIter[0]
        #TEC,CPhase=Sol.x.reshape((2,na))

        return Sol.x.ravel()#TEC,CPhase
    
        # TEC-=TEC[0]
        # CPhase-=CPhase[0]
        # #GThis=np.abs(self.G.T)*TECToZ(TEC.reshape((-1,1)),CPhase.reshape((-1,1)),self.CentralFreqs.reshape((1,-1)))
        # #T.timeit("done")
        
        # return GThis.T,TEC,CPhase

    def X2G(self,X):
        Mode=self.Mode
        na=1
        if "CPhase" in Mode:
            TEC,CPhase=X.reshape((len(Mode),na))
        else:
            TEC,=X.reshape((len(Mode),na))
            CPhase=np.zeros((1,na),np.float32)
        G=self.GIn.reshape((-1,1))
        GOut=np.abs(G).T*TECToZ(TEC.reshape((-1,1)),CPhase.reshape((-1,1)),self.CentralFreqs.reshape((1,-1)))
        return GOut.ravel()


    

class ClassFitTEC_Cross():
    def __init__(self,gains,nu,Tol=5e-2,Incr=1,
                 #Mode=["TEC","CPhase"],
                 Mode=["TEC"]):
        self.nf,self.na=gains.shape
        self.Mode=Mode
        self.LMode=len(Mode)
        self.G=gains.copy()
        Norm(self.G)
        self.Tol=Tol

        Ga=np.abs(self.G)
        Ga[Ga==0]=1
        self.G/=Ga
        self.CentralFreqs=self.nu=nu
        self.NFreq=nu.size
        na=self.na
        self.nbl=(na**2-na)//2
        self.CurrentX=None
        log.print("Number of Antennas: %i"%self.na)
        log.print("Number of Freqs:    %i"%nu.size)
        log.print("Number of Points:   %i"%(nu.size*self.na**2))

        self.Y=np.array([(self.G[iFreq].reshape((-1,1))*self.G[iFreq].conj().reshape((1,-1))).ravel() for iFreq in range(self.NFreq)]).ravel()
        self.nu_Y=np.array([(self.nu[iFreq]*np.ones((self.na,self.na)).ravel()) for iFreq in range(self.NFreq)]).ravel()
        self.A0=np.array([(np.mgrid[0:na:1,0:na:1][0]).ravel() for iFreq in range(self.NFreq)]).ravel()
        self.A1=np.array([(np.mgrid[0:na:1,0:na:1][1]).ravel() for iFreq in range(self.NFreq)]).ravel()
        

        self.Incr=Incr
        Mask=np.where(self.A1>self.A0)[0]
        self.Mask=Mask
        self.Y=self.Y[Mask][::self.Incr]
        self.nu_Y=self.nu_Y[Mask][::self.Incr]
        self.A0=self.A0[Mask][::self.Incr]
        self.A1=self.A1[Mask][::self.Incr]




        self.x0=None
        self.indA0=[np.where(self.A0==iAnt)[0] for iAnt in range(na)]
        self.indA1=[np.where(self.A1==iAnt)[0] for iAnt in range(na)]



    def doFit(self,NIter=100):
        if self.x0 is None and self.CurrentX is None:
            self.CurrentX=np.zeros((self.LMode*self.na,),np.float32)+1e-10
            #self.CurrentX=np.random.randn(2*self.na)

        self.Current_iIter=0
        for iIter in range(NIter):
            self.doLMIter()
            #self.Plot()
            self.Current_iIter=iIter
            if self.Diff<self.Tol:
                log.print("Convergence in %i steps"%(iIter+1))
                break

        return self.CurrentX

    def GiveGPredict(self,X):
        t=X[0:self.na].reshape((1,-1))
        c=None
        if "CPhase" in self.Mode:
            c=X[self.na:].reshape((1,-1))
        z=TECToZ(t,c,self.nu.reshape((-1,1)))
        return z

    def doLMIter(self):
        T=ClassTimeIt.ClassTimeIt()
        T.disable()
        #J,H=
        self.giveJacobianHessian()
        T.timeit("J, H")
        z=self.GiveGPredict(self.CurrentX)
        Y=np.array([(z[iFreq].reshape((-1,1))*z[iFreq].conj().reshape((1,-1))).ravel() for iFreq in range(self.NFreq)]).ravel()

        r=self.Y-Y[self.Mask][::self.Incr]

        v=self.JHy(r)
        H=self.DiagJHJ()
        T.timeit("diff")

        Hinv=killMS.Array.ModLinAlg.invSVD(H)
        T.timeit("inv")
        
        X = self.CurrentX + np.real(np.dot(Hinv,v.reshape((-1,1))).ravel())
        

        xx=self.CurrentX.copy()
        xx[xx==0]=1e-6
        self.Diff=np.max(np.abs((X-xx)/xx))



        z0=self.GiveGPredict(self.CurrentX)
        Norm(z0)
        self.CurrentX=X
        z=self.GiveGPredict(self.CurrentX)
        Norm(z)
        self.Diff=np.max(np.abs(np.angle(z*z0.conj())))
        #print self.Diff
        
        return 

        # HinvJH=np.dot(scipy.sparse.coo_matrix(Hinv),J.T.conj())
        # T.timeit("HinvJH")
        # HinvJHy=np.dot(HinvJH,scipy.sparse.coo_matrix(r.reshape((-1,1))))
        # T.timeit("HinvJHy")

        # self.CurrentX+=np.real(np.array(HinvJHy.todense())).ravel()
        # T.timeit("X")
        # #self.CurrentX+=np.real(Dot(Hinv,J.T.conj(),r.reshape((-1,1))).ravel())
        
    def setX0(self,x0):
        self.CurrentX=x0

    def findX0(self):
        NTEC=101
        NConstPhase=51
        TECGridAmp=0.1
        if self.LMode==2:
            TECGrid,CPhase=np.mgrid[-TECGridAmp:TECGridAmp:NTEC*1j,-np.pi:np.pi:NConstPhase*1j]
            Z=TECToZ(TECGrid.reshape((-1,1)),CPhase.reshape((-1,1)),self.CentralFreqs.reshape((1,-1)))
        elif self.LMode==1:
            TECGrid,CPhase=np.mgrid[-TECGridAmp:TECGridAmp:NTEC*1j],None
            Z=TECToZ(TECGrid.reshape((-1,1)),CPhase,self.CentralFreqs.reshape((1,-1)))

        self.Z=Z
        self.TECGrid,self.CPhase=TECGrid,CPhase

        self.CurrentX=np.zeros((self.LMode,self.na),np.float32)

        for iAnt in range(self.na):
            g=self.G[:,iAnt]
            g0=g/np.abs(g)
            W=np.ones(g0.shape,np.float32)
            W[g==1.]=0
            Z=self.Z
            for iTry in range(5):
                R=(g0.reshape((1,-1))-Z)*W.reshape((1,-1))
                Chi2=np.sum(np.abs(R)**2,axis=1)
                iTec=np.argmin(Chi2)
                rBest=R[iTec]
                if np.max(np.abs(rBest))==0: break
                Sig=np.sum(np.abs(rBest*W))/np.sum(W)
                ind=np.where(np.abs(rBest*W)>5.*Sig)[0]
                if ind.size==0: break
                W[ind]=0
            self.CurrentX[0,iAnt]=self.TECGrid.ravel()[iTec]
            if "CPhase" in self.Mode:
                self.CurrentX[1,iAnt]=self.CPhase.ravel()[iTec]
        self.CurrentX=self.CurrentX.ravel()
        
    def Plot(self):
        z=self.GiveGPredict(self.CurrentX)
        Norm(z)

        import pylab
        pylab.clf()
        pylab.plot(self.nu,np.angle(self.G),color="black")
        pylab.plot(self.nu,np.angle(z),color="gray")
        pylab.draw()
        pylab.show(False)
        pylab.pause(0.1)

    def JHy(self,y):
        T=ClassTimeIt.ClassTimeIt("JHy")
        T.disable()
        v=np.zeros((self.LMode,self.na),np.complex64)
        for iAnt in range(self.na):
            v[0,iAnt]+=np.sum(self.J_TEC[self.indA0[iAnt]].conj()*y[self.indA0[iAnt]])
            v[0,iAnt]+=np.sum(-self.J_TEC[self.indA1[iAnt]].conj()*y[self.indA1[iAnt]])

            if "CPhase" in self.Mode:
                v[1,iAnt]+=np.sum(self.J_Phase[self.indA0[iAnt]].conj()*y[self.indA0[iAnt]])
                v[1,iAnt]+=np.sum(-self.J_Phase[self.indA1[iAnt]].conj()*y[self.indA1[iAnt]])

        v=v.ravel()
        T.timeit("Prod")
        # PSparse=np.array(np.dot(self.J.T.conj(),scipy.sparse.coo_matrix(y.reshape((-1,1)))).todense()).ravel()
        # T.timeit("PSparse")
        return v

    def DiagJHJ(self):
        T=ClassTimeIt.ClassTimeIt("JHy")
        T.disable()
        H=np.zeros((self.LMode,self.na,self.LMode,self.na),np.complex64)
        for iAnt in range(self.na):
            #Jt[self.indA0[iAnt],iAnt]=J_TEC[self.indA0[iAnt]]
            #Jt[self.indA1[iAnt],iAnt]=-J_TEC[self.indA1[iAnt]]
            H[0,iAnt,0,iAnt]+=np.sum(self.J_TEC[self.indA0[iAnt]].conj()*self.J_TEC[self.indA0[iAnt]])
            H[0,iAnt,0,iAnt]+=np.sum(self.J_TEC[self.indA1[iAnt]].conj()*self.J_TEC[self.indA1[iAnt]])

            if "CPhase" in self.Mode:
                H[0,iAnt,1,iAnt]+=np.sum(self.J_TEC[self.indA0[iAnt]].conj()*self.J_Phase[self.indA0[iAnt]])
                H[0,iAnt,1,iAnt]+=np.sum(self.J_TEC[self.indA1[iAnt]].conj()*self.J_Phase[self.indA1[iAnt]])
                
                H[1,iAnt,0,iAnt]+=np.sum(self.J_Phase[self.indA0[iAnt]].conj()*self.J_TEC[self.indA0[iAnt]])
                H[1,iAnt,0,iAnt]+=np.sum(self.J_Phase[self.indA1[iAnt]].conj()*self.J_TEC[self.indA1[iAnt]])
            
                H[1,iAnt,1,iAnt]+=np.sum(self.J_Phase[self.indA0[iAnt]].conj()*self.J_Phase[self.indA0[iAnt]])
                H[1,iAnt,1,iAnt]+=np.sum(self.J_Phase[self.indA1[iAnt]].conj()*self.J_Phase[self.indA1[iAnt]])

        T.timeit("H")

        H=H.reshape((self.LMode*self.na,self.LMode*self.na))
        return H
        A=np.log10(np.abs(self.H))
        
        B=np.log10(np.abs(H))
        vmin,vmax=A.min(),A.max()
        import pylab
        pylab.clf()
        pylab.subplot(1,2,1)
        pylab.imshow(A,interpolation="nearest",vmin=vmin,vmax=vmax)
        pylab.colorbar()
        pylab.subplot(1,2,2)
        pylab.imshow(B,interpolation="nearest",vmin=vmin,vmax=vmax)
        pylab.colorbar()
        pylab.draw()
        pylab.show(False)

        stop
        return H


    def giveJacobianHessian(self):
        T=ClassTimeIt.ClassTimeIt("J")
        J=np.zeros((self.Y.size,self.na*2),np.complex64)
        Jt=J[:,0:self.na]
        Jc=J[:,self.na:]
        
        TEC=self.CurrentX[0:self.na]
        dTEC=TEC[self.A0]-TEC[self.A1]
        if "CPhase" in self.Mode:
            ConstPhase=self.CurrentX[self.na:]
            dConstPhase=ConstPhase[self.A0]-ConstPhase[self.A1]
        else:
            dConstPhase=0
            
        Phase=K/self.nu_Y*dTEC+dConstPhase
        Z=np.exp(1j*Phase)

        self.J_TEC=J_TEC=1j*K/self.nu_Y*Z
        self.J_Phase=J_Phase=1j*Z
        # print(J_TEC,self.J_Phase)
        
        return
        T.timeit("first")
        for iAnt in range(self.na):
            Jt[self.indA0[iAnt],iAnt]=J_TEC[self.indA0[iAnt]]
            Jt[self.indA1[iAnt],iAnt]=-J_TEC[self.indA1[iAnt]]

            Jc[self.indA0[iAnt],iAnt]=J_Phase[self.indA0[iAnt]]
            Jc[self.indA1[iAnt],iAnt]=-J_Phase[self.indA1[iAnt]]

        T.timeit("build")
        self.J=J
        self.Jsp=Jsp=scipy.sparse.coo_matrix(J)
        T.timeit("sp")
        # import pylab
        # pylab.clf()
        # pylab.subplot(1,2,1)
        # pylab.imshow(Jt.real,interpolation="nearest",aspect="auto")
        # pylab.subplot(1,2,2)
        # pylab.imshow(Jc.real,interpolation="nearest",aspect="auto")
        # pylab.draw()
        # pylab.show(False)
        # stop

        #print np.count_nonzero(J)/float(J.size)

        T.timeit("prod")
        H=np.array(np.dot(Jsp.T.conj(),Jsp).todense())
        self.H=H
        T.timeit("Hsp")
        
        return J,H


