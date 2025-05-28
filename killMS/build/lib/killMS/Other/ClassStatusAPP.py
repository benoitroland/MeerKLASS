import multiprocessing
from DDFacet.Array import shared_dict
from DDFacet.Other import AsyncProcessPool as APP
from DDFacet.Other import Multiprocessing
import time
import numpy as np
NCPU=int(multiprocessing.cpu_count())

class ClassStatusAPP():
    def __init__(self,path=""):
        self.DicoWorkers=shared_dict.SharedDict("%sDicoWorkers"%path,reset=0)
        self.nTotWorkers=NCPU

        
    def initWorkersInfo(self):
        Names=sorted(list(set(self.getWorkersNames())))
        self.nTotWorkers=len(Names)+1
        self.DicoNameToID={}
        for iName,Name in enumerate(Names):
            self.DicoNameToID[Name]=iName+1
            
        #self.DicoWorkers["WorkerCurrentTask"]=np.zeros((self.nTotWorkers,),dtype=("|S200"))
        
        
    def setStatus(self,Str):
        #print("[%s] set %s"%(multiprocessing.current_process().name,Str))
        self.DicoWorkers.reload()
        self.DicoWorkers[multiprocessing.current_process().name]=Str

    def Print(self):
        self.DicoWorkers.reload()
        print(list(self.DicoWorkers.keys()))
        for k in list(self.DicoWorkers.keys()):
            print(k,self.DicoWorkers[k])

    def getWorkersNames(self):
        for iWorker in range(NCPU*2):
            APP.APP.runJob("_workerName:%i"%(iWorker),
                           self._workerName,
                           args=())#,serial=True)
        APP.APP.runJob("_workerName:%i"%(iWorker+1),
                       self._workerName,
                       args=(),
                       io=0)#,serial=True)
        Names=APP.APP.awaitJobResults("_workerName*")
        return Names
    
    def _workerName(self):
        N=multiprocessing.current_process().name
        time.sleep(0.2)
        return N


CS_APP=ClassStatusAPP()
        
import os

def getListDirs():
    files = os.listdir("/dev/shm")
    files.sort(key=lambda x: os.path.getctime("/dev/shm/%s"%x))
    print(files)
    return files

def startMain():
    CS=ClassStatusAPP(path=getListDirs()[-1])
    APP.APP.init(ncpu=NCPU, affinity=0, num_io_processes=1)
    APP.APP.registerJobHandlers(CS)
    APP.APP.startWorkers()
    APP.APP.awaitWorkerStart()
    CS.initWorkersInfo()
    CS.setStatus("a+=%s"%multiprocessing.current_process().name)
    APP.APP.shutdown()
    Multiprocessing.cleanupShm()
    
    
# if __name__=="__main__":
#     path=getListDirs()[-2]
#     stop
#     CS=ClassStatusAPP(path=)
    
#     while True:
#         CS.Print()
#         print(CS.DicoWorkers.path)
#         time.sleep(0.1)
