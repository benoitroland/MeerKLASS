#!/usr/bin/env python
# Upload to archive directory
# To store on Leiden systems in /disks/paradata/shimwell/LoTSS-DR2/archive

# version 2 -- simple rsync wrapper with timeout

# *** This script is obsolute and should not be used
# *** directly. rclone_upload.py makes use of some of its functions
# *** 

from __future__ import print_function
from builtins import range
from subprocess import call
import os
import glob
from time import sleep
from surveys_db import use_database,SurveysDB,tag_field

def myglob(g,workdir):
    f=glob.glob(workdir+'/'+g)
    return [os.path.basename(file) for file in f]

def images(rootname,workdir):
    l=[rootname+'.'+f+'.fits' for f in ['app.restored','int.restored','int.model','int.residual']]
    n=1
    while os.path.isfile(workdir+'/'+rootname+'.mask%02i.fits' % n):
        n+=1
    if n>1:
        l.append(rootname+'.mask%02i.fits' % (n-1))
    l.append(rootname+'.DicoModel')
    l.append(rootname+'.tessel.reg')
    return l

def shiftimages(rootname):
    return [rootname+f+'.fits' for f in ['_shift.app.facetRestored','_shift.int.facetRestored']]

def do_rsync(name,basedir,f):
    other=False
    if use_database():
        with SurveysDB(readonly=True) as sdb:
            idd=sdb.get_field(name)
        if not (idd['lotss_field']>0) or idd['proprietary_date'] is not None:
            other=True
            
    workdir=basedir+'/'+name

    if os.environ['DDF_PIPELINE_CLUSTER']!='paracluster':
        target=os.environ['DDF_PIPELINE_LEIDENUSER']+'@ssh.strw.leidenuniv.nl:'
    else:
        target=''

    f=[fle for fle in f if os.path.exists(workdir+'/'+fle)]
        
    while True:
        s='cd '+workdir+'; rsync -avz --progress --safe-links --perms --chmod=ugo+rX --partial --timeout=60 '+' '.join(f)+' '+target+'/disks/paradata/shimwell/LoTSS-DR2/archive/'+name
        if other:
            s=s.replace('DR2/archive','DR2/archive_other')
        print('Running command:',s)
        retval=call(s,shell=True)
        if retval==0:
            break
        print('Non-zero return value',retval)
        if retval!=30:
            raise RuntimeError('rsync failed unexpectedly')
        sleep(10)

def do_delete_keep(name,basedir):
    target=os.environ['DDF_PIPELINE_LEIDENUSER']+'@ssh.strw.leidenuniv.nl'
    if not os.path.isdir(basedir+'/'+name+'/KEEP'):
        raise RuntimeError('KEEP directory does not exist')
    else:
        command='cd /disks/paradata/shimwell/LoTSS-DR2/archive/'+name+'; rm -rf summary.txt logs '
        g=glob.glob(basedir+'/'+name+'/KEEP/*')
        fl=[os.path.basename(f) for f in g]
        command+=' '.join(fl)
        print(command)
        os.system('ssh '+target+' "'+command+'"')
        
def do_upload_compressed(name,basedir):
    workdir=basedir+'/'+name
    f=myglob('*.archive',workdir)
    f+=images('image_full_ampphase_di_m.NS',workdir)
    f+=images('image_full_low_m',workdir)

    do_rsync(name,basedir,f)
    with SurveysDB() as sdb:
        idd=sdb.get_field(name)
        idd['archive_version']=2
        sdb.set_field(idd)

def do_upload_vlow(name,basedir):
    workdir=basedir+'/'+name
    f=myglob('*.fz',workdir)
    do_rsync(name,basedir,f)
    with SurveysDB() as sdb:
        idd=sdb.get_field(name)
        idd['vlow_reprocess']="Uploaded"
        sdb.set_field(idd)

        
def do_upload(name,basedir,skipstokes=False):

    workdir=basedir+'/'+name

    f=['summary.txt']
    f+=['logs','SOLSDIR']
    f+=['astromap.fits','facet-offset.txt']
    f+=myglob('*-fit_state.pickle',workdir)
    f+=['image_dirin_SSD_m.npy.ClusterCat.npy']
    f+=myglob('DynSpecs*.tgz',workdir)
    f+=myglob('*.png',workdir)
    f+=myglob('*mslist*txt',workdir)
    f+=myglob('DDS*smoothed*.npz',workdir)
    f+=myglob('DDS*full_slow*.npz',workdir)
    f+=myglob('*crossmatch-results*',workdir)
    f+=myglob('*crossmatch-*.fits',workdir)
    f+=images('image_full_ampphase_di_m.NS',workdir)
    f+=images('image_full_low_m',workdir)
    f+=shiftimages('image_full_ampphase_di_m.NS')
    for i in range(3):
        f+=shiftimages('image_full_ampphase_di_m.NS_Band%i' %i)
    if not skipstokes:
        f+=myglob('image_full_low_stokesV.dirty.*',workdir)
        f+=myglob('image_full_low_stokesV.SmoothNorm.fits',workdir)
        f+=myglob('image_full_low_QU.cube.*',workdir)
        f+=myglob('image_full_vlow_QU.cube.*',workdir)
    f+=myglob('*.archive',workdir)

    do_rsync(name,basedir,f)
    compressed_done=(len(myglob('*.archive',workdir))>0)
    with SurveysDB() as sdb:
        idd=sdb.get_field(name)
        if idd['proprietary_date'] is not None:
            idd['status']='Proprietary'
        else:
            idd['status']='Archived'
        tag_field(sdb,idd,workdir=workdir)
        sdb.set_field(idd)
    
if __name__=='__main__':
    import sys
    do_upload(sys.argv[1],sys.argv[2])
