#!/usr/bin/python

# clone of monitor.py to do vlow monitoring

from __future__ import print_function
from time import sleep
import datetime
from surveys_db import SurveysDB,get_next
import os
import threading
from make_vlow_stokesI import do_download,update_status
from make_vlow_stokesI_sub import do_upload
import glob

cluster=os.environ['DDF_PIPELINE_CLUSTER']
download_thread=None
download_name=None
tidy_thread=None
upload_name=None
upload_thread=None
basedir='/beegfs/car/mjh/vlow'
totallimit=20

os.chdir(basedir)

while True:

    with SurveysDB(readonly=True) as sdb:
        sdb.cur.execute('select * from fields where vlow_image is not NULL order by end_date desc')
        result=sdb.cur.fetchall()
        sdb.cur.execute('select * from fields where vlow_image="Required" order by decl desc')
        result2=sdb.cur.fetchall()
        if len(result2)>0:
            nextfield=result2[0]['id']
        else:
            nextfield=None

    d={}
    for r in result:
        status=r['vlow_image']
        if status in d:
            d[status]=d[status]+1
        else:
            d[status]=1
    print('\n\n-----------------------------------------------\n\n')
    print('DDF-pipeline reprocess status on cluster',cluster)
    print(datetime.datetime.now())
    print()
    failed=0
    for k in sorted(d.keys()):
        print("%-20s : %i" % (k,d[k]))
        if 'ailed' in k:
            failed+=d[k]

    print()
    ksum=len(glob.glob(basedir+'/*'))-failed
    if ksum<0: ksum=0
    print(ksum,'live files out of',totallimit)
    print('Next field to work on is',nextfield)
            
    if download_thread is not None:
        print('Download thread is running (%s)' % download_name)
    if upload_thread is not None:
        print('Upload thread is running (%s)' % upload_name)
        
    if download_thread is not None and not download_thread.isAlive():
        print('Download thread seems to have terminated')
        # run the qsub here!
        download_thread=None

    if upload_thread is not None and not upload_thread.isAlive():
        print('Upload thread seems to have terminated')
        upload_thread=None

    if ksum<totallimit and nextfield is not None and download_thread is None:
        download_name=nextfield
        print('We need to download a new file (%s)!' % download_name)
        download_thread=threading.Thread(target=do_download, args=(download_name,basedir))
        download_thread.start()

    if upload_thread is None and 'WSCLEAN complete' in d:
        for r in result:
            if r['vlow_image']=='WSCLEAN complete':
                upload_name=r['id']
                print('We need to upload a new file (%s)!' % upload_name)
                upload_thread=threading.Thread(target=do_upload, args=(upload_name,basedir))
                upload_thread.start()
                break
        
    if 'Downloaded' in d:
        for r in result:
            if r['vlow_image']=='Downloaded':
                field=r['id']
                command="qsub -v FIELD=%s -N vlow-%s /home/mjh/pipeline-master/ddf-pipeline/torque/vlow_image.qsub" % (field, field)
                if os.system(command)==0:
                    update_status(field,"Queued")

    if 'DDF complete' in d:
        for r in result:
            if r['vlow_image']=='DDF complete':
                field=r['id']
                command="qsub -v FIELD=%s -N subtract-%s /home/mjh/pipeline-py3/ddf-pipeline/torque/vlow_subtract.qsub" % (field, field)
                if os.system(command)==0:
                    update_status(field,"Subtract queued")

    if 'Subtracted' in d:
        for r in result:
            if r['vlow_image']=='Subtracted':
                field=r['id']
                command="qsub -v FIELD=%s -N subimage-%s /home/mjh/pipeline-py3/ddf-pipeline/torque/vlow_image_subtracted.qsub" % (field, field)
                if os.system(command)==0:
                    update_status(field,"WSCLEAN queued")

    if 'Uploaded' in d:
        for r in result:
            if r['vlow_image']=='Uploaded':
                field=r['id']
                print('Tidying uploaded directory for',field)
                target='/data/lofar/DR2/fields/'+field
                for f in ['image_full_vlow_nocut_m.app.restored.fits','image_full_vlow_nocut_m.int.restored.fits','WSCLEAN_low-MFS-image.fits','WSCLEAN_low-MFS-image-int.fits']:
                    command='cp '+basedir+'/'+field+'/'+f+' '+target
                    print('running',command)
                    os.system(command)
                command='rm -r '+basedir+'/'+field
                os.system(command)
                update_status(field,'Archived')

    print('\n\n-----------------------------------------------\n\n')
        
    sleep(60)
