#!/usr/bin/env python

from __future__ import print_function
from surveys_db import SurveysDB
import glob
import os
import sys

# remove archived data! Yikes!

with SurveysDB() as sdb:
    sdb.cur.execute('select * from fields where (status="Verified" or status="Archived" or status="Proprietary") and archive_version=4 and clustername="Herts" and location like "%beegfs%"')
    results=sdb.cur.fetchall()

dryrun=True
try:
    if sys.argv[1]=="delete":
        dryrun=False
except:
    pass

if len(sys.argv)>2:
    selected=sys.argv[2:]
else:
    selected=None

if selected:
    for i in range(len(selected)):
        if selected[i]=='.':
            selected[i]=os.getcwd().split('/')[-1]
    
deleted=[]
for r in results:
    location=r['location']
    if not location:
        continue
    if not os.path.isdir(location):
        continue
    print(r['id'],r['location'])
    g=len(glob.glob(location+'/*.archive'))
    g2=len(glob.glob(location+'/*.fz'))
    print(r['id'], g, g2)
    if selected is not None and r['id'] not in selected:
        print('Failed selection, skipping')
        continue
    if g>3 or selected is not None:
        print('Deleting files from',r['id'],'at',location)
        print('In DR2 dir:',len(glob.glob('/data/lofar/DR2/fields/'+r['id']+'/*NS*int*.fits')),len(glob.glob('/data/lofar/DR2/fields/'+r['id']+'/*.fz')))
        if dryrun:
            print('Only not really')
        else:
            os.system('rm -r '+location+'/*')
            os.system('rmdir '+location)
            #os.system('rm ~/ddfp-'+r['id']+'.*')
            os.system('rm /beegfs/car/mjh/torque/ddfp-'+r['id']+'.*')
        rr=r
        rr['location']=''
        deleted.append(rr)
    else:
        print('Not removing',r['id'])

if not dryrun:
    with SurveysDB() as sdb:
        for r in deleted:
            sdb.set_field(r)
