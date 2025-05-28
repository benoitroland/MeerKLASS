#!/bin/bash
#. /home/cyril.tasse/.script_init/Python

########################################################################
########################### killMS
########################################################################
if (($# < 1))
then
    echo "Specify a path"
    return
fi

export DEVDIR=$1
if [ ! -d "$DEVDIR" ];
then
    echo $DEVDIR" directory does not exist."
    return
fi


export VE_FOLDER=$1
export OPENBLAS_NUM_THREADS=1
export OPENBLAS_MAX_THREADS=1
#export TMPDIR=/data/cyril.tasse/'tmp'
#mkdir -p $TMPDIR
export PYTHONPATH_FIRST=1
export NUMEXPR_MAX_THREADS=96
export PYTHONHASHSEED=0

export THIS_USER
THIS_USER=$USER
#export DDF_PIPELINE_CATALOGS=/home/cyril.tasse/CATALOGS
export KILLMS_DIR=$VE_FOLDER
export DDFACET_DIR=$KILLMS_DIR
export DDF_DIR=$DDFACET_DIR
echo -e Source directory for killMS: $BLEU $KILLMS_DIR $NORMAL
export PYTHONPATH=$KILLMS_DIR:$PYTHONPATH
export PYTHONPATH=$KILLMS_DIR/DDFacet:$PYTHONPATH
export PYTHONPATH=$KILLMS_DIR/killMS:$PYTHONPATH
export PYTHONPATH=$KILLMS_DIR/FindCluster:$PYTHONPATH
export PYTHONPATH=$KILLMS_DIR/DynSpecMS:$PYTHONPATH
export PYTHONPATH=$KILLMS_DIR/lotss-query:$PYTHONPATH
export PYTHONPATH=$KILLMS_DIR/ddf-pipeline:$PYTHONPATH
export PYTHONPATH=$KILLMS_DIR/ddf-pipeline/utils:$PYTHONPATH
export PYTHONPATH=$KILLMS_DIR/ddf-pipeline/scripts:$PYTHONPATH
#export PYTHONPATH=/home/cyril.tasse/MCMCounts:$PYTHONPATH
#export PATH=/home/cyril.tasse/MCMCounts:$PATH
export PATH=/cep/lofar/ds9:$PATH

#export RCLONE_CONFIG_DIR=/home/cyril.tasse/MACARON

export LD_LIBRARY_PATH=$KILLMS_DIR/DDFacet/DDFacet/cbuild:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=/home/cyril.tasse/VE_nancep4/sfft/build:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$VE_FOLDER/sources/LOFARBeam/cbuild
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$VE_FOLDER/lib


export PATH=$KILLMS_DIR/killMS/killMS:$PATH
export PATH=$KILLMS_DIR/DDFacet/DDFacet:$PATH
export PATH=$KILLMS_DIR/drawMS:$PATH
export PATH=$KILLMS_DIR/AnalyseDynSpecMS:$PATH
export PATH=$KILLMS_DIR/DDFacet/SkyModel:$PATH
export PATH=$KILLMS_DIR/DynSpecMS:$PATH
export PATH=$KILLMS_DIR/FindCluster:$PATH
export PATH=$KILLMS_DIR/DIDAP:$PATH
export PATH=$KILLMS_DIR/ddf-pipeline/scripts:$PATH

#export FINDCLUSTER_DATA_DIR=/data/cyril.tasse/DataDeepFields

export PS1="\[\e[1;92m\][src@${KILLMS_DIR}]\[\e[m\] \u@\h:\w\$ "

echo "DDF exec: "$(which DDF.py)
echo "kMS exec: "$(which kMS.py)

#exec "$@"
