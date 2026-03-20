#!/bin/bash


G=/code/
export PATH=${PATH}:${G}/niftyreg_install/bin/:${G}/c3d/bin/

export NIFTYREG_INSTALL=${G}/niftyreg_install

export NIFTYMIC_ITK_DIR=${G}/ITK_NiftyMIC-build

. /etc/fsl/fsl.sh
. /code/venv/bin/activate

echo $PATH

exec "$@"
