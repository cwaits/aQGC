if [[ `basename $PWD` != "aQGC" ]]; then echo "Execute from aQGC dir"; exit; fi

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${PWD}/lib
export prodBase=$PWD

if [[ $HOSTNAME == "login.snowmass21.io" ]]; then
    #important
    module load python/2.7.15
    module load py-numpy/1.15.2-py2.7
    module load py-six/1.11.0-py2.7
    . /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.22.08/x86_64-centos7-gcc48-opt/bin/thisroot.sh
    
    #convenient
    module load emacs
else
    which root >> /dev/null
    if [[ $? -ne 0 ]]; then
	echo "Please ensure ROOT is available and retry"
    fi
fi
