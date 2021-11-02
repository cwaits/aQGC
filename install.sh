if [[ `basename $PWD` != "aQGC" ]]; then echo "Execute from aQGC dir"; exit; fi

baseDir=$PWD
source setup.sh

#MadGraph
wget https://launchpad.net/mg5amcnlo/3.0/3.1.x/+download/MG5_aMC_v3.1.1.tar.gz
tar -xzvf MG5_aMC_v3.1.1.tar.gz
if [[ $? -ne 0 ]]; then
    echo "ERROR getting madgraph"
    exit
fi
rm MG5_aMC_v3.1.1.tar.gz

#Pythia
wget https://pythia.org/download/pythia83/pythia8306.tgz
tar -xzvf pythia8306.tgz
rm pythia8306.tgz
cd pythia8306
./configure --prefix=$baseDir
make -j
make install
if [[ $? -ne 0 ]]; then
    echo "ERROR compiling pythia"
    exit
fi
cd ..

#Delphes
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${baseDir}/lib    
export PYTHIA8=$baseDir
git clone https://github.com/delphes/delphes.git
cd delphes
make HAS_PYTHIA8=true -j -I${baseDir}/include/Pythia8
if [[ $? -ne 0 ]]; then
    echo "ERROR compiling delphes"
    exit
fi

#echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${baseDir}/lib" >> ~/.bash_profile
