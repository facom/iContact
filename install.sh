#!/bin/bash
DIR=$(pwd)
ARCH="32"
#ARCH="64"

########################################
#INSTALL GSL
########################################
echo "Installing GSL..."
if [ ! -e util/lib/libgsl.a ];then
    cd util
    tar zxf gsl.tar.gz
    cd gsl
    ./configure --prefix=$DIR/util && make && make install
    cd ..
    #rm -rf gsl
    echo "Done."
else
    echo "Already installed."
fi
cd $DIR

########################################
#INSTALL CSPICE
########################################
echo "Installing SPICE..."
if [ ! -e util/lib/cspice.a ];then
    cd util
    cp cspice$ARCH.tar.Z cspice.tar.Z
    uncompress cspice.tar.Z
    tar xf cspice.tar
    cd cspice
    ./makeall.csh
    cp include/* ../include
    cp lib/* ../lib
    cd ..
    rm -rf cspice.tar 
    #rm -rf cspice
    echo "Done."
else
    echo "Already installed."
fi
cd $DIR

