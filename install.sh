#!/bin/bash

mkdir External_Libraries
cd External_Libraries
wget https://github.com/USCiLab/cereal/archive/v1.1.1.tar.gz
tar -xvf v1.1.1.tar.gz
mv cereal-1.1.1 cereal
wget http://bitbucket.org/eigen/eigen/get/3.2.4.tar.gz
tar -xvf 3.2.4.tar.gz
mv eigen-eigen-10219c95fe65 eigen
git clone git@github.com:open-source-parsers/jsoncpp.git
cd jsoncpp
git checkout tags/1.6.2
cmake .
make
python amalgamate.py
cp dist/libjsoncpp.so ../../bin/
cd ..
wget http://www.openscenegraph.org/downloads/stable_releases/OpenSceneGraph-3.0/source/OpenSceneGraph-3.0.0.zip
unzip OpenSceneGraph-3.0.0.zip
mv OpenSceneGraph-3.0.0 osg
cd osg
./configure
make
cd ../..
scons
