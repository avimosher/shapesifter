#!/bin/bash

if [ ! -d External_Libraries ]
then
    mkdir External_Libraries
fi
cd External_Libraries

if [ ! -d cereal ]
then
    wget https://github.com/USCiLab/cereal/archive/v1.1.1.tar.gz
    tar -xvf v1.1.2.tar.gz
    mv cereal-1.1.2 cereal
fi

if [ ! -d eigen ]
then
    wget http://bitbucket.org/eigen/eigen/get/3.2.4.tar.gz
    tar -xvf 3.2.4.tar.gz
    mv eigen-eigen-10219c95fe65 eigen
fi

if [ ! -d jsoncpp ]
then
    git clone git@github.com:open-source-parsers/jsoncpp.git
fi

cd jsoncpp
git checkout tags/1.6.2
cmake .
make
python amalgamate.py
g++ -o dist/jsoncpp.os -c -fPIC -Idist dist/jsoncpp.cpp
g++ -o ../../bin/libjsoncpp.so -shared dist/jsoncpp.os
cd ..

if [ ! -d osg ]
then
    wget http://www.openscenegraph.org/downloads/stable_releases/OpenSceneGraph-3.0/source/OpenSceneGraph-3.0.0.zip
    unzip OpenSceneGraph-3.0.0.zip
    mv OpenSceneGraph-3.0.0 osg
fi

cd osg
./configure
make
cd ..

cd ..
scons
