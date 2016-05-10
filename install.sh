#!/bin/bash

if [ ! -d External_Libraries ]
then
    mkdir External_Libraries
fi
cd External_Libraries

if [ ! -d cereal ]
then
    wget https://github.com/USCiLab/cereal/archive/v1.1.2.tar.gz
    tar -xvf v1.1.2.tar.gz
    mv cereal-1.1.2 cereal
    patch cereal/include/types/polymorphic.hpp patches/cereal/polymorphic.patch
    patch cereal/include/details/polymorphic_impl.hpp patches/cereal/polymorphic_impl.patch
fi

if [ ! -d eigen ]
then
    hg clone https://avimosher@bitbucket.org/eigen/eigen
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

if [ ! -d catch ]
then
    mkdir catch
    cd catch
    wget https://raw.github.com/philsquared/Catch/master/single_include/catch.hpp
    cd ..
fi

if [ ! -d ESBTL ]
then
    wget http://downloads.sourceforge.net/project/esbtl/ESBTL-1.0-beta01.tar.bz2
    tar -xvf ESBTL-1.0-beta01.tar.bz2
    mv ESBTL-1.0-beta01 ESBTL
    patch ESBTL/include/ESBTL/line_reader.h patches/ESBTL/line_reader.h.patch
fi

cd ..
scons
