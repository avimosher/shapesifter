#!/bin/bash

platform='unknown'
unamestr=`uname`
if [[ "$unamestr" == 'Linux' ]];then
    platform='Linux'
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:Linux
elif [[ "$unamestr" == 'Darwin' ]];then
    platform='OSX'
    export DYLD_FALLBACK_LIBRARY_PATH=$DYLD_FALLBACK_LIBRARY_PATH:OSX
fi

args=("$*")
script=$(cat)
echo "$script" | ./expand_scene.sh | ${platform}/simulate
