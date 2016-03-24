#!/bin/bash

args=("$*")
script=$(cat)
echo "$script" | ./expand_scene.sh | ./viewer $args
