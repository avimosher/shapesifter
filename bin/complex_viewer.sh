#!/bin/bash

args=("$@")
script=$(cat)
echo "$script" | ./expand_scene.py | ./viewer $args
