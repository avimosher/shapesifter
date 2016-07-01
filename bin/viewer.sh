#!/bin/bash

args=("$*")
script=$(cat)
echo "$script" | ./expand_scene.sh | ./viewer/node_modules/electron-prebuilt/cli.js ./viewer/main.js
