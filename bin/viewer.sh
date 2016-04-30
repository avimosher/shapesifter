#!/bin/bash

args=("$*")
script=$(cat)
echo "$script" | ./expand_scene.sh | ./viewer-js/node_modules/electron-prebuilt/cli.js ./viewer-js/main.js
