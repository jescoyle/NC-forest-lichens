#! /bin/bash

for file in *.tif; do
    mv "$file" "${file//Yadkin2_T/Yadkin1_T}"
done
