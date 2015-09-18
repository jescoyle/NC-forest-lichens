#! /bin/bash

for file in *.bmp; do
    mv "$file" "${file//L_bw_Moments/}"
done

