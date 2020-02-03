#!/bin/bash
#test

# 4 steps

## 1. Download GOES data from Worldview
# write python script using Hauke's script

## 2. (if necessary) break gif in separated png files
# e.g. convert movies/nasa-worldview-2020-01-24T10_00Z-to-2020-01-24T16_00Z.gif movies/temp/frame.png

## 3. Add grid, circle and dropsonde point on each image
# call python script

## 4. Merge processed images into movie
# e.g. convert -loop 0 -delay 20 figures/forGIF/frame_0??.png movies/movie.gif

