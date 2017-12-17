#!/usr/bin/env bash

ffmpeg -r 30 -s 1024x1024 -i data/%d.pgm -vcodec libx264 -crf 15 -pix_fmt yuv420p out.mp4
