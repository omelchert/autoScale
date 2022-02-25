# FILE:  scalingScript.sh
# TASK:  invoke scaling analysis program for different intervals xr
# USAGE: bash scalingScript.sh
for MAX in  1.5 1.25 1.0 0.75 ; 
do
 for MIN in -2.25 -2.0 -1.75 -1.5 -1.25 -1. -0.75 ;
 do
  python autoScale.py -f inputFiles.dat -o scaled_L512_256_128.out \
  	        -xc 0.5927 -a 0.75 -b 0.104 -showS -xr ${MIN} ${MAX}
 done
done
