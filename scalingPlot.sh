# FILE:  scalingPlot.sh
# TASK:  create scaling plot to inspect data collapse
# USAGE: bash scalingPlot.sh <xc> <a> <b>
P='./ORDER_PARAMETER' 	# path to input files
gnuplot -persist << EOF

set key samplen 1. left					# customize key
set xl "(x-xc) L^a"; set yl "y L^b"			# set x,y labels 
xc=$1; a=$2; b=$3					# set scaling parameters
set label 1 "xc=$1\na =$2\nb =$3" at graph 0.7,0.2 	# list scaling parameters
sx(x,L)=(x-xc)*L**a; sy(y,L)=y*L**b			# def scaling assumption

p "$P/orderParam_L16.dat"  u (sx(\$1,16)): (sy(\$2,16))  w lp t "L=16"\
, "$P/orderParam_L32.dat"  u (sx(\$1,32)): (sy(\$2,32))  w lp t "  32"\
, "$P/orderParam_L64.dat"  u (sx(\$1,64)): (sy(\$2,64))  w lp t "  64"\
, "$P/orderParam_L128.dat" u (sx(\$1,128)):(sy(\$2,128)) w lp t " 128"\
, "$P/orderParam_L256.dat" u (sx(\$1,256)):(sy(\$2,256)) w lp t " 256"\
, "$P/orderParam_L512.dat" u (sx(\$1,512)):(sy(\$2,512)) w lp t " 512"
EOF
