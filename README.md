# autoScale 

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

This archive contains `autoScale.py`, a standalone python tool for performing an
automated finite-size scaling (FSS) analysis. FSS is done in a model
independent manner, based on the minimization of an objective function
measuring the quality of the data collapse, as detailed in 

>  J. Houdayer, A. K. Hartmann, "Low-temperature behavior of two-dimensional Gaussian Ising spin glasses", Phys. Rev. B 70, 014418 (2004), DOI: [10.1103/PhysRevB.70.014418](https://doi.org/10.1103/PhysRevB.70.014418)

The archive further includes raw data, allowing a user to perfom an exemplary
scaling analysis for the  two-dimensional site percolatin problem. 


### Prerequisites

autoScale has been under development since 2007. It was implemented under
Python 2.3.4 (and later versions of Python 2) and it imports only basic Python
modules that come with any python installation. 

## Included materials

The repository contains: 

```
autoScale/
├── LICENSE.md
├── README.md
├── ORDER_PARAMETER
├── autoScale.py
├── autoScale_guide.pdf
├── main_fig_FSS.py
├── fig_FSS.png
├── inputFiles.dat
├── scaled_L512_256_128.out
├── error.out
├── scalingPlot.sh
└── scalingScript.sh
```

Subfolder `ORDER_PARAMETER/` contains the raw data for the order parameter,
i.e.\ the relative size of the largest cluster, for the two-dimensional
site-percolation problem on two-dimensional lattices.

`autoScale.py` is the stanalone python script that performs the FSS analysis.

`autoScale_guide.pdf` is a user guide.

`inputFiles.dat` is a file serving as an easy interface between `autoScale.py`
and the raw data.

`scalingPlot.sh` is a bash script that allows to display the scaled data using [gnuplot](http://gnuplot.sourceforge.net).

`scalingScript.sh` is a bash script that runs a sequence of FSS analyses on the provided data.

`main_fig_FSS.py` is a python script that generates the figure below (it requires the functionality of numpy and pythons matplotlib).

`scaled_L512_256_128.out` contains the results produced by `scalingScript.sh`, and `err.out` contains the error analysis.

The repository further contains
* `LICENSE`, a license file.
* `Readme.md`, this file.

For a more detailed description of functions, defined in the above modules,
their parameters and return values we refer to the example cases and
documentation provided within the code.

## Exemplary results for two-dimensional site percolation 

![alt text](https://github.com/omelchert/autoScale/blob/main/fig_FSS.png)

### Brief explanation of the above figure

Illustration of the data collapse technique for the order parameter of the site
percolation problem in 2d. (a) unscaled data for different system sizes close
to the critical point (each data point represents an average over 12800
individual occupied/empty configurations). (b) data collapse after the scaling
parameter optimization using `autoScale.py`.

### Reporting the results obtained using autoScale

Performing a FSS analysis for the enclodes order parameter data for the
two-dimensional site percolation problem within the region [-1.5:1.0] of the
scaled control parameter is done in the following way:

```
python autoScale.py -f inputFiles.dat -xc 0.5927 -a 0.75 -b 0.104 -xr -1.5 1.
dx = [-1.500000:1.000000]  xc = 0.592708  a = 0.747566  b = 0.104459  S = 1.076691 
```

Upon termination, `autoScale.py` lists the parameter for which the best data collapse
was achieved. In this case, the quality of the data collapse is approximately S=1.08.
The quality function S measures the mean-squared distance of the rescaled
quantities at finite system size to their master curve in units of the standard
error, similar to a chi-square-test.
Error bars for each of the parameter can be obtained by invoking the additional flag
`-getErrors`, resulting in 

```
python autoScale.py -f inputFiles.dat -xc 0.5927 -a 0.75 -b 0.104 -xr -1.5 1. -getError
# S+1 error analysis yields:
# Scaling analysis restricted to
  xr = [-1.500000 : 1.000000]
# <scalePar>  <-Err>  <+Err>
  xc = 0.592708 0.000075 0.000075
   a = 0.747566 0.005424 0.005419
   b = 0.104459 0.000971 0.000964
```

For a given parameter, say, parameter xc, the other parameters (a and b) are
kept fixed, and the value of xc is first increased and then decreased until the
quality S increases its value to S+1. The larger of the two xc-values is then
quoted as the error on xc. This is referred to as a S+1 analysis.


## Further information

For further information about program options and file-formats look up
`autoScale.py` directly, or see the users guide enclosed as
`autoScale_guide.pdf`, also available under 

> O. Melchert, "autoScale.py - A program for automatic finite-size scaling analyses: A user's guide,"  [arXiv:0910.5403](https://arxiv.org/abs/0910.5403) (2009).

## Availability of the software

autoScale is provided as a (system-)local software tool. There is no need to install it once you
got a local
[clone](https://help.github.com/en/github/creating-cloning-and-archiving-repositories/cloning-a-repository)
of the repository, e.g. via

``$ git clone https://github.com/omelchert/autoScale.git``

## License 

This project is licensed under the MIT License - see the
[LICENSE.md](LICENSE.md) file for details.

