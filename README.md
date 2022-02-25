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

![alt text](https://github.com/omelchert/py-ecpn/blob/main/results/fig_01/fig_FSS.png)

### Brief explanation of the above figure

Illustration of the data collapse technique for the order parameter of the site
percolation problem in 2d. (a) unscaled data for different system sizes close
to the critical point (each data point represents an average over 12800
individual occupied/empty configurations). (b) data collapse after the scaling
parameter optimization using `autoScale.py`.


## Further information

For further information about program options and file-formats look up
`autoScale.py` directly, or see the users guide under 

> O. Melchert, "autoScale.py - A program for automatic finite-size scaling analyses: A user's guide," arXiv (2009), see [here](https://arxiv.org/abs/0910.5403)

## Availability of the software

autoScale is provided as a (system-)local software tool. There is no need to install it once you
got a local
[clone](https://help.github.com/en/github/creating-cloning-and-archiving-repositories/cloning-a-repository)
of the repository, e.g. via

``$ git clone https://github.com/omelchert/autoScale``

## License 

This project is licensed under the MIT License - see the
[LICENSE.md](LICENSE.md) file for details.

