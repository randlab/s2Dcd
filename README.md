# README

This is the README file of the `s2Dcd`, a Python package that allows
to obtain 3D multiple-point statistics (MPS) simulations by using only
2D training images (TIs).  More detais about the methodology can be
found in the paper by *Comunian, Renard and Straubhaar*, DOI:
[10.1016/j.cageo.2011.07.009](http://dx.doi.org/10.1016/j.cageo.2011.07.009).


## Requirements

This version of `s2Dcd` heavily relies on the `geone` module. You can
find more info about it and installation instructions at [this
link](https://github.com/randlab/geone).

## Quick start installation

> **NOTE:** To run the `s2Dcd`, a MPS simulation engine like the
[DeeSse](http://www.randlab.org/research/deesse/) is needed! If you
have installed the module `geone` then that should be already
included.

Clone or download the `s2Dcd` package on a local directory, by using for example:
```
git clone git@github.com:randlab/s2Dcd.git
```
Then inside the dowloaded directory

```
  pip install .
```


## Quick start examples

### Animation
For an animation that illustrated how the `s2Dcd` approach works, check [this link](./logo/animation.md).

### Simple example

Have a look at [this
link](examples/01_Strebelle/s2Dcd_run-ex01.ipynb) for a commented
Jupyter notebook. You can also find the same file as Python script in the `examples/01_Strebelle` folder.







