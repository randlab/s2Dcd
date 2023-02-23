![s2Dcd logo](./logo/logo_thom.png "s2Dcd logo")


[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7006646.svg)](https://doi.org/10.5281/zenodo.7006646)



# README

This is the README file of the `s2Dcd`, a Python package that allows
to obtain 3D multiple-point statistics (MPS) simulations by using only
2D training images (TIs).  More details about the methodology can be
found in the paper by *Comunian, Renard and Straubhaar*, DOI:
[10.1016/j.cageo.2011.07.009](http://dx.doi.org/10.1016/j.cageo.2011.07.009).


## Requirements

This version of `s2Dcd` heavily relies on the `geone` module. You can
find more info about it and installation instructions at the link
[https://github.com/randlab/geone](https://github.com/randlab/geone).

## Quick start installation

> **NOTE:** To run the `s2Dcd`, a MPS simulation engine like the
[DeeSse](http://www.randlab.org/research/deesse/) is needed! If you
have installed the module `geone` then that should be already
included. Note however that the version of the DeeSse included in the `geone` repository has
**some limitations**. If you want to use the full functionalities of the package, please ask for a license.

Clone or download the `s2Dcd` package on a local directory, by using for example:
```
git clone git@github.com:randlab/s2Dcd.git
```
or
```
git clone git@bitbucket.org:alecomunian/s2dcd.git
```
(actually, the two repositories should contain exactly the same version of the code.)

Then inside the downloaded directory
```
pip install .
```

If the installation worked properly, then you should be able to perform an
`import s2Dcd` from a Python console without any error/warning.

## Examples

### Animation

For an animation that illustrates how the `s2Dcd` approach works, check [this link](./logo/animation.md).

### Simple example

Have a look at
[examples/01_Strebelle/s2Dcd_run-ex01.ipynb](examples/01_Strebelle/s2Dcd_run-ex01.ipynb)
for a commented Jupyter notebook. You can also find the same file as
Python script in
[examples/01_Strebelle/s2Dcd_run-ex01.py](examples/01_Strebelle/s2Dcd_run-ex01.py).

## More info

### Maintainers

At the moment, the code is maintained by
[A.Comunian](https://sites.unimi.it/alecomunian). Don't hesitate to
contact him if you have some suggestions of questions about the
`s2Dcd`.

### Source

The source file of the `s2Dcd` package is available both on *GitHub* and *Bitbucket* at the following links:

- [https://github.com/randlab/s2Dcd](https://github.com/randlab/s2Dcd)
- [https://bitbucket.org/alecomunian/s2dcd](https://bitbucket.org/alecomunian/s2dcd)








