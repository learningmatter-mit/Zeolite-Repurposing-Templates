# Code for: "Repurposing Templates for Zeolite Synthesis from Simulations and Data Mining"

This repository has the code and data to reproduce all plots from the following paper:

D. Schwalbe-Koda et al. Repurposing Templates for Zeolite Synthesis from Simulations and Data Mining. Under Review (2022).

The code used to reproduce the plots is available at the [code](code/) folder. The data is available at the Materials Data Facility under the DOI [10.18126/c5z9-zej7](https://doi.org/10.18126/c5z9-zej7), and is loaded automatically wtihin the Jupyter Notebooks.

## Usage

The Jupyter notebooks require the following packages to run correctly:

```
pandas
numpy
matplotlib
```

These packages can be installed using the [environment.yml](environment.yml) file in this repo and `conda`:

```
conda env create -f environment.yml
```

## Related links

 - [VOID](https://github.com/learningmatter-mit/VOID): code to automate docking of molecules in zeolites.
 - [GULPy](https://github.com/learningmatter-mit/gulpy): Python interface to the GULP simulation software.
 - [Phase competition data](https://github.com/learningmatter-mit/Zeolite-Phase-Competition): repository containing more details on the binding energy data

## Citing

The bibtex citation for this paper is the following:

```
@article{schwalbe2022repurposing,
  title={Repurposing Templates for Zeolite Synthesis from Simulations and Data Mining},
  author={Schwalbe-Koda, Daniel and Santiago-Reyes, Omar A. and Corma, Avelino and Rom{\'a}n-Leshkov, Yuriy and Moliner, Manuel and G{\'o}mez-Bombarelli, Rafael},
  journal={Under Review},
  volume={},
  pages={},
  year={2022},
}
```
