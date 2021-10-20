Local expectation calculation using hold-out kriging
====================================================

K Cowtan, 7/2021

This directory contains python code for calculation of a local
expectation for every station on a given period using hold-out
kriging. There are three programs:

* `calc_expectation.py`    Calculate the local expectations
* `plot_expectation.py`    Plot the local expectations
* `plot_neighbours.py`     Make maps and html pages for navigating plots
* `plot_changepoints.py`   Make cumulative sum plots for deviation from expected
* `calc_homogenization.py` Proof-of-concept homogenization program

Usage:
------

`python3 calc_expectation.py`

`python3 calc_expectation.py -filter=01`

`python3 calc_expectation.py -years=1750,2021`


calc_expectations.py
--------------------

Read data from `../DATA/df_temp.pkl`, process and write a file in an
extention of the same format to `./df_temp_expect.pkl`
36 extra columns are added:
* n1-n12 contain monthly normals - subtract these from the observations
* e1-e12 contain local expectations
* s1-s12 contain standard deviations for the local expectation

Command line arguments:
* -filter=prefix: Specify a subset of stations by the start of the code,
  e.g. 0 for Europe. 01 runs in minutes, 0 runs in ~ 1 hour. 
* -years=start,end: Specify a the range of years to process.
  Default 1780,2020
  

plot_expectations.py
--------------------

Plot the expectations from `./df_temp_expect.pkl` in the
`graphs/` directory.

plot_expectations.py
--------------------

Plot the station neigbourhoods from `./df_temp_expectation.pkl` in the
`graphs/` directory, and make html pages in the `html/` directory.

plot_changepoints.py
--------------------

Plot the station changepoints vs expectation from `./df_temp_expect.pkl`.

calc_homogenization.py
--------------------

The command line options are the same as calc_expectation.py
Output goes to `./df_temp_homog.pkl`
