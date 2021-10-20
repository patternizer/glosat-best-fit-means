# glosat-homogenisation

Python codebase for GloSAT algorithms to homogenise land surface air temperature observations. Part of ongoing work for the [GloSAT](https://www.glosat.org) project: www.glosat.org 

## Contents

Two complementary approaches for homogenising instrumental LSAT timeseries are under development:

* `best_fit_means/` - python codeabse for the Best Fit Means algorithm ( Tim Osborn )
* `local_expectation_krig/` - python codebase for the Local Expectation Kriging algorithm ( Kevin Cowtan )

The format of the NetCDF output files is being designed to align with the input file format required by the EUSTACE processing system ( Colin Morice )

## Instructions for use

The first step is to clone the latest glosat-homogenisation code and step into the installed Github directory: 

    $ git clone https://github.com/patternizer/glosat-homogenisation.git
    $ cd glosat-homogenisation

Then create a DATA/ directory and copy to it the latest version of the pickled pandas dataframe GloSAT archive file: df_temp.pkl.

### Using Standard Python

The code is designed to run in an environment using Miniconda3-latest-Linux-x86_64.

## License

The code is distributed under terms and conditions of the [Open Government License](http://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/).

## Contact information

* [Michael Taylor](michael.a.taylor@uea.ac.uk)

