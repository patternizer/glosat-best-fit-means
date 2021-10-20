![image](https://github.com/patternizer/glosat-homogenisation/blob/main/MODEL-1-monthly-x1r-SE1r-725092(boston_city_wso)-744920(bho).png)
![image](https://github.com/patternizer/glosat-homogenisation/blob/main/MODEL-1-fit-725092(boston_city_wso)-744920(bho).png)

# Best Fit Means Algorithm

Python codebase for the Best Fit Means homogenisation algorithm being developed with Professor Tim Osborn for land air surface temperature station data. Part of ongoing work for the [GloSAT](https://www.glosat.org) project: www.glosat.org 

## Contents

* `baseline-estimator-model-1.py` - python script for Model 1A (uncorrelated standard errors) and Model 1B (modeling out the correlation) to estimate the baseline normal from single neighbouring station timeseries within a lasso radius and with the constraint that each monthly normal has at least 15 years of values
* `baseline-estimator-model-2a.py` - python script for Model 2A using the filtered neighbouring station ensemble mean in the segment and reference baseline region
* `baseline-estimator-model-2b.py` - python script for Model 2B optimising on the sampling error of the mean of core stations
* `baseline-estimator-model-2c.py` - python script for Model 2C maximizing co-located neighbour overlaps and segment length
* `baseline-estimator-model-3.py` - python script for Model 3 to solve the system of linear equations using core neighbours (in progress)


## Instructions for use

The first step is to clone the latest glosat-homogenisation code and step into the installed Github directory: 

    $ git clone https://github.com/patternizer/glosat-homogenisation.git
    $ cd glosat-homogenisation

Then create a DATA/ directory and copy to it the pickled dataframe GloSAT archive file df_temp.pkl.

### Using Standard Python

The code is designed to run in an environment using Miniconda3-latest-Linux-x86_64.

    $ python baseline-estimator-model-1.py
    $ python baseline-estimator-model-2a.py
    $ python baseline-estimator-model-2b.py
    $ python baseline-estimator-model-2c.py

## License

The code is distributed under terms and conditions of the [Open Government License](http://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/).

## Contact information

* [Michael Taylor](michael.a.taylor@uea.ac.uk)


