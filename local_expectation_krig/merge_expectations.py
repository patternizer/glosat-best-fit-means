#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#------------------------------------------------------------------------------
# PROGRAM: merge_expectations.py
#------------------------------------------------------------------------------
# Verion 0.1
# 8 August, 2021
# Michael Taylor
# https://patternizer.github.io
# patternizer AT gmail DOT com
# michael DOT a DOT taylor AT uea DOT ac DOT uk
#------------------------------------------------------------------------------

import os
import glob
import pandas as pd

pkllist = glob.glob("*.pkl")
df_temp_expect = pd.concat([pd.read_pickle(pkllist[i], compression='bz2') for i in range(len(pkllist))])
df_temp_expect.to_pickle('df_temp_expect.pkl', compression='bz2')


