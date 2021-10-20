import sys, os, numpy, pandas


# MAIN PROGRAM


# read station data
dsrc = pandas.read_pickle("df_temp_homog.pkl", compression='bz2')
#pandas.options.display.max_columns = 50
print(dsrc.columns)
print(dsrc.describe(include="all"))


# extract 1 row per station
dsub = dsrc.drop_duplicates( subset="stationcode" )
print(dsub.columns)
print(dsub.describe(include="all"))

# extract column info
codes = dsub.loc[:,"stationcode"].values

# and full data
dcodes = dsrc.loc[:,"stationcode"].values
dyears = dsrc.loc[:,"year"].values
dtemps = dsrc.loc[:,["1","2","3","4","5","6","7","8","9","10","11","12"]].values
dnorms = dsrc.loc[:,["n1","n2","n3","n4","n5","n6","n7","n8","n9","n10","n11","n12"]].values

# write
with open( "df_temp_homog.dat", "w" ) as f:
  for i in range(dcodes.shape[0]):
    f.write( "{:11s}{:4d}TAVG".format(dcodes[i],dyears[i]) )
    for j in range(12):
      if not pandas.isna(dtemps[i,j]):
        f.write("{:5d}   ".format(int(round(100*(dtemps[i,j]-dnorms[i,j])))))
      else:
        f.write("-9999   ")
    f.write("\n")
