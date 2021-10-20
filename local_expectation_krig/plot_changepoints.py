import sys, os, numpy, pandas
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# moving average helper function
def moving_average(x, w):
  """
  Calculate moving average of a vector
  
  Parameters:
    x (vector of float): the vector to be smoothed
    w (int): the number of samples over which to average
    
  Returns:
    (vector of float): smoothed vector, which is shorter than the input vector
  """
  return numpy.convolve(x, numpy.ones(w), 'valid') / w


# Change point detection using Taylor (2000)
def changepoints( dorig, **opts ):
  """
  Change point detection using cumulative sum plots: Taylor (2000)
  
  Parameters:
    dorig (vector of float): original data with no seasonal cycle, e.g.
      difference between obs and local expectation. No missing values.
    opts (dictionary): additional options, including:
      "nsig": sigma threshold for accepting a changepoint
      "nbuf": minimum number of months between changepoints
      "errorfn": a function which gives the expected sigma of the cusum plot
    
  Returns:
    (list of float): list of indices of changepoints
  """  
  nsig, nbuf, errf = opts["nsig"], opts["nbuf"], opts["stderror"]
  if dorig.shape[0] < 3*nbuf: return []
  dnorm = dorig - numpy.mean(dorig)
  cusum = numpy.cumsum(dnorm)
  custd = errf(dorig)
  ratio = numpy.absolute(cusum)/numpy.maximum(custd,1.0e-30)
  i = numpy.argmax(ratio[nbuf:-nbuf])+nbuf
  if ratio[i] <= nsig: return []
  return changepoints(dorig[:i],**opts)+[i]+[x+i for x in changepoints(dorig[
i:],**opts)]


# calculate changepoints on data with breaks
  """
  Change point detection wrapper to allow missing data and return a
  vector of flags
  
  Parameters:
    dorig (vector of float): original data with no seasonal cycle, e.g.
      difference between obs and local expectation. No missing values.
    opts (dictionary): additional options for changepoints function
    
  Returns:
    (vector of unit8): vector of station fragment flags
  """  
def changemissing(dnorm,**opts):
  mask = ~numpy.isnan(dnorm)
  diff = dnorm[mask]
  chg = changepoints(diff,**opts)
  index = numpy.arange( dnorm.shape[0] )[mask]
  flags = numpy.full( dnorm.shape, 0, numpy.uint8 )
  for i in chg:
    flags[index[i]:] += 1
  return flags


# function defining the standard uncertainty profile for the cusum plot
def errorfn(data):
  """
  Standard uncertainty function for the cumulative sum plot. This is
  determined globally based on the frequency and distribution of
  inhomogeneities.
  
  Parameters:
    data (vector of float): Initial data (only the length is used, values ignored).
    
  Returns:
    (vector of float): vector of standard uncertainties
  """  
  svalue = 0.425  # std of station shift
  smnths = 300    # mean months between station shifts
  n = len(data)
  rt = smnths/n   # weighting function for short segments
  w = (rt+0.25)/(rt+0.75)  # (maybe rt/(rt+0.5)
  x = numpy.linspace( 1.0/n, 1.0-1.0/n, n )
  return svalue*w*n*x*(1.0-x)



# MAIN PROGRAM

# command line arguments
year0,year1 = None,None
stationfilter = None

for arg in sys.argv[1:]:
  if arg.split("=")[0] == "-years":    # year calc
    year0,year1 = [int(x) for x in arg.split("=")[1].split(",")]
  if arg.split("=")[0] == "-filter":    # year calc
    stationfilter = arg.split("=")[1]

# read station data
dsrc = pandas.read_pickle("df_temp_expect.pkl", compression='bz2')
#pandas.options.display.max_columns = 50
print(dsrc.columns)
print(dsrc.describe(include="all"))

# filter if required
if stationfilter:
  stnmask = dsrc["stationcode"].str.startswith(stationfilter)
else:
  stnmask = numpy.full( [dsrc.shape[0]], True )
dflt = dsrc[stnmask]

# extract 1 row per station
dsub = dflt.drop_duplicates( subset="stationcode" )
print(dsub.columns)
print(dsub.describe(include="all"))

# extract column info
codes = dsub.loc[:,"stationcode"].values
names = dsub.loc[:,"stationname"].values
ctrys = dsub.loc[:,"stationcountry"].values

# and full data
dcodes = dflt.loc[:,"stationcode"].values
dyears = dflt.loc[:,"year"].values
dtemps = dflt.loc[:,["1","2","3","4","5","6","7","8","9","10","11","12"]].values
dnorms = dflt.loc[:,["n1","n2","n3","n4","n5","n6","n7","n8","n9","n10","n11","n12"]].values
dexpts = dflt.loc[:,["e1","e2","e3","e4","e5","e6","e7","e8","e9","e10","e11","e12"]].values
dstdvs = dflt.loc[:,["s1","s2","s3","s4","s5","s6","s7","s8","s9","s10","s11","s12"]].values

# years and months
if year0 == None or year1 == None:
  year0,year1 = numpy.min(dyears), numpy.max(dyears)
years  = list(range(year0,year1+1))
months = list(range(1,13))

# make data table
nstn = len(codes)
nmon = 12*len(years)
datat  = numpy.full( [nmon,nstn], numpy.nan )
datan  = numpy.full( [nmon,nstn], numpy.nan )
datae  = numpy.full( [nmon,nstn], numpy.nan )
datas  = numpy.full( [nmon,nstn], numpy.nan )
dates = numpy.full( [nmon,2], 0.0 )
for j in range(nmon):
  y = years[j//12]
  m = months[j%12]
  dates[j,0] = y
  dates[j,1] = m

# fill table by station and month
for r in range(dcodes.shape[0]):
  if year0 <= dyears[r] <= year1:
    i = 12*(dyears[r]-year0)
    j = numpy.argmax( dcodes[r] == codes )
    for m in range(12):
      datat[i+m,j] = dtemps[r,m]
      datan[i+m,j] = dnorms[r,m]
      datae[i+m,j] = dexpts[r,m]
      datas[i+m,j] = dstdvs[r,m]

# plot each station in turn
datef = dates[:,0] + dates[:,1]/12 - 1/24
for s in range(datat.shape[1]):
  # find start and end dates
  try:
    dateobs = datef[~numpy.isnan(datat[:,s])]
    date0, date1 = round(dateobs[0]-7.5,10), round(dateobs[-1]+7.5,10)
  except:
    date0, date1 = dates[0,0], dates[-1,0]

  # calculate differece series and cumulative sum plot
  x = datef
  diff = datat[:,s]-datan[:,s]-datae[:,s]
  mask = ~numpy.isnan(diff)
  dmask = diff[mask]
  print(s,codes[s],dmask.shape)
  if dmask.shape[0] > 12:
    dmask -= numpy.mean(dmask)
    cusum = numpy.cumsum(dmask)
    cu = numpy.full( diff.shape, numpy.nan )
    ce = numpy.full( diff.shape, numpy.nan )
    cu[mask] = cusum
    ce[mask] = errorfn(dmask)
    chgf = changemissing(diff,nsig=2,nbuf=60,stderror=errorfn)
    chg = numpy.nonzero(numpy.diff(chgf)>0.5)[0]
    plt.figure(figsize=(6,2.5))
    plt.fill_between( x, -2*ce, 2*ce, color='b', alpha=0.1 )
    plt.plot( x, cu )
    plt.axhline(0,ls=":",c='k')
    for i in chg: plt.axvline(x=x[i],ls="-",c="g")
    ymax = max( numpy.nanmax(numpy.abs(cu)), numpy.nanmax( 2*ce ) )
    plt.ylim(-1.05*ymax,1.05*ymax)
    plt.xlim([date0,date1])
    plt.subplots_adjust(left=0.08,bottom=0.10,right=0.96,top=0.98)
    plt.savefig( "graphs/cusum_{:s}_obs.svg".format(codes[s]) )
    plt.close()

    # SAVE: cumsum to CSV

    df = pandas.DataFrame({'cu':cu, 'ce':ce}, index=x)
    df.to_csv( "graphs/cusum_{:s}_obs.csv".format(codes[s]) )
