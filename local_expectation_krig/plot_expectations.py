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

  # plot covering whole data
  f, (ax0, ax1) = plt.subplots(2, 1, gridspec_kw={'height_ratios':[2,1],'hspace':0},sharex=True,figsize=(6,4))
  l0, = ax0.plot(moving_average(datef,12),moving_average(datat[:,s]-datan[:,s],12),'k-',lw=0.8)
  l1, = ax0.plot(moving_average(datef,12),moving_average(datae[:,s]           ,12),'b-',lw=0.5)
  ax0.fill_between(moving_average(datef,12),moving_average(datae[:,s]-2*datas[:,s],12),moving_average(datae[:,s]+2*datas[:,s],12),color='b',alpha=0.2,lw=0)
  ax0.legend([l0,l1],["Observations","Local expectation"],loc="lower center",ncol=2)
  ax1.plot(datef,datat[:,s]-datan[:,s]-datae[:,s],'k-',lw=0.5)
  ax1.fill_between(datef,-2*datas[:,s],2*datas[:,s],color='b',alpha=0.2,lw=0)
  ax0.axhline(0,color='k',lw=0.5,alpha=0.2)
  ax1.axhline(0,color='k',lw=0.5,alpha=0.2)
  plt.subplots_adjust(left=0.08,bottom=0.08,right=0.96,top=0.96)
  plt.savefig( "graphs/station_{:s}_full.svg".format(codes[s]) )
  plt.close()

  # plot covering observations data
  f, (ax0, ax1) = plt.subplots(2, 1, gridspec_kw={'height_ratios':[2,1],'hspace':0},sharex=True,figsize=(6,4))
  l0, = ax0.plot(moving_average(datef,12),moving_average(datat[:,s]-datan[:,s],12),'k-',lw=0.8)
  l1, = ax0.plot(moving_average(datef,12),moving_average(datae[:,s]           ,12),'b-',lw=0.5)
  ax0.fill_between(moving_average(datef,12),moving_average(datae[:,s]-2*datas[:,s],12),moving_average(datae[:,s]+2*datas[:,s],12),color='b',alpha=0.2,lw=0)
  ax0.legend([l0,l1],["Observations","Local expectation"],loc="lower center",ncol=2)
  ax1.plot(datef,datat[:,s]-datan[:,s]-datae[:,s],'k-',lw=0.5)
  ax1.fill_between(datef,-2*datas[:,s],2*datas[:,s],color='b',alpha=0.2,lw=0)
  ax0.axhline(0,color='k',lw=0.5,alpha=0.2)
  ax1.axhline(0,color='k',lw=0.5,alpha=0.2)
  plt.xlim([date0,date1])
  plt.subplots_adjust(left=0.08,bottom=0.08,right=0.96,top=0.96)
  plt.savefig( "graphs/station_{:s}_obs.svg".format(codes[s]) )
  plt.close()
