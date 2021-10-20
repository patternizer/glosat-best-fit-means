import sys, math, numpy, pandas, ruptures


# prepare distance matrix from lat/lon in degrees
def prepare_dists(lats, lons):
  """
  Prepare distance matrix from vectors of lat/lon in degrees assuming
  spherical earth
  
  Parameters:
    lats (vector of float): latitudes
    lons (vector of float): latitudes
  
  Returns:
    (matrix of float): distance matrix in km
  """
  las = numpy.radians(lats)
  lns = numpy.radians(lons)
  dists = numpy.zeros([las.size,las.size])
  for i in range(lns.size):
    dists[i,:] = 6371.0*numpy.arccos( numpy.minimum( (numpy.sin(las[i])*numpy.sin(las) + numpy.cos(las[i])*numpy.cos(las)*numpy.cos(lns[i]-lns) ), 1.0 ) )
  return dists


# interpolate a list of locations using ordinary krigging
def interpolate( obs, cov, tor=0.0 ):
  """
  Interpolate values at locations described by the given covariance matrix
  using ordinary kriging with errors.
  
  Parameters:
    obs (vector of float): observations (some of which may be missing)
    cov (matrix of float): covariances or correlation matrix
    tor (float): error parameter, see https://hal.archives-ouvertes.fr/cel-02285439v2/document
  
  Returns:
    (vector of float): infilled values
  """
  # set up matrices
  data = obs.flatten()
  obsflag = numpy.logical_not( numpy.isnan(data) )
  a = cov[obsflag,:][:,obsflag] + (tor**2)*numpy.identity(numpy.count_nonzero(obsflag))
  b = cov[obsflag,:]
  a = numpy.vstack( [
        numpy.hstack( [ a, numpy.ones([a.shape[0],1]) ] ),
        numpy.hstack( [ numpy.ones([1,a.shape[1]]), numpy.zeros([1,1]) ] )
    ] )
  b = numpy.vstack( [ b,numpy.ones([1,b.shape[1]]) ] )
  # solve for basis function weigths
  try:
    x = numpy.linalg.solve( a, b )
  except:
    x = numpy.dot( numpy.linalg.pinv(a), b )
  # calculate temperatures and store
  xs = numpy.zeros_like( cov )
  xs[obsflag,:] = x[:-1,:]
  ys = obs.copy()
  ys[ numpy.isnan(ys) ] = 0.0
  return numpy.dot( ys, xs )


# interpolate a list of locations using ordinary krigging with holdout
def interpolatex( obs, cov, tor=0.0 ):
  """
  Interpolate values at locations described by the given covariance matrix
  using ordinary kriging with errors and approximate station holdout,
  i.e. each location only contributes to the estimates at other locations.
  
  Parameters:
    obs (vector of float): observations (some of which may be missing)
    cov (matrix of float): covariances or correlation matrix
    tor (float): error parameter, see https://hal.archives-ouvertes.fr/cel-02285439v2/document
  
  Returns:
    (vector of float),(vector of float): infilled values, variance estimates
  """
  # set up matrices
  data = obs.flatten()
  obsflag = numpy.logical_not( numpy.isnan(data) )
  a = cov[obsflag,:][:,obsflag] + (tor**2)*numpy.identity(numpy.count_nonzero(obsflag))
  b = cov[obsflag,:]
  a = numpy.vstack( [
        numpy.hstack( [ a, numpy.ones([a.shape[0],1]) ] ),
        numpy.hstack( [ numpy.ones([1,a.shape[1]]), numpy.zeros([1,1]) ] )
    ] )
  b = numpy.vstack( [ b,numpy.ones([1,b.shape[1]]) ] )
  # solve for basis function weigths .. numpy.dot( numpy.linalg.pinv(a), b )
  try:
    x = numpy.linalg.solve( a, b )
  except:
    x = numpy.dot( numpy.linalg.pinv(a), b )
  # calculate temperatures and store
  xs = numpy.zeros_like( cov )
  xs[obsflag,:] = x[:-1,:]
  numpy.fill_diagonal( xs, 0.0 )
  xs = xs / numpy.sum( xs, axis=0 )
  ys = obs.copy()
  ys[ numpy.isnan(ys) ] = 0.0
  ys = numpy.dot( ys, xs )
  vs = numpy.diagonal( cov ) - numpy.diagonal( numpy.dot( cov, xs ) )
  return ys, vs


# Change point detection using Killick 2012, Truong 2020
def changepoints( dorig, **opts ):
  """
  Change point detection using PELT (Killick 2012),
  implemented in the ruptures package (Truong 2020)
  
  Parameters:
    dorig (vector of float): original data with no seasonal cycle, e.g.
      difference between obs and local expectation. No missing values.
    opts (dictionary): additional options, including:
      "nbuf": minimum number of months between changepoints
    
  Returns:
    (list of float): list of indices of changepoints
  """  
  min_size = opts["nbuf"]
  penalty_value = 10
  algo = ruptures.KernelCPD(kernel="linear",min_size=min_size).fit(dorig)
  result = algo.predict(pen=penalty_value)
  return result[:-1]


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
  chg = []
  if diff.shape[0] > 2*opts["nbuf"]: chg = changepoints(diff,**opts)
  index = numpy.arange( dnorm.shape[0] )[mask]
  flags = numpy.full( dnorm.shape, 0, numpy.uint8 )
  for i in chg:
    flags[index[i]:] += 1
  return flags


# OLS fit for norms using seasonal cycle and fragment means
def fitnorms( y, flags, nfourier=0 ):
  """
  Fit station fragment norms station fragment using mean and seasonal cycle 
  
  Parameters:
    y (vector of float): Station temperature series
    flags (vector of uint8): Station fragment flags
    nfourier (int): Number of fourier orders used to fit annual cycle changes

  Returns:
    (vector of float): vector of norms by month
  """
  frags = numpy.unique(flags)
  nfrags = frags.shape[0]
  npar = 12 + (nfrags-1)*(1+2*nfourier)
  x = numpy.zeros([y.shape[0],npar])
  # annual cycle
  p = 0
  for m in range(12):
    x[m::12,p] = 1.0
    p += 1
  # constant shifts for fragments
  for f in range(nfrags-1):
    x[:,p] = 1.0*(flags==f)
    p += 1
  # cosine shifts for fragments
  for n in range(nfourier):
    for f in range(nfrags-1):
      x[:,p  ] = _cost[:,n]*(flags==f)
      x[:,p+1] = _sint[:,n]*(flags==f)
      p += 2
  xm = x[ ~numpy.isnan(y), : ]
  ym = y[ ~numpy.isnan(y) ]
  if ym.shape[0] < xm.shape[1]: return numpy.zeros_like(y)
  p = numpy.linalg.lstsq(xm,ym,rcond=None)[0]
  return numpy.dot(x,p)



# MAIN PROGRAM
def main():
  # command line arguments
  year0,year1 = 1780,2020
  base0,base1 = 1961,1990
  stationfilter = None
  nfourier = 0

  for arg in sys.argv[1:]:
    if arg.split("=")[0] == "-years":    # year calc
      year0,year1 = [int(x) for x in arg.split("=")[1].split(",")]
    if arg.split("=")[0] == "-bases":    # year calc
      base0,base1 = [int(x) for x in arg.split("=")[1].split(",")]
    if arg.split("=")[0] == "-filter":    # year calc
      stationfilter = arg.split("=")[1]
    if arg.split("=")[0] == "-fourier":    # year calc
      nfourier = int(arg.split("=")[1])

  # years and months
  bases  = list(range(base0,base1+1))
  years  = list(range(year0,year1+1))
  months = list(range(1,13))

  # read station data
  dsrc = pandas.read_pickle("../DATA/df_temp.pkl", compression='bz2')

  # filter if required
  if stationfilter:
    stnmask = dsrc["stationcode"].str.startswith(stationfilter)
  else:
    stnmask = numpy.full( [dsrc.shape[0]], True )
  stnmask = numpy.logical_and( stnmask, dsrc["stationlat"].notna() )
  stnmask = numpy.logical_and( stnmask, dsrc["stationlon"].notna() )
  dflt = dsrc[stnmask]

  # extract 1 row per station
  dsub = dflt.drop_duplicates( subset="stationcode" )
  print(dsub.columns)
  print(dsub.describe())

  # extract column info
  codes = dsub.loc[:,"stationcode"].values
  names = dsub.loc[:,"stationname"].values
  ctrys = dsub.loc[:,"stationcountry"].values
  lats  = dsub.loc[:,"stationlat"].values
  lons  = dsub.loc[:,"stationlon"].values
  dists = prepare_dists( lats, lons )

  # and full data
  dcodes = dflt.loc[:,"stationcode"].values
  dyears = dflt.loc[:,"year"].values
  dtemps = dflt.loc[:,["1","2","3","4","5","6","7","8","9","10","11","12"]].values

  # make data table
  nstn = len(codes)
  nmon = 12*len(years)
  data  = numpy.full( [nmon,nstn], numpy.nan )
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
        data[i+m,j] = dtemps[r,m]

  # set up covariance data
  cov = numpy.exp( -dists/900.0 )

  # tabulate fourier coefficients
  global _cost, _sint
  _cost = numpy.zeros([nmon,nfourier])
  _sint = numpy.zeros([nmon,nfourier])
  for n in range(nfourier):
    _cost[:,n] = numpy.cos( numpy.arange(nmon) * numpy.pi * (n+1) / 6.0 )
    _sint[:,n] = numpy.sin( numpy.arange(nmon) * numpy.pi * (n+1) / 6.0 )

  # baseline period mask
  basemask = numpy.logical_and( dates[:,0] >= base0, dates[:,0] < base1 )

  # set initial baselines
  print(dates[:,0],year0)
  norms = numpy.full( [nmon,nstn], numpy.nan )
  mnths = numpy.arange(nmon)%12
  print(mnths)
  # set baselines
  for s in range(nstn):
    for m in range(12):
      mmask = (mnths==m)
      norms[mmask,s] = numpy.nanmean(data[mmask,s])

  # set initial fragment flags
  flags = numpy.full( [nmon,nstn], 0, numpy.uint8 )

  # now iteratively fit baselines based on all available data
  ncycle = 20
  for cycle in range(ncycle):

    # calculate updated station anomalies using current baselines
    dnorm = data - norms

    # infill from updated station anomalies
    fill = numpy.full( [nmon,nstn], numpy.nan )
    var  = numpy.full( [nmon,nstn], numpy.nan )
    for j in range(nmon):
      fill[j,:],var[j,:] = interpolatex( dnorm[j,:], cov, 0.1 )
      if j%120 == 0: print( cycle, dates[j,0] )

    # apply empirical scale to variances
    d2 = (dnorm-fill)**2
    vmsk = numpy.logical_and( ~numpy.isnan(d2), var>0.0 )
    print( numpy.count_nonzero( ~numpy.isnan( var ) ), numpy.count_nonzero( vmsk ), var.shape )
    var *= numpy.mean( d2[vmsk] ) / numpy.mean( var[vmsk] )
    print( numpy.count_nonzero( ~numpy.isnan( var ) ), numpy.count_nonzero( vmsk ), var.shape )

    # diagnostics
    print( cycle, numpy.std(fill), numpy.nanstd(fill), flush=True )

    # update flags
    if cycle>10:
      for s in range(nstn):
        diff = dnorm[:,s]-fill[:,s]
        flags[:,s] = changemissing( diff, nbuf=60 )

    # normalize filled data on baseline period
    fillm = fill[basemask,:]
    for m in range(12):
      fill[m::12,:] -= numpy.mean( fillm[m::12,:], axis=0 )

    # update station baselines to match filled data
    print( numpy.unique(flags) )
    for s in range(nstn):
      norms[:,s] = fitnorms( data[:,s]-fill[:,s], flags[:,s], nfourier=nfourier )


  # DATA OUTPUT
  # We need to add missing years (rows) to the source dataframe, then add
  # extra columns for the additional information
  # Make output data frame
  nyr = len(years)
  ddst = pandas.DataFrame( columns=["year","stationcode"] )
  # now populate from source data frame
  ddst["year"] = numpy.tile( years, nstn )
  ddst["stationcode"] = numpy.repeat( codes, nyr )

  print(dflt.shape, ddst.shape)

  # make new columns for source dataframe
  # we need to do this in two steps because updating the dataframe is very slow
  print("Making data\n",flush=True)
  nvals = numpy.full( [ddst.shape[0],12], numpy.nan )
  evals = numpy.full( [ddst.shape[0],12], numpy.nan )
  svals = numpy.full( [ddst.shape[0],12], numpy.nan )
  for r in range(len(ddst)):
    year = ddst.at[ r, "year" ]
    code = ddst.at[ r, "stationcode" ]
    j = numpy.argmax( code == codes )
    if year0 <= year <= year1 and codes[j] == code:
      i = 12*(year-year0)
      nvals[r,:] = norms[i:i+12,j]
      evals[r,:] = fill[i:i+12,j]
      svals[r,:] = numpy.sqrt(var[i:i+12,j])
    if r%1000 == 0: print(" row ",r,flush=True)

  # add the new columns to the source data frame
  print("Copying data\n",flush=True)
  ncols = ["n1","n2","n3","n4","n5","n6","n7","n8","n9","n10","n11","n12"]
  ecols = ["e1","e2","e3","e4","e5","e6","e7","e8","e9","e10","e11","e12"]
  scols = ["s1","s2","s3","s4","s5","s6","s7","s8","s9","s10","s11","s12"]
  for m in range(12):
    ddst[ncols[m]] = nvals[:,m]
  for m in range(12):
    ddst[ecols[m]] = evals[:,m]
  for m in range(12):
    ddst[scols[m]] = svals[:,m]
  print( numpy.nanmean(ddst.loc[:,"n1"]), numpy.nanmean(ddst.loc[:,"e1"]), numpy.nanmean(ddst.loc[:,"s1"]) )

  # now join with the original data
  djoin = pandas.merge( dflt, ddst, how="outer", on=["year","stationcode"] )

  # save the data
  print("Saving data\n",flush=True)
  if stationfilter:
    djoin.to_pickle( "df_temp_homog_{:s}.pkl".format(stationfilter.zfill(2)), compression="bz2" )
  else:
    djoin.to_pickle( "df_temp_homog.pkl", compression="bz2" )



# Main program launcher
if __name__ == '__main__':
    main()

