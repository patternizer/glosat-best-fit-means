import sys, numpy, pandas
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


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


# MAIN PROGRAM

stationfilter = None
for arg in sys.argv[1:]:
  if arg.split("=")[0] == "-filter":    # year calc
    stationfilter = arg.split("=")[1]

# read station data
dsrc = pandas.read_pickle("df_temp_expect.pkl", compression='bz2')

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
year0 = dsub.loc[:,"stationfirstyear"].values
year1 = dsub.loc[:,"stationlastyear"].values
lats = dsub.loc[:,"stationlat"].values
lons = dsub.loc[:,"stationlon"].values
dists = prepare_dists( lats, lons )


# loop over stations
nnear = 20
for i in range( len(codes) ):
  # create web page
  hfile = "html/station.{:s}.html".format(codes[i])
  with open( hfile, "w" ) as file:
    # get the "nnear" nearest stations
    idist = numpy.argsort( dists[:,i] )[:nnear]
    # html station header
    file.write("<html>\n<head><title>{:s}</title></head>\n<body>\n<h1>{:s} [{:s}]</h1>\n".format(codes[i],names[i],ctrys[i]))
    file.write("<img src='../graphs/station_{:s}_full.svg'/>\n".format(codes[i]))
    file.write("<div style='clear:both'></div>\n")
    file.write("<img src='../graphs/station_{:s}_obs.svg'/>\n".format(codes[i]))
    file.write("<div style='clear:both'></div>\n")
    file.write("<img src='../graphs/cusum_{:s}_obs.svg'/>\n".format(codes[i]))
    file.write("<div style='clear:both'></div>\n")
    file.write("<img src='../graphs/homogenization_{:s}_obs.svg'/>\n".format(codes[i]))
    file.write("<div style='clear:both'></div>\n")
    file.write("<img src='../graphs/neighbours.{:s}.png' style='float:left'/>".format(codes[i]))
    # write table of neighbours
    file.write("<table><tr><th>Neighbour</th><th>Name</th><th>Country</th><th>Distance</th><th>Lon/Lat</th><th>Years</th></tr>\n")
    for j in idist:
      file.write("<tr><td><a href='station.{:s}.html'>{:s}</a></td><td>{:s}</td><td>{:s}</td><td>{:d}</td><td>{:6.1f},{:6.1f}</td><td>{:d}-{:d}</td></tr>\n".format(codes[j],codes[j],names[j],ctrys[j],int(dists[j,i]),lons[j],lats[j],int(year0[j]),int(year1[j])))
    file.write("</table>\n")
    # html station footer
    file.write("</body></html>\n")

    # find neighbourhood size
    clon, clat = lons[i], lats[i]
    dlon = numpy.max( numpy.abs( lons[idist] - clon ) )
    dlat = numpy.max( numpy.abs( lats[idist] - clat ) )
    print(codes[i],names[i],lons[i],lats[i],dlon,dlat)

    # plot neighbourhood
    plt.figure( figsize=[4,6] )

    # plot 1st panel shows where in the world we are
    proj = ccrs.PlateCarree()
    gs = matplotlib.gridspec.GridSpec(2, 1, height_ratios=[1,2])
    gs.update(hspace=0.02)

    ax = plt.subplot(gs[0],projection=proj)
    ax.set_global()
    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.OCEAN)
    ax.plot( lons[idist], lats[idist], 'ro', transform=proj )

    # plot 2nd panel shows station neighbourhood
    proj = ccrs.PlateCarree()
    ax = plt.subplot(gs[1],projection=proj)

    # do coordinate transformation to get square plot (except at poles)
    LL = proj.transform_point(clon-dlon, clat-dlat, ccrs.Geodetic())
    UR = proj.transform_point(clon+dlon, clat+dlat, ccrs.Geodetic())
    EW = UR[0] - LL[0]
    SN = UR[1] - LL[1]
    # get side of the square extent (in map units, usually meters)
    side = max(EW, SN)    # larger value is in effect
    mid_x, mid_y = LL[0]+EW/2.0, LL[1]+SN/2.0  # center location
    extent = [mid_x-side/2.0, mid_x+side/2.0, mid_y-side/2.0, mid_y+side/2.0]
    extent[2] = max(extent[2],-90)
    extent[3] = min(extent[3], 90)

    # this sets square extent
    try:
      ax.set_extent(extent, crs=proj)
    except:
      ax.set_extent([-180,180,-90,90], crs=proj)
    
    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.OCEAN)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAKES, alpha=0.5)
    ax.add_feature(cfeature.RIVERS)
    ax.plot( lons[idist], lats[idist], 'ro' )
    for k in idist:
      ax.annotate( codes[k], [lons[k],lats[k]], fontsize=6, alpha=0.5 )
    
    plt.subplots_adjust(left=0.01,right=0.99,top=0.99,bottom=0.01)
    plt.savefig("graphs/neighbours.{:s}.png".format(codes[i]))
    plt.close()
