import sys, pandas, numpy, pickle

months = range(1,13)

stations = {}
# read ghcn data
for line in open(sys.argv[1].split(".")[0]+".inv"):
  id = line[0:11]
  lati = float( line[12:20] )
  lngi = float( line[21:30] )
  stations[id] = {"lat":lati,"lon":lngi,"data":{}}
# read and store the temperature data
for line in open(sys.argv[1]):
  id = line[0:11]
  if id in stations:
    year = int( line[11:15] )
    stations[id]["data"][year] = {}
    for month in months:
      temp = int( line[11+8*month:16+8*month] )
      flag = line[16+8*month:19+8*month]
      if temp != -9999 and flag[0:2] == "  ":
        stations[id]["data"][year][month] = 0.01*temp


cols = ['year', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12',
       'stationcode', 'stationlat', 'stationlon', 'stationelevation',
       'stationname', 'stationcountry', 'stationfirstyear', 'stationlastyear',
       'stationsource', 'stationfirstreliable']

nrow = 0
for id in stations:
  nrow += len(stations[id]["data"])

df = pandas.DataFrame(index=range(nrow),columns=cols,dtype=float)
df["year"] = 0
df.astype({"year":"int"})
df.astype({"stationcode":"string","stationname":"string","stationcountry":"string"})
df["stationcode"] = ""
df["stationname"] = ""
df["stationcountry"] = ""
print(df.shape)
print(df.dtypes)

irow = 0
for id in sorted(stations):
  for year in sorted(stations[id]["data"]):
    df.at[irow,"year"] = year
    df.at[irow,"stationcode"] = id
    df.at[irow,"stationlat"] = stations[id]["lat"]
    df.at[irow,"stationlon"] = stations[id]["lon"]
    df.at[irow,"stationelevation"] = 1.0
    df.at[irow,"stationname"] = "station"
    df.at[irow,"stationcountry"] = "country"
    df.at[irow,"stationfirstyear"] = 1900
    df.at[irow,"stationlastyear"] = 1999
    df.at[irow,"stationsource"] = 0
    df.at[irow,"stationfirstreliable"] = 0
    for m in months:
      if m in stations[id]["data"][year]:
        df.at[irow,cols[m]] = stations[id]["data"][year][m]
    irow += 1

print(df.dtypes)
print()
print(df.isna().sum())
print()

pandas.set_option('display.expand_frame_repr', False)
print(df)

df.to_pickle( sys.argv[1].split(".")[0]+".pkl", compression="bz2" )
