import sys, math, numpy


# MAIN PROGRAM
years = list(range(1900,2000))
datef = []
for y in years:
  for m in range(1,13):
    datef.append(y+m/12.0-1/24.0)

# read ghcn inventory
stations1 = {}
tnull = numpy.empty([len(datef)])
tnull.fill( numpy.nan )
for line in open(sys.argv[1]):
  stations1[line[0:11]] = { "lati":float(line[12:20]), "lngi":float(line[21:30]), "data":tnull.copy() }
# read and store the temperature data
for line in open(sys.argv[2]):
  id = line[0:11]
  if id in stations1:
    year = int( line[11:15] )
    if year in years:
      for m in range(12):
        temp = int( line[19+8*m:24+8*m] )
        flag = line[24+8*m:37+8*m]
        if temp != -9999 and flag[0:2] == "  ":
          stations1[id]["data"][12*(year-years[0])+m] = 0.01*temp

# read ghcn inventory
stations2 = {}
tnull = numpy.empty([len(datef)])
tnull.fill( numpy.nan )
for line in open(sys.argv[1]):
  stations2[line[0:11]] = { "lati":float(line[12:20]), "lngi":float(line[21:30]), "data":tnull.copy() }
# read and store the temperature data
for line in open(sys.argv[3]):
  id = line[0:11]
  if id in stations1:
    year = int( line[11:15] )
    if year in years:
      for m in range(12):
        temp = int( line[19+8*m:24+8*m] )
        flag = line[24+8*m:37+8*m]
        if temp != -9999 and flag[0:2] == "  ":
          stations2[id]["data"][12*(year-years[0])+m] = 0.01*temp

for key1 in stations1:
  t1 = numpy.nanmean(stations1[key1]["data"].reshape(-1,12),axis=1)
  t2 = numpy.nanmean(stations2[key1]["data"].reshape(-1,12),axis=1)
  print(key1,numpy.nanstd(t1-t2))

