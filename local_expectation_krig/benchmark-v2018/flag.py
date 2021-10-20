import sys, numpy

l1 = open( sys.argv[1] ).readlines()
l2 = open( sys.argv[2] ).readlines()
s  = numpy.full([len(l1),12],-1)
t1 = numpy.full([len(l1),12],numpy.nan)
t2 = numpy.full([len(l2),12],numpy.nan)
for y in range(len(l1)):
  for m in range(12):
    s[y,m]  = int(l1[y][:11])
    t1[y,m] = int( l1[y][19+8*m:24+8*m] )
    t2[y,m] = int( l2[y][19+8*m:24+8*m] )
    if t1[y,m] < -9000: t1[y,m] = numpy.nan
    if t2[y,m] < -9000: t2[y,m] = numpy.nan

d = (t2-t1).flatten()
d[12:] = d[12:]-d[:-12]
d[:12] = 0

s = s.flatten()
sf = numpy.zeros_like(s)
sf[12:] = (s[12:]==s[:-12])

f = numpy.logical_and( ~numpy.isnan(d), numpy.abs(d)>0.5 )
f = numpy.logical_and( f, sf )

for i in range(f.shape[0]-23):
  if f[i] > 0:
    f[i+1:i+23] = 0

f = numpy.reshape( f, t1.shape )
for y in range(len(l1)):
  s = l1[y]
  for m in range(12):
    if f[y,m]:
      s = s[:19+8*m]+"{:5d}".format(20000+int(s[19+8*m:24+8*m]))+s[24+8*m:]
  print(s,end="")
