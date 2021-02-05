import scanninglaw.times as times
from scanninglaw.source import Source

dr2_sl = times.dr2_sl(version='cog')
s5_hvs1 = Source('22h54m51.68s',
                 '-51d11m44.19s',
           photometry={'gaia_g':16.02},
           frame='icrs')
scans = dr2_sl(s5_hvs1, 
               return_fractions=True, fov=1)
print('times FoV1:', scans['times'][0])
print('fracs FoV1:', scans['fractions'][0])

>> times FoV1: [[1714.17326, ..., 2325.42742]]
>> fracs FoV1: [[0.94576, ..., 0.95739]]