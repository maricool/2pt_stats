# This script reads in npy binary format from compute_log_weight.py
# and saves it at the precision and ascii format expected by the
# COSEBIs C code. --Ami Choi
import numpy as np
import mpmath as mp

# Change top two lines if necessary
nmode = 20
tmin, tmax = 0.5, 300.0

norms = 'tmin%0.1f_tmax%0.1f_Nmax%d_normalisations_tlog.npy'%(tmin,tmax,nmode)
roots = 'tmin%0.1f_tmax%0.1f_Nmax%d_roots.npy'%(tmin,tmax,nmode)
out_norms = 'Normalization_%0.2f-%0.2f.table'%(tmin,tmax)
out_roots = 'Root_%0.2f-%0.2f.table'%(tmin,tmax)

############

np.set_printoptions(precision=50)
mp.pretty = True

data = np.load(roots,allow_pickle=True)

# Output roots file
out_roots = open(out_roots,'w')
for idx in range(nmode):
    out_roots.write(
        '%d\t'%(idx+1)+'%0.50f\t'*len(data[idx])%tuple(data[idx])+'\n')

data = np.load(norms,allow_pickle=True)
#print data[0]
print("%0.50f"%data[0])

# Output text file
out_norms = open(out_norms,'w')
for idx in range(nmode):
    out_norms.write('%d\t%0.50f\n'%(idx+1,data[idx]))

### End of script ###
