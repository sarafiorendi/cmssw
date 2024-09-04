import uproot
import pandas as pd
import numpy as np
import awkward as ak
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('f_digi', help='file path')
parser.add_argument('f_redigi', help='file path')
parser.add_argument('--noplot', default='HLT_*,L1_*,Flag_*', help='coma-separated list of names not to plot, default HLT_*,L1_*')
parser.add_argument('--plot-only-failing', dest='plot_fail_only', action='store_true')
args = parser.parse_args()

import fnmatch
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
from pdb import set_trace
eps = 10**-7

output_folder = 'validation_results'
if not os.path.isdir(output_folder):
  os.makedirs(output_folder)

# probably logging would be better
logfile = open(output_folder + '/validation_log.html', 'w')
logfile.write('''<html>
<body>
<pre>
''')

def to_html(txt):
  return txt.replace('<', '&lt;').replace('>', '&gt;')

var_nbins = {
 'eventNumber': 3000,  
 'detId': 30,  
 'clusterCenter': 200,  
 'clusterZ': 30,  
 'clusterGlobalY': 30,  
 'clusterSize': 30,  
 'clusterGlobalZ': 100,  
 'clusterGlobalX': 100,  
 'clusterR': 30,  
 'clusterLocalX': 300,
 'clusterLocalY': 300,  
}


color_code = {
  'green' : '\033[1;32m %s \033[0m',
  'red' : '\033[1;31m %s \033[0m',
  'orange' : '\033[1;35m %s \033[0m', # Could not find it, use purple
  'black' : '%s',
}
def log(txt, color = 'black'):
  print(color_code[color] % txt)
  logfile.write('<code style="color: %s">%s</code>\n' % (color, to_html(txt)))


class NanoFrame(object):
  def __init__(self, infile):
    self.uf = uproot.open(infile)
    self.tt = self.uf['ClusterTree']

  def __getitem__(self, key):
    return self.tt.tree[key].arrays()

  def keys(self):
    return list(self.tt.keys())

def byval_validation(v1, v2):
#   if not np.isfinite(v1).all() or not np.isfinite(v2).all():
#     v1 = v1[np.isfinite(v1)]
#     v2 = v2[np.isfinite(v2)]

  try:
    if v1.type == 'bool' : ## or np.issubdtype(v1.type, np.integer):
      return np.array_equal(v1, v2)
    else:
      ret_val = ((np.abs(v1 - v2) / (abs(v1) + eps)) < 0.001)
      return ak.all(ret_val)
  except ValueError:
    return False
    
    

noplot = args.noplot.split(',')

def stat_validation(v1, v2, name = '', val_valid = False, nbins = 20):
#   if not np.isfinite(v1).all() or not np.isfinite(v2).all():
#     log(name + '--> CONTAINS INFs/NANs!', 'orange')
#     v1 = v1[np.isfinite(v1)]
#     v2 = v2[np.isfinite(v2)]

  if v1.ndim == 0 and v2.ndim == 0:
    return True
  elif v1.ndim == 0 or v2.ndim == 0:
    return False

  v1 = v1.to_numpy()
  v2 = v2.to_numpy()

  M = max(v1.max(), v2.max())
  m = min(v1.min(), v2.min())
  m = m * 0.9 if m > 0 else m * 1.2
  M = M * 1.2 if M > 0 else M * 0.9
  if 'int' in str(v1.dtype):
    m = int(m) - 1
    M = int(M) + 1
    nbins = min(M - m, nbins*2)
  plt.clf()
  h1, _, _ = plt.hist(v1, range = (m,M), bins = nbins, label = 'digi', histtype = 'step')
  h2, _, _ = plt.hist(v2, range = (m,M), bins = nbins, label = 're-digi', histtype = 'step')
  plt.xlabel(name)
  ret_val = (h1 == h2).all()
  plt.legend(loc='best')
  skip = any(fnmatch.fnmatch(name, i) for i in noplot)
  skip = skip or (args.plot_fail_only and ret_val and val_valid)
  if not skip:
    plt.savefig(output_folder + '/%s.png' % name)
  plt.clf()
  return ret_val.all()

def plot_branch(vals, name = '', nbins = 20):
  if not np.isfinite(vals).all():
    log(name + '--> CONTAINS INFs/NANs!', 'orange')
    vals = vals[np.isfinite(vals)]

  if vals.shape[0] == 0:
    return 

  M = vals.max()
  m = vals.min()
  m = m * 0.9 if m > 0 else m * 1.2
  M = M * 1.2 if M > 0 else M * 0.9
  if 'int' in str(vals.dtype):
    m = int(m) - 1
    M = int(M) + 1
    nbins = min(M - m, nbins*2)
  
  plt.clf()
  plt.hist(vals, range = (m,M), bins = nbins, label = 're-digi', histtype = 'step')
  plt.legend(loc='best')
  plt.savefig(output_folder + '/%s.png' % name)
  plt.clf()
  return 


old = NanoFrame(args.f_digi)
new = NanoFrame(args.f_redigi)

#
# Size Checks
#
def writer(pct):
  if pct < 5: return ''
  else: return '%.1f%%' % pct

#
# Branch checks
#
old_k = set(old.keys())
new_k = set(new.keys())
intersection = old_k.intersection(new_k)

log('Branch diff:')
for branch in (new_k - old_k):
  v_new = new[branch]
  if hasattr(v_new, 'content'):
    v_new = v_new.content
  plot_branch(v_new, branch)
  log(' '.join(['+', branch]), 'green')

for branch in (old_k - new_k):
  log(' '.join(['-', branch]), 'red')

log('\n\n')


for branch in sorted(intersection):
  v_old = old[branch]
  v_new = new[branch]

  if hasattr(v_old, branch):
    v_old = getattr(v_old, branch)
    v_new = getattr(v_new, branch)

  val_valid  = byval_validation(v_old, v_new)
  stat_valid = stat_validation(v_old, v_new, branch, val_valid, nbins = var_nbins[branch] if branch in var_nbins.keys() else 20)

  if val_valid and stat_valid:
    log(' '.join([branch, '--> OK!']), 'green')
  elif stat_valid:
    log(' '.join([branch, '--> FAILS BY VALUE CHECK ONLY!']), 'orange')
  else:
    log(' '.join([branch, '--> FAILS ALL CHECKS!']), 'red')

logfile.write('''
</pre>
</body>
</html>
''')
