import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
# matplotlib.use("pgf")
# matplotlib.rcParams.update({
#     "pgf.texsystem": "pdflatex",
#     'font.family': 'serif',
#     'text.usetex': True,
#     'pgf.rcfonts': False,
# })

import sys

workdir = sys.argv[1]

dfen=pd.read_csv(workdir+"/valen.dat", header=None)
dfen*=10.22/64
# dfen=dfen.iloc[lambda x: x.index%100==0]
dfet=pd.read_csv(workdir+"/valet.dat", header=None)
dfet*=10.22/64
# dfet=dfet.iloc[lambda x: x.index%100==0]
dfencum=pd.read_csv(workdir+"/valencum.dat", header=None)
dfencum*=10.22/64
# dfencum=dfencum.iloc[lambda x: x.index%100==0]
dfetcum=pd.read_csv(workdir+"/valetcum.dat", header=None)
dfetcum*=10.22/64
# dfetcum=dfetcum.iloc[lambda x: x.index%100==0]

fig, ax = plt.subplots()
ax.plot(dfen, linewidth=0.5, linestyle='None', color='blue', alpha=0.3, label='en', marker='.', markersize=1)
ax.plot(np.array(dfet.index)*10,dfet, linewidth=0.5, linestyle='None', color='red', alpha=0.3, label='et', marker='.', markersize=1)
ax.plot(dfencum, linewidth=0.8, color='blue', label='encum')
ax.plot(np.array(dfetcum.index)*10, dfetcum, linewidth=0.8, color='red', label='etcum')
ax.legend()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set(position='zero', linestyle='None')
ax.spines['left'].set(position='zero', linestyle='None')
# ax.annotate("", xy=(0, 0), xytext=(520000,0), arrowprops=dict(arrowstyle="<|-", linewidth=0.3))
# ax.annotate("", xy=(0,-13), xytext=(0,10), arrowprops=dict(arrowstyle="<|-", linewidth=0.3))
plt.show()
