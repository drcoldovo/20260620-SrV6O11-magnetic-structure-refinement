import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd
from pathlib import Path


# ========================== Paths ==========================
# get the project root based on this file’s location
try:
    ROOT = Path(__file__).resolve().parent.parent  # if file inside src/
except NameError:
    ROOT = Path.cwd()  # console fallback

DATA_DIR = ROOT / "data"
DATA_5K = DATA_DIR / "NPD_5K.dat"
DATA_60K = DATA_DIR / "NPD_60K.dat"
DATA_100K = DATA_DIR / "NPD_100K.dat"
DATA_300K = DATA_DIR / "NPD_300K.dat"


# ========================== Settings ==========================
font = {'size': 22}
matplotlib.rc('font', **font)
plt.rcParams.update({
    "text.usetex": True,})
plt.rc('text.latex', preamble=r'\usepackage[cm]{sfmath}'
                              r'\usepackage{xfrac}')


# ========================== Load Data ==========================
df_5K = pd.read_csv(DATA_5K, delimiter='\s+', names=["d", "I", "err"])
df_60K = pd.read_csv(DATA_60K, delimiter='\s+', names=["d", "I", "err"])
df_100K = pd.read_csv(DATA_100K, delimiter='\s+', names=["d", "I", "err"])
df_300K = pd.read_csv(DATA_300K, delimiter='\s+', names=["d", "I", "err"])


# ========================== Plot ==========================
fig = plt.figure(figsize=(20, 8))
ax = plt.axes([0.10, 0.2, 0.87, 0.6])

# ax.plot(df_5K["d"], df_5K["I"], label="5 K", lw=1, color='#1240AB')
# ax.plot(df_60K["d"], df_60K["I"] + 0.5, label="60 K", lw=1, color='#FFAA00')
# ax.plot(df_100K["d"], df_100K["I"] + 1, label="100 K", lw=1, color='#FF0000')

ax.plot(df_100K["d"], df_5K["I"] - df_100K["I"], label="5 K - 100 K", lw=1, color='#1240AB')
ax.plot(df_100K["d"], df_60K["I"] - df_100K["I"], label="60 K - 100 K", lw=1, color='#FFAA00')

# Axes formatting
ax.set_xlim(2.1, 2.2)
# ax.set_xticks([0, 1/3, 1/2], [r'$\Gamma$', r'K', r'M'])
# ax.set_ylim(-0.4, 0.4)
# ax.set_yticks([0, 0.5, 1])

ax.set_xlabel(r'$d$-spacing ($\mathrm{\AA}$)', fontsize=22)
ax.set_ylabel(r'Intensity (arb. u.)', fontsize=22)
# ax.set_title(title_str, fontsize=22)

ax.tick_params(direction='in', top=True, right=True)
ax.set_axisbelow(False)

plt.legend(frameon=False, loc=2, fontsize=18)
# plt.savefig(ROOT / 'results/fig_magnetization/dMdT_2D_2K_JHU.svg')
plt.show()
