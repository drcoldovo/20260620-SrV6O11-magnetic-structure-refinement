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
list_x = np.loadtxt(DATA_x, delimiter=',')
list_y = np.loadtxt(DATA_y, delimiter=',')
list_c = np.loadtxt(DATA_c, delimiter=',')


# ========================== Plot ==========================
fig = plt.figure(figsize=(6, 5))
ax = plt.axes([0.20, 0.2, 0.67, 0.6])


# Axes formatting
# ax.set_xlim(x_min, x_max)
ax.set_xticks([0, 1/3, 1/2], [r'$\Gamma$', r'K', r'M'])
# ax.set_ylim(y_min, y_max)
# ax.set_yticks([0, 0.5, 1])

# ax.set_xlabel(r'($hh$)', fontsize=22)
# ax.set_ylabel(r'$\hbar\omega$ (meV)', fontsize=22)
# ax.set_title(title_str, fontsize=22)

ax.tick_params(direction='in', top=True, right=True)
ax.set_axisbelow(False)

# plt.legend(frameon=False, loc=2, fontsize=18)
# plt.savefig(ROOT / 'results/fig_magnetization/dMdT_2D_2K_JHU.svg')
plt.show()
