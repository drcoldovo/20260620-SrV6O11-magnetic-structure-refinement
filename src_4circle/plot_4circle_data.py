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

DATA_DIR = ROOT / "data/HFIR-HB-1A/exp1100/Datafiles"


# ========================== Settings ==========================
font = {'size': 22}
matplotlib.rc('font', **font)
plt.rcParams.update({
    "text.usetex": True,})
plt.rc('text.latex', preamble=r'\usepackage[cm]{sfmath}'
                              r'\usepackage{xfrac}')


# ========================== Load Data ==========================
filen_4K = DATA_DIR / "HB1A_exp1100_scan0205.dat"
# Extract column names from the "#   Pt. ..." line.
with open(filen_4K, "r", encoding="utf-8") as f:
    for line in f:
        if line.lstrip().startswith("#"):
            fields = line.lstrip()[1:].split()
            if fields and fields[0] == "Pt.":
                columns = fields
                break
    else:
        raise ValueError("Column header beginning with '# Pt.' not found.")
# Read the whitespace-separated numerical data.
df_4K = pd.read_csv(filen_4K, sep=r"\s+", comment="#", header=None, names=columns,)

filen_100K = DATA_DIR / "HB1A_exp1100_scan0373.dat"
# Extract column names from the "#   Pt. ..." line.
with open(filen_100K, "r", encoding="utf-8") as f:
    for line in f:
        if line.lstrip().startswith("#"):
            fields = line.lstrip()[1:].split()
            if fields and fields[0] == "Pt.":
                columns = fields
                break
    else:
        raise ValueError("Column header beginning with '# Pt.' not found.")
# Read the whitespace-separated numerical data.
df_100K = pd.read_csv(filen_100K, sep=r"\s+", comment="#", header=None, names=columns,)


# ========================== Plot ==========================
fig = plt.figure(figsize=(5, 5))
ax = plt.axes([0.22, 0.2, 0.67, 0.6])

ax.errorbar(x=df_4K["omega"], y=df_4K["detector"]/df_4K["mcu"]*60, yerr=np.sqrt(df_4K["detector"]/df_4K["mcu"]*60),
            capsize=2, fmt='o', color='k', mfc='k', label=r'4 K')
ax.errorbar(x=df_100K["omega"], y=df_100K["detector"]/df_100K["mcu"]*60, yerr=np.sqrt(df_100K["detector"]/df_100K["mcu"]*60),
            capsize=2, fmt='o', color='k', mfc='none', label=r'100 K')

# Axes formatting
# ax.set_xlim(2.1, 2.2)
# ax.set_xticks([0, 1/3, 1/2], [r'$\Gamma$', r'K', r'M'])
# ax.set_ylim(-0.4, 0.4)
# ax.set_yticks([0, 0.5, 1])

ax.set_xlabel(r'omega (deg)', fontsize=22)
ax.set_ylabel(r'Intensity (arb. u.)', fontsize=22)
ax.set_title(r"Rocking Scan at (111)", fontsize=22)

ax.tick_params(direction='in', top=True, right=True)
ax.set_axisbelow(False)

plt.legend(frameon=False, loc=2, fontsize=14)
# plt.savefig(ROOT / 'results/fig_magnetization/dMdT_2D_2K_JHU.svg')
plt.show()

