# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np

from mantid.geometry import CrystalStructure, ReflectionGenerator, ReflectionConditionFilter
from mantid.kernel import V3D

import datetime

from mantid import config
config['Q.convention'] = 'Inelastic'

# --- sample parameters

cif_file = "/HFIR/HB1A/IPTS-33477/shared/hexaferrite.cif"

UB_str = "-0.0240786678434794 -0.00237374586730765 0.0228962739971058 0.196116873592074 0.0929225205042786 0.00285588062808851 0.00725628947373798 0.174510843798493 -0.00120924113281822"

lamda = 2.379
d_min = 0.7

k_vec = [0, 0, 0.0]
output_macro = '/HFIR/HB1A/IPTS-33477/shared/hexaferrite_scanlist_Xtal_Nuclear.macro'
title = "(Nuc_@h) @(k) @(l)) rocking, T = @(coldtip.position) K"

k_vec = [0, 0, 1.5]
output_macro = '/HFIR/HB1A/IPTS-33477/shared/hexaferrite_scanlist_Xtal_CM.macro'
title = "CM_@(temp.position)K_"

k_vec = [0, 0, 0.673]
output_macro = '/HFIR/HB1A/IPTS-33477/shared/hexaferrite_scanlist_Xtal_ICM_plus.macro'
title = "ICM_plus (@(h) @(k) @(l)) rocking, T = @(coldtip.position) K"


k_vec = [0, 0, -0.673]
output_macro = '/HFIR/HB1A/IPTS-33477/shared/hexaferrite_scanlist_Xtal_ICM_minus.macro'
title = "ICM_minus (@(h) @(k) @(l)) rocking, T = @(coldtip.position) K"
 

# k_vec = [0, 0, 0.673]
# output_macro = '/HFIR/HB1A/IPTS-33477/shared/hexaferrite_scanlist_Xtal_ICM_extra.macro'
# title = "ICM_extra (@(h) @(k) @(l)) rocking, T = @(coldtip.position) K"

# --- 

srange = 1.2
step = 0.12

tmax = 4 # seconds
tmin = 1

# --- goniometer settigs

theta2_min = -70 # deg
theta2_max = -2 # deg

chi_min = -50 # deg
chi_max = 60 # deg

al_width = 1 # deg

theta2_rate = 40 / 45
chi_rate = 94 / 90
phi_rate = 93 / 360
pt_overhead = 40 / 21
pt_overhead = 1.5

### ---

def angular_distance(setting1, setting2):
    return np.sqrt(sum((a - b) ** 2 for a, b in zip(setting1, setting2)))

def greedy_sort(settings):
    n = len(settings)
    visited = np.zeros(n, dtype=bool)
    path = []
    current_index = 0
    path.append(current_index)
    visited[current_index] = True

    distance_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            distance_matrix[i, j] = angular_distance(settings[i], settings[j])

    for _ in range(n - 1):
        next_index = np.argmin([distance_matrix[current_index, j] if not visited[j] else np.inf for j in range(n)])
        visited[next_index] = True
        path.append(next_index)
        current_index = next_index
    return [i for i in path]

CreateSingleValuedWorkspace(OutputWorkspace='sample')

LoadCIF(Workspace='sample', InputFile=cif_file)

generator = ReflectionGenerator(mtd['sample'].sample().getCrystalStructure())

hkls = generator.getHKLsUsingFilter(d_min, 1e6, ReflectionConditionFilter.Centering)

k = V3D(*k_vec)

hkls = [V3D(*(hkl+k)) for hkl in hkls] #+ [V3D(*(hkl-k)) for hkl in hkls]

hkls = list(set(hkls))

ds = np.array(generator.getDValues(hkls))
F2s = np.array(generator.getFsSquared(hkls))

hkls = np.array(hkls)

UB = np.array(UB_str.split()).astype(float).reshape(3, 3)

UB[0, :] *= -1
UB[2, :] *= -1

SetUB(Workspace='sample', UB=UB)

B = mtd['sample'].sample().getOrientedLattice().getB()

a = mtd['sample'].sample().getOrientedLattice().a()
b = mtd['sample'].sample().getOrientedLattice().b()
c = mtd['sample'].sample().getOrientedLattice().c()
alpha = mtd['sample'].sample().getOrientedLattice().alpha()
beta = mtd['sample'].sample().getOrientedLattice().beta()
gamma = mtd['sample'].sample().getOrientedLattice().gamma()


newphi = np.zeros_like(ds)
newchi = np.zeros_like(ds)

Qc = np.einsum('ij,kj->ik', B, hkls)

q = np.sqrt(np.sum(Qc**2, axis=0))

theta =  np.arcsin(lamda * q / 2)
theta2 = theta * 2 * 180 / np.pi

omega = theta2 / 2

hphi = np.einsum('ij,kj->ik', UB, hkls)

for i in range(len(ds)):
    newphi[i] = np.arctan(hphi[1, i] / hphi[0, i]) * 180 / np.pi
    if newphi[i] > 0:
        if hphi[0, i] < 0:
            newphi[i] -= 180
    else:
        if hphi[0, i] < 0:
            newphi[i] += 180

    newchi[i] = np.arctan2(hphi[2, i], np.sqrt(hphi[0, i]**2 + hphi[1, i]**2)) * 180 / np.pi


theta2 *= -1
omega *= -1

mask = (theta2 >= theta2_min) & (theta2 < theta2_max) \
     & (newchi > chi_min) & (newchi < chi_max) \
     & (np.abs(theta2 - 61) > al_width) & (np.abs(theta2 - 72) > al_width)

hkls = hkls[mask]
ds = ds[mask]
F2s = F2s[mask]
theta2 = theta2[mask]
omega = omega[mask]
newphi = newphi[mask]
newchi = newchi[mask]

if np.linalg.norm(k) > 0:
    F2s *= 0

# n = 20

# sort = np.argsort(F2s)[::-1][:n]
# hkls = hkls[sort][:n]
# ds = ds[sort][:n]
# F2s = F2s[sort][:n]
# theta2 = theta2[sort][:n]
# omega = omega[sort][:n]
# newphi = newphi[sort][:n]
# newchi = newchi[sort][:n]

sort = greedy_sort(np.column_stack([theta2 * theta2_rate, newchi * chi_rate,  newphi * phi_rate]))

F2_min = np.min(F2s)
F2_max = np.max(F2s)

F2_diff = F2_max - F2_min

print("a={:4f} b={:4f} c={:4f} alpha={:4f} beta={:4f} gamma={:4f}".format(a, b, c, alpha, beta, gamma))

totaltime = 0 

counts = 0

prev_phi = 0
prev_chi = 0
prev_theta2 = 0

with open(output_macro, 'w') as f:
    for i in sort:
        h, k, l = hkls[i]

        f.write('br {:.3f} {:.3f} {:.3f}\n'.format(h, k, l))

        f.write('scantitle "{}"\n'.format(title))

        scale = np.arcsinh(F2s[i]) - np.arcsinh(F2_min)
        scale /= np.arcsinh(F2_max) - np.arcsinh(F2_min) + 1e-12
        t = tmax - round(scale * (tmax - tmin))

        f.write('preset mcu {:2.0f}\n'.format(t))

        f.write('scanrel omega {:.2f} -{:.2f} {:.2f}\n\n'.format(srange, srange, step))

        n = (round(2 * srange / step) + 1) * srange

        br_time = np.abs([(theta2[i] - prev_theta2) * theta2_rate,
                          (newchi[i] - prev_chi) * chi_rate, 
                          (newphi[i] - prev_phi) * phi_rate])

        totaltime += t * n 
        totaltime += pt_overhead * n
        totaltime += np.max(br_time)

        prev_phi = newphi[i]
        prev_chi = newchi[i]
        prev_theta2 = theta2[i]

        counts += 1

print("total peaks = {}".format(counts))

print("total time = {}".format(datetime.timedelta(seconds=int(totaltime))))
