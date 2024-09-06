#!/usr/bin/python3

# Generates a PDF with plots of all summary curves from a reference
# case and a 'new' simulation.
# sys.argv[1] = Base name for reference case (ie no extension)
# sys.argv[2] = Base name for new simulation case (ie no extension)
# sys.argv[3] = 'Pretty' name of the test. Used for filename of output file.

from datetime import datetime, timedelta
from opm.io.ecl import ESmry
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.dates import DateFormatter
from matplotlib.backends.backend_pdf import PdfPages

ref_file = ESmry(sys.argv[1] + '.SMSPEC')
sim_file = ESmry(sys.argv[2] + '.SMSPEC')
test_name = sys.argv[3]
case_name = os.path.basename(sys.argv[1])

ref_time = ref_file.dates()
sim_time = sim_file.dates()

ref_time_in_secs = [(v - ref_time[0]).total_seconds() for v in ref_time]
sim_time_in_secs = [(v - sim_time[0]).total_seconds() for v in sim_time]

plt.rcParams['font.size'] = 8

# Attempt at sorting the graphs in descending eyeball norm order
deviation={}
for r in ref_file.keys():
    if r == 'TIME' or r == 'YEARS':
        continue
    try:
        ref = ref_file[r]
        sim = sim_file[r]
    except:
        continue

    if len(ref) == 0 and len(sim) == 0:
        continue

    if not (any(ref) or any(sim)):
        continue

    # Resample new simulation on reference points
    sim_resampled = np.interp(ref_time_in_secs, sim_time_in_secs, sim)

    dev = np.abs(sim_resampled - ref)
    deviation[r] = np.linalg.norm(dev, 2) / np.linalg.norm(ref, 2)

p = PdfPages(f'{test_name}.pdf')
for r in sorted(deviation, key = lambda x: deviation[x], reverse=True):
    try:
        ref = ref_file[r]
        sim = sim_file[r]
    except:
        continue

    fig, ax = plt.subplots()
    ax.plot(ref_time, ref, linestyle='dashed', linewidth=0.5, marker='o', markersize=1.0)
    ax.plot(sim_time, sim, linewidth=0.5, marker='x', markersize=1.0)
    ax.legend(['Reference', 'New simulation'])
    plt.title(r)
    u = ref_file.units(r)
    if u:
        plt.ylabel(u)
    myFmt = DateFormatter("%Y-%b")
    ax.xaxis.set_major_formatter(myFmt)
    ax.xaxis.set_major_locator(plt.MaxNLocator(20))
    plt.grid()
    fig.autofmt_xdate()
    fig.savefig(p, format='pdf')
    plt.close(fig)

p.close()
