#!/usr/bin/python
import TTCF
import FortranDriver 
import numpy as np
import matplotlib.pyplot as plt

md = FortranDriver.MDInterface()
md.npart = 500
#md.eps = 0.213
## Solid
#md.temp = 0.2
#md.density=1.5
## Liquid
#md.temp = 1.0
#md.density=0.6
## Triple-point
md.temp = 0.8
md.density = 0.8
## Gas
#md.temp = 4.0
#md.density=0.3

# Equilibriate
t_eq = 500
md.run(t_eq)
#md.fieldstrength = 3.0
#md.toggle_nemd()

# Get the RDF
rdf_dict = md.compute_rdf()

fig, ax = plt.subplots()

ax.set_xlabel("r (LJ units)")
ax.set_ylabel("g(r)")
ax.set_title(f"Radial distribution for LJ near triple-point")

ax.plot(rdf_dict['r'], rdf_dict['rdf'])
plt.hlines(1, 0, 2.5, linestyles={'dashed'})
plt.show()
