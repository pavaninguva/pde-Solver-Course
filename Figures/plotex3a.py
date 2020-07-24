import glob as glob
import re
import matplotlib.pyplot as plt
from matplotlib import rc
import pandas as pd

rc('font', family='serif')
rc('xtick', labelsize='medium')
rc('ytick', labelsize='medium')
rc("axes", titlesize="medium")
rc("legend", fontsize="medium")

fig = plt.figure(figsize=(6,6))

initial = pd.read_csv("./ex3a_initial.csv")
initial_x = initial["x"]
initial_y = initial["$c$"]

initial_plot = plt.plot(initial_x, initial_y, label ="Initial Concentration Profile")

advection = pd.read_csv("./ex3a_advection.csv")
advection_x = advection["x"]
advection_y = advection["$c$"]

advection_plot = plt.plot(advection_x, advection_y, label="Convection-only")

diffusion = pd.read_csv("./ex3a_diffusion.csv")
diffusion_x = diffusion["x"]
diffusion_y = diffusion["$c$"]

diffusion_plot = plt.plot(diffusion_x, diffusion_y, label="Diffusion-only")

convdiffusion = pd.read_csv("./ex3a_condif.csv")
convdiffusion_x = convdiffusion["x"]
convdiffusion_y = convdiffusion["$c$"]

convdiffusion_plot = plt.plot(convdiffusion_x, convdiffusion_y, label="Convection and Diffusion")

plt.xlabel("Distance")
plt.ylabel("Concentration")
plt.ylim(0.0, 0.6)
plt.legend(loc="upper right")

plt.savefig("ex3a.png", dpi=600)