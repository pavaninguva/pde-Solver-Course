import matplotlib.pyplot as plt
from matplotlib import rc
import pandas as pd

rc('font', family='serif')
rc('xtick', labelsize='medium')
rc('ytick', labelsize='medium')
rc("axes", titlesize="medium")
rc("legend", fontsize="medium")


fig = plt.figure(figsize=(6,6))

gibbs = pd.read_csv("./ex4gibbs.csv")
print (gibbs)

xvals = gibbs["Time"]
yvals = gibbs["Gibbs"]

plt.plot(xvals, yvals)

plt.xlabel(r"$\tilde{t}$")
plt.ylabel(r"$G_{System}$")

plt.savefig("ex4_gibbs.png", dpi=600)