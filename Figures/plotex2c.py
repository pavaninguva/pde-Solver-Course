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


# Generate list of the relevant csv files: 
csv_list = glob.glob("ex2c*.csv")

# Setup plot
fig = plt.figure(figsize=(6, 6))

for csv in csv_list:
    df = pd.read_csv(csv)
    xaxis = df["x"]
    yaxis = df["$c$"]
    fulltime = float(re.search("ex2c_(.*?).csv", csv).group(1))
    time = round(fulltime, 1)
    label_string = "t = %s" %time
    plt.plot(xaxis, yaxis, label=label_string)


plt.xlabel("Distance")
plt.ylabel("Concentration")
plt.legend()

plt.savefig("ex2c.png", dpi=600)
