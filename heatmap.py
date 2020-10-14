import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


np.random.seed(0)
uniform_data = np.random.rand(10, 12)
print(uniform_data)

data = np.loadtxt("hp_steps.csv", delimiter=",", skiprows=1)

plt.figure()
sns.heatmap(data)
plt.savefig('heatmap.png')
#print(data)
#plt.imshow(data, aspect="auto", interpolation = "none")
#plt.colorbar()
#plt.show()