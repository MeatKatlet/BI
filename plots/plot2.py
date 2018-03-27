import matplotlib
matplotlib.use('TkAgg')
import numpy as np
import matplotlib.pyplot as plt

x = np.arange(0, 5, 1)

x2 = np.arange(0, 5, 1)
#stdev = np.array([1, 19])

#stdev = np.full((1, 20), 7)


y = np.sin(x)*10
y2 = np.sin(x)*10

#plt.plot(x, y)


#plt.boxplot(data)
#medianprops = dict(linewidth=0, color='red')
#cap=dict(linewidth=5)
#whiskerprops=dict(linewidth=5)
#plt.boxplot([[1,10],[2,5]], widths=0, positions=[5,15], showbox=True,showfliers=False,showmeans=False, medianprops=medianprops)

#plt.show()


plt.errorbar(x, y, yerr=[7,3,5,7,1], color='red', ls='--', marker='o', capsize=5, capthick=1, ecolor='black')
plt.errorbar(x2, y2, yerr=[3,3,3,3,3], color='red', ls='--', marker='o', capsize=5, capthick=1, ecolor='black')

plt.show()

