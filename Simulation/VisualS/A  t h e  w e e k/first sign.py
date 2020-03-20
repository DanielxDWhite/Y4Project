import matplotlib.pyplot as plt
import numpy as np

Off = [9.79,9.16,8.94,9.28,9.64,9.5,9.22,9.10,8.94,9.37,8.53,8.95,7.26,8.33,9.5,9.28,7.93,9.01,9.11]
TOff = np.arange(len(Off))


On = [10.1,10.4,10,10.5,10.5,10.2,9.26,9.49,10.1,10,9.52,9.23]
TOn = np.arange(len(On))

plt.scatter(TOff,Off,c='k',s=300)
plt.scatter(TOn,On,c='r',s=300)

Str = 'off = {}\non = {}\nmode p 200 = 6.9mW\nBase detuning = 60-70MHz'.format(round(np.std(Off),4),round(np.std(On),4))
plt.text(0,7,Str,fontsize=14)

plt.errorbar(TOff, Off,yerr=np.std(Off))
plt.errorbar(TOn, On,yerr=np.std(On))

plt.title('The First Sign of Slowing!',size=14)
plt.xlabel('"Time"',size=14)
plt.ylabel('# Atoms x10^5',size=14)
plt.plot(TOff,Off)
plt.plot(TOn,On)

plt.show()