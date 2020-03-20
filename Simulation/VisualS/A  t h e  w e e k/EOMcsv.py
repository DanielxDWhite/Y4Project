import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

fig = plt.figure()
ax1 = plt.subplot2grid((6,1), (0,0), colspan=1)
ax2 = plt.subplot2grid((6,1), (1,0),sharex=ax1 )
ax3 = plt.subplot2grid((6,1), (2,0),sharex=ax1 )
ax4 = plt.subplot2grid((6,1), (3,0),sharex=ax1 )
ax5 = plt.subplot2grid((6,1), (4,0),sharex=ax1 )
ax6 = plt.subplot2grid((6,1), (5,0),sharex=ax1 )

df1 = pd.read_excel('EOM_1.xlsx', skiprows=1)
H1,B1 = [],[]
for index,row in df1.iterrows():
    B1.append(row[0])
    H1.append(row[1])
ax1.fill_between(B1,H1)

df2 = pd.read_excel('EOM_2.xlsx', skiprows=1)
H2,B2 = [],[]
for index,row in df2.iterrows():
    B2.append(row[0])
    H2.append(row[1])
ax2.fill_between(B2,H2)

df3 = pd.read_excel('EOM_3.xlsx', skiprows=1)
H3,B3 = [],[]
for index,row in df3.iterrows():
    B3.append(row[0])
    H3.append(row[1])
ax3.fill_between(B3,H3)

df4 = pd.read_excel('EOM_4.xlsx', skiprows=1)
H4,B4 = [],[]
for index,row in df4.iterrows():
    B4.append(row[0])
    H4.append(row[1])
ax4.fill_between(B4,H4)

df5 = pd.read_excel('EOM_5.xlsx', skiprows=1)
H5,B5 = [],[]
for index,row in df5.iterrows():
    B5.append(row[0])
    H5.append(row[1])
ax5.fill_between(np.multiply(B5,2),H5)

df6 = pd.read_excel('EOM_6.xlsx', skiprows=1)
H6,B6 = [],[]
for index,row in df6.iterrows():
    B6.append(row[0])
    H6.append(row[1])
ax6.fill_between(np.multiply(B6,2),H6)


plt.title('EOM for csv',size=19)
plt.xlabel('Freq',size=17)
plt.ylabel('Intensity',size=17)

plt.show()