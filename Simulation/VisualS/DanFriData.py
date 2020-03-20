'''
Daniel's Friday Data
'''
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
'''
Data Order:  1.220 4.220 4.250 1.250 1.280 4.280
'''

D = pd.read_excel('allDansdata.xlsx')

#print(D[[24]])
an,af,ac = [],[],[]
bn,bf,bc = [],[],[]
cn,cf,cc = [],[],[]
dn,df,dc = [],[],[]
en,ef,ec = [],[],[]
fn,ff,fc = [],[],[]
for index,row in D[[24]].iterrows():
    an.append(row[0])
for index,row in D[[25]].iterrows():
    af.append(row[0])
for index,row in D[[26]].iterrows():
    ac.append(row[0])
for index,row in D[[27]].iterrows():
    bn.append(row[0])
for index,row in D[[28]].iterrows():
    bf.append(row[0])
for index,row in D[[29]].iterrows():
    bc.append(row[0])
for index,row in D[[30]].iterrows():
    cn.append(row[0])
for index,row in D[[31]].iterrows():
    cf.append(row[0])
for index,row in D[[32]].iterrows():
    cc.append(row[0])
for index,row in D[[33]].iterrows():
    dn.append(row[0])
for index,row in D[[34]].iterrows():
    df.append(row[0])
for index,row in D[[35]].iterrows():
    dc.append(row[0])
for index,row in D[[36]].iterrows():
    en.append(row[0])
for index,row in D[[37]].iterrows():
    ef.append(row[0])
for index,row in D[[38]].iterrows():
    ec.append(row[0])
for index,row in D[[39]].iterrows():
    fn.append(row[0])
for index,row in D[[40]].iterrows():
    ff.append(row[0])
for index,row in D[[41]].iterrows():
    fc.append(row[0])

an =np.asarray(an)
af =np.asarray(af)
ac =np.asarray(ac)
bn =np.asarray(bn)
bf =np.asarray(bf)
bc =np.asarray(bc)
cn =np.asarray(cn)
cf =np.asarray(cf)
cc =np.asarray(cc)
dn =np.asarray(dn)
df =np.asarray(df)
dc =np.asarray(dc)
en =np.asarray(en)
ef =np.asarray(ef)
ec =np.asarray(ec)
fn =np.asarray(fn)
ff =np.asarray(ff)
fc =np.asarray(fc)

ac = ac[~ np.isnan(ac)]
af = af[~ np.isnan(af)]
an = an[~ np.isnan(an)]
bc = bc[~ np.isnan(bc)]
bf = bf[~ np.isnan(bf)]
bn = bn[~ np.isnan(bn)]
cc = cc[~ np.isnan(cc)]
cf = cf[~ np.isnan(cf)]
cn = cn[~ np.isnan(cn)]
dc = dc[~ np.isnan(dc)]
df = df[~ np.isnan(df)]
dn = dn[~ np.isnan(dn)]
ec = ec[~ np.isnan(ec)]
ef = ef[~ np.isnan(ef)]
en = en[~ np.isnan(en)]
fc = fc[~ np.isnan(fc)]
ff = ff[~ np.isnan(ff)]
fn = fn[~ np.isnan(fn)]

D1 = [an,af]
D2 = [bn,bf]
D3 = [cn,cf]
D4 = [dn,df]
D5 = [en,ef]
D6 = [fn,ff]

'''# # # # # # #'''# # # # # # # ''' # # # # # # # ''' # # # # 

fig = plt.figure()
ax1 = plt.subplot2grid((2,1), (0,0), colspan=1)
ax2 = plt.subplot2grid((2,1), (1,0),sharex=ax1 )
'''
Data Order:  1.220 4.220 4.250 1.250 1.280 4.280
'''
ax1.axvline(115,c='pink',linewidth=10)
ax2.axvline(115,c='pink',linewidth=10)

ax1.scatter(D1[0],D1[1],c='r',s=250)
ax1.scatter(D4[0],10*D4[1],c='g',s=250)
ax1.scatter(D5[0],10*D5[1],c='b',s=250)

ax2.scatter(D2[0],D2[1],c='r',s=250)
ax2.scatter(D3[0],D3[1],c='g',s=250)
ax2.scatter(D6[0],10*D6[1],c='b',s=250)


plt.show()


