'''   James's Data & method
'''

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
file1=('MarchResults Sheet1.csv')
data = pd.read_csv(file1, skiprows=44)
Mode, Base, EOM, Beta, N=data.iloc[:,0], data.iloc[:,1],data.iloc[:,3],data.iloc[:,2],data.iloc[:,7]

print(N)