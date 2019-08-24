from matplotlib import pyplot as plt
import numpy as np
Y = [[3000]*8]*4
X = [[2290,1552,1030,808,795,804,823,1050],[1154,845,1254,828,858,1319,9999,9999],
     [897,890,1332,978,9999,9999,9999,9999],[1034,1154,9999,9999,9999,9999,9999,9999]]
for i in range(4):
    for j in range(8):
        X[i][j]=795/X[i][j]

ax=plt.gca()  

ax.set_xticklabels( ('2', '4', '8', '16', '32', '64', '128', '256', '512'))
ax.set_yticklabels( ('0', '1', '2', '4', '8'))
plt.imshow(X,interpolation='nearest', cmap=plt.cm.hot)  
plt.colorbar(shrink=0.5)
#a.set_ticklabels(('11','22','33','44','55','66','77','88'))
plt.show() 
