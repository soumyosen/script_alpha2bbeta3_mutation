import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

fig, axs = plt.subplots(nrows=3, ncols=3, figsize=(9,7))
fig.subplots_adjust(hspace =.2, wspace=.2)
axs = axs.ravel()

for i in range(9):
     df = pd.read_csv("contact_%s.dat" % i)
     #print(df)
     x_axis = df.columns[1:].tolist()
     y_axis = df['Residues'].tolist()
     #print(x_axis)
     #print(y_axis)
     xTickMarks = [" "]
     yTickMarks = [" "]
     for j in x_axis:
         xTickMarks.append(j)
     
     for j in y_axis:
         yTickMarks.append(j)
     
     xTickMarks.append(" ")
     yTickMarks.append(" ")
     
     #print(xTickMarks)
     #print(yTickMarks)
     
     
     y = [yTickMarks.index(i) for i in y_axis]
     y_len = len(y)
     x = [xTickMarks.index(i) for i in x_axis]
     x0 = np.repeat(x[0],y_len)
     x1 = np.repeat(x[1],y_len)
     
     colors = {1:'b', 0:'r'}
     
     #print(x)
     #print(y)
    
     axs[i].scatter(x0, y, s=25, c=[colors[r] for r in df[x_axis[0]]])
     axs[i].scatter(x1, y, s=25, c=[colors[r] for r in df[x_axis[1]]])
     axs[i].set_xlim([0, len(xTickMarks)-1])
     axs[i].set_ylim([0, len(yTickMarks)-1])
     axs[i].set_xticks(range(len(xTickMarks)))
     axs[i].set_yticks(range(len(yTickMarks)))
     axs[i].set_xticklabels(xTickMarks)
     axs[i].set_yticklabels(yTickMarks)


axs[0].scatter(0.2, 4.6, s=25, c='b')
axs[0].scatter(1.5, 4.6, s=25, c='r')
axs[0].text(0.3, 4.5, ': present', color='blue', fontsize=12)
axs[0].text(1.6, 4.5, ': absent', color='red', fontsize=12)


#circle1 = plt.Circle((2.5, 4.5), 0.05, color='r')
#circle2 = plt.Circle((2.5, 4.0), 0.05, color='blue')
#axs[0].set_aspect(1)
#axs[0].add_patch(circle1)
#axs[0].add_patch(circle2)

plt.tight_layout()
#plt.show()
plt.savefig("interaction.png", dpi=400)
