#!/usr/bin/env python
import argparse
import numpy as np
import matplotlib.pyplot as plt


relative = 0
edges_ref = []
extrapolated=[-91.36401158,
-92.16210385,
-93.62562205,
-95.71532926,
-98.77681228]
pt2=[-88.58433806,
-89.44906798,
-90.38841926,
-91.7283117,
-93.47335424]
tpsci=[-72.12045094,
-74.44468695,
-75.74580553,
-76.3689305,
-77.49977338]
methods = ['TPSCI', 'PT2', 'Extrapolated']
J_0_1=[tpsci[0], pt2[0], extrapolated[0]]
J_1_2=[tpsci[1], pt2[1], extrapolated[1]]
J_2_3=[tpsci[2], pt2[2], extrapolated[2]]
J_3_4=[tpsci[3], pt2[3], extrapolated[3]]
J_4_5=[tpsci[4], pt2[4], extrapolated[4]]
J_1=[J_0_1[0], J_1_2[0], J_2_3[0], J_3_4[0], J_4_5[0]]
J_2=[J_0_1[1], J_1_2[1], J_2_3[1], J_3_4[1], J_4_5[1]]
J_3=[J_0_1[2], J_1_2[2], J_2_3[2], J_3_4[2], J_4_5[2]]
S=["J(S0,S1)","J(S1,S2)","J(S2,S3)","J(S3,S4)","J(S4,S5)"]
print(J_0_1)
J=[J_0_1, J_1_2, J_2_3, J_3_4, J_4_5]
colors = ['red', 'blue', 'green', 'orange', 'purple']
print(colors)
plt.figure(figsize=(5, 5))  # Set figure size
bar_width = 0.3  # Set the desired width of the bar
plt.bar(S, J_1, color=colors[0], width=bar_width, label='TPSCI')
plt.bar(S, J_2, color=colors[1], width=bar_width, label='TPSCI+PT2')
plt.bar(S, J_3, color=colors[2], width=bar_width, label='EXTRAPOLATED')

plt.xlabel('Methods', fontsize=20)  # Set font size
plt.ylabel('values of J ($cm^{-1}$)', fontsize=20)  # Set font size
plt.gca().invert_yaxis()

plt.legend()  # Add legend

plt.show()
plt.xlabel('Methods', fontsize=20)  # Set font size
plt.ylabel('values of J ($cm^{-1}$)', fontsize=20)  # Set font size

# Add text annotations to the bars
# for i, v in enumerate(J_0_1):
#     plt.text(i, v, f"{v:.2f}", ha='center', va='bottom')

# Add negative sign to y-axis labeling
plt.gca().invert_yaxis()




# plt.title('Optimization Time for Different Methods')
# plt.savefig("plot_bar_ph.pdf", dpi=600, bbox_inches='tight')
plt.show()







 




    