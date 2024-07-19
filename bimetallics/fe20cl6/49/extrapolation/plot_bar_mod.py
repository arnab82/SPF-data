import numpy as np
import matplotlib.pyplot as plt

# Define the data
extrapolated = [-91.36401158, -92.16210385, -93.62562205, -95.71532926, -98.77681228]
pt2 = [-88.58433806, -89.44906798, -90.38841926, -91.7283117, -93.47335424]
tpsci = [-72.12045094, -74.44468695, -75.74580553, -76.3689305, -77.49977338]
methods = ['TPSCI', 'TPSCI+PT2', 'Extrapolated']
methods_2 = ['TPSCI', 'TPSCI+PT2']
J_values = [tpsci, pt2, extrapolated]
J_values_2 = [tpsci, pt2]
S = ["J(S0,S1)", "J(S1,S2)", "J(S2,S3)", "J(S3,S4)", "J(S4,S5)"]

# Print the values
for i, s in enumerate(S):
    print(f"{s}:")
    print(f"  TPSCI: {tpsci[i]}")
    print(f"  TPSCI+PT2: {pt2[i]}")
    print(f"  Extrapolated: {extrapolated[i]}")
    print()
J_0_1=[tpsci[0], pt2[0], extrapolated[0]]
J_1_2=[tpsci[1], pt2[1], extrapolated[1]]
J_2_3=[tpsci[2], pt2[2], extrapolated[2]]
J_3_4=[tpsci[3], pt2[3], extrapolated[3]]
J_4_5=[tpsci[4], pt2[4], extrapolated[4]]
J_1=[J_0_1[0], J_1_2[0], J_2_3[0], J_3_4[0], J_4_5[0]]
J_2=[J_0_1[1], J_1_2[1], J_2_3[1], J_3_4[1], J_4_5[1]]
J_3=[J_0_1[2], J_1_2[2], J_2_3[2], J_3_4[2], J_4_5[2]]
# Set colors for the bars
colors = ['red', 'blue', 'green']

# Set figure size
plt.figure(figsize=(10, 6))

# Set the desired width of the bar
bar_width = 0.3  

# Plotting the bars
# for i, method in enumerate(methods):
#     plt.bar(np.arange(len(S)) + i * bar_width, J_values[i], color=colors[i], width=bar_width, label=method)
for i, method in enumerate(methods_2):
    plt.bar(np.arange(len(S)) + i * bar_width, J_values_2[i], color=colors[i], width=bar_width, label=method)
# plt.bar(S, J_1, color=colors[0], width=bar_width, label='TPSCI')
# Set labels and legend
plt.xlabel('Methods and States', fontsize=14)
plt.ylabel('Values of J ($cm^{-1}$)', fontsize=14)
plt.xticks(np.arange(len(S)) + bar_width, S)
plt.ylim(-120,0)
plt.gca().invert_yaxis()
plt.legend(loc='upper left')  # Fixing the legend position to top left
plt.axhline(y=-117, color='gray', linestyle='-', alpha=0.9,linewidth=2.0)
plt.axhline(y=-81, color='violet', linestyle='--', alpha=0.9,linewidth=2.0)
# Show plot
plt.tight_layout()
plt.savefig("plot_bar_J_tpsci_pt2.pdf", dpi=600, bbox_inches='tight')
plt.show()
