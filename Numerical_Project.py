import math
import numpy as np
import matplotlib.pyplot as plt

# Define Given Parameters
h_LE = 200  # Convection Coefficient for leading edge (W/m^2 K)
h_TE = 100  # Convection Coefficient for trailing edge (W/m^2 K)
h_top = 150  # Convection Coefficient for top edge (W/m^2 K)
h_bot = h_top  # Convection Coefficient for bottom edge (W/m^2 K)
T_inf = 800  # Ambient temperature surrounding blade (deg C)
k = 11.7  # Thermal Conductivity of Inconcel (W/m K)

# Define Blade Dimensions
length = 0.075  # (m)
width = 0.025  # (m)

# Define mesh size
n = 25  # grid size
x = np.linspace(0, length, n)  # Location of nodes along x-dir
y = np.linspace(-width / 2, width / 2, n)  # Location of nodes along y-dir
dx = abs(x[1] - x[2])  # x-distance between nodes (m)
dy = abs(y[2] - y[1])  # y-distance between nodes (m)

# Simplify Parameters
r = dx / dy  # ratio of x/y distance between nodes
Bi_xtop = h_top * dx / k  # Biot number for the top boundary
Bi_xbot = h_bot * dx / k  # Biot number for the bottom boundary
Bi_yLE = h_LE * dy / k  # Biot number for the leading edge boundary
Bi_yTE = h_TE * dy / k  # Biot number for the trailing edge boundary

print("Solving for Energy Sink Distribution along the Cross-Section of the Blade...")


# Energy Sink Distribution Parameters
x_o = np.multiply([12, 28, 43, 55, 67, 67], 1e-3)
y_o = np.multiply([0, 0, 0, 0, 5, -5], 1e-3)
A_n = -15 * 1e5
sigma_xn = np.multiply([5, 4.5, 4, 4, 3, 3], 1e-3)
sigma_yn = np.multiply([3, 2.5, 2, 2, 1, 1], 1e-3)
theta = np.pi / 2
a_n = np.zeros(6)
b_n = np.zeros(6)
c_n = np.zeros(6)
for i in range(6):
    a_n[i] = np.divide((math.cos(theta) ** 2), (2 * sigma_xn[i] ** 2)) + np.divide((math.sin(theta) ** 2), (2 * sigma_yn[i] ** 2))
    b_n[i] = -np.divide((math.sin(2 * theta)), (4 * sigma_xn[i] ** 2)) + np.divide((math.sin(2 * theta)), (4 * sigma_yn[i] ** 2))
    c_n[i] = np.divide((math.sin(theta) ** 2), (2 * sigma_xn[i] ** 2)) + np.divide((math.cos(theta) ** 2), (2 * sigma_yn[i] ** 2))

# Energy Sink Coordinates
x_egen = np.linspace(0, length, n)
y_egen = np.linspace(-width / 2, width / 2, n)
XX, YY = np.meshgrid(x_egen, y_egen)
YY = np.flipud(YY)

# Energy Sink Function
e_gen = lambda x, y: A_n * np.exp(-(np.multiply(a_n[1], (np.abs(x - x_o[1]) ** 2)) +
                                    np.multiply((2 * b_n[1]), np.abs(x - x_o[1]), np.abs(y - y_o[1])) +
                                    np.multiply(c_n[1], (np.abs(y-y_o[1]) **2))))

E = e_gen(XX, YY)       # initialize energy sink matrix

# Plot energy sink within blade (3a)
plt.contourf(XX, YY, E)
plt.show()
