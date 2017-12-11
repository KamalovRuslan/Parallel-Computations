#!/usr/bin/python

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np

with open("result.txt", 'r') as f:
	data = f.read()
data = data.split('\n')
solve = []
for row in data:
	row = row.split(' ')[:-1]
	row = [float(elem) for elem in row]
	solve.append(np.array(row))
solve = np.array(solve)
dim = solve.shape[0]
x = np.linspace(0.0, 2.0, solve.shape[0])
y = np.linspace(0.0, 2.0, solve.shape[1])
xx, yy = np.meshgrid(x, y, sparse=True)

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("u(x, y)")

surf = ax.plot_surface(xx, yy, solve, cmap=cm.coolwarm, linewidth=0, antialiased=False)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
fig.colorbar(surf, shrink=0.5, aspect=5)
ax.view_init(30, 30)

fig.savefig('images/dim_{}.eps'.format(dim))
plt.close(fig)
