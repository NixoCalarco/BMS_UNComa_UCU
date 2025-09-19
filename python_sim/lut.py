import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator

emf = np.loadtxt("./archivos/EMF_flipdim.txt", dtype=float)

emf_x = emf[:,0]
emf_y = emf[:,1]

dV = emf[1:, 0] - emf[:-1, 0]
dX = emf[1:, 1] - emf[:-1, 1]
Df = dV / dX
Df = np.column_stack((emf[1:, 1], Df))

x_to_y = PchipInterpolator(emf_x, emf_y)

res = 1000
x_grid = np.linspace(emf_x.min(), emf_x.max(), res)
x_to_y_grid = x_to_y(x_grid)

df_interp = PchipInterpolator(Df[:,0], Df[:,1])
df_grid = np.linspace(Df[:,0].min(), Df[:,0].max(), res)
df_values = df_interp(df_grid)

np.savetxt("lut.csv", np.column_stack([x_grid, x_to_y_grid]), delimiter=",")
np.savetxt("lut_df.csv", np.column_stack([df_grid, df_values]), delimiter=",")

plt.figure()
plt.plot(Df[:, 0], Df[:, 1])
plt.grid(True)
plt.show()
