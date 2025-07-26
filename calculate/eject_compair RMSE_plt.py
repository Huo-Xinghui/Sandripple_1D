import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.ndimage import gaussian_filter
from scipy.ndimage import zoom

mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['text.usetex'] = True
label_size = 40
ticks_size = 30

fig = plt.figure(1, figsize=(16, 12))
data = np.loadtxt('eject_compair_rmse_final.txt', delimiter=',')
read_x = data[:, 0]
read_y = data[:, 1]
read_z = data[:, 2]
mesh_x = np.unique(read_x)
mesh_y = np.unique(read_y)
mesh_z = read_z.reshape(len(mesh_y), len(mesh_x))
mesh_z_smooth = gaussian_filter(mesh_z, sigma=1.0)
read_z_smooth = mesh_z_smooth.flatten()
mesh_z_fine = zoom(mesh_z_smooth, (10, 10), order=3)
#plt.imshow(mesh_z_smooth, cmap='jet', extent=(mesh_x.min(), mesh_x.max(), mesh_y.min(), mesh_y.max()), origin='lower', aspect='auto')
plt.contourf(mesh_x, mesh_y, mesh_z_smooth, levels=np.linspace(0, 4, 1000), cmap='jet', extend='both')
cbar = plt.colorbar(ticks=np.linspace(0, 4, 9))
cbar.set_label('$\Delta^{\mathrm{norm}}_{\mathrm{ave}}$', rotation=90, labelpad=ticks_size, fontsize=label_size)
cbar.ax.tick_params(labelsize=ticks_size)
plt.xlabel('$\\varepsilon$', fontsize=label_size)
plt.ylabel('$\\nu$', fontsize=label_size)
plt.tick_params(labelsize=ticks_size)

# 添加 inset axes
inset_ax = fig.add_axes([0.52, 0.65, 0.2, 0.2])  # [left, bottom, width, height]，可根据需要调整

xmin = 1e-4
xmax = 6e-4
xnum = 100000
dx = xmax / xnum
ddx = dx * 0.5
mu = -8.30271
sigma = 0.25778
x = np.linspace(xmin+dx, xmax, xnum)
x = x - ddx
y = (1 / (x * sigma * np.sqrt(2 * np.pi))) * np.exp(-(np.log(x) - mu) ** 2 / (2 * sigma ** 2))
y = y / y.sum()
cum_y = np.cumsum(y)
x_plot = x*1e6
inset_ax.plot(x_plot, cum_y, 'b-', label='Fited CDF', linewidth=4)  # 线宽为4
inset_ax.set_xlim(x_plot.min(), x_plot.max())
inset_ax.set_xticks([x_plot.min(), x_plot.max()])
inset_ax.set_yticks([cum_y.min(), cum_y.max()])
x_data = np.array([0.000150, 0.000180, 0.000212, 0.000250, 0.000300, 0.000355, 0.000425, 0.000500])
x_data_plot = x_data * 1e6
p_data = np.array([0.00785, 0.08593, 0.28992, 0.52017, 0.76283, 0.91208, 0.97593, 1.0])
inset_ax.plot(x_data_plot, p_data, 'ro', label='Wi89', markersize=15)  # 点大小为12
inset_ax.set_xlabel('$d (\mu \mathrm{m})$', fontsize=ticks_size, labelpad=-10)
inset_ax.set_ylabel('CDF', fontsize=ticks_size, labelpad=-10)
inset_ax.legend(fontsize=ticks_size, loc='upper left', bbox_to_anchor=(-1.3, 1.1))
inset_ax.tick_params(labelsize=ticks_size, direction='in', top=True, right=True)
#inset_ax.xaxis.set_label_position('top')    # x轴标签放到上方
#inset_ax.xaxis.tick_top()                   # x轴刻度也放到上方（可选）
#inset_ax.yaxis.set_label_position('right')  # y轴标签放到右侧
#inset_ax.yaxis.tick_right()                 # y轴刻度也放到右侧（可选）

plt.show()