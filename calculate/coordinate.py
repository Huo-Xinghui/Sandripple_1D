import matplotlib.pyplot as plt
import numpy as np

# 定义向量的起点和方向
x = [0.68, -0.72, -0.13]
y = [0.22, 0.03, 0.98]
z = [-0.7, -0.7, 0.18]

# 创建一个新的图形
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# 绘制向量箭头
ax.quiver(0, 0, 0, x[0], y[0], z[0], color='r', label='Vector 1')
ax.quiver(0, 0, 0, x[1], y[1], z[1], color='g', label='Vector 2')
ax.quiver(0, 0, 0, x[2], y[2], z[2], color='b', label='Vector 3')

# 设置坐标轴标签
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

# 设置图例
ax.legend()

# 设置标题
ax.set_title('3D Vector Arrows')

# 显示图形
plt.show()