# -*-*-*- Author: EkalHxH -*-*-*-
# version: 1.0 (2024-11-30)
#
#      ___           ___           ___           ___       ___           __            ___
#     /  /\         /  /\         /  /\         /  /\     /  /\         |  |\         /  /\
#    /  /::\       /  /:/        /  /::\       /  /:/    /  /:/         |  |:|       /  /:/
#   /  /:/\:\     /  /:/        /  /:/\:\     /  /:/    /  /:/          |  |:|      /  /:/
#  /  /::\ \:\   /  /::\____   /  /::\ \:\   /  /:/    /  /::\ ___      |__|:|__   /  /::\ ___
# /__/:/\:\ \:\ /__/:/\:::::\ /__/:/\:\_\:\ /__/:/    /__/:/\:\  /\ ____/__/::::\ /__/:/\:\  /\
# \  \:\ \:\_\/ \__\/~|:|~~~~ \__\/  \:\/:/ \  \:\    \__\/  \:\/:/ \__\::::/~~~~ \__\/  \:\/:/
#  \  \:\ \:\      |  |:|          \__\::/   \  \:\        \__\::/     |~~|:|          \__\::/
#   \  \:\_\/      |  |:|          /  /:/     \  \:\       /  /:/      |  |:|          /  /:/
#    \  \:\        |__|:|         /__/:/       \  \:\     /__/:/       |__|:|         /__/:/
#     \__\/         \__\|         \__\/         \__\/     \__\/         \__\|         \__\/

# ********************************************************************************************************

# 用于读取并输出床面内部的粒径分布情况

# ********************************************************************************************************

# 导入必要的库
import os # 用于文件操作
import re # 用于正则表达式
import numpy as np # 用于数值计算
import pandas as pd # 用于数据处理
import matplotlib.pyplot as plt # 用于绘图
from mpl_toolkits.axes_grid1.inset_locator import inset_axes # 用于绘制放大图
from typing import Dict # 用于类型提示

pd.set_option('future.no_silent_downcasting', True) # 设置pandas不显示警告

# 定义一个函数来读取任意文件
def read_file(file_path):
	if not os.path.exists(file_path):
		print(f"文件 {file_path} 不存在。")
		exit()

	with open(file_path, 'r', encoding='utf-8') as file:
		line = file.readlines()
		return line

# 读取多个文件内容, 并将其存储到一个字典中
def read_inside_file(folder_path, start_file, end_file):
	time_step_data_dict: Dict[int, pd.DataFrame] = {}
	for i in range(start_file, end_file+1):
		file_path = os.path.join(folder_path, f"InsideData_{i}.plt")
		results = read_file(file_path)

		current_time = None
		for line in results:
			if "variables" in line:
				continue
			elif "zone" in line:
				if current_time is not None:
					current_time_data = pd.DataFrame(current_time_data_list, columns=["x", "y", "z", "d"])
					time_step_data_dict[current_time] = current_time_data
				current_time_list = [float(num) for num in re.findall(r'\d+\.?\d*', line)]
				current_time = round(current_time_list[0])
				current_time_data_list = []
			elif "i=" not in line:
				data_values = [float(num) for num in line.split()]
				if len(data_values) == 4:
					current_time_data_list.append(data_values)
		current_time_data = pd.DataFrame(current_time_data_list, columns=["x", "y", "z", "d"])
		time_step_data_dict[current_time] = current_time_data

	return time_step_data_dict

# 主程序
if __name__ == "__main__":
	# 判断操作系统
	sys_OS = "w" # "w" for windows, "l" for linux
	if sys_OS == "l":
		linux_flag = True
	elif sys_OS == "w":
		linux_flag = False
	else:
		print("Invalid input!")
		exit()

# ----------------------------------------------------------------------------------------
	# 定义常数
	interval = 120 #源文件时间间隔
	file_interval = 240 #两个源文件之间的时间间隔
	start = interval*29 #求平均的起始时间，若average_on_time为False，则表示所显示的瞬时时刻
	end = interval*30 #求平均的终止时间，若average_on_time为False，则无效
	plot_case = 1 #对照case_dict的键值
	plot_section = 1 #从0开始，在average_on_y为False时表示所显示的截面在y方向的位置
	average_on_time = False #是否对时间求平均
	average_on_y = True #是否对y方向求平均
	remove_zero = True #是否在求平均前将0值替换为None，以保障平均后不出现不真实表面
	bed_profile = False #是否绘制床面廓线
	correct_d = True #是否对d值进行修正
	plot_type = 1 #绘图方式 0: pcolormesh, 1: contourf
# ----------------------------------------------------------------------------------------

	# 定义文件路径
	if linux_flag:
		working_dir = "/home/ekalhxh/ripple/coll11"
	else:
		working_dir = "E:/Data/Sandripples1DFluid/ripple/coll11"

	# 定义文件名字典
	case_dict = {
		1: "uStar040_300log80_0_2650_3600",
		2: "uStar050_300log80_0_2650_3600",
		3: "uStar060_300log80_0_2650_3600",
		4: "uStar050_400log80_0_2650_3600",
		5: "uStar050_400log100_0_2650_3600",
		6: "uStar050_400log120_0_2650_3600",
	}

	folder_name = case_dict[plot_case] #工作目录名

	# 从目录名中提取粒径信息
	parts = folder_name.split("_")
	dia_name = parts[1]
	if "and" in dia_name:
		dia_name_list = dia_name.split("and")
		dia1 = float(dia_name_list[0])/1e6
		dia2 = float(dia_name_list[1])/1e6
	elif "stdd" in dia_name:
		dia_name_list = dia_name.split("stdd")
		dia1 = (float(dia_name_list[0])-1*float(dia_name_list[1]))/1e6
		dia2 = (float(dia_name_list[0])+1*float(dia_name_list[1]))/1e6
	elif "log" in dia_name:
		dia_name_list = dia_name.split("log")
		dia1 = (float(dia_name_list[0])-1*float(dia_name_list[1]))/1e6
		dia2 = (float(dia_name_list[0])+1*float(dia_name_list[1]))/1e6
	dia = (dia1+dia2)/2

	# 根据起止时间确定读取的Inside文件序号
	last_time = int(parts[4])
	real_end = min(end, last_time)
	start_file_num = start // file_interval
	end_file_num = real_end // file_interval

	folder_path = f"{working_dir}/{folder_name}/Inside" #Inside文件所在文件夹路径
	time_step_inside_dict = read_inside_file(folder_path, start_file_num, end_file_num) #读取Inside文件

	y_section_time_step_dict: Dict[int, pd.DataFrame] = {} #存储各时刻y方向各截面的字典
	y_ave_section_time_step_dict: Dict[int, pd.DataFrame] = {} #存储各时刻y方向平均截面的字典
	for t, current_inside_df in time_step_inside_dict.items():
		# 若remove_zero为True，则将各截面0值替换为None，为之后求平均做准备
		new_current_inside_df = current_inside_df.copy()
		if remove_zero:
			new_current_inside_df['d'] = new_current_inside_df['d'].replace(0, None)

		# 将截面数据存储到字典中，并得到输出的截面数据plot_df
		if average_on_time:
			if t>=start and t<=end:
				if average_on_y:
					ave_section = new_current_inside_df.groupby(['x', 'z']).mean().reset_index()
					y_ave_section_time_step_dict[t] = ave_section
					section_key = None
				else:
					y_sections = new_current_inside_df.groupby('y')
					sections_list = list(y_sections)
					section_key, section_df = sections_list[plot_section]
					y_section_time_step_dict[t] = section_df
		else:
			if t == start:
				if average_on_y:
					plot_df = new_current_inside_df.groupby(['x', 'z']).mean().reset_index()
					section_key = None
				else:
					y_sections = new_current_inside_df.groupby('y')
					sections_list = list(y_sections)
					section_key, plot_df = sections_list[plot_section]
	if average_on_time:
		if average_on_y:
			plot_df = pd.concat(y_ave_section_time_step_dict.values()).groupby(['x', 'z']).mean().reset_index()
		else:
			plot_df = pd.concat(y_section_time_step_dict.values()).groupby(['x', 'z']).mean().reset_index()
	plot_df['d'] = plot_df['d'].replace(np.nan, 0).infer_objects(copy=False) #将None值替换回0

	# 对d值进行修正
	if correct_d:
		# 获取唯一的x和z值
		x_unique = plot_df['x'].unique()
		z_unique = plot_df['z'].unique()
		# 对x_unique和z_unique中的值由小到大排序
		x_unique.sort()
		z_unique.sort()

		# 将 DataFrame 转换为 NumPy 数组
		x_values = plot_df['x'].values
		z_values = plot_df['z'].values
		d_values = plot_df['d'].values

		# 遍历所有x值
		for i, x in enumerate(x_unique):
			if i == 0 or i == len(x_unique) - 1:
				continue
			# 获取当前x值对应的所有行索引
			x_indices = np.where(x_values == x)[0]

			# 遍历所有z值
			for j, z in enumerate(z_unique):
				# 获取当前节点的d值
				current_indices = x_indices[z_values[x_indices] == z]
				if len(current_indices) == 0:
					continue
				current_d = d_values[current_indices[0]]
				if current_d == 0:
					continue

				# 获取当前节点左右邻居的x值
				left_x = x_unique[i - 1]
				right_x = x_unique[i + 1]

				# 获取左右邻居的d值
				left_indices = np.where((x_values == left_x) & (z_values == z))[0]
				left_d = d_values[left_indices[0]] if len(left_indices) > 0 else None

				right_indices = np.where((x_values == right_x) & (z_values == z))[0]
				right_d = d_values[right_indices[0]] if len(right_indices) > 0 else None

				# 获取当前节点下面邻居的z值
				if j <= 1:
					below1_z = z_unique[j]
					below2_z = below1_z
				else:
					below1_z = z_unique[j - 1]
					below2_z = z_unique[j - 2]

				# 获取下面邻居的d值
				below1_indices = np.where((x_values == x) & (z_values == below1_z))[0]
				below1_d = d_values[below1_indices[0]] if len(below1_indices) > 0 else None

				below2_indices = np.where((x_values == x) & (z_values == below2_z))[0]
				below2_d = d_values[below2_indices[0]] if len(below2_indices) > 0 else None

				# 检查左右邻居的d值是否相等
				if left_d == right_d and current_d > dia:
					d_values[current_indices[0]] = left_d

				# 检查下面邻居的d值是否相等
				if below1_d == below2_d and current_d > dia:
					d_values[current_indices[0]] = below1_d

		# 将修改后的d值更新回DataFrame
		plot_df['d'] = d_values

	# 提取床面廓线
	nonzero_df = plot_df[plot_df['d'] != 0] #过滤出d不为0的网格节点
	idx = nonzero_df.groupby('x')['z'].idxmax() #找到每列x中z最大的行索引
	# 提取床面廓线的x和z坐标
	x_bed = nonzero_df.loc[idx, 'x'].tolist()
	z_bed = nonzero_df.loc[idx, 'z'].tolist()

	# 创建网格
	x_unique = np.sort(plot_df['x'].unique())
	z_unique = np.sort(plot_df['z'].unique())
	X, Z = np.meshgrid(x_unique, z_unique)

	# 创建数据透视表并转换为numpy数组
	d_grid = plot_df.pivot_table(index='z', columns='x', values='d', aggfunc='mean').values

	# 将d_grid中的0值进行掩码处理
	d_grid_masked = np.ma.masked_equal(d_grid, 0)

	# 创建图形和轴
	fig, ax = plt.subplots(figsize=(15, 5))

	if plot_type == 0:
		# 使用pcolormesh绘制伪彩色图
		c = ax.pcolormesh(X, Z, d_grid_masked, shading='auto', cmap='viridis')
	elif plot_type == 1:
		# 使用contourf绘制等高线图
		c = ax.contourf(X, Z, d_grid_masked, cmap='viridis')

	# 添加颜色条
	cbar = fig.colorbar(c, ax=ax, orientation='horizontal', fraction=0.1, pad=0.2)
	cbar.ax.ticklabel_format(style='scientific', scilimits=(0,0))

	# 绘制床面廓线，并将其覆盖在图形上
	if bed_profile:
		ax.plot(x_bed, z_bed, color='black', linewidth=1.5)

	# 设置主图的显示范围
	ax.set_ylim(0.015, 0.035)

	# 设置标题和标签
	if average_on_time:
		if average_on_y:
			ax.set_title(f"Average InsideData at y=All, t={start}-{real_end}")
		else:
			ax.set_title(f"Average InsideData at y={section_key}, t={start}-{real_end}")
	else:
		if average_on_y:
			ax.set_title(f"InsideData at y=All, t={start}")
		else:
			ax.set_title(f"InsideData at y={section_key}, t={start}")
	ax.set_xlabel("x")
	ax.set_ylabel("z")

	# 设置坐标轴纵横比和图形大小
	ax.set_aspect(3, 'box')

	# 插入放大图,并设置放大区域的颜色范围
	axins = inset_axes(ax, width="50%", height="30%", loc='upper right')
	if plot_type == 0:
		c_inset = axins.pcolormesh(X, Z, d_grid_masked, shading='auto', cmap='viridis')
	elif plot_type == 1:
		c_inset = axins.contourf(X, Z, d_grid_masked, cmap='viridis')
		c_inset = axins.contour(X, Z, d_grid_masked, colors='black', linewidths=0.5)
	axins.set_xlim(0.2, 0.3)
	axins.set_ylim(0.0175, 0.025)
	axins.set_xticklabels('')
	axins.set_yticklabels('')
	axins.set_aspect(1, 'box')

	# 设置颜色映射范围
	c.set_clim(dia1, dia2)
	c_inset.set_clim(dia1, dia2) #设置放大区域的颜色范围

	# 添加放大区域的边框
	ax.indicate_inset_zoom(axins, edgecolor="black")

	plt.show()