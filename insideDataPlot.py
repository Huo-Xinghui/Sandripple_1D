# 导入必要的库
import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from typing import Dict

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
	# 定义常数
	interval = 120
	file_interval = 240
	start = interval*20
	end = interval*23
	plot_case = 18 #对照case_dict的键值
	plot_section = 1 #从0开始
	average_on_time = False
	average_on_y = True
	plot_type = 0 #0: pcolormesh, 1: contourf
	# 定义文件路径
	working_dir = "/home/ekalhxh/ripple/coll"
	#working_dir = "E:/Data/Sandripples1DFluid/ripple/coll"
	# 定义文件名字典
	case_dict = {
		1: "uStar040_150and350_0_2650_3600",
		2: "uStar045_150and350_0_2650_3600",
		3: "uStar050_150and350_0_2650_3600",
		4: "uStar055_150and350_0_2650_3600",
		5: "uStar060_150and350_0_2650_3600",
		6: "uStar050_150and450_0_2650_3600",
		7: "uStar050_150and550_0_2650_3600",
		8: "uStar050_200and400_0_2650_3600",
		9: "uStar050_250and350_0_2650_3600",
		10: "uStar050_300stdd5_0_2650_3600",
		11: "uStar050_300stdd10_0_2650_3600",
		12: "uStar050_300stdd20_0_2650_3600",
		13: "uStar050_300stdd50_0_2650_3600",
		14: "uStar035_300stdd100_0_2650_3600",
		15: "uStar040_300stdd100_0_2650_3600",
		16: "uStar045_300stdd100_0_2650_3600",
		17: "uStar050_300stdd100_0_2650_3600",
		18: "uStar055_300stdd100_0_2650_3600",
		19: "uStar060_300stdd100_0_2650_3600",
		20: "uStar065_300stdd100_0_2650_3600",
	}
	folder_name = case_dict[plot_case]
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
	last_time = int(parts[4])
	real_end = min(end, last_time)
	start_file_num = start // file_interval
	end_file_num = real_end // file_interval
	folder_path = f"{working_dir}/{folder_name}/Inside"
	time_step_inside_dict = read_inside_file(folder_path, start_file_num, end_file_num)
	sections_file = os.path.join(working_dir, f"{folder_name}/sections.dat")
	ave_section_file = os.path.join(working_dir, f"{folder_name}/ave_section.dat")
	if os.path.exists(sections_file):
		os.remove(sections_file)
	if os.path.exists(ave_section_file):
		os.remove(ave_section_file)
	y_section_time_step_dict: Dict[int, pd.DataFrame] = {}
	y_ave_section_time_step_dict: Dict[int, pd.DataFrame] = {}
	for t, new_current_inside_df in time_step_inside_dict.items():
		new_current_inside_df = new_current_inside_df.copy()
		new_current_inside_df['d'] = new_current_inside_df['d'].replace(0, None)
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
	plot_df['d'] = plot_df['d'].replace(np.nan, 0)
	#print(f"y={section_key}\n")
	#print(plot_df)

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
	#ax.set_aspect(3, 'box')

	# 插入放大图
	axins = inset_axes(ax, width="50%", height="30%", loc='upper right')
	if plot_type == 0:
		axins.pcolormesh(X, Z, d_grid_masked, shading='auto', cmap='viridis')
	elif plot_type == 1:
		axins.contourf(X, Z, d_grid_masked, cmap='viridis')
		axins.contour(X, Z, d_grid_masked, colors='black', linewidths=0.5)
	axins.set_xlim(0.17, 0.26)
	axins.set_ylim(0.0175, 0.0225)
	axins.set_xticklabels('')
	axins.set_yticklabels('')

	# 添加放大区域的边框
	ax.indicate_inset_zoom(axins, edgecolor="black")

	# 设置颜色映射范围
	#c.set_clim(dia1, dia2)

	plt.show()