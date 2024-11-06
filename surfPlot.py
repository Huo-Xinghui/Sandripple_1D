# 导入必要的库
import os # 用于文件操作
import re # 用于正则表达式
import matplotlib.pyplot as plt # 用于绘图
import numpy as np # 用于数值计算
from dataclasses import dataclass # 用于定义数据类
from typing import List # 用于定义列表
from typing import Dict # 用于定义字典

# 定义一个数据类，用于存储文件中的数据
@dataclass
class file_data:
	data: List[float]

# 定义一个函数来读取任意文件
def read_file(file_path):
	if not os.path.exists(file_path):
		print(f"文件 {file_path} 不存在。")
		exit()

	with open(file_path, 'r', encoding='utf-8') as file:
		line = file.readlines()
		return line

# 读取廓线文件内容, 并将其存储到一个字典中
def read_data_file1(folder_path, file_name):
	file_path = os.path.join(folder_path, file_name) # 拼接文件路径
	all_lines = read_file(file_path) # 读取文件内容

	time_step_data_dict: Dict[int, List[file_data]] = {} # 定义一个字典，用于存储数据
	for line in all_lines:
		if "vs" in line:
			nums_in_first_line = [int(num) for num in re.findall(r'\d+\.?\d*', line)]
			current_time = nums_in_first_line[0]
			point_num = nums_in_first_line[1]
			data_list = []
			point_count = 0
		else:
			point_count += 1
			columns = line.split()
			columns_float = [float(column) for column in columns]
			current_data = file_data(
				data=columns_float
			)
			data_list.append(current_data)
			if point_count == point_num:
				time_step_data_dict[current_time] = data_list
	return time_step_data_dict

# 读取时间序列文件内容, 并将其存储到一个列表中
def read_data_file2(folder_path, file_name):
	file_path = os.path.join(folder_path, file_name) # 拼接文件路径
	all_lines = read_file(file_path)

	data_list = [] # 定义一个列表，用于存储数据
	for line in all_lines:
		if "vs" not in line:
			columns = line.split()
			columns_float = [float(column) for column in columns]
			current_data = file_data(
				data=columns_float
			)
			data_list.append(current_data)
	return data_list

# 主程序
if __name__ == "__main__":
	# 判断操作系统
	sys_OS = "l" # "w" for windows, "l" for linux
	if sys_OS == "l":
		linux_flag = True
	elif sys_OS == "w":
		linux_flag = False
	else:
		print("Invalid input!")
		exit()

	# 固定参数
	nu = 1.51e-5 # 运动粘度
	interval = 30 # 源文件时间间隔
	# 控制参数
	output_num = 0 # 出图类型：0为床面廓线，1为相关性，2为波长, 3为波高, 4为波速, 5为床面粒径分布
	start = 30 # 起始时间
	end = 2400 # 结束时间
	average_start = 30 # 开始计算平均值的时间
	average_end = 900 # 结束计算平均值的时间
	profile_offset = 2e-4 # 廓线图的纵向偏移
	corr_offset = 1e-8 # 相关性图的纵向偏移
	diameter_offset = 2e-5 # 廓线图的纵向偏移
	case_num = 30 # 算例号

	# 定义算例字典
	case_num_dict = {
		#0: 0,
		#1: 1,
		2: 2,
		3: 3,
		4: 4,
		5: 5,
		6: 6,
		#7: 7,
		#8: 8,
		#9: 9,
		#10: 10,
		#11: 11,
		#12: 12,
		#13: 13,
		#14: 14,
		#15: 15,
		#16: 16,
		17: 17,
		18: 18,
		19: 19,
		20: 20,
		21: 21,
		#22: 22,
		#23: 23,
		#24: 24,
		#25: 25,
		26: 26,
		27: 27,
		28: 28,
		29: 29,
		30: 30,
		#31: 31,
		#32: 32,
		#33: 33,
		#34: 34,
		#35: 35,
		#36: 36,
		#37: 37,
		#38: 38,
		#39: 39,
		#40: 40,
		#41: 41,
		#42: 42,
		#43: 43,
	}

	# 定义文件路径
	if linux_flag:
		working_dir = "/home/ekalhxh/ripple/coll"
	else:
		working_dir = "E:/Data/Sandripples1DFluid/ripple/coll"

	# 定义文件名字典
	case_dict = {
		0: "uStar030_250_0_2650_3600",
		1: "uStar035_250_0_2650_3600",
		2: "uStar040_250_0_2650_3600",
		3: "uStar045_250_0_2650_3600",
		4: "uStar050_250_0_2650_3600",
		5: "uStar055_250_0_2650_3600",
		6: "uStar060_250_0_2650_3600",
		7: "uStar065_250_0_2650_3600",
		8: "uStar045_200_0_2650_3600",
		9: "uStar045_250_1_2000_3600",
		10: "uStar045_250_1_3000_3600",
		11: "uStar045_250_1_4000_3600",
		12: "uStar045_250_2_2650_3600",
		13: "uStar045_250_3_2650_3600",
		14: "uStar045_300_0_2650_3600",
		15: "uStar045_350_0_2650_3600",
		16: "uStar045_400_0_2650_3600",
		17: "uStar040_150and350_0_2650_3600",
		18: "uStar045_150and350_0_2650_3600",
		19: "uStar050_150and350_0_2650_3600",
		20: "uStar055_150and350_0_2650_3600",
		21: "uStar060_150and350_0_2650_3600",
		22: "uStar050_150and450_0_2650_3600",
		23: "uStar050_150and550_0_2650_3600",
		24: "uStar050_200and400_0_2650_3600",
		25: "uStar050_250and350_0_2650_3600",
		26: "uStar050_300stdd5_0_2650_3600",
		27: "uStar050_300stdd10_0_2650_3600",
		28: "uStar050_300stdd20_0_2650_3600",
		29: "uStar050_300stdd50_0_2650_3600",
		30: "uStar035_300_0_2650_3600",
		31: "uStar040_300_0_2650_3600",
		32: "uStar045_300_0_2650_3600",
		33: "uStar050_300_0_2650_3600",
		34: "uStar055_300_0_2650_3600",
		35: "uStar060_300_0_2650_3600",
		36: "uStar035_300stdd100_0_2650_3600",
		37: "uStar040_300stdd100_0_2650_3600",
		38: "uStar045_300stdd100_0_2650_3600",
		#39: "uStar050_300stdd100_0_2650_2774",
		40: "uStar055_300stdd100_0_2650_3600",
		41: "uStar060_300stdd100_0_2650_3600",
		42: "uStar065_300stdd100_0_2650_3600",
	}
	output_file_dict = {
		0: "surfProfile.dat",
		1: "autoCorr.dat",
		2: "surfTopology.dat",
		3: "surfTopology.dat",
		4: "surfTopology.dat",
		5: "surfDiameter.dat",
	}

	output_file = output_file_dict[output_num]
	if output_num == 0 or output_num == 1 or output_num == 5:
		folder_name = case_dict[case_num]
		folder_path = f"{working_dir}/{folder_name}"
		all_plot_data = read_data_file1(folder_path, output_file)
		offset = 0
		# 设置图像的宽度和高度
		plt.figure(figsize=(10, 3))  # 宽度为10，高度为3
		for t in range(start, end+1, interval):
			current_plot_data = all_plot_data[t]
			x = [point.data[0] for point in current_plot_data]
			y = [point.data[1]+offset for point in current_plot_data]
			if output_num == 0:
				offset += profile_offset
			elif output_num == 1:
				offset += corr_offset
			elif output_num == 5:
				offset += diameter_offset
			plt.plot(x, y, color='black', linewidth=0.5)  # 设置线的颜色为黑色，线的粗细为0.5

		plt.xlabel('x')
		plt.ylabel('t')
		plt.gca().yaxis.set_ticks([])  # 移除纵坐标刻度
		plt.xlim(min(x), max(x))  # 设置横坐标范围为数据的最小值和最大值
		plt.show()
	elif output_num >= 2 and output_num <= 4:
		y_list = []
		y_max = 0
		y_min = 1000
		for i, case in case_num_dict.items():
			folder_name = case_dict[case]
			folder_path = f"{working_dir}/{folder_name}"
			parts = folder_name.split("_")
			ustar_name = parts[0]
			ustar_str = ustar_name[6:]
			ustar_int = int(ustar_str)
			ustar = ustar_int / 100
			dia_name = parts[1]
			if "and" in dia_name:
				dia_name_list = dia_name.split("and")
			else:
				dia_name_list = [dia_name, dia_name]
			dia1 = float(dia_name_list[0])/1e6
			dia2 = float(dia_name_list[1])/1e6
			dia =  (dia1 + dia2) / 2
			rho_name = parts[2]
			if rho_name == "0":
				rho = 1.263
			else:
				rho = float(rho_name)
			rhoP = float(parts[3])
			s = rhoP/rho
			g_hat = 9.8 * (1 - 1/s)
			Sh = rho*ustar**2/(rhoP*g_hat*dia)
			Ga = (s*g_hat*dia**3)**0.5/nu
			all_plot_data = read_data_file2(folder_path, output_file)
			x = [point.data[0] for point in all_plot_data]
			if output_num == 2:
				y = [point.data[1] for point in all_plot_data]
			elif output_num == 3:
				y = [point.data[2] for point in all_plot_data]
			elif output_num == 4:
				y = [point.data[3] for point in all_plot_data]
			i_start = start // interval - 1
			i_end = end // interval - 1
			y_max = max(y_max, max(y[i_start:i_end]))
			y_min = min(y_min, min(y[i_start:i_end]))
			# Apply 3-point box filter to y
			y_filtered = np.convolve(y, np.ones(3)/3, mode='same')
			y = y_filtered
			i_start = average_start // interval - 1
			i_end = average_end // interval - 1
			average_y = np.mean(y[i_start:i_end])
			y_list.append((ustar, dia1, dia2, dia, Sh, Ga, average_y))
			plt.figure(1)
			if dia1 == dia2:
				plt.plot(x, y, label=f"$u*={ustar}$, $d={dia}$")
			else:
				plt.plot(x, y, label=f"$u*={ustar}$, $d_1={dia1}$, $d_2={dia2}$")
		plt.xlabel('$t$')
		if output_num == 2:
			plt.ylabel('$\lambda$')
			plt.ylim(0, y_max)
		elif output_num == 3:
			plt.yscale('log')
			plt.ylabel('$A$')
			plt.ylim(y_min, y_max)
		elif output_num == 4:
			plt.ylabel('$c$')
			plt.ylim(0, y_max)
		plt.legend(loc='upper left')
		plt.xlim(start, end)
		plt.show()

		plt.figure(2)
		plt.plot([item[0] for item in y_list], [item[6] for item in y_list], 'o')
		plt.xlabel('$u*$')
		if output_num == 2:
			plt.ylabel('$\lambda$')
		elif output_num == 3:
			plt.ylabel('$A$')
		elif output_num == 4:
			plt.ylabel('$c$')
		plt.show()