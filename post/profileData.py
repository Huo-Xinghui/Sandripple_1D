# -*-*-*- Author: EkalHxH -*-*-*-
# version: 1.0 (2024-12-13)
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
"""用于统计流场廓线相关数据"""
# ********************************************************************************************************

# 导入必要的库
import math # 导入数学库
import os # 导入操作系统库
import re # 导入正则表达式库
import pandas as pd # 导入pandas库
import numpy as np # 导入numpy库
import matplotlib.pyplot as plt # 导入matplotlib库
from dataclasses import dataclass # 导入dataclass库
from dataclasses import field # 导入field库
from typing import List # 导入List类型
from typing import Dict # 导入Dict类型
from tqdm import tqdm # 导入进度条库

def read_file(file_path: str) -> list[str]:
	"""读取文件"""
	if not os.path.exists(file_path):
		print(f"文件 {file_path} 不存在。")
		exit()

	with open(file_path, 'r', encoding='utf-8') as file:
		data = file.readlines()
		return data

def find_indices_with_keyword(data: List[str], keyword: str) -> List[int]:
	return [i for i, line in enumerate(data) if keyword in line]

def read_profile_file(folder_path: str) -> list[np.ndarray]:
	"""读取廓线文件"""
	file_path = os.path.join(folder_path, "Profile.plt")
	results = read_file(file_path)
	indices_with_zone = find_indices_with_keyword(results, "zone")

	time_step_data = list[np.ndarray]()
	for i in range(len(indices_with_zone)):
		index1 = indices_with_zone[i] + 1
		index2 = indices_with_zone[i+1] - 1 if i < len(indices_with_zone) - 1 else len(results) - 1
		lines = results[index1:index2]
		current_data = np.ndarray(shape=(len(lines), 6), dtype=float)
		for j, line in enumerate(lines):
			values = [float(num) for num in line.split()]
			for k, value in enumerate(values):
				current_data[j, k] = value
		time_step_data.append(current_data)
	return time_step_data

def read_d_file(folder_path: str) -> list[np.ndarray]:
	"""读取粒径文件"""
	file_path = os.path.join(folder_path, "d_in_air.dat")
	results = read_file(file_path)
	index_with_variables = find_indices_with_keyword(results, "variables")
	dm_list = []
	d50_list = []
	d90_list = []
	index1 = index_with_variables[0] + 1
	lines = results[index1:]
	for line in lines:
		values = [float(num) for num in line.split()]
		dm_list.append(values[1])
		d50_list.append(values[2])
		d90_list.append(values[3])
	return dm_list, d50_list, d90_list

def define_file_path(OS_name: str) -> str:
	"""定义文件路径"""
	if OS_name == 'l':
		working_dir = "/home/ekalhxh/ripple/coll11"
	elif OS_name == 'w':
		#working_dir = "E:/Data/Sandripples1DFluid/ripple/coll15"
		working_dir = "E:/Data/Q_on_flat_bed"
	else:
		print("Invalid OS name.")
		exit()
	return working_dir

def extract_case_parameters(case_name: str, nu) -> Dict[str, float]:
	"""提取案例参数"""
	parts = case_name.split("_")
	u_star = float(parts[0][6:])/100
	dia_name = parts[1]
	if "and" in dia_name:
		dia_name_list = dia_name.split("and")
		dia1 = float(dia_name_list[0])/1e6
		dia2 = float(dia_name_list[1])/1e6
		dia =  (dia1 + dia2) / 2
		stdd = dia2 - dia1
	elif "stdd" in dia_name:
		dia_name_list = dia_name.split("stdd")
		dia = float(dia_name_list[0])/1e6
		stdd = float(dia_name_list[1])/1e6
		dia1 = dia - 3*stdd
		dia2 = dia + 3*stdd
	elif "log" in dia_name:
		dia_name_list = dia_name.split("log")
		dia = float(dia_name_list[0])/1e6
		stdd = float(dia_name_list[1])/1e6
		dia1 = dia - 3*stdd
		dia2 = dia + 3*stdd
	else:
		dia = float(dia_name)/1e6
		dia1 = dia
		dia2 = dia
		stdd = 0.0
	rho_name = parts[2]
	if rho_name == "0":
		rho = 1.263
	else:
		rho = float(rho_name)
	rho_p = float(parts[3])
	last_time = float(parts[4])
	s = rho_p/rho
	g_hat = 9.8 * (1 - 1/s)
	Sh = rho*u_star**2/(rho_p*g_hat*dia)
	Ga = (s*g_hat*dia**3)**0.5/nu
	parameters = {
		"u_star": u_star,
		"dia": dia,
		"dia1": dia1,
		"dia2": dia2,
		"stdd": stdd,
		"rho": rho,
		"rho_p": rho_p,
		"last_time": last_time,
		"Sh": Sh,
		"Ga": Ga,
		"s": s,
		"g_hat": g_hat
	}
	return parameters

def xy_label(temporal: bool, output: int) -> tuple[str, str]:
	"""定义坐标轴标签"""
	if temporal:
		x_label = "Time"
	else:
		x_label = "Height"
	if output == 0:
		y_label = "u (m/s)"
	elif output == 1:
		y_label = "tau_p (Pa)"
	elif output == 2:
		y_label = "tau_f (Pa)"
	elif output == 3:
		y_label = "F_p (N)"
	elif output == 4:
		y_label = "phi_p"
	elif output == 5:
		y_label = "Sh"
	else:
		print("Invalid output type.")
		exit()
	return x_label, y_label

if __name__ == '__main__':

	#-----------------用户输入区-------------------------
	OS_name = 'w' # 输入操作系统名称，'w'代表Windows，'l'代表Linux
	nu = 1.51e-5 # 运动粘度
	interval = 60 # 源文件输出时间间隔
	start = 120 # 起始时间
	end = 300 # 终止时间
	temporal = True # 是否输出时间序列数据
	semi_log = True # 是否使用半对数坐标
	output = 2 # 输出数据类型，0代表u, 1代表tau_p, 2代表tau_f, 3代表F_p, 4代表phi_p, 5代表Sh
	location = 0 # 输出数据所在高度，0代表床面，-1代表顶面, valid when temporal is True
	kapa = 0.42
	rhof = 1.263
	rhop = 2650
	s = rhop/rhof
	g = 9.8 * (1 - 1/s)
	#---------------------------------------------------

	# 定义文件名字典
	#case_dict = {
		#0: "uStar040_300_0_2650_300",
		#1: "uStar050_300_0_2650_300",
		#2: "uStar060_300_0_2650_300",
		#3: "uStar040_400_0_2650_300",
		#4: "uStar050_400_0_2650_300",
		#5: "uStar060_400_0_2650_300",
		#6: "uStar040_296log366_0_2650_300",
		#7: "uStar050_296log366_0_2650_300",
		#8: "uStar060_296log366_0_2650_300",
		#9: "uStar050_260log537_0_2650_300",
		#10: "uStar055_260log537_0_2650_300",
		#11: "uStar060_260log537_0_2650_300",
		#12: "uStar050_268log491_0_2650_300",
		#13: "uStar055_268log491_0_2650_300",
		#14: "uStar060_268log491_0_2650_300",
		#15: "uStar040_352log522_0_2650_300",
		#16: "uStar050_352log522_0_2650_300",
		#17: "uStar060_352log522_0_2650_300",
		#18: "uStar040_279log368_0_2650_300",
		#19: "uStar050_279log368_0_2650_300",
		#20: "uStar060_279log368_0_2650_300",
	#}
	case_dict = {
		#0: "uStar030_300log50_0_2650_300",
		#1: "uStar040_300log50_0_2650_300",
		#2: "uStar050_300log50_0_2650_300",
		#3: "uStar060_300log50_0_2650_300",
		#16: "uStar035_300log50_0_2650_300",
		#17: "uStar045_300log50_0_2650_300",
		#18: "uStar055_300log50_0_2650_300",
		#19: "uStar065_300log50_0_2650_300",
		#4: "uStar030_300log100_0_2650_300",
		#5: "uStar040_300log100_0_2650_300",
		#6: "uStar050_300log100_0_2650_300",
		#7: "uStar060_300log100_0_2650_300",
		#20: "uStar035_300log100_0_2650_300",
		#21: "uStar045_300log100_0_2650_300",
		#22: "uStar055_300log100_0_2650_300",
		#23: "uStar065_300log100_0_2650_300",
		#8: "uStar030_300log200_0_2650_300",
		#9: "uStar040_300log200_0_2650_300",
		#10: "uStar050_300log200_0_2650_300",
		#11: "uStar060_300log200_0_2650_300",
		#24: "uStar035_300log200_0_2650_300",
		#25: "uStar045_300log200_0_2650_300",
		#26: "uStar055_300log200_0_2650_300",
		#27: "uStar065_300log200_0_2650_300",
		12: "uStar030_300log300_0_2650_300",
		13: "uStar040_300log300_0_2650_300",
		14: "uStar050_300log300_0_2650_300",
		15: "uStar060_300log300_0_2650_300",
		28: "uStar035_300log300_0_2650_300",
		29: "uStar045_300log300_0_2650_300",
		30: "uStar055_300log300_0_2650_300",
		31: "uStar065_300log300_0_2650_300",
		#32: "uStar040_430log100_0_2650_300",
		#33: "uStar050_430log100_0_2650_300",
		#34: "uStar060_430log100_0_2650_300",
		#35: "uStar030_167log100_0_2650_300",
		#36: "uStar040_167log100_0_2650_300",
		#37: "uStar050_167log100_0_2650_300",
		#38: "uStar030_269log100_0_2650_300",
		#39: "uStar040_269log100_0_2650_300",
		#40: "uStar050_269log100_0_2650_300",
		#41: "uStar030_321log100_0_2650_300",
		#42: "uStar040_321log100_0_2650_300",
		#43: "uStar050_321log100_0_2650_300",
		#44: "uStar030_240log50_0_2650_300",
		#45: "uStar035_240log50_0_2650_300",
		#46: "uStar040_240log50_0_2650_300",
		#47: "uStar045_240log50_0_2650_300",
		#48: "uStar050_240log50_0_2650_300",
		#49: "uStar055_240log50_0_2650_300",
		#50: "uStar035_269log100_0_2650_300",
		#51: "uStar045_269log100_0_2650_300",
		#52: "uStar055_269log100_0_2650_300",
		#53: "uStar040_400log50_0_2650_300",
		#54: "uStar050_400log50_0_2650_300",
		#55: "uStar060_400log50_0_2650_300",
		#56: "uStar030_250log25_0_2650_300",
		#57: "uStar040_250log25_0_2650_300",
		#58: "uStar050_250log25_0_2650_300",
		#59: "uStar060_250log25_0_2650_300",
		#60: "uStar030_271log121_0_2650_300",
		#61: "uStar040_271log121_0_2650_300",
		#62: "uStar050_271log121_0_2650_300",
		#63: "uStar060_271log121_0_2650_300",
		#64: "uStar030_317log252_0_2650_300",
		#65: "uStar040_317log252_0_2650_300",
		#66: "uStar050_317log252_0_2650_300",
		#67: "uStar060_317log252_0_2650_300",
		#68: "uStar030_347log537_0_2650_300",
		#69: "uStar040_347log537_0_2650_300",
		#70: "uStar050_347log537_0_2650_300",
		#71: "uStar060_347log537_0_2650_300"
	}

	# 定义文件路径
	working_dir = define_file_path(OS_name)

	# 读取廓线数据
	time_averaged_data = []
	profile_list = []
	dm_case_list = []
	for case in tqdm(case_dict.values(), desc="Processing", unit="case"):
		case_folder = f"{working_dir}/{case}/Field"
		case_data = read_profile_file(case_folder)
		case_parameters = extract_case_parameters(case, nu)
		case_folder = f"{working_dir}/{case}"
		dm_list, d50_list, d90_list = read_d_file(case_folder)
		dm_case_list.append(np.mean(dm_list))
		rho_p = case_parameters['rho_p']
		g_hat = case_parameters['g_hat']
		dia = case_parameters['dia']
		Sh_str = f"{case_parameters['Sh']:.2f}"
		s_str = f"{case_parameters['s']:.2f}"
		Ga_str = f"{case_parameters['Ga']:.2f}"
		stdd_str = f"{case_parameters['stdd']*1e6:.0f}"
		#laber_str = f"Sh={Sh_str}, s={s_str}, Ga={Ga_str}"
		laber_str = f"stdd = {stdd_str} μm"
		if temporal:
			times = [i for i in range(start, end+interval, interval)]
			points = []
			for t in times:
				time_step = t
				time_index = (time_step - start) // interval
				current_data = case_data[time_index]
				current_dm = dm_list[time_index]
				if output < 5:
					point = current_data[location, output + 1]
				elif output == 5:
					tau_f = current_data[location, 3]
					point = tau_f/(rho_p*g_hat*current_dm)
					#point = np.sqrt(tau_f/1.263)
				else:
					print("Invalid output type.")
					exit()
				points.append(point)
			plt.plot(times, points, label=laber_str)
			time_averaged_data.append(np.mean(points))
		else:
			profiles = []
			for t in range(start, end+interval, interval):
				time_step = t
				time_index = (time_step - start) // interval
				current_data = case_data[time_index]
				heights = [current_data[i, 0] for i in range(len(current_data))]
				if output < 5:
					current_profile = [current_data[i, output + 1] for i in range(len(current_data))]
				else:
					print("Invalid output type.")
					exit()
				profiles.append(current_profile)
			average_profile = np.mean(profiles, axis=0)
			profile_list.append(average_profile)
			if semi_log:
				plt.semilogx(heights, average_profile, label=laber_str)
			else:
				plt.plot(heights, average_profile, label=laber_str)
			d = case_parameters['dia']

	# 确定Bagnold node
	cross_x = []
	cross_y = []
	n = len(profile_list)
	for i in range(n):
		for j in range(i+1, n):
			# 两条曲线的纵坐标
			y1 = np.array(profile_list[i])
			y2 = np.array(profile_list[j])
			# 遍历相邻点，找交点
			for k in range(len(heights)-1):
				# 判断是否跨越
				if (y1[k]-y2[k])*(y1[k+1]-y2[k+1]) < 0:
            	    # 线性插值计算交点
					x0, x1 = heights[k], heights[k+1]
					y10, y11 = y1[k], y1[k+1]
					y20, y21 = y2[k], y2[k+1]
					# 交点横坐标
					t = (y20 - y10) / ((y11 - y10) - (y21 - y20))
					x_cross = x0 + t * (x1 - x0)
					y_cross = y10 + t * (y11 - y10)
					cross_x.append(x_cross)
					cross_y.append(y_cross)

	if not temporal:
		case_averaged_dm = np.mean(dm_case_list)
		# 计算平均交点坐标
		z_f = np.mean(cross_x)
		u_f = np.mean(cross_y)
		print("平均交点横坐标：", z_f)
		print("平均交点纵坐标：", u_f)
		H = heights[-1]
		z_0 = d/30
		alpha1 = kapa/np.log(H/z_0)
		alpha2 = kapa/np.log(H/z_f)
		u_d = alpha1*alpha2*u_f/(alpha2 - alpha1)
		S_d = rhof*u_d**2/(rho_p*g*case_averaged_dm)
		print("u_d =", u_d, "m/s")
		print("S_d =", S_d)

	x_str, y_str = xy_label(temporal, output)
	plt.xlabel(x_str)
	plt.ylabel(y_str)
	plt.legend()
	plt.show()
	if temporal:
		time_averaged_data = np.array(time_averaged_data)
		case_averaged_data = np.mean(time_averaged_data)
		case_averaged_stdd = np.std(time_averaged_data)
		print(f"Time averaged data: {time_averaged_data}")
		print(f"Case averaged data: {case_averaged_data}")
		print(f"Case std: {case_averaged_stdd}")