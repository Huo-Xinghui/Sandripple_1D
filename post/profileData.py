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

def define_file_path(OS_name: str) -> str:
	"""定义文件路径"""
	if OS_name == 'l':
		working_dir = "/home/ekalhxh/ripple/coll11"
	elif OS_name == 'w':
		working_dir = "E:/Data/Sandripples1DFluid/ripple/coll11"
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
	start = 60 # 起始时间
	end = 600 # 终止时间
	temporal = True # 是否输出时间序列数据
	semi_log = True # 是否使用半对数坐标
	output = 5 # 输出数据类型，0代表u, 1代表tau_p, 2代表tau_f, 3代表F_p, 4代表phi_p, 5代表Sh
	location = 0 # 输出数据所在高度，0代表床面，-1代表顶面, valid when temporal is True
	#---------------------------------------------------

	# 定义文件名字典
	case_dict = {
		0: "uStar035_300_0_2650_300",
		1: "uStar040_300_0_2650_300",
		2: "uStar045_300_0_2650_300",
		3: "uStar050_300_0_2650_300",
		4: "uStar055_300_0_2650_300",
		5: "uStar060_300_0_2650_300",
		6: "uStar065_300_0_2650_300",
		7: "uStar030_200_0_2650_300",
		8: "uStar040_200_0_2650_300",
		9: "uStar050_200_0_2650_300",
		10: "uStar060_200_0_2650_300",
		11: "uStar040_400_0_2650_300",
		12: "uStar050_400_0_2650_300",
		13: "uStar060_400_0_2650_300",
		14: "uStar050_500_0_2650_300",
		15: "uStar050_600_0_2650_300",
		16: "uStar030_300_2_2650_300",
		17: "uStar040_300_2_2650_300",
		18: "uStar050_300_2_2650_300",
		19: "uStar060_300_2_2650_300",
		20: "uStar050_300_3_2650_300",
		21: "uStar050_300_0_500_180",
		22: "uStar050_300_0_1000_300",
		23: "uStar050_300_0_2000_300",
		24: "uStar050_300_0_3000_300",
		25: "uStar050_300_0_4000_300",
	}

	# 定义文件路径
	working_dir = define_file_path(OS_name)

	# 读取廓线数据
	time_averaged_data = []
	for case in tqdm(case_dict.values(), desc="Processing", unit="case"):
		case_folder = f"{working_dir}/{case}/Field"
		case_data = read_profile_file(case_folder)
		case_parameters = extract_case_parameters(case, nu)
		rho_p = case_parameters['rho_p']
		g_hat = case_parameters['g_hat']
		dia = case_parameters['dia']
		Sh_str = f"{case_parameters['Sh']:.2f}"
		s_str = f"{case_parameters['s']:.2f}"
		Ga_str = f"{case_parameters['Ga']:.2f}"
		laber_str = f"Sh={Sh_str}, s={s_str}, Ga={Ga_str}"
		if temporal:
			times = [i for i in range(start, end+interval, interval)]
			points = []
			for t in times:
				time_step = t
				time_index = (time_step - start) // interval
				current_data = case_data[time_index]
				if output < 5:
					point = current_data[location, output + 1]
				elif output == 5:
					tau_f = current_data[location, 3]
					point = tau_f/(rho_p*g_hat*dia)
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
			if semi_log:
				plt.semilogy(heights, average_profile, label=laber_str)
			else:
				plt.plot(heights, average_profile, label=laber_str)
	x_str, y_str = xy_label(temporal, output)
	plt.xlabel(x_str)
	plt.ylabel(y_str)
	plt.legend()
	plt.show()
	time_averaged_data = np.array(time_averaged_data)
	case_averaged_data = np.mean(time_averaged_data)
	print(f"Time averaged data: {time_averaged_data}")
	print(f"Case averaged data: {case_averaged_data}")