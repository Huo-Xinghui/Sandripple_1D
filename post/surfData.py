# -*-*-*- Author: EkalHxH -*-*-*-
# version: 1.0 (2024-11-5)
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
# 用于处理床面数据，输出处理后的数据到surfProfile.dat, surfDiameter.dat, autoCorr.dat, surfTopology.dat

# cases.dat 用于记录每个case的参数
# surfProfile.dat 用于记录每个时间步的床面廓线数据
# surfDiameter.dat 用于记录每个时间步的床面粒径数据
# autoCorr.dat 用于记录每个时间步的床面自相关数据
# surfTopology.dat 用于记录每个时间步的波长、振幅、波速数据

# 用法: 设置start和end变量来控制处理的时间范围，working_dir变量来控制工作目录，case_dict字典来控制处理的case
# ********************************************************************************************************

# 导入必要的库
import os # 用于文件操作
import re # 用于正则表达式
import numpy as np # 用于数值计算
from dataclasses import dataclass # 用于定义数据类
from typing import List # 用于类型提示
from typing import Dict # 用于类型提示
from scipy.signal import find_peaks # 用于寻找峰值
from tqdm import tqdm # 用于显示进度条

# 定义一个数据类来存储网格节点数据
@dataclass
class grid_node_data:
	x: float
	y: float
	z: float
	d: float

# 定义一个函数来读取任意文件
def read_file(file_path):
	if not os.path.exists(file_path):
		print(f"文件 {file_path} 不存在。")
		exit()

	with open(file_path, 'r', encoding='utf-8') as file:
		line = file.readlines()
		return line

# 读取多个文件内容, 并将其存储到一个字典中
def read_surface_file(folder_path, start_file, end_file):
	time_step_data_dict: Dict[int, List[List[grid_node_data]]] = {}
	for i in range(start_file, end_file+1):
		file_path = os.path.join(folder_path, f"SurfaceData_{i}.plt")
		results = read_file(file_path)

		for line in results:
			if "variables" in line:
				continue
			elif "zone" in line:
				current_time_list = [float(num) for num in re.findall(r'\d+\.?\d*', line)]
				current_time = round(current_time_list[0])
			elif "i=" in line:
				ij_range = [int(num) for num in re.findall(r'\d+', line)]
				i_range = ij_range[0]
				j_range = ij_range[1]
				current_surface = [[None for _ in range(i_range)] for _ in range(j_range)]
				current_i = 0
				current_j = 0
			else:
				columns = line.split()
				current_grid_node = grid_node_data(
					x=float(columns[0]),
					y=float(columns[1]),
					z=float(columns[2]),
					d=float(columns[3])
				)
				current_surface[current_j][current_i] = current_grid_node
				current_i += 1
				if current_i >= i_range:
					current_i = 0
					current_j += 1
					if current_j >= j_range:
						time_step_data_dict[current_time] = current_surface
	return time_step_data_dict

# 主程序
if __name__ == "__main__":
	# 操作系统
	sys_OS = "w" # "w" for windows, "l" for linux
	if sys_OS == "l":
		linux_flag = True
	elif sys_OS == "w":
		linux_flag = False
	else:
		print("Invalid input!")
		exit()

	# 定义常数
	nu = 1.51e-5 # 运动粘度
	interval = 30 # 源文件输出时间间隔
	dt = 60 # 计算波速的时间间隔
	file_interval = 240 # 两个文件之间的时间间隔
	start = 30 # 开始时间
	end = 600 # 结束时间

	# 定义文件路径
	if linux_flag:
		working_dir = "/home/ekalhxh/ripple/coll12"
	else:
		working_dir = "E:/Data/Sandripples1DFluid/ripple/coll12"

	# 删除已存在的cases.dat文件
	cases_file = os.path.join(working_dir, "cases.dat")
	if os.path.exists(cases_file):
		os.remove(cases_file)

	# 定义文件名字典
	case_dict = {
		1: "uStar040_300_0_2650_600",
		2: "uStar045_300_0_2650_600",
		3: "uStar050_300_0_2650_600",
		4: "uStar055_300_0_2650_600",
		5: "uStar060_300_0_2650_600",
		6: "uStar050_200_0_2650_600",
		7: "uStar050_400_0_2650_600",
		8: "uStar050_500_0_2650_600",
		9: "uStar050_300_2_2650_600",
		10: "uStar050_300_0_1000_600",
		11: "uStar050_300_0_4000_600",
	}

	# 遍历字典中的每个case
	for i, folder_name in tqdm(case_dict.items(), desc="Processing", unit="case"):
		# 从文件名中提取参数
		parts = folder_name.split("_")
		ustar_name = parts[0]
		ustar_str = ustar_name[6:]
		ustar_int = int(ustar_str)
		ustar = ustar_int / 100
		dia_name = parts[1]
		if "and" in dia_name:
			dia_name_list = dia_name.split("and")
			dia1 = float(dia_name_list[0])/1e6
			dia2 = float(dia_name_list[1])/1e6
			dia =  (dia1 + dia2) / 2
		elif "stdd" in dia_name:
			dia_name_list = dia_name.split("stdd")
			dia1 = float(dia_name_list[0])/1e6
			dia2 = float(dia_name_list[1])/1e6
			dia =  dia1
		else:
			dia1 = float(dia_name)/1e6
			dia2 = dia1
			dia =  dia1
		rho_name = parts[2]
		if rho_name == "0":
			rho = 1.263
		else:
			rho = float(rho_name)
		rhoP = float(parts[3])

		# 计算s, Sh和Ga
		s = rhoP/rho
		g_hat = 9.8 * (1 - 1/s)
		Sh = rho*ustar**2/(rhoP*g_hat*dia)
		Ga = (s*g_hat*dia**3)**0.5/nu

		# 将各参数写入cases.dat文件
		with open(cases_file, 'a') as file:
			file.write(f"{folder_name}\n")
			file.write(f"u* = {ustar:.2f}, dia1={dia1:.6f}, dia2={dia2:.6f}, dia={dia:.6f}\n")
			file.write(f"rho={rho:.4f}, rhoP={rhoP:.1f}\n")
			file.write(f"s={s:.4f}, Sh={Sh:.4f}, Ga={Ga:.4f}\n")
			file.write("\n")

		# 读取床面数据
		last_time = int(parts[4])
		real_end = min(end, last_time)
		start_file_num = start // file_interval
		end_file_num = real_end // file_interval
		folder_path = f"{working_dir}/{folder_name}/Surface"
		time_step_data = read_surface_file(folder_path, start_file_num, end_file_num)

		# 删除已存在的surfProfile.dat, surfDiameter.dat, autoCorr.dat, surfTopology.dat文件
		surface_profile_file = os.path.join(working_dir, f"{folder_name}/surfProfile.dat")
		surface_diameter_file = os.path.join(working_dir, f"{folder_name}/surfDiameter.dat")
		auto_corr_file = os.path.join(working_dir, f"{folder_name}/autoCorr.dat")
		surface_topology_file = os.path.join(working_dir, f"{folder_name}/surfTopology.dat")
		if os.path.exists(surface_profile_file):
			os.remove(surface_profile_file)
		if os.path.exists(surface_diameter_file):
			os.remove(surface_diameter_file)
		if os.path.exists(auto_corr_file):
			os.remove(auto_corr_file)
		if os.path.exists(surface_topology_file):
			os.remove(surface_topology_file)

		# 初始化
		wavelengths = [] # 波长
		amplitudes = [] # 振幅
		time = [] # 时间
		time_step_profile: Dict[int, List[float]] = {} # 床面高度廓线
		time_step_profile_d: Dict[int, List[float]] = {} # 床面粒径廓线
		lag = 0

		# 遍历每个时间步
		for t, surface in time_step_data.items():
			if t>=start and t<=end:
				time.append(t)
				# 计算床面高度、粒径的平均廓线
				profile_x = np.array([node.x for node in surface[0]])
				surface_z = np.array([[node.z for node in row] for row in surface])
				surface_d = np.array([[node.d for node in row] for row in surface])
				profile_z = np.mean(surface_z, axis=0)
				profile_d = np.mean(surface_d, axis=0)
				average_z = np.mean(profile_z)
				profile_z = profile_z - average_z
				time_step_profile[t] = profile_z
				time_step_profile_d[t] = profile_d

				# 将床面廓线、粒径写入文件
				with open(surface_profile_file, 'a') as file:
					file.write(f"x vs z at {t} with {len(profile_x)} points\n")
					for i in range(len(profile_x)):
						file.write(f"{profile_x[i]} {profile_z[i]}\n")
				with open(surface_diameter_file, 'a') as file:
					file.write(f"x vs d at {t} with {len(profile_x)} points\n")
					for i in range(len(profile_x)):
						file.write(f"{profile_x[i]} {profile_d[i]}\n")

				# 计算床面自相关
				auto_corr_x = profile_x
				auto_corr = np.correlate(profile_z, profile_z, mode='full')
				auto_corr_z = auto_corr[len(auto_corr)//2:]

				# 将床面自相关写入文件
				with open(auto_corr_file, 'a') as file:
					file.write(f"x vs corr at {t} with {len(auto_corr_x)} points\n")
					for i in range(len(auto_corr_x)):
						file.write(f"{auto_corr_x[i]} {auto_corr_z[i]}\n")

				# 计算波长、振幅
				peaks, _ = find_peaks(auto_corr_z, height=0)
				amplitude = 2*np.sqrt(2*auto_corr_z[0])
				amplitudes.append(amplitude)
				if len(peaks) > 0:
					lag = peaks[0]
				displacement = auto_corr_x[lag]
				wavelengths.append(displacement)

		# 计算波速
		velocities = [] # 波速
		di = dt // interval # 计算波速的键值间隔
		lag = 0
		for i in range(len(time)-di):
			t1 = time[i]
			t2 = time[i+di]
			z1 = time_step_profile[t1]
			z2 = time_step_profile[t2]
			cross_corr = np.correlate(z1, z2, mode='full')
			cross_start = len(z1) - 1
			cross_corr_z = cross_corr[cross_start:]
			peaks, _ = find_peaks(cross_corr_z, height=0)
			if len(peaks) > 0:
				lag = peaks[0]
			displacement = profile_x[lag]
			time_diff = t2 - t1
			velocity = displacement / time_diff
			velocities.append(velocity)
		# 确保波速和时间长度相同, 如果不同则用最后一个值填充
		while len(velocities) < len(time):
			velocities.append(velocities[-1])

		# 将波长、振幅、波速写入文件
		with open(surface_topology_file, 'w') as file:
			file.write(f"lambda amplitude velocity vs t with {len(time)} points\n")
			for i in range(len(time)):
				file.write(f"{time[i]} {wavelengths[i]} {amplitudes[i]} {velocities[i]}\n")