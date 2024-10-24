# 导入必要的库
import os
import re
import matplotlib.pyplot as plt
import numpy as np
from dataclasses import dataclass
from typing import List
from typing import Dict
from scipy.signal import find_peaks
from tqdm import tqdm

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
	all_lines = []
	time_step_data_dict: Dict[int, List[List[grid_node_data]]] = {}
	for i in range(start_file, end_file+1):
		file_path = os.path.join(folder_path, f"SurfaceData_{i}.plt")
		results = read_file(file_path)
		all_lines.extend(results)

		for line in all_lines:
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

def auto_correlation(data1):
	mean_data1 = np.mean(data1)
	data2 = data1
	corr_list = []
	for _ in range(len(data1)):
		mean_data12 = np.mean([data1[i] * data2[i] for i in range(len(data1))])
		corr = mean_data12 - mean_data1 * mean_data1
		corr_list.append(corr)
		data2_head = data2[-1]
		data2 = np.insert(data2[:-1], 0, data2_head)
	return corr_list

# 主程序
if __name__ == "__main__":
	# 定义常数
	nu = 1.51e-5
	interval = 30
	dt = 60 # 计算波速的时间间隔
	file_interval = 240
	user_input = False
	if user_input:
		# 用户输入
		start = int(input("Please input start time (s), default=30: ") or 60)
		end = int(input("Please input end time (s), default=1200: ") or 120)
	else:
		start = 30
		end = 3600
	# 定义文件路径
	working_dir = "/home/ekalhxh/ripple/coll"
	cases_file = os.path.join(working_dir, "cases.dat")
	if os.path.exists(cases_file):
		os.remove(cases_file)
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
		30: "uStar050_300stdd100_0_2650_2774",
	}
	for i, folder_name in tqdm(case_dict.items(), desc="Processing", unit="case"):
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
		s = rhoP/rho
		g_hat = 9.8 * (1 - 1/s)
		Sh = rho*ustar**2/(rhoP*g_hat*dia)
		Ga = (s*g_hat*dia**3)**0.5/nu
		with open(cases_file, 'a') as file:
			file.write(f"{folder_name}\n")
			file.write(f"u* = {ustar:.2f}, dia1={dia1:.4f}, dia2={dia2:.4f}, dia={dia:.4f}\n")
			file.write(f"rho={rho:.4f}, rhoP={rhoP:.1f}\n")
			file.write(f"s={s:.4f}, Sh={Sh:.4f}, Ga={Ga:.4f}\n")
			file.write("\n")

		last_time = int(parts[4])
		real_end = min(end, last_time)
		start_file_num = start // file_interval
		end_file_num = real_end // file_interval
		folder_path = f"{working_dir}/{folder_name}/Surface"
		time_step_data = read_surface_file(folder_path, start_file_num, end_file_num)
		proflie_t = []
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
		wavelength_t = []
		amplitude_t = []
		migration_t = []
		time = []
		time_step_profile: Dict[int, List[float]] = {}
		time_step_profile_d: Dict[int, List[float]] = {}
		lag = 0
		for t, surface in time_step_data.items():
			if t>=start and t<=end:
				profile_x = np.array([node.x for node in surface[0]])
				surface_z = np.array([[node.z for node in row] for row in surface])
				surface_d = np.array([[node.d for node in row] for row in surface])
				profile_z = np.mean(surface_z, axis=0)
				profile_d = np.mean(surface_d, axis=0)
				average_z = np.mean(profile_z)
				profile_z = profile_z - average_z
				time_step_profile[t] = profile_z
				time_step_profile_d[t] = profile_d
				with open(surface_profile_file, 'a') as file:
					file.write(f"x vs z at {t} with {len(profile_x)} points\n")
					for i in range(len(profile_x)):
						file.write(f"{profile_x[i]} {profile_z[i]}\n")
				with open(surface_diameter_file, 'a') as file:
					file.write(f"x vs d at {t} with {len(profile_x)} points\n")
					for i in range(len(profile_x)):
						file.write(f"{profile_x[i]} {profile_d[i]}\n")
				auto_corr_x = profile_x
				auto_corr = np.correlate(profile_z, profile_z, mode='full')
				auto_corr_z = auto_corr[len(auto_corr)//2:]
				#auto_corr_z = auto_correlation(profile_z)
				with open(auto_corr_file, 'a') as file:
					file.write(f"x vs corr at {t} with {len(auto_corr_x)} points\n")
					for i in range(len(auto_corr_x)):
						file.write(f"{auto_corr_x[i]} {auto_corr_z[i]}\n")
				peaks, _ = find_peaks(auto_corr_z, height=0)
				if len(peaks) > 0:
					lag = peaks[0]
				displacement = auto_corr_x[lag]
				amplitude = 2*np.sqrt(2*auto_corr_z[0])
				time.append(t)
				wavelength_t.append(displacement)
				amplitude_t.append(amplitude)
		times = list(time_step_profile.keys())
		velocities = []
		di = dt // interval
		lag = 0
		for i in range(len(times)-di):
			t1 = times[i]
			t2 = times[i+di]
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
		# Ensure velocities has the same length as times by appending the last element
		while len(velocities) < len(times):
			velocities.append(velocities[-1])
		with open(surface_topology_file, 'w') as file:
			file.write(f"lambda amplitude velocity vs t with {len(time)} points\n")
			for i in range(len(time)):
				file.write(f"{time[i]} {wavelength_t[i]} {amplitude_t[i]} {velocities[i]}\n")