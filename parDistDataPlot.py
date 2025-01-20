# -*-*-*- Author: EkalHxH -*-*-*-
# version: 1.0 (2024-12-12)
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
"""用于统计颗粒分布"""
# ********************************************************************************************************

# 导入必要的库
import math # 导入数学库
import os # 导入操作系统库
import re # 导入正则表达式库
import matplotlib.pyplot as plt # 导入绘图库
import numpy as np # 导入numpy库
from dataclasses import dataclass # 导入dataclass库
from dataclasses import field # 导入field库
from typing import List # 导入List类型
from typing import Dict # 导入Dict类型
from tqdm import tqdm # 导入进度条库

@dataclass
class ParticleData:
	loc: List[float] = field(default_factory=lambda: [0.0, 0.0, 0.0])
	vel: List[float] = field(default_factory=lambda: [0.0, 0.0, 0.0])
	vMag: float = 0.0
	h: float = 0.0
	hMax: float = 0.0
	dia: float = 0.0
	rho: float = 0.0
	survL: float = 0.0
	survT: float = 0.0
	kEng: float = 0.0
	collNum: int = 0

	def calculate_mass(self) -> float:
		"""计算质量"""
		return math.pi*self.dia**3/6*self.rho

	def calculate_momentum(self) -> List[float]:
		"""计算动量"""
		mass = self.calculate_mass()
		return [mass*vel for vel in self.vel]

	def calculate_energy(self) -> float:
		"""计算能量"""
		mass = self.calculate_mass()
		return 0.5*mass*self.vMag**2 + mass*9.8*self.h

@dataclass
class TimeStepData:
	particles: List[ParticleData] = field(default_factory=list)
	time: int = 0
	count: int = 0

	def parse_particle_data(self, lines: List[str], time: int, coll_flag: str):
		for line in lines:
			columns = line.split()
			try:
				x = float(columns[0])
			except ValueError:
				x = 0.0
			try:
				y = float(columns[1])
			except ValueError:
				y = 0.0
			try:
				z = float(columns[2])
			except ValueError:
				z = 0.0
			loc=[x, y, z]
			try:
				u = float(columns[3])
			except ValueError:
				u = 0.0
			try:
				v = float(columns[4])
			except ValueError:
				v = 0.0
			try:
				w = float(columns[5])
			except ValueError:
				w = 0.0
			vel=[u, v, w]
			try:
				vMag = float(columns[6])
			except ValueError:
				vMag = 0.0
			try:
				h = float(columns[7])
			except ValueError:
				h = 0.0
			try:
				hMax = float(columns[8])
			except ValueError:
				hMax = 0.0
			try:
				dia = float(columns[9])
			except ValueError:
				dia = 0.0
			try:
				survL = float(columns[10])
			except ValueError:
				survL = 0.0
			try:
				survT = float(columns[11])
			except ValueError:
				survT = 0.0
			try:
				collNum = int(float(columns[12]))
			except ValueError:
				collNum = 0

			if h > 0.0003*5 or w > 0 or survL <0.0:
			#if h < 0.0 or vMag <= 0.0:
				continue
			if coll_flag == 'c' and collNum == 0:
				continue
			elif coll_flag == 'nc' and collNum != 0:
				continue

			particle = ParticleData(loc=loc, vel=vel, vMag=vMag, h=h, hMax=hMax, dia=dia, survL=survL, survT=survT, collNum=collNum)
			self.particles.append(particle)
			self.count += 1
		self.time = time

# 定义一个函数来读取任意文件
def read_file(file_path):
	if not os.path.exists(file_path):
		print(f"文件 {file_path} 不存在。")
		return None

	with open(file_path, 'r', encoding='utf-8') as file:
		data = file.readlines()
		return data

def find_indices_with_keyword(lines: List[str], keyword: str) -> List[int]:
	return [i for i, line in enumerate(lines) if keyword in line]

def extract_times(lines: List[str], indices: List[int]) -> List[int]:
	times = []
	for index in indices:
		line = lines[index]
		time = [float(num) for num in re.findall(r'\d+\.?\d*', line)]
		times.extend(time)
	times = [round(time) for time in times]
	return times

# 定义一个函数来读取多个颗粒文件内容
def read_particle_file(folder_path: str, end_num: int, coll_flag: str) -> List[TimeStepData]:
	all_data = []
	for i in range(end_num):
		file_path = os.path.join(folder_path, f"ParticleData_{i}.plt")
		results = read_file(file_path)
		if results:
			all_data.extend(results)

	indices_with_zone = find_indices_with_keyword(all_data, "zone")
	indices_with_variables = find_indices_with_keyword(all_data, "variables")

	extracted_times = extract_times(all_data, indices_with_zone)

	time_step_data = []
	for i in range(len(indices_with_zone)):
		time = extracted_times[i]
		index1 = indices_with_zone[i] + 1
		index2 = indices_with_zone[i+1] - 1 if i < len(indices_with_zone) - 1 else len(all_data) - 1
		if index2 in indices_with_variables:
			index2 -= 1
		current_step_data = TimeStepData()
		current_step_data.parse_particle_data(all_data[index1:index2], time, coll_flag)
		time_step_data.append(current_step_data)
	return time_step_data

# 定义一个函数来将start_time至end_time之间的颗粒数据存入新变量
def gather_particle_data(time_step_data: List[TimeStepData], start_time: int, end_time: int) -> List[ParticleData]:
	particle_data = []
	for current_data in time_step_data:
		if current_data.time >= start_time and current_data.time <= end_time:
			particle_data.extend(current_data.particles)
	return particle_data

# 统计List[ParticleData]中某个物理量的分布，bin_num为直方图的柱数
def calculate_distribution(particle_data: List[ParticleData], key: str, bin_num: int) -> Dict[str, List[float]]:
	"""key: 'h', 'hMax', 'dia', 'vMag', 'survL', 'survT', 'kEng', 'collNum'"""
	data = [getattr(particle, key) for particle in particle_data]
	min_data = min(data)
	max_data = max(data)
	bin_width = (max_data - min_data) / bin_num
	bins = [min_data + bin_width*i for i in range(bin_num+1)]
	bin_mid = [(bins[i] + bins[i+1]) / 2 for i in range(bin_num)]
	counts = [0 for _ in range(bin_num)]
	for value in data:
		for i in range(bin_num):
			if bins[i] <= value < bins[i+1]:
				counts[i] += 1
				break
	total_count = sum(counts)
	pdf = [count / total_count / bin_width for count in counts]
	distribution = {key: bin_mid, "n": counts, "pdf": pdf}
	return distribution

def write_file_head(file_path, dia_name, uStar, dia1, dia2, rho, rhoP, s, g_hat, Sh, Ga, variables):
	with open(file_path, 'w') as file:
		if "and" in dia_name:
			file.write(f"""uStar = {uStar} dia1 = {dia1} dia2 = {dia2} rho = {rho} rhoP = {rhoP}\n
s = {s} g_hat = {g_hat} Sh = {Sh} Ga = {Ga}\n
variables = {variables}\n""")
		elif "stdd" in dia_name:
			file.write(f"""uStar = {uStar} dia = {dia1} stdd = {dia2} rho = {rho} rhoP = {rhoP}\n
s = {s} g_hat = {g_hat} Sh = {Sh} Ga = {Ga}\n
variables = {variables}\n""")
		else:
			file.write(f"""uStar = {uStar} dia = {dia1} rho = {rho} rhoP = {rhoP}\n
s = {s} g_hat = {g_hat} Sh = {Sh} Ga = {Ga}\n
variables = {variables}\n""")

def write_file_data(file_path, data):
	with open(file_path, 'a') as file:
		file.write(data)

# 主程序
if __name__ == "__main__":
	# 操作系统
	sys_OS = "w" # "w" for windows, "l" for linux
	if sys_OS == "l":
		linux_flag = True
	elif sys_OS == "w":
		linux_flag = False
	else:
		print("Invalid OS.")
		exit()

	# 定义常量
	nu = 1.51e-5 # 运动粘度
	interval = 60 # 源文件输出时间间隔
	file_interval = 240 # 两个文件之间的时间间隔
	start_time = 60 # 起始时间
	end_time = 300 # 终止时间
	bin_num = 50 # 直方图的柱数
	output_flag = 'survL' # 'h', 'hMax', 'dia', 'vMag', 'survL', 'survT', 'kEng', 'collNum'
	plot_type_flag = 'box' # 'dist' for distribution, 'box' for box plot
	y_label_flag = 'pdf' # 'n' for number, 'pdf' for probability density
	if start_time % interval != 0 or end_time % interval != 0:
		print("Invalid start/end time.")
		exit()

	# 定义文件路径
	if linux_flag:
		working_dir = "/home/ekalhxh/ripple/coll1"
	else:
		working_dir = "E:/Data/Sandripples1DFluid/ripple/coll"

	# 定义文件名字典, key为case编号，value为文件夹名字和碰撞信息的标志, 'c' for collision, 'nc' for no collision
	case_dict = {
		0: ["uStar035_300_0_2650_3600", 'nc'],
		1: ["uStar040_300_0_2650_3600", 'nc'],
		2: ["uStar045_300_0_2650_3600", 'nc'],
		3: ["uStar050_300_0_2650_3600", 'nc'],
		4: ["uStar055_300_0_2650_3600", 'nc'],
		5: ["uStar060_300_0_2650_3600", 'nc'],
		6: ["uStar065_300_0_2650_3600", 'nc']
	}

	# 定义横坐标字典
	x_label_dict = {
		'h': 'hop height (m)',
		'hMax': 'max hop height (m)',
		'dia': 'diameter (m)',
		'vMag': 'velocity magnitude (m/s)',
		'survL': 'survival length (m)',
		'survT': 'survival time (s)',
		'kEng': 'kinetic energy (J)',
		'collNum': 'collision number'
	}

	plt.figure(figsize=(10, 6))
	x_plot = []
	y_plot = []
	for case_value in tqdm(case_dict.values(), desc="Processing", unit="case"):
		# 从文件夹名字中提取参数
		folder_name = case_value[0]
		coll_flag = case_value[1]
		parts = folder_name.split("_")
		uStar = int(parts[0][6:])/100
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

		# 计算s，g_hat，Sh，Ga
		s = rhoP/rho
		g_hat = 9.8 * (1 - 1/s)
		Sh = rho*uStar**2/(rhoP*g_hat*dia)
		Ga = (s*g_hat*dia**3)**0.5/nu

		# 读取颗粒数据
		last_time = int(parts[4])
		end_file_num = last_time // file_interval
		folder_path = f"{working_dir}/{folder_name}/Particle"
		time_step_data = read_particle_file(folder_path, end_file_num, coll_flag)

		for current_data in time_step_data:
			for particle in current_data.particles:
				particle.rho = rhoP
				particle.kEng = particle.calculate_energy()

		# 提取时间段内的颗粒数据
		analyzed_particles = gather_particle_data(time_step_data, start_time, end_time)
		if plot_type_flag == 'dist':
			distribution_dict = calculate_distribution(analyzed_particles, output_flag, bin_num)

			# 文件路径
			hopheight_dist_file = os.path.join(working_dir, f"{folder_name}/hopheight_dist.dat")
			max_hopheight_dist_file = os.path.join(working_dir, f"{folder_name}/max_hopheight_dist.dat")
			d_in_air_dist_file = os.path.join(working_dir, f"{folder_name}/d_in_air_dist.dat")
			vel_dist_file = os.path.join(working_dir, f"{folder_name}/vel_dist.dat")
			hoplength_dist_file = os.path.join(working_dir, f"{folder_name}/hoplength_dist.dat")
			hoptime_dist_file = os.path.join(working_dir, f"{folder_name}/hoptime_dist.dat")
			k_energy_dist_file = os.path.join(working_dir, f"{folder_name}/k_energy_dist.dat")
			coll_num_dist_file = os.path.join(working_dir, f"{folder_name}/coll_num_dist.dat")

			file_dict = {
				'h': hopheight_dist_file,
				'hMax': max_hopheight_dist_file,
				'dia': d_in_air_dist_file,
				'vMag': vel_dist_file,
				'survL': hoplength_dist_file,
				'survT': hoptime_dist_file,
				'kEng': k_energy_dist_file,
				'collNum': coll_num_dist_file
			}

			# 删除已有文件, 并写入文件头和数据
			if os.path.exists(file_dict[output_flag]):
				os.remove(file_dict[output_flag])
			write_file_head(file_dict[output_flag], dia_name, uStar, dia1, dia2, rho, rhoP, s, g_hat, Sh, Ga, "bin pdf")
			for i in range(bin_num):
				write_file_data(file_dict[output_flag], f"{distribution_dict[output_flag][i]} {distribution_dict[y_label_flag][i]}\n")

			x_plot = distribution_dict[output_flag]
			y_plot = distribution_dict[y_label_flag]
			plt.loglog(x_plot, y_plot, 'o', label=f"$u^*={uStar}$")
		elif plot_type_flag == 'box':
			analyzed_data = [getattr(particle, output_flag) for particle in analyzed_particles]
			#analyzed_data = [math.log(data) for data in analyzed_data]
			analyzed_data = np.array([data**0.1 for data in analyzed_data])
			plt.boxplot(analyzed_data, positions=[uStar], widths=0.03, showfliers=True)
			x_plot.append(uStar)
			y_plot.append(np.mean(analyzed_data))
		else:
			print("Invalid plot type.")
			exit()
	if plot_type_flag == 'dist':
		plt.legend()
		plt.xlabel(x_label_dict[output_flag])
		plt.ylabel(y_label_flag)
	elif plot_type_flag == 'box':
		plt.xlabel("$u^*$")
		plt.ylabel(x_label_dict[output_flag])
		plt.xlim(0.3, 0.7)
		plt.plot(x_plot, y_plot, 'ro-')
	plt.show()
