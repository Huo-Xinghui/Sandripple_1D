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
"""用于绘制颗粒相关统计量的图像"""
# ********************************************************************************************************

# 导入必要的库
import math # 导入数学库
import os # 导入操作系统库
import re # 导入正则表达式库
import matplotlib.pyplot as plt # 导入绘图库
from dataclasses import dataclass # 导入dataclass库
from dataclasses import field # 导入field库
from typing import List # 导入List类型
from typing import Dict # 导入Dict类型

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
	collNum: int = 0

	def calculate_momentum(self) -> List[float]:
		"""计算动量"""
		return [self.vel[i]*math.pi*self.dia**3/6*self.rho for i in range(len(self.vel))]

	def __post_init__(self):
		"""初始化后计算速度的模"""
		self.vMag = math.sqrt(sum([self.vel[i]**2 for i in range(len(self.vel))]))

@dataclass
class TimeStepData:
	time: int = 0
	count: int = 0
	particles: List[ParticleData]

	def add_particle(self, particle: ParticleData):
		"""添加粒子"""
		self.particles.append(particle)
		self.count += 1

	def get_particle_count(self) -> int:
		"""返回粒子数量"""
		return len(self.particles)

	def calculate_total_momentum(self) -> List[float]:
		"""计算总动量"""
		total_momentum = [0.0, 0.0, 0.0]
		for particle in self.particles:
			particle_momentum = particle.calculate_momentum()
			for i in range(len(total_momentum)):
				total_momentum[i] += particle_momentum[i]
		return total_momentum

@dataclass
class ParticleNumData:
	time: List[int] = field(default_factory=list)
	particleNum: List[int] = field(default_factory=list)
	iterNum: List[int] = field(default_factory=list)

	def add_data(self, time: int, particleNum: int, iterNum: int):
		"""添加数据"""
		self.time.append(time)
		self.particleNum.append(particleNum)
		self.iterNum.append(iterNum)

	def plot_num_vs_time(self, start_time=None, end_time=None):
		"""绘制颗粒数随时间变化的图像"""
		if start_time is not None and end_time is not None:
			time_filtered = [t for t in self.time if start_time <= t <= end_time]
			particleNum_filtered = [self.particleNum[i] for i, t in enumerate(self.time) if start_time <= t <= end_time]
		else:
			time_filtered = self.time
			particleNum_filtered = self.particleNum

		plt.plot(time_filtered, particleNum_filtered)
		plt.xlabel('Time')
		plt.ylabel('Number')
		plt.title('Number vs Time')
		plt.show()

	def plot_num_vs_iter(self, start_iter=None, end_iter=None):
		"""绘制颗粒数随迭代次数变化的图像"""
		if start_iter is not None and end_iter is not None:
			iter_filtered = [i for i in self.iterNum if start_iter <= i <= end_iter]
			particleNum_filtered = [self.particleNum[i] for i, t in enumerate(self.iterNum) if start_iter <= t <= end_iter]
		else:
			iter_filtered = self.iterNum
			particleNum_filtered = self.particleNum

		plt.plot(iter_filtered, particleNum_filtered)
		plt.xlabel('Iteration')
		plt.ylabel('Number')
		plt.title('Number vs Iteration')
		plt.show()

# 定义一个函数来读取任意文件
def read_file(file_path):
	if not os.path.exists(file_path):
		print(f"文件 {file_path} 不存在。")
		return None

	with open(file_path, 'r', encoding='utf-8') as file:
		data = file.readlines()
		return data

# 定义一个函数来读取颗粒数
def read_num_file(folder_path):
	file_path = os.path.join(folder_path, "ParticleNum.plt")
	results = read_file(file_path)
	if results is not None:
		particle_num_data = ParticleNumData()
		for line in results[1:]:
			columns = line.split()
			time = round(float(columns[0]))
			particle_num = int(columns[1])
			iter_num = int(columns[2])
			particle_num_data.add_data(time, particle_num, iter_num)
		return particle_num_data
	return None

def find_indices_with_keyword(data: List[str], keyword: str) -> List[int]:
	return [i for i, line in enumerate(data) if keyword in line]

def extract_times(data: List[str], indices: List[int]) -> List[int]:
	times = []
	for index in indices:
		line = data[index]
		numbers = [float(num) for num in re.findall(r'\d+\.?\d*', line)]
		times.extend(numbers)
	times = [round(num) for num in times]
	return times

def parse_particle_data(lines: List[str]) -> List[ParticleData]:
	particles = []
	for line in lines:
		columns = line.split()
		particle = ParticleData(
			loc=[float(columns[0]), float(columns[1]), float(columns[2])],
			vel=[float(columns[3]), float(columns[4]), float(columns[5])],
			vMag=float(columns[6]),
			h=float(columns[7]),
			hMax=float(columns[8]),
			dia=float(columns[9]),
			survL=float(columns[10]),
			survT=float(columns[11]),
			collNum=int(columns[12])
		)
		particles.append(particle)
	return particles

# 定义一个函数来读取多个颗粒文件内容
def read_data_file(folder_path: str, file_num: int) -> Dict[int, TimeStepData]:
	all_data = []
	for i in range(file_num):
		file_path = os.path.join(folder_path, f"ParticleData_{i}.plt")
		results = read_file(file_path)
		if results:
			all_data.extend(results)

	indices_with_zone = find_indices_with_keyword(all_data, "zone")
	indices_with_variables = find_indices_with_keyword(all_data, "variables")

	extracted_time = extract_times(all_data, indices_with_zone)

	time_step_data: Dict[int, TimeStepData] = {}

	for i in range(len(indices_with_zone)):
		index1 = indices_with_zone[i] + 1
		index2 = indices_with_zone[i+1] - 1 if i < len(indices_with_zone) - 1 else len(all_data) - 1
		if index2 in indices_with_variables:
			index2 -= 1
		particles = parse_particle_data(all_data[index1:index2+1])
		extracted_num = len(particles)
		time_step_data[extracted_time[i]] = TimeStepData(time=extracted_time[i], count=extracted_num, particles=particles)

	return time_step_data

def plot_Q_vs_time(time, flux_x, output_Q_t):
	if output_Q_t:
		plt.plot(time, flux_x)
		plt.xlabel('Time')
		plt.ylabel('Q')
		plt.title('Q vs Time')
		plt.show()
	else:
		return

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

	# 定义常量
	xMax = 0.5 # 计算域x方向长度
	yMax = 0.01 # 计算域y方向长度
	nu = 1.51e-5 # 运动粘度
	interval = 60 # 源文件输出时间间隔
	file_interval = 240 # 两个文件之间的时间间隔
	start = 30 # 开始时间
	end = 3600 # 结束时间
	Sh_d = 0.00298 # 冲击临界
	area = xMax*yMax # 计算域底面积

	# 定义文件路径
	if linux_flag:
		working_dir = "/home/ekalhxh/ripple/coll1"
	else:
		working_dir = "E:/Data/Sandripples1DFluid/ripple/coll1"

	# 定义文件名字典
	case_dict = {
		#0: "uStar030_250_0_2650_3600",
		#1: "uStar035_250_0_2650_3600",
		#2: "uStar040_250_0_2650_3600",
		#3: "uStar045_250_0_2650_3600",
		#4: "uStar050_250_0_2650_3600",
		#5: "uStar055_250_0_2650_3600",
		#6: "uStar060_250_0_2650_3600",
		#7: "uStar065_250_0_2650_3600",
		#8: "uStar045_200_0_2650_3600",
		#9: "uStar045_250_1_2000_3600",
		#10: "uStar045_250_1_3000_3600",
		#11: "uStar045_250_1_4000_3600",
		#12: "uStar045_250_2_2650_3600",
		#13: "uStar045_250_3_2650_3600",
		#14: "uStar045_300_0_2650_3600",
		#15: "uStar045_350_0_2650_3600",
		#16: "uStar045_400_0_2650_3600",
		#17: "uStar040_150and350_0_2650_3600",
		#18: "uStar045_150and350_0_2650_3600",
		#19: "uStar050_150and350_0_2650_3600",
		#20: "uStar055_150and350_0_2650_3600",
		#21: "uStar060_150and350_0_2650_3600",
		#22: "uStar050_150and450_0_2650_3600",
		#23: "uStar050_150and550_0_2650_3600",
		#24: "uStar050_200and400_0_2650_3600",
		#25: "uStar050_250and350_0_2650_3600",
		#26: "uStar050_300stdd5_0_2650_3600",
		#27: "uStar050_300stdd10_0_2650_3600",
		#28: "uStar050_300stdd20_0_2650_3600",
		#29: "uStar050_300stdd50_0_2650_3600",
		30: "uStar035_300_0_2650_3600",
		31: "uStar040_300_0_2650_3600",
		32: "uStar045_300_0_2650_3600",
		33: "uStar050_300_0_2650_3600",
		34: "uStar055_300_0_2650_3600",
		35: "uStar060_300_0_2650_3600",
		36: "uStar065_300_0_2650_3600",
		#37: "uStar035_300stdd100_0_2650_3600",
		#38: "uStar040_300stdd100_0_2650_3600",
		#39: "uStar045_300stdd100_0_2650_3600",
		#40: "uStar050_300stdd100_0_2650_3600",
		#41: "uStar055_300stdd100_0_2650_3600",
		#42: "uStar060_300stdd100_0_2650_3600",
		#43: "uStar065_300stdd100_0_2650_3600",
	}

		z_num = end // 60
		recalculate = input("Recalculate the data (y/n), default=no: ")
		average_data_file = os.path.join(working_dir, "averageData.dat")
		if recalculate.lower() == "y":
			if os.path.exists(average_data_file):
				os.remove(average_data_file)
			if not midair_collision:
				my_dict = {
					0: "uStar030_250_0_2650_3600_nc",
					1: "uStar035_200_0_2650_120_nc",
					2: "uStar035_250_0_2650_3600_nc",
					3: "uStar035_250_1_2000_120_nc",
					4: "uStar035_250_1_3000_120_nc",
					5: "uStar035_250_1_4000_120_nc",
					6: "uStar035_250_2_2650_120_nc",
					7: "uStar035_250_3_2650_120_nc",
					8: "uStar035_300_0_2650_120_nc",
					9: "uStar035_350_0_2650_120_nc",
					10: "uStar035_400_0_2650_120_nc",
					11: "uStar040_250_0_2650_3600_nc",
					12: "uStar045_200_0_2650_340_nc",
					13: "uStar045_250_0_2650_3600_nc",
					14: "uStar045_250_1_2000_120_nc",
					15: "uStar045_250_1_3000_120_nc",
					16: "uStar045_250_1_4000_120_nc",
					17: "uStar045_250_2_2650_337_nc",
					18: "uStar045_250_3_2650_120_nc",
					19: "uStar045_300_0_2650_439_nc",
					20: "uStar045_350_0_2650_120_nc",
					21: "uStar045_400_0_2650_677_nc",
					22: "uStar050_200_1_1000_120_nc",
					23: "uStar050_250_0_2650_3600_nc",
					24: "uStar050_250_1_1000_120_nc",
					25: "uStar050_250_3_2650_120_nc",
					26: "uStar050_300_1_1000_120_nc",
					27: "uStar055_200_0_2650_120_nc",
					28: "uStar055_200_1_1000_120_nc",
					29: "uStar055_250_0_2650_3600_nc",
					30: "uStar055_250_1_1000_120_nc",
					31: "uStar055_250_1_2000_120_nc",
					32: "uStar055_250_1_3000_120_nc",
					33: "uStar055_250_1_4000_120_nc",
					34: "uStar055_250_2_2650_120_nc",
					35: "uStar055_250_3_2650_120_nc",
					36: "uStar055_300_0_2650_120_nc",
					37: "uStar055_300_1_1000_120_nc",
					38: "uStar055_350_0_2650_120_nc",
					39: "uStar055_400_0_2650_120_nc",
					40: "uStar060_200_0_2650_120_nc",
					41: "uStar060_250_0_2650_3600_nc",
					42: "uStar060_250_1_1000_120_nc",
					43: "uStar060_250_1_2000_120_nc",
					44: "uStar060_250_1_3000_120_nc",
					45: "uStar060_250_1_4000_120_nc",
					46: "uStar060_250_2_2650_120_nc",
					47: "uStar060_300_1_1000_120_nc",
					48: "uStar065_250_0_2650_3600_nc",
				}
			elif Bernoulli:
				my_dict = {
					0: "uStar035_150and350_0_2650_120",
					1: "uStar035_150and450_0_2650_120",
					2: "uStar045_150and350_0_2650_3600",
					3: "uStar045_150and450_0_2650_3600",
					4: "uStar055_150and350_0_2650_120",
					5: "uStar035_250and450_0_2650_120",
					6: "uStar035_250and550_0_2650_120",
					7: "uStar045_250and450_0_2650_3600",
					8: "uStar045_250and550_0_2650_3600",
					9: "uStar055_250and450_0_2650_120",
					10: "uStar055_250and550_0_2650_120",
					11: "uStar035_350and550_0_2650_120",
					12: "uStar045_350and550_0_2650_3600",
					13: "uStar055_350and550_0_2650_120",
				}
			else:
				my_dict = {
					0: "uStar030_250_0_2650_3600",
					1: "uStar035_200_0_2650_120",
					2: "uStar035_250_0_2650_3600",
					3: "uStar035_250_1_2000_120",
					4: "uStar035_250_1_3000_120",
					5: "uStar035_250_1_4000_120",
					6: "uStar035_250_2_2650_120",
					7: "uStar035_250_3_2650_120",
					8: "uStar035_300_0_2650_120",
					9: "uStar035_350_0_2650_120",
					10: "uStar035_400_0_2650_120",
					11: "uStar040_250_0_2650_3600",
					12: "uStar045_200_0_2650_3600",
					13: "uStar045_250_0_2650_3600",
					14: "uStar045_250_1_2000_3600",
					15: "uStar045_250_1_3000_3600",
					16: "uStar045_250_1_4000_3600",
					17: "uStar045_250_2_2650_3600",
					18: "uStar045_250_3_2650_3600",
					19: "uStar045_300_0_2650_3600",
					20: "uStar045_350_0_2650_3600",
					21: "uStar045_400_0_2650_3600",
					22: "uStar050_200_1_1000_120",
					23: "uStar050_250_0_2650_3600",
					24: "uStar050_250_1_1000_120",
					25: "uStar050_250_3_2650_120",
					26: "uStar050_300_1_1000_120",
					27: "uStar055_200_0_2650_120",
					28: "uStar055_200_1_1000_120",
					29: "uStar055_250_0_2650_3600",
					30: "uStar055_250_1_1000_120",
					31: "uStar055_250_1_2000_120",
					32: "uStar055_250_1_3000_120",
					33: "uStar055_250_1_4000_120",
					34: "uStar055_250_2_2650_120",
					35: "uStar055_250_3_2650_120",
					36: "uStar055_300_0_2650_120",
					37: "uStar055_300_1_1000_120",
					38: "uStar055_350_0_2650_120",
					39: "uStar055_400_0_2650_120",
					40: "uStar060_200_0_2650_120",
					41: "uStar060_250_0_2650_3600",
					42: "uStar060_250_1_1000_120",
					43: "uStar060_250_1_2000_120",
					44: "uStar060_250_1_3000_120",
					45: "uStar060_250_1_4000_120",
					46: "uStar060_250_2_2650_120",
					47: "uStar060_300_1_1000_120",
					48: "uStar065_250_0_2650_3600",
				}
			uStar_list = []
			dia1_list = []
			dia2_list = []
			dia_ave_list = []
			rho_list = []
			rhoP_list = []
			s_list = []
			Q_list = []
			Q_star_list = []
			M_list = []
			M_star_list = []
			Sh_list = []
			Ga_list = []
			d_list = []
			for i, value in my_dict.items():
				parts = value.split("_")
				uStar_name = parts[0]
				uStar_str = uStar_name[6:]
				uStar_int = int(uStar_str)
				uStar = uStar_int / 100
				uStar_list.append(uStar)
				dia_name = parts[1]
				if "and" in dia_name:
					dia1_name, dia2_name = dia_name.split("and")
					dia1 = float(dia1_name)/1e6
					dia2 = float(dia2_name)/1e6
				else:
					dia1 = float(dia_name)/1e6
					dia2 = dia1
				dia_ave = (dia1 + dia2) / 2
				dia1_list.append(dia1)
				dia2_list.append(dia2)
				dia_ave_list.append(dia_ave)
				rho_name = parts[2]
				if rho_name == "0":
					rho = 1.263
				else:
					rho = float(rho_name)
				rho_list.append(rho)
				rhoP_name = parts[3]
				rhoP = float(rhoP_name)
				rhoP_list.append(rhoP)
				s = rhoP/rho
				s_list.append(s)
				g_hat = 9.8 * (1 - 1/s)
				Sh = rho*uStar**2/(rhoP*g_hat*dia_ave)
				Sh_list.append(Sh)
				Ga = (s*g_hat*dia_ave**3)**0.5/nu
				Ga_list.append(Ga)
				last_time = int(parts[4])
				zone_num = last_time // 60
				file_num = zone_num // 4 + 1
				folder_name = value
				folder_path = f"{working_dir}/{folder_name}/Particle"
				results = read_data_file(folder_path, file_num)
				time_step_data = results
				flux_x_list = []
				flux_x_star_list = []
				mass_list = []
				mass_star_list = []
				d_ave_list = []
				for current_time in range(start, end+1, 60):
					if z_num <= zone_num:
						current_data = time_step_data[current_time]
						d_ave = sum([particle.dia for particle in current_data.particles]) / current_data.number
						m_sum = sum([math.pi*particle.dia**3/6*rhoP for particle in current_data.particles])
						x_momentum_sum = sum([particle.vel[0]*(math.pi*particle.dia**3)/6*rhoP for particle in current_data.particles])
						mass = m_sum / area
						mass_star = mass / (rhoP * dia_ave)
						flux_x = x_momentum_sum / area
						flux_x_star = flux_x / (rhoP * dia_ave * (s*g_hat*dia_ave)**0.5)
					else:
						d_ave = 0
						mass = 0
						mass_star = 0
						flux_x = 0
						flux_x_star = 0
					mass_list.append(mass)
					mass_star_list.append(mass_star)
					flux_x_list.append(flux_x)
					flux_x_star_list.append(flux_x_star)
					d_ave_list.append(d_ave)
				ave_M = sum(mass_list) / len(mass_list)
				ave_M_star = sum(mass_star_list) / len(mass_star_list)
				ave_Q = sum(flux_x_list) / len(flux_x_list)
				ave_Q_star = sum(flux_x_star_list) / len(flux_x_star_list)
				ave_d = sum(d_ave_list) / len(d_ave_list)
				Q_list.append(ave_Q)
				Q_star_list.append(ave_Q_star)
				M_list.append(ave_M)
				M_star_list.append(ave_M_star)
				d_list.append(ave_d)
			with open(average_data_file, 'w') as file:
				file.write("uStar dia1 dia2 d_ave rho rhoP Sh Ga s Q Q* M M*\n")
				for i in range(len(uStar_list)):
					file.write(f"{uStar_list[i]} {dia1_list[i]} {dia2_list[i]} {d_list[i]} {rho_list[i]} {rhoP_list[i]} \
{Sh_list[i]} {Ga_list[i]} {s_list[i]} {Q_list[i]} {Q_star_list[i]} {M_list[i]} {M_star_list[i]}\n")
		else:
			if os.path.exists(average_data_file):
				with open(average_data_file, 'r') as file:
					lines = file.readlines()
				uStar_list = []
				dia1_list = []
				dia2_list = []
				d_list = []
				rho_list = []
				rhoP_list = []
				Sh_list = []
				Ga_list = []
				s_list = []
				Q_list = []
				Q_star_list = []
				M_list = []
				M_star_list = []
				for line in lines[1:]:
					columns = line.split()
					uStar_list.append(float(columns[0]))
					dia1_list.append(float(columns[1]))
					dia2_list.append(float(columns[2]))
					d_list.append(float(columns[3]))
					rho_list.append(float(columns[4]))
					rhoP_list.append(float(columns[5]))
					Sh_list.append(float(columns[6]))
					Ga_list.append(float(columns[7]))
					s_list.append(float(columns[8]))
					Q_list.append(float(columns[9]))
					Q_star_list.append(float(columns[10]))
					M_list.append(float(columns[11]))
					M_star_list.append(float(columns[12]))
			else:
				print(f"文件 {average_data_file} 不存在。")
				exit()
		plt.plot(Sh_list, Q_star_list, 'ro')
		plt.xlabel('Sheilds number')
		plt.ylabel('Q')
		plt.title(f"Sheilds number vs Q* during {start}s to {end}s")
		plt.show()
	else:
		file_num = int(input("Please input the number of files, default=16: ") or 16)
		uStar = float(input("Please input u star (m/s): "))
		dia1 = float(input("Please input particle diameter (um), default=250: ") or 250)
		rho = float(input("Please input fluid density (kg/m^3), default=1.263: ") or 1.263)
		rhoP = float(input("Please input particle density (kg/m^3), default=2650: ") or 2650)
		specified_time = int(input("Please input specified time (s), default=60: ") or 60)
		Bernoulli = False
		Bernoulli_input = input("Bernoulli distribution (y/n), default=no: ")
		if Bernoulli_input.lower() == "y":
			Bernoulli = True
			dia2 = float(input("Please input the second particle diameter (um): "))
			if dia2 <= dia1:
				print("The second particle diameter should be larger than the first one.")
				exit()
		my_uStar_dict = {
        	0.3: "uStar030",
        	0.35: "uStar035",
        	0.4: "uStar040",
        	0.45: "uStar045",
        	0.5: "uStar050",
        	0.55: "uStar055",
        	0.6: "uStar060",
        	0.65: "uStar065",
        	0.7: "uStar070"
    	}
		uStar_name = my_uStar_dict[uStar]
		if rho == 1.263:
			rho_name = 0
		else:
			rho_name = int(rho)
		if Bernoulli:
			dia_name = f"{int(dia1)}and{int(dia2)}"
		else:
			dia_name = int(dia1)
		rhoP_name = int(rhoP)
		folder_name = f"{uStar_name}_{dia_name}_{rho_name}_{rhoP_name}"
		folder_path = f"{working_dir}/{folder_name}/Particle"
		dia1 = dia1 * 1e-6
		results1 = read_num_file(folder_path)
		results2 = read_data_file(folder_path, file_num)
		time, particleNum, iterNum = results1
		time_step_data = results2
		flux_x_list = []
		time_list = []
		for current_time in range(60, 3601, 60):
			current_data = time_step_data[current_time]
			x_momentum_sum = sum([particle.vel[0]*(math.pi*particle.dia**3)/6*rhoP for particle in current_data.particles])
			flux_x = x_momentum_sum / area
			flux_x_list.append(flux_x)
			current_time = current_data.time
			time_list.append(current_time)
		plot_Q_vs_time(time_list, flux_x_list, plot_Q_t)
		plot_Num_vs_time(time, particleNum, plot_N_t)