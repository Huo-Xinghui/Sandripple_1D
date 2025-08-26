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
"""用于处理颗粒数据，计算颗粒的通量、携带量、能量等信息"""
# ********************************************************************************************************

# 导入必要的库
import math # 导入数学库
import os # 导入操作系统库
import re # 导入正则表达式库
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
	particles: List[ParticleData]
	time: int = 0
	count: int = 0

	def calculate_total_energy(self) -> Dict[str, float]:
		"""计算总能量"""
		total_energy = 0.0
		coll_energy = 0.0
		non_coll_energy = 0.0

		for particle in self.particles:
			energy = particle.calculate_energy()
			total_energy += energy
			if particle.collNum > 0:
				coll_energy += energy
			else:
				non_coll_energy += energy

		return {
			"total": total_energy,
			"coll": coll_energy,
			"non_coll": non_coll_energy
		}

	def calculate_dia_in_air(self) -> Dict[str, float]:
		"""计算空中颗粒平均直径"""
		dia_sum = 0.0
		coll_dia_sum = 0.0
		non_coll_dia_sum = 0.0
		coll_count = 0
		non_coll_count = 0

		for particle in self.particles:
			dia_sum += particle.dia
			if particle.collNum > 0:
				coll_dia_sum += particle.dia
				coll_count += 1
			else:
				non_coll_dia_sum += particle.dia
				non_coll_count += 1
		return {
			"total": dia_sum / self.count,
			"coll": coll_dia_sum / coll_count,
			"non_coll": non_coll_dia_sum / non_coll_count
		}

	def calculate_mass_flux(self) -> Dict[str, List[float]]:
		"""计算质量通量"""
		total_momentum_sum = [0.0, 0.0, 0.0]
		coll_momentum_sum = [0.0, 0.0, 0.0]
		non_coll_momentum_sum = [0.0, 0.0, 0.0]

		for particle in self.particles:
			momentum = particle.calculate_momentum()
			for i in range(3):
				total_momentum_sum[i] += momentum[i]
				if particle.collNum > 0:
					coll_momentum_sum[i] += momentum[i]
				else:
					non_coll_momentum_sum[i] += momentum[i]

		return {
			"total": [m / area for m in total_momentum_sum],
			"coll": [m / area for m in coll_momentum_sum],
			"non_coll": [m / area for m in non_coll_momentum_sum]
		}

	def calculete_carrying_capacity(self) -> Dict[str, float]:
		"""计算颗粒携带量"""
		total_mass_sum = 0.0
		coll_mass_sum = 0.0
		non_coll_mass_sum = 0.0

		for particle in self.particles:
			mass = particle.calculate_mass()
			total_mass_sum += mass
			if particle.collNum > 0:
				coll_mass_sum += mass
			else:
				non_coll_mass_sum += mass
		return {
			"total": total_mass_sum / area,
			"coll": coll_mass_sum / area,
			"non_coll": non_coll_mass_sum / area
		}

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

# 定义一个函数来读取任意文件
def read_file(file_path) -> List[str]:
	"""读取文件"""
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

		particle = ParticleData(loc=loc, vel=vel, vMag=vMag, h=h, hMax=hMax, dia=dia, survL=survL, survT=survT, collNum=collNum)
		particles.append(particle)
	return particles

# 定义一个函数来读取多个颗粒文件内容
def read_particle_file(folder_path: str, end_num: int) -> List[TimeStepData]:
	all_data = []
	for i in range(end_num):
		file_path = os.path.join(folder_path, f"ParticleData_{i}.plt")
		results = read_file(file_path)
		if results:
			all_data.extend(results)

	indices_with_zone = find_indices_with_keyword(all_data, "zone")
	indices_with_variables = find_indices_with_keyword(all_data, "variables")

	extracted_time = extract_times(all_data, indices_with_zone)

	time_step_data = []
	for i in range(len(indices_with_zone)):
		index1 = indices_with_zone[i] + 1
		index2 = indices_with_zone[i+1] - 1 if i < len(indices_with_zone) - 1 else len(all_data) - 1
		if index2 in indices_with_variables:
			index2 -= 1
		particles = parse_particle_data(all_data[index1:index2+1])
		extracted_num = len(particles)
		time_step_data.append(TimeStepData(time=extracted_time[i], count=extracted_num, particles=particles))

	return time_step_data

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
		elif "log" in dia_name:
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
		print("Invalid input!")
		exit()

	# 定义常量
	xMax = 0.1 # 计算域x方向长度
	yMax = 0.01 # 计算域y方向长度
	nu = 1.51e-5 # 运动粘度
	interval = 60 # 源文件输出时间间隔
	file_interval = 240 # 两个文件之间的时间间隔
	Sh_d = 0.00298 # 冲击临界
	area = xMax*yMax # 计算域底面积

	# 定义文件路径
	if linux_flag:
		working_dir = "/home/ekalhxh/ripple/coll13"
	else:
		working_dir = "E:/Data/Q_on_ero_bed"

	# 定义文件名字典
	case_dict = {
		0: "uStar030_300log50_0_2650_300",
		1: "uStar040_300log50_0_2650_300",
		2: "uStar050_300log50_0_2650_300",
		#3: "uStar060_300log50_0_2650_300",
		#4: "uStar030_300log100_0_2650_300",
		#5: "uStar040_300log100_0_2650_300",
		#6: "uStar050_300log100_0_2650_300",
		#7: "uStar060_300log100_0_2650_300",
		#8: "uStar030_300log200_0_2650_300",
		#9: "uStar040_300log200_0_2650_300",
		#10: "uStar050_300log200_0_2650_300",
		#11: "uStar060_300log200_0_2650_300",
		12: "uStar030_300log300_0_2650_300",
		13: "uStar040_300log300_0_2650_300",
		14: "uStar050_300log300_0_2650_300",
		#15: "uStar060_300log300_0_2650_300",
		16: "uStar035_300log50_0_2650_300",
		17: "uStar045_300log50_0_2650_300",
		18: "uStar055_300log50_0_2650_300",
		#19: "uStar065_300log50_0_2650_300",
		#20: "uStar035_300log100_0_2650_300",
		#21: "uStar045_300log100_0_2650_300",
		#22: "uStar055_300log100_0_2650_300",
		#23: "uStar065_300log100_0_2650_300",
		#24: "uStar035_300log200_0_2650_300",
		#25: "uStar045_300log200_0_2650_300",
		#26: "uStar055_300log200_0_2650_300",
		#27: "uStar065_300log200_0_2650_300",
		28: "uStar035_300log300_0_2650_300",
		29: "uStar045_300log300_0_2650_300",
		30: "uStar055_300log300_0_2650_300",
		#31: "uStar065_300log300_0_2650_300",
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
		#52: "uStar055_269log100_0_2650_300"
	}

	for folder_name in tqdm(case_dict.values(), desc="Processing", unit="case"):
		# 从文件夹名字中提取参数
		parts = folder_name.split("_")
		uStar = float(parts[0][6:])/100
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
		elif "log" in dia_name:
			dia_name_list = dia_name.split("log")
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
		remainder = last_time % file_interval
		if remainder != 0:
			end_file_num += 1
		folder_path = f"{working_dir}/{folder_name}/Particle"
		time_step_data = read_particle_file(folder_path, end_file_num)
		particle_num_data = read_num_file(folder_path)

		for current_data in time_step_data:
			for particle in current_data.particles:
				particle.rho = rhoP

		# 文件路径
		d_in_air_file = os.path.join(working_dir, f"{folder_name}/d_in_air.dat")
		mass_flux_x_file = os.path.join(working_dir, f"{folder_name}/mass_flux_x.dat")
		mass_flux_y_file = os.path.join(working_dir, f"{folder_name}/mass_flux_y.dat")
		mass_flux_z_file = os.path.join(working_dir, f"{folder_name}/mass_flux_z.dat")
		carrying_capacity_file = os.path.join(working_dir, f"{folder_name}/carrying_capacity.dat")
		total_energy_file = os.path.join(working_dir, f"{folder_name}/total_energy.dat")
		particle_num_file = os.path.join(working_dir, f"{folder_name}/particle_num.dat")

		file_dict = {
			1: d_in_air_file,
			2: mass_flux_x_file,
			3: mass_flux_y_file,
			4: mass_flux_z_file,
			5: carrying_capacity_file,
			6: total_energy_file,
			7: particle_num_file
		}

		# 删除已有文件
		for file_path in file_dict.values():
			if os.path.exists(file_path):
				os.remove(file_path)

		# 写入文件头
		write_file_head(d_in_air_file, dia_name, uStar, dia1, dia2, rho, rhoP, s, g_hat, Sh, Ga, "time d d_c d_nc")
		write_file_head(mass_flux_x_file, dia_name, uStar, dia1, dia2, rho, rhoP, s, g_hat, Sh, Ga, "time Q Q_c Q_nc Q_star Q_star_c Q_star_nc")
		write_file_head(mass_flux_y_file, dia_name, uStar, dia1, dia2, rho, rhoP, s, g_hat, Sh, Ga, "time Q Q_c Q_nc Q_star Q_star_c Q_star_nc")
		write_file_head(mass_flux_z_file, dia_name, uStar, dia1, dia2, rho, rhoP, s, g_hat, Sh, Ga, "time Q Q_c Q_nc Q_star Q_star_c Q_star_nc")
		write_file_head(carrying_capacity_file, dia_name, uStar, dia1, dia2, rho, rhoP, s, g_hat, Sh, Ga, "time M M_c M_nc M_star M_star_c M_star_nc")
		write_file_head(total_energy_file, dia_name, uStar, dia1, dia2, rho, rhoP, s, g_hat, Sh, Ga, "time E E_c E_nc")
		write_file_head(particle_num_file, dia_name, uStar, dia1, dia2, rho, rhoP, s, g_hat, Sh, Ga, "time Num iteration")

		with open(particle_num_file, 'a') as file:
			for i in range(len(particle_num_data.time)):
				file.write(f"{particle_num_data.time[i]} {particle_num_data.particleNum[i]} {particle_num_data.iterNum[i]}\n")

		for current_data in time_step_data:
			d_in_air = current_data.calculate_dia_in_air()
			mass_flux = current_data.calculate_mass_flux()
			carrying_capacity = current_data.calculete_carrying_capacity()
			total_energy = current_data.calculate_total_energy()

			mass_flux_star = {
		        "total": [m / (rhoP * d_in_air['total'] * (s * g_hat * d_in_air['total']) ** 0.5) for m in mass_flux['total']],
		        "coll": [m / (rhoP * d_in_air['coll'] * (s * g_hat * d_in_air['coll']) ** 0.5) for m in mass_flux['coll']],
		        "non_coll": [m / (rhoP * d_in_air['non_coll'] * (s * g_hat * d_in_air['non_coll']) ** 0.5) for m in mass_flux['non_coll']]
		    }

			carrying_capacity_star = {
		        "total": carrying_capacity['total'] / (rhoP * d_in_air['total']),
				"coll": carrying_capacity['coll'] / (rhoP * d_in_air['coll']),
				"non_coll": carrying_capacity['non_coll'] / (rhoP * d_in_air['non_coll'])
		    }

			write_file_data(d_in_air_file, f"{current_data.time} {d_in_air['total']} {d_in_air['coll']} {d_in_air['non_coll']}\n")
			write_file_data(mass_flux_x_file, f"{current_data.time} {mass_flux['total'][0]} {mass_flux['coll'][0]} {mass_flux['non_coll'][0]} "
            								  f"{mass_flux_star['total'][0]} {mass_flux_star['coll'][0]} {mass_flux_star['non_coll'][0]}\n")
			write_file_data(mass_flux_y_file, f"{current_data.time} {mass_flux['total'][1]} {mass_flux['coll'][1]} {mass_flux['non_coll'][1]} "
											  f"{mass_flux_star['total'][1]} {mass_flux_star['coll'][1]} {mass_flux_star['non_coll'][1]}\n")
			write_file_data(mass_flux_z_file, f"{current_data.time} {mass_flux['total'][2]} {mass_flux['coll'][2]} {mass_flux['non_coll'][2]} "
											  f"{mass_flux_star['total'][2]} {mass_flux_star['coll'][2]} {mass_flux_star['non_coll'][2]}\n")
			write_file_data(carrying_capacity_file, f"{current_data.time} {carrying_capacity['total']} {carrying_capacity['coll']} {carrying_capacity['non_coll']} "
												  	f"{carrying_capacity_star['total']} {carrying_capacity_star['coll']} {carrying_capacity_star['non_coll']}\n")
			write_file_data(total_energy_file, f"{current_data.time} {total_energy['total']} {total_energy['coll']} {total_energy['non_coll']}\n")