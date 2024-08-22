# 导入必要的库
import math
import os
import re
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import List
from typing import Dict

@dataclass
class ParticleData:
	loc: List[float]
	vel: List[float]
	vMag: float
	h: float
	hMax: float
	dia: float
	survL: float
	survT: float
	collNum: int

@dataclass
class TimeStepData:
	time: int
	number: int
	particles: List[ParticleData]

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
		time, particleNum, iterNum = [], [], []
		for line in results[1:]:
			columns = line.split()
			time.append(round(float(columns[0])))
			particleNum.append(int(columns[1]))
			iterNum.append(int(columns[2]))
		return time, particleNum, iterNum

# 定义一个函数来读取多个颗粒文件内容
def read_data_file(folder_path, file_num):
	all_data = []
	for i in range(file_num):
		file_path = os.path.join(folder_path, f"ParticleData_{i}.plt")
		results = read_file(file_path)
		if results is not None:
			all_data.extend(results)

	indices_with_zone = find_indices_with_zone(all_data)
	indices_with_variables = find_indices_with_variables(all_data)

	extracted_time = []
	for index in indices_with_zone:
		line = all_data[index]
		numbers = [float(num) for num in re.findall(r'\d+\.?\d*', line)]
		extracted_time.extend(numbers)
	extracted_time = [round(num) for num in extracted_time]

	time_step_data: Dict[float, TimeStepData] = {}

	for i in range(len(indices_with_zone)):
		index1 = indices_with_zone[i] + 1
		if i == len(indices_with_zone) - 1:
			index2 = len(all_data) - 1
		else:
			index2 = indices_with_zone[i+1] - 1
		if index2 in indices_with_variables:
			index2 -= 1
		particles = []
		for line in all_data[index1:index2+1]:
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
		extracted_num = len(particles)
		time_step_data[extracted_time[i]] = TimeStepData(time=extracted_time[i], number=extracted_num, particles=particles)

	return time_step_data

def find_indices_with_zone(all_data):
	indices = []
	for i, line in enumerate(all_data):
		if "zone" in line:
			indices.append(i)
	return indices

def find_indices_with_variables(all_data):
	indices = []
	for i, line in enumerate(all_data):
		if "variables" in line:
			indices.append(i)
	return indices

def plot_Num_vs_time(time, particleNum, output_N_t):
	if output_N_t:
		plt.plot(time, particleNum)
		plt.xlabel('Time')
		plt.ylabel('Number')
		plt.title('Number vs Time')
		plt.show()
	else:
		return

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
	# 定义一些常量
	working_dir = "E:/Data/Sandripples1DFluid"
	file_num = 16
	xMax = 0.5
	yMax = 0.01
	plot_N_t = False
	plot_Q_t = False
	output_ave = True
	area = xMax*yMax
	# 用户输入
	if output_ave:
		my_dict = {
			0: "uStar03_250_0_2650",
			1: "uStar035_250_0_2650",
			2: "uStar04_250_0_2650",
			3: "uStar045_250_0_2650",
			4: "uStar05_250_0_2650",
			5: "uStar055_250_0_2650",
			6: "uStar06_250_0_2650",
			7: "uStar065_250_0_2650",
			8: "uStar045_200_0_2650",
			9: "uStar045_300_0_2650",
			10: "uStar045_350_0_2650",
			11: "uStar045_400_0_2650",
			12: "uStar045_250_2_2650",
			13: "uStar045_250_3_2650",
			14: "uStar045_250_1_2000",
			15: "uStar045_250_1_3000",
			16: "uStar045_250_1_4000",
			17: "uStar045_150and350_0_2650",
			18: "uStar045_250and450_0_2650"
		}
		uStar_list = []
		dia1_list = []
		dia2_list = []
		rho_list = []
		rhoP_list = []
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
				dia2 = None
			dia1_list.append(dia1)
			dia2_list.append(dia2)
			rho_name = parts[2]
			if rho_name == "0":
				rho = 1.263
			else:
				rho = float(rho_name)
			rho_list.append(rho)
			rhoP_name = parts[3]
			rhoP = float(rhoP_name)
			rhoP_list.append(rhoP)
		#ave_Q = sum(flux_x_list) / len(flux_x_list)
		#print(f"Average Q: {ave_Q}")
		plt.plot(uStar_list, rho_list, 'ro')
		plt.xlabel('uStar')
		plt.ylabel('rho')
		plt.title('uStar vs rho')
		plt.show()
	else:
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
		if results1 is not None:
			time, particleNum, iterNum = results1
			plot_Num_vs_time(time, particleNum, plot_N_t)
		if results2 is not None:
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