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
		time_step_data[extracted_time[i]] = TimeStepData(number=extracted_num, particles=particles)

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

# 主程序
if __name__ == "__main__":
	# 定义文件路径
	working_dir = "E:/Data/Sandripples1DFluid/noMidairColl"
	# 定义常数
	xMax = 0.5
	yMax = 0.01
	zcell = 0.0005
	nu = 1.51e-5
	Sh_d = 0.00298
	area = xMax*yMax
	user_input = False
	if user_input:
		# 用户输入
		start = int(input("Please input start time (s), default=60: ") or 60)
		end = int(input("Please input end time (s), default=120: ") or 120)
		z_num = end // 60
		recalculate = input("Recalculate the data (y/n), default=no: ")
	else:
		start = 60
		end = 120
		z_num = end // 60
		recalculate = "y"
	imp_data_file = os.path.join(working_dir, "impDataNC.dat")
	imp_particle_file = os.path.join(working_dir, "impParticleNC.dat")
	if os.path.exists(imp_particle_file):
		os.remove(imp_particle_file)
	if recalculate.lower() == "y":
		if os.path.exists(imp_data_file):
			os.remove(imp_data_file)
		# 定义文件名字典
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
			#49: "uStar050_250_0_2650_120_nc"
		}
		uStar_list = []
		dia_list = []
		rho_list = []
		rhoP_list = []
		s_list = []
		Sh_list = []
		Ga_list = []
		hMax_list = []
		survT_list = []
		thT_list = []
		T_list = []
		survL_list = []
		di_list = []
		Phi_list = []
		hMax_tot_list = []
		survT_tot_list = []
		thT_tot_list = []
		survL_tot_list = []
		collNum_tot_list = []
		s_tot_list = []
		Sh_tot_list = []
		Ga_tot_list = []
		for i, value in my_dict.items():
			parts = value.split("_")
			uStar_name = parts[0]
			uStar_str = uStar_name[6:]
			uStar_int = int(uStar_str)
			uStar = uStar_int / 100
			uStar_list.append(uStar)
			dia_name = parts[1]
			dia1 = float(dia_name)/1e6
			dia_list.append(dia1)
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
			Sh = rho*uStar**2/(rhoP*g_hat*dia1)
			Sh_list.append(Sh)
			Ga = (s*g_hat*dia1**3)**0.5/nu
			Ga_list.append(Ga)
			last_time = int(parts[4])
			zone_num = last_time // 60
			file_num = zone_num // 4 + 1
			folder_name = value
			folder_path = f"{working_dir}/{folder_name}/Particle"
			results = read_data_file(folder_path, file_num)
			time_step_data = results
			hMax_ave_list = []
			survT_ave_list = []
			survL_ave_list = []
			d_ave_list = []
			Phi_ave_list = []
			thT_ave_list = []
			for current_time in range(start, end+1, 60):
				if z_num <= zone_num:
					current_data = time_step_data[current_time]
					impact_particle = [particle for particle in current_data.particles if particle.h <= zcell and particle.vel[2] < 0]
					len_impact_particle = len(impact_particle)
					hm_list = [particle.hMax for particle in impact_particle]
					st_list = [particle.survT for particle in impact_particle]
					t_list = [2*math.sqrt(2*particle.hMax/9.8) for particle in impact_particle]
					sl_list = [particle.survL for particle in impact_particle]
					n_list = [particle.collNum for particle in impact_particle]
					hMax_tot_list.extend(hm_list)
					survT_tot_list.extend(st_list)
					thT_tot_list.extend(t_list)
					survL_tot_list.extend(sl_list)
					collNum_tot_list.extend(n_list)
					s_tot_list.extend([s]*len_impact_particle)
					Sh_tot_list.extend([Sh]*len_impact_particle)
					Ga_tot_list.extend([Ga]*len_impact_particle)
					hMax_ave = sum(hm_list) / len_impact_particle
					survT_ave = sum(st_list) / len_impact_particle
					thT_ave = sum(t_list) / len_impact_particle
					survL_ave = sum(sl_list) / len_impact_particle
					d_ave = sum([particle.dia for particle in impact_particle]) / len_impact_particle
					Phi_ave = sum([-particle.vel[2]*(math.pi*particle.dia**3)/6*rhoP for particle in impact_particle]) / (area*zcell)
				else:
					hMax_ave = 0
					survT_ave = 0
					survL_ave = 0
					d_ave = 0
					Phi_ave = 0
				hMax_ave_list.append(hMax_ave)
				survT_ave_list.append(survT_ave)
				thT_ave_list.append(thT_ave)
				survL_ave_list.append(survL_ave)
				d_ave_list.append(d_ave)
				Phi_ave_list.append(Phi_ave)
			ave_hMax = sum(hMax_ave_list) / len(hMax_ave_list)
			ave_survT = sum(survT_ave_list) / len(survT_ave_list)
			ave_thT = sum(thT_ave_list) / len(thT_ave_list)
			ave_survL = sum(survL_ave_list) / len(survL_ave_list)
			ave_d = sum(d_ave_list) / len(d_ave_list)
			ave_Phi = sum(Phi_ave_list) / len(Phi_ave_list)
			hMax_list.append(ave_hMax)
			survT_list.append(ave_survT)
			thT_list.append(ave_thT)
			survL_list.append(ave_survL)
			di_list.append(ave_d)
			Phi_list.append(ave_Phi)
		with open(imp_data_file, 'w') as file:
			file.write("uStar dia Sh Ga s hMax survT thSurvT survL Phi\n")
			for i in range(len(uStar_list)):
				file.write(f"{uStar_list[i]} {dia_list[i]} {Sh_list[i]} {Ga_list[i]} {s_list[i]} {hMax_list[i]} {survT_list[i]} {thT_list[i]} {survL_list[i]} {Phi_list[i]}\n")
		with open(imp_particle_file, 'w') as file:
			file.write("s Sh Ga hMax survT thSurvT survL collNum\n")
			for i in range(len(hMax_tot_list)):
				file.write(f"{s_tot_list[i]} {Sh_tot_list[i]} {Ga_tot_list[i]} {hMax_tot_list[i]} {survT_tot_list[i]} {thT_tot_list[i]} {survL_tot_list[i]} {collNum_tot_list[i]}\n")
	else:
		if os.path.exists(imp_data_file):
			with open(imp_data_file, 'r') as file:
				lines = file.readlines()
			uStar_list = []
			dia_list = []
			Sh_list = []
			Ga_list = []
			s_list = []
			hMax_list = []
			survT_list = []
			thT_list = []
			survL_list = []
			Phi_list = []
			for line in lines[1:]:
				columns = line.split()
				uStar_list.append(float(columns[0]))
				dia_list.append(float(columns[1]))
				Sh_list.append(float(columns[2]))
				Ga_list.append(float(columns[3]))
				s_list.append(float(columns[4]))
				hMax_list.append(float(columns[5]))
				survT_list.append(float(columns[6]))
				thT_list.append(float(columns[7]))
				survL_list.append(float(columns[8]))
				Phi_list.append(float(columns[9]))
		else:
			print(f"文件 {imp_data_file} 不存在。")
			exit()
	plt.plot(hMax_tot_list, survT_tot_list, 'ro')
	plt.plot(hMax_tot_list, thT_tot_list, 'b*')
	plt.xlabel('hMax')
	plt.ylabel('survT')
	plt.title(f"hMax vs survT during {start}s to {end}s")
	plt.show()