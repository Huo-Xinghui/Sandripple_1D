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
import os # 导入操作系统库
import re # 导入正则表达式库
import matplotlib.pyplot as plt # 导入绘图库

def read_file(file_path) -> list:
	"""读取文件"""
	if not os.path.exists(file_path):
		print(f"文件 {file_path} 不存在。")
		exit()

	with open(file_path, 'r', encoding='utf-8') as file:
		data = file.readlines()
		return data

def read_par_file(file_path) -> tuple:
	"""读取颗粒相关的.dat文件"""
	lines = read_file(file_path)
	parameters = {}
	variables = []
	line_start = 0

	def extract_value(pattern, line, parameter_name):
		"""提取参数名称后面的值"""
		match = re.search(pattern, line)
		if match:
			parameters[parameter_name] = float(match.group(1))
		else:
			raise ValueError(f"Error in reading {parameter_name}.")

	# 读取文件中的各参数值并储存在带条目的parameters列表中
	for i, line in enumerate(lines):
		if "uStar" in line:
			extract_value(r"uStar\s*=\s*([\d.eE+-]+)", line, "u_star")
			if re.search(r"dia\s*=\s*([\d.eE+-]+)", line):
				extract_value(r"dia\s*=\s*([\d.eE+-]+)", line, "d")
				parameters["d1"] = None
				parameters["d2"] = None
				parameters["delta_d"] = None
				if re.search(r"stdd\s*=\s*([\d.eE+-]+)", line):
					extract_value(r"stdd\s*=\s*([\d.eE+-]+)", line, "stdd")
				else:
					parameters["stdd"] = None
			elif re.search(r"dia1\s*=\s*([\d.eE+-]+)", line) and re.search(r"dia2\s*=\s*([\d.eE+-]+)", line):
				extract_value(r"dia1\s*=\s*([\d.eE+-]+)", line, "d1")
				extract_value(r"dia2\s*=\s*([\d.eE+-]+)", line, "d2")
				parameters["delta_d"] = parameters["d2"] - parameters["d1"]
				parameters["d"] = parameters["d1"] + parameters["delta_d"] / 2
				parameters["stdd"] = None
			else:
				raise ValueError("Error in reading particle diameter.")
			extract_value(r"rho\s*=\s*([\d.eE+-]+)", line, "rho")
			extract_value(r"rhoP\s*=\s*([\d.eE+-]+)", line, "rhoP")
		elif "g_hat" in line:
			extract_value(r"s\s*=\s*([\d.eE+-]+)", line, "s")
			extract_value(r"g_hat\s*=\s*([\d.eE+-]+)", line, "g_hat")
			extract_value(r"Sh\s*=\s*([\d.eE+-]+)", line, "Sh")
			extract_value(r"Ga\s*=\s*([\d.eE+-]+)", line, "Ga")
		elif "variable" in line:
			variables = line.split()
			line_start = i + 1
			break

	return parameters, variables, lines[line_start:] # lines中存储数据所在的行

def get_entrys(x_axis, y_axis, collision, nondim) -> tuple:
	"""获取x轴和y轴的条目"""

	# x轴条目字典, key的格式：(x_axis, nondim)
	x_entrys = {
		(0, True): "t_star",
		(0, False): "t",
		(1, False): "u_star",
		(2, True): "Sh",
		(3, True): "Ga",
		(4, True): "s",
		(5, False): "d",
		(6, False): "stdd",
		(7, False): "delta_d"
	}
	# y轴条目字典，key的格式：(y_axis, collision, nondim)
	y_entrys = {
		(0, 0, True): "num",
		(0, 0, False): "num",
		(1, 0, True): "Q_star",
		(1, 0, False): "Q",
		(1, 1, True): "Q_c_star",
		(1, 1, False): "Q_c",
		(1, 2, True): "Q_nc_star",
		(1, 2, False): "Q_nc",
		(2, 0, True): "E_star",
		(2, 0, False): "E",
		(2, 1, True): "E_c_star",
		(2, 1, False): "E_c",
		(2, 2, True): "E_nc_star",
		(2, 2, False): "E_nc",
		(3, 0, False): "d",
		(3, 1, False): "d_c",
		(3, 2, False): "d_nc",
		(4, 0, True): "M_star",
		(4, 0, False): "M",
		(4, 1, True): "M_c_star",
		(4, 1, False): "M_c",
		(4, 2, True): "M_nc_star",
		(4, 2, False): "M_nc"
	}

	# 获取x轴和y轴的条目
	x_entry = x_entrys.get((x_axis, nondim))
	y_entry = y_entrys.get((y_axis, collision, nondim))

	# 检查x轴和y轴的条目是否为空
	if x_entry is None or y_entry is None:
		raise ValueError("Error in getting entrys.")

	return x_entry, y_entry

def process_data(lines, y_axis, parameters) -> list:
	"""将数据存储到带条目的列表中"""
	data = []
	d = parameters["d"] # 用于无量纲化计算
	g_hat = parameters["g_hat"] # 用于无量纲化计算

	def creat_entry(values, keys):
		"""创建条目并存储数据"""
		entry = {
			"t": float(values[0]),
			"t_star": float(values[0]) / (d/g_hat)**0.5
		}
		entry.update({key: float(values[i+1]) for i, key in enumerate(keys)})
		return entry

	# 列表条目字典
	y_axis_keys = {
		0: ["num", "iteration"],
		1: ["Q", "Q_c", "Q_nc", "Q_star", "Q_c_star", "Q_nc_star"],
		2: ["E", "E_c", "E_nc"],
		3: ["d", "d_c", "d_nc"],
		4: ["M", "M_c", "M_nc", "M_star", "M_c_star", "M_nc_star"]
	}

	# 获取当前y_axis条件下数据列表的各条目名称
	keys = y_axis_keys.get(y_axis)
	if keys is None:
		raise ValueError(f"Invalid y_axis: {y_axis}")

	# 将数据存储到带条目的列表中
	for line in lines:
		values = line.split()
		entry = creat_entry(values, keys)
		data.append(entry)

	return data

def smooth_data(data) -> list:
	"""平滑处理数据(三点盒式滤波)"""
	smoothed_data = []
	for i in range(1, len(data)-1):
		smoothed_entry = {}
		for key in data[i].keys():
			if key == "t" or key == "t_star":
				smoothed_entry[key] = data[i][key]
			else:
				smoothed_entry[key] = (data[i-1][key] + data[i][key] + data[i+1][key]) / 3
		smoothed_data.append(smoothed_entry)
	return smoothed_data

if __name__ == '__main__':
	# 控制参数
	x_axis = 2 # x轴类型：0为time，1为u*, 2为Sh, 3为Ga, 4为s, 5为d，6为stdd, 7为delta_d
	y_axis = 1 # y轴类型：0为空中颗粒数，1为颗粒通量，2为颗粒总能量，3为空中颗粒粒径，4为空中颗粒承载量
	direction = 0 # 统计量方向：0为x方向，1为y方向，2为z方向
	collision = 0 # 输出何种颗粒的统计量：0为全部颗粒，1为碰撞颗粒，2为未碰撞颗粒
	nondim = True # 是否无量纲化
	smooth = True # 是否平滑处理(三点盒式滤波)
	legend = 2 # 图例类型：0为u*，1为d，2为Sh，3为Ga，4为s，5为stdd，6为delta_d
	start = 60 # 统计量平均值的开始时间或者变量随时间变化的开始时间
	end = 120 # 统计量平均值的结束时间或者变量随时间变化的结束时间

	# 操作系统
	sys_OS = "w" # "w" for Windows, "l" for Linux
	if sys_OS == "w":
		# Windows
		working_dir = "E:/Data/Sandripples1DFluid/ripple/coll6"
	elif sys_OS == "l":
		# Linux
		working_dir = "home/ekalhxh/ripple/coll"
	else:
		print("Invalid OS system!")
		exit()

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
		30: "uStar035_300_0_2650_1500",
		31: "uStar040_300_0_2650_1500",
		32: "uStar045_300_0_2650_1500",
		33: "uStar050_300_0_2650_1500",
		34: "uStar055_300_0_2650_1500",
		35: "uStar060_300_0_2650_1500",
		#36: "uStar065_300_0_2650_3600",
		#37: "uStar035_300stdd100_0_2650_3600",
		#38: "uStar040_300stdd100_0_2650_3600",
		#39: "uStar045_300stdd100_0_2650_3600",
		#40: "uStar050_300stdd100_0_2650_3600",
		#41: "uStar055_300stdd100_0_2650_3600",
		#42: "uStar060_300stdd100_0_2650_3600",
		#43: "uStar065_300stdd100_0_2650_3600",
	}

	label_dict = {
		"t": r"$t \, (\mathrm{s})$",
		"t_star": r"$t^\ast$",
		"u_star": r"$u^\ast \, (\mathrm{m/s})$",
		"Sh": r"$Sh$",
		"Ga": r"$Ga$",
		"s": r"$s$",
		"stdd": r"$\sigma_d \, (\mathrm{m})$",
		"delta_d": r"$\Delta d \, (\mathrm{m})$",
		"num": r"$N$",
		"Q": r"$Q \, (\mathrm{kg/ms})$",
		"Q_c": r"$Q_c \, (\mathrm{kg/ms})$",
		"Q_nc": r"$Q_{nc} \, (\mathrm{kg/ms})$",
		"Q_star": r"$Q^\ast$",
		"Q_c_star": r"${Q_c}^\ast$",
		"Q_nc_star": r"${Q_{nc}}^\ast$",
		"E": r"$E \, (\mathrm{J})$",
		"E_c": r"$E_c \, (\mathrm{J})$",
		"E_nc": r"$E_{nc} \, (\mathrm{J})$",
		"E_star": r"$E^\ast$",
		"E_c_star": r"${E_c}^\ast$",
		"E_nc_star": r"${E_{nc}}^\ast$",
		"d": r"$d \, (\mathrm{m})$",
		"d_c": r"$d_c \, (\mathrm{m})$",
		"d_nc": r"$d_{nc} \, (\mathrm{m})$",
		"M": r"$M \, (\mathrm{kg/m^2})$",
		"M_c": r"$M_c \, (\mathrm{kg/m^2})$",
		"M_nc": r"$M_{nc} \, (\mathrm{kg/m^2})$",
		"M_star": r"$M^\ast$",
		"M_c_star": r"${M_c}^\ast$",
		"M_nc_star": r"${M_{nc}}^\ast$",
	}

	# 定义读取文件名字典
	output_file_dict = {
		"Q_x": "mass_flux_x.dat",
		"Q_y": "mass_flux_y.dat",
		"Q_z": "mass_flux_z.dat",
		"M": "carrying_capacity.dat",
		"d": "d_in_air.dat",
		"E": "total_energy.dat",
		"N": "particle_num.dat",
	}

	output_list = ["N", "Q_", "E", "d", "M"]
	direction_list = ["x", "y", "z"]
	legend_list = ["u*", "d", "Sh", "Ga", "s", "stdd", "delta_d"]

	output_file_key = output_list[y_axis]
	dir_file_key = direction_list[direction]
	if output_file_key == "Q_":
		output_file_key = "Q_" + dir_file_key

	x_entry, y_entry = get_entrys(x_axis, y_axis, collision, nondim)
	x_label = label_dict[x_entry]
	y_label = label_dict[y_entry]

	plot_x_points = []
	plot_y_points = []
	plt.figure()
	for case_name in case_dict.values():
		read_file_path = os.path.join(working_dir, case_name, output_file_dict[output_file_key])
		tuple_val = read_par_file(read_file_path)
		parameters = tuple_val[0]
		lines = tuple_val[2]
		data = process_data(lines, y_axis, parameters)
		if smooth:
			data = smooth_data(data)
		if x_entry == "t" or x_entry == "t_star":
			plot_x_list = [entry[x_entry] for entry in data if entry["t"] >= start and entry["t"] <= end]
			plot_y_list = [entry[y_entry] for entry in data if entry["t"] >= start and entry["t"] <= end]
			label_str = f"{legend_list[legend]} = {parameters[legend_list[legend]]}"
			plt.scatter(plot_x_list, plot_y_list, label=label_str)
		else:
			plot_x_point = parameters[x_entry]
			plot_y_list = [entry[y_entry] for entry in data if entry["t"] >= start and entry["t"] <= end]
			plot_y_point = sum(plot_y_list) / len(plot_y_list)
			plot_x_points.append(plot_x_point)
			plot_y_points.append(plot_y_point)
	if x_entry != "t" and x_entry != "t_star":
		plt.scatter(plot_x_points, plot_y_points, label="average")
	else:
		plt.legend()
	plt.xlabel(x_label)
	plt.ylabel(y_label)
	save_file_path = os.path.join(working_dir, f"{y_entry} vs {x_entry}.eps")
	plt.savefig(save_file_path)
	plt.show()