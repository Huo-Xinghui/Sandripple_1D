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
import numpy as np

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
	x_axis = 1 # x轴类型：0为time，1为u*, 2为Sh, 3为Ga, 4为s, 5为d，6为stdd, 7为delta_d
	y_axis = 3 # y轴类型：0为空中颗粒数，1为颗粒通量，2为颗粒总能量，3为空中颗粒粒径，4为空中颗粒承载量
	direction = 0 # 统计量方向：0为x方向，1为y方向，2为z方向
	collision = 0 # 输出何种颗粒的统计量：0为全部颗粒，1为碰撞颗粒，2为未碰撞颗粒
	nondim = False # 是否无量纲化
	smooth = False # 是否平滑处理(三点盒式滤波)
	fit = True # 是否进行拟合
	legend = 5 # 图例类型：0为u*，1为d，2为Sh，3为Ga，4为s，5为stdd，6为delta_d
	start = 120 # 统计量平均值的开始时间或者变量随时间变化的开始时间
	end = 300 # 统计量平均值的结束时间或者变量随时间变化的结束时间

	# 操作系统
	sys_OS = "w" # "w" for Windows, "l" for Linux
	if sys_OS == "w":
		# Windows
		working_dir = "E:/Data/Q_on_flat_bed"
		#working_dir = "E:/Data/Sandripples1DFluid/ripple/coll13"
	elif sys_OS == "l":
		# Linux
		working_dir = "/home/ekalhxh/ripple/coll13"
	else:
		print("Invalid OS system!")
		exit()

	# 定义文件名字典
	#case_dict = {
	#	0: "uStar040_300_0_2650_300",
	#	1: "uStar045_300_0_2650_300",
	#	2: "uStar050_300_0_2650_300",
	#	3: "uStar055_300_0_2650_300",
	#	4: "uStar060_300_0_2650_300",
	#	5: "uStar065_300_0_2650_300",
	#	6: "uStar040_400_0_2650_300",
	#	7: "uStar045_400_0_2650_300",
	#	8: "uStar050_400_0_2650_300",
	#	9: "uStar055_400_0_2650_300",
	#	10: "uStar060_400_0_2650_300",
	#	11: "uStar065_400_0_2650_300",
	#	12: "uStar040_296log366_0_2650_300",
	#	13: "uStar050_296log366_0_2650_300",
	#	14: "uStar060_296log366_0_2650_300",
	#	15: "uStar050_260log537_0_2650_300",
	#	16: "uStar055_260log537_0_2650_300",
	#	17: "uStar060_260log537_0_2650_300",
	#	18: "uStar050_268log491_0_2650_300",
	#	19: "uStar055_268log491_0_2650_300",
	#	20: "uStar060_268log491_0_2650_300",
	#	21: "uStar040_352log522_0_2650_300",
	#	22: "uStar050_352log522_0_2650_300",
	#	23: "uStar060_352log522_0_2650_300",
	#	24: "uStar040_279log368_0_2650_300",
	#	25: "uStar050_279log368_0_2650_300",
	#	26: "uStar060_279log368_0_2650_300",
	#}
	case_dict = {
		0: "uStar030_300log50_0_2650_300",
		1: "uStar040_300log50_0_2650_300",
		2: "uStar050_300log50_0_2650_300",
		3: "uStar060_300log50_0_2650_300",
		4: "uStar030_300log100_0_2650_300",
		5: "uStar040_300log100_0_2650_300",
		6: "uStar050_300log100_0_2650_300",
		7: "uStar060_300log100_0_2650_300",
		8: "uStar030_300log200_0_2650_300",
		9: "uStar040_300log200_0_2650_300",
		10: "uStar050_300log200_0_2650_300",
		11: "uStar060_300log200_0_2650_300",
		12: "uStar030_300log300_0_2650_300",
		13: "uStar040_300log300_0_2650_300",
		14: "uStar050_300log300_0_2650_300",
		15: "uStar060_300log300_0_2650_300",
		16: "uStar035_300log50_0_2650_300",
		17: "uStar045_300log50_0_2650_300",
		18: "uStar055_300log50_0_2650_300",
		19: "uStar065_300log50_0_2650_300",
		20: "uStar035_300log100_0_2650_300",
		21: "uStar045_300log100_0_2650_300",
		22: "uStar055_300log100_0_2650_300",
		23: "uStar065_300log100_0_2650_300",
		24: "uStar035_300log200_0_2650_300",
		25: "uStar045_300log200_0_2650_300",
		26: "uStar055_300log200_0_2650_300",
		27: "uStar065_300log200_0_2650_300",
		28: "uStar035_300log300_0_2650_300",
		29: "uStar045_300log300_0_2650_300",
		30: "uStar055_300log300_0_2650_300",
		31: "uStar065_300log300_0_2650_300",
		32: "uStar040_430log100_0_2650_300",
		33: "uStar050_430log100_0_2650_300",
		34: "uStar060_430log100_0_2650_300",
		35: "uStar030_167log100_0_2650_300",
		36: "uStar040_167log100_0_2650_300",
		37: "uStar050_167log100_0_2650_300",
		38: "uStar030_269log100_0_2650_300",
		39: "uStar040_269log100_0_2650_300",
		40: "uStar050_269log100_0_2650_300",
		41: "uStar030_321log100_0_2650_300",
		42: "uStar040_321log100_0_2650_300",
		43: "uStar050_321log100_0_2650_300",
		44: "uStar030_240log50_0_2650_300",
		45: "uStar035_240log50_0_2650_300",
		46: "uStar040_240log50_0_2650_300",
		47: "uStar045_240log50_0_2650_300",
		48: "uStar050_240log50_0_2650_300",
		49: "uStar055_240log50_0_2650_300",
		50: "uStar035_269log100_0_2650_300",
		51: "uStar045_269log100_0_2650_300",
		52: "uStar055_269log100_0_2650_300",
		53: "uStar040_400log50_0_2650_300",
		54: "uStar050_400log50_0_2650_300",
		55: "uStar060_400log50_0_2650_300",
		56: "uStar030_250log25_0_2650_300",
		57: "uStar040_250log25_0_2650_300",
		58: "uStar050_250log25_0_2650_300",
		59: "uStar060_250log25_0_2650_300",
		60: "uStar030_271log121_0_2650_300",
		61: "uStar040_271log121_0_2650_300",
		62: "uStar050_271log121_0_2650_300",
		63: "uStar060_271log121_0_2650_300",
		64: "uStar030_317log252_0_2650_300",
		65: "uStar040_317log252_0_2650_300",
		66: "uStar050_317log252_0_2650_300",
		67: "uStar060_317log252_0_2650_300",
		68: "uStar030_347log537_0_2650_300",
		69: "uStar040_347log537_0_2650_300",
		70: "uStar050_347log537_0_2650_300",
		71: "uStar060_347log537_0_2650_300",
		72: "uStar030_290log97_0_2650_300",
		73: "uStar040_290log97_0_2650_300",
		74: "uStar050_290log97_0_2650_300",
		75: "uStar060_290log97_0_2650_300",
		76: "uStar030_197log65_0_2650_300",
		77: "uStar040_197log65_0_2650_300",
		78: "uStar050_197log65_0_2650_300",
		79: "uStar060_197log65_0_2650_300",
	}

	# Creyssels et al. (2009)的数据
	Cr09_Sh = [0.012, 0.022, 0.035, 0.05, 0.068, 0.098]
	Cr09_delta_Sh = [2*0.05*Sh for Sh in Cr09_Sh]
	Cr09_Q_star = [0.00396, 0.00984, 0.0163, 0.02464, 0.0376, 0.06144]
	Cr09_delta_Q_star = [0.05*Q_star for Q_star in Cr09_Q_star]
	Cr09_data = {
		"Sh": Cr09_Sh,
		"delta_Sh": Cr09_delta_Sh,
		"Q_star": Cr09_Q_star,
		"delta_Q_star": Cr09_delta_Q_star
	}

	# 定义参照量字典
	others_case_dict = {
		"Cr09": Cr09_data
	}

	# 参照字典中的条目数
	other_case_num = len(others_case_dict)

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
	if other_case_num > 0 and x_axis == 2 and y_axis == 1 and nondim:
		for case_name, case_data in others_case_dict.items():
			plot_x_list = case_data["Sh"]
			plot_x_error = case_data["delta_Sh"]
			plot_y_list = case_data["Q_star"]
			plot_y_error = case_data["delta_Q_star"]
			plt.errorbar(plot_x_list, plot_y_list, xerr=plot_x_error, yerr=plot_y_error, fmt='ro', label=case_name)
	if x_entry != "t" and x_entry != "t_star":
		plt.scatter(plot_x_points, plot_y_points)
		if (x_axis == 2 and fit) and (y_axis == 1 or y_axis == 4):
			# 线性拟合
			coeffs = np.polyfit(plot_x_points, plot_y_points, 1)
			fit_line = np.poly1d(coeffs)
			# 计算y=0时的x值（x轴截距）
			x_intercept = -coeffs[1] / coeffs[0]
			# 扩展x轴范围到x轴截距
			x_fit = np.linspace(min(min(plot_x_points), x_intercept), max(max(plot_x_points), x_intercept), 10000)
			# 在x轴截距处添加标注
			plt.annotate(f'x = {x_intercept:.6f}',
				xy=(x_intercept, 0),
				xytext=(x_intercept, max(plot_y_points)/4),
				ha='center',
				va='bottom',
				arrowprops=dict(arrowstyle='->'))
			plt.plot(x_fit, fit_line(x_fit), '--', label=f'fit line: $y = {coeffs[0]:.6f}x {coeffs[1]:+.6f}$')
	plt.legend()
	plt.xlabel(x_label)
	plt.ylabel(y_label)
	save_file_path = os.path.join(working_dir, f"{y_entry} vs {x_entry}.eps")
	plt.savefig(save_file_path)
	plt.show()