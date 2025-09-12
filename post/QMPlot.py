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
"""用于绘制颗粒x,z方向通量的图像"""
# ********************************************************************************************************

# 导入必要的库
import os # 导入操作系统库
import re # 导入正则表达式库
import matplotlib.pyplot as plt # 导入绘图库
import numpy as np
import matplotlib as mpl

# Set up matplotlib parameters
mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = r'\usepackage{amssymb}'

def read_file(file_path) -> list:
	"""读取文件"""
	if not os.path.exists(file_path):
		print(f"文件 {file_path} 不存在。")
		exit()

	with open(file_path, 'r', encoding='utf-8') as file:
		data = file.readlines()
		return data

def read_par_file(file_path) -> tuple:
	"""
	读取.dat文件.
	dict: parameters={"d", "stdd", "s", "Sh", "Ga"...}.
	list: variables (行中每个元素的名称) = ["time", "Q", "Q_c", "Q_nc"...].
	list: lines (存储数据所在的行) = [<line1>, <line2>, ...]
	"""
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

def get_entrys(y_axis, nondim) -> tuple:
	"""获取x轴和y轴的条目"""

	# x轴条目字典, key的格式：(x_axis, nondim)
	x_entrys = {
		False: "u_star",
		True: "Sh",
	}
	# y轴条目字典，key的格式：(y_axis, nondim)
	y_entrys = {
		(0, True): "Q_star",
		(0, False): "Q",
		(1, True): "M_star",
		(1, False): "M",
		(2, False): "E",
		(2, True): "E",
	}

	# 获取x轴和y轴的条目
	x_entry = x_entrys.get(nondim)
	y_entry = y_entrys.get((y_axis, nondim))

	# 检查x轴和y轴的条目是否为空
	if x_entry is None or y_entry is None:
		raise ValueError("Error in getting entrys.")

	return x_entry, y_entry

def process_data(lines, y_axis) -> list:
	"""将数据存储到data列表中, data中的元素为dict, dict的key包含{t, tstar, Q或M..}"""
	data = []

	def creat_entry(values, keys):
		"""创建条目并存储数据"""
		entry = {
			"t": float(values[0]),
		}
		entry.update({key: float(values[i+1]) for i, key in enumerate(keys)})
		return entry

	# 列表条目字典
	y_axis_keys = {
		0: ["Q", "Q_c", "Q_nc", "Q_star", "Q_c_star", "Q_nc_star"],
		1: ["M", "M_c", "M_nc", "M_star", "M_c_star", "M_nc_star"],
		2: ["E", "E_c", "E_nc"],
		3: ["d", "d50", "d90"]
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

def rslt_plot(ax, rslts_dict, dictkey, legend_str, fitline, arrow, lims, sizes, style):
	x = rslts_dict[dictkey]["x"]
	y = rslts_dict[dictkey]["y"]
	x_max = lims["x_max"]
	ax.loglog(x, y, style, markersize=sizes["marker"], label=legend_str)
	if fitline:
		coeffs = rslts_dict[dictkey]["fit"]
		fit_line = np.poly1d(coeffs)
		# 计算y=0时的x值（x轴截距）
		x_intercept = -coeffs[1] / coeffs[0]
		# 扩展x轴范围到x轴截距
		x_fit = np.linspace(min(min(x), x_intercept), max(x_max, x_intercept), 10000)
		# 在x轴截距处添加标注
		if arrow:
			plt.annotate(f'x = {x_intercept:.6f}',
				xy=(x_intercept, 0),
				xytext=(x_intercept, max(y)/4),
				ha='center',
				va='bottom',
				arrowprops=dict(arrowstyle='->'))
		ax.plot(x_fit, fit_line(x_fit), '--', color='gray', linewidth=sizes["line"])

# 控制参数
ampf = 1/0.6 # amplification factor
label_size = 15*ampf
ticks_size = 10*ampf
marker_size = 10*ampf
line_width = 1.5
marker_width = 2
sizes = {
	"marker": marker_size,
	"line": line_width
}
rhof = 1.263
rhop = 2650
s = rhop / rhof
g = 9.81 * (1.0 - 1.0 / s)
y_axis = 1 # y轴类型：0为颗粒通量，1为空中颗粒承载量, 2为颗粒总能量
direction = 0 # 统计量方向：0为x方向，1为z方向
nondim = True # 是否无量纲化
nondim_d = 0 # 无量纲化所使用的粒径: 0为dm, 1为d50, 2为d90, 3为平均空中粒径, 4为空中d50, 5为空中d90
comparison = True # 是否与他人行比较, only Q*
start = 60 # 统计量平均值的开始时间
end = 300 # 统计量平均值的结束时间
use_lims = False # 是否使用坐标轴限制
lims = {
	"x_min": 0,
	"x_max": 0.14,
	"y_min": 0,
	"y_max": 0.1,
}

# 操作系统
sys_OS = "w" # "w" for Windows, "l" for Linux
if sys_OS == "w":
	# Windows
	#working_dir = "E:/Data/Q_on_flat_bed"
	#working_dir = "E:/Data/Q_on_flat_bed_Eeff"
	#working_dir = "E:/Data/mono_in_air"
	working_dir = "E:/Data/mono_on_bed"
	#working_dir = "E:/Data/Sandripples1DFluid/ripple/coll13"
elif sys_OS == "l":
	# Linux
	working_dir = "/home/ekalhxh/ripple/coll13"
else:
	print("Invalid OS system!")
	exit()

# 定义文件名字典
case_dict1 = {
	0: "uStar030_300log50_0_2650_300",
	2: "uStar040_300log50_0_2650_300",
	4: "uStar050_300log50_0_2650_300",
	6: "uStar060_300log50_0_2650_300",
	#1: "uStar035_300log50_0_2650_300",
	#3: "uStar045_300log50_0_2650_300",
	#5: "uStar055_300log50_0_2650_300",
	#7: "uStar065_300log50_0_2650_300",
}
case_dict2 = {
	0: "uStar030_300log100_0_2650_300",
	2: "uStar040_300log100_0_2650_300",
	4: "uStar050_300log100_0_2650_300",
	6: "uStar060_300log100_0_2650_300",
	#1: "uStar035_300log100_0_2650_300",
	#3: "uStar045_300log100_0_2650_300",
	#5: "uStar055_300log100_0_2650_300",
	#7: "uStar065_300log100_0_2650_300",
}
case_dict3 = {
	0: "uStar030_300log200_0_2650_300",
	2: "uStar040_300log200_0_2650_300",
	4: "uStar050_300log200_0_2650_300",
	6: "uStar060_300log200_0_2650_300",
	#1: "uStar035_300log200_0_2650_300",
	#3: "uStar045_300log200_0_2650_300",
	#5: "uStar055_300log200_0_2650_300",
	#7: "uStar065_300log200_0_2650_300",
}
case_dict4 = {
	0: "uStar030_300log300_0_2650_300",
	2: "uStar040_300log300_0_2650_300",
	4: "uStar050_300log300_0_2650_300",
	6: "uStar060_300log300_0_2650_300",
	#1: "uStar035_300log300_0_2650_300",
	#3: "uStar045_300log300_0_2650_300",
	#5: "uStar055_300log300_0_2650_300",
	#7: "uStar065_300log300_0_2650_300",
}

case_dict5 = {
	0: "uStar040_430log100_0_2650_300",
	1: "uStar050_430log100_0_2650_300",
	2: "uStar060_430log100_0_2650_300",
}

case_dict6 = {
	0: "uStar030_167log100_0_2650_300",
	1: "uStar040_167log100_0_2650_300",
	2: "uStar050_167log100_0_2650_300",
}

case_dict7 = {
	0: "uStar030_269log100_0_2650_300",
	1: "uStar040_269log100_0_2650_300",
	2: "uStar050_269log100_0_2650_300",
	3: "uStar035_269log100_0_2650_300",
	4: "uStar045_269log100_0_2650_300",
	5: "uStar055_269log100_0_2650_300",
}

case_dict8 = {
	0: "uStar030_321log100_0_2650_300",
	1: "uStar040_321log100_0_2650_300",
	2: "uStar050_321log100_0_2650_300",
}

case_dict9 = {
	0: "uStar030_240log50_0_2650_300",
	1: "uStar040_240log50_0_2650_300",
	2: "uStar050_240log50_0_2650_300",
	3: "uStar035_240log50_0_2650_300",
	4: "uStar045_240log50_0_2650_300",
	5: "uStar055_240log50_0_2650_300",
}

case_dict10 = {
	1: "uStar040_400log50_0_2650_300",
	2: "uStar050_400log50_0_2650_300",
	3: "uStar060_400log50_0_2650_300"
}

case_dict11 = {
	1: "uStar030_250log25_0_2650_300",
	2: "uStar040_250log25_0_2650_300",
	3: "uStar050_250log25_0_2650_300",
	4: "uStar060_250log25_0_2650_300",
	#5: "uStar035_250log25_0_2650_300",
	#6: "uStar045_250log25_0_2650_300",
	#7: "uStar055_250log25_0_2650_300",
	#8: "uStar065_250log25_0_2650_300"
}

case_dict12 = {
	1: "uStar030_271log121_0_2650_300",
	2: "uStar040_271log121_0_2650_300",
	3: "uStar050_271log121_0_2650_300",
	4: "uStar060_271log121_0_2650_300",
	#5: "uStar035_271log121_0_2650_300",
	#6: "uStar045_271log121_0_2650_300",
	#7: "uStar055_271log121_0_2650_300",
	#8: "uStar065_271log121_0_2650_300"
}

case_dict13 = {
	1: "uStar030_317log252_0_2650_300",
	2: "uStar040_317log252_0_2650_300",
	3: "uStar050_317log252_0_2650_300",
	4: "uStar060_317log252_0_2650_300",
	#5: "uStar035_317log252_0_2650_300",
	#6: "uStar045_317log252_0_2650_300",
	#7: "uStar055_317log252_0_2650_300",
	#8: "uStar065_317log252_0_2650_300"
}

case_dict14 = {
	1: "uStar030_347log537_0_2650_300",
	2: "uStar040_347log537_0_2650_300",
	3: "uStar050_347log537_0_2650_300",
	4: "uStar060_347log537_0_2650_300",
	#5: "uStar035_347log537_0_2650_300",
	#6: "uStar045_347log537_0_2650_300",
	#7: "uStar055_347log537_0_2650_300",
	#8: "uStar065_347log537_0_2650_300",
}

case_dict15 = {
	1: "uStar030_290log97_0_2650_300",
	2: "uStar040_290log97_0_2650_300",
	3: "uStar050_290log97_0_2650_300",
	4: "uStar060_290log97_0_2650_300"
}

case_dict16 = {
	1: "uStar030_197log65_0_2650_300",
	2: "uStar040_197log65_0_2650_300",
	3: "uStar050_197log65_0_2650_300",
	4: "uStar060_197log65_0_2650_300"
}

case_dict17 = {
	1: "uStar030_150log50_0_2650_300",
	2: "uStar040_150log50_0_2650_300",
	3: "uStar050_150log50_0_2650_300",
	4: "uStar060_150log50_0_2650_300",
}
case_dict18 = {
	1: "uStar030_150log100_0_2650_300",
	2: "uStar040_150log100_0_2650_300",
	3: "uStar050_150log100_0_2650_300",
	4: "uStar060_150log100_0_2650_300",
}
case_dict19 = {
	1: "uStar030_150log200_0_2650_300",
	2: "uStar040_150log200_0_2650_300",
	3: "uStar050_150log200_0_2650_300",
	4: "uStar060_150log200_0_2650_300",
}
case_dict20 = {
	1: "uStar030_150log300_0_2650_300",
	2: "uStar040_150log300_0_2650_300",
	3: "uStar050_150log300_0_2650_300",
	4: "uStar060_150log300_0_2650_300",
}

cases_dict = {
	"d300stdd50": case_dict1,
	"d300stdd100": case_dict2,
	"d300stdd200": case_dict3,
	"d300stdd300": case_dict4,
	#"d430stdd100": case_dict5,
	#"d167stdd100": case_dict6,
	#"d269stdd100": case_dict7,
	#"d321stdd100": case_dict8,
	#"d240stdd50": case_dict9,
	#"d400stdd50": case_dict10,
	"d250stdd25": case_dict11,
	"d271stdd121": case_dict12,
	"d317stdd252": case_dict13,
	"d347stdd537": case_dict14,
	#"d290stdd97": case_dict15,
	#"d197stdd65": case_dict16,
	#"d150stdd50": case_dict17,
	#"d150stdd100": case_dict18,
	#"d150stdd200": case_dict19,
	#"d150stdd300": case_dict20,
}

d_dict ={
	"d300stdd50": {"d50": 2.96e-4, "d90": 3.66e-4},
	"d300stdd100": {"d50": 2.85e-4, "d90": 4.32e-4},
	"d300stdd200": {"d50": 2.61e-4, "d90": 5.37e-4},
	"d300stdd300": {"d50": 2.49e-4, "d90": 5.89e-4},
	"d430stdd100": {"d50": 4.09e-4, "d90": 6.18e-4},
	"d167stdd100": {"d50": 1.57e-4, "d90": 2.33e-4},
	"d269stdd100": {"d50": 2.68e-4, "d90": 3.61e-4},
	"d321stdd100": {"d50": 3.01e-4, "d90": 4.44e-4},
	"d240stdd50": {"d50": 2.35e-4, "d90": 3.06e-4},
	"d400stdd50": {"d50": 3.97e-4, "d90": 4.66e-4},
	"d250stdd25": {"d50": 2.48e-4, "d90": 2.82e-4},
	"d271stdd121": {"d50": 2.50e-4, "d90": 4.15e-4},
	"d317stdd252": {"d50": 2.65e-4, "d90": 5.84e-4},
	"d347stdd537": {"d50": 2.81e-4, "d90": 6.85e-4},
	"d290stdd97": {"d50": 2.75e-4, "d90": 4.16e-4},
	"d197stdd65": {"d50": 1.86e-4, "d90": 2.80e-4},
	"d150stdd50": {"d50": 2.96e-4, "d90": 3.66e-4},
	"d150stdd100": {"d50": 2.85e-4, "d90": 4.32e-4},
	"d150stdd200": {"d50": 2.61e-4, "d90": 5.37e-4},
	"d150stdd300": {"d50": 2.49e-4, "d90": 5.89e-4},
}

# Creyssels et al. (2009)的数据
d_cr09 = 2.42e-4
rhop_cr09 = 2500
rhof_cr09 = 1.2365
s_cr09 = rhop_cr09 / rhof_cr09
g_cr09 = 9.81*(1.0 - 1.0/s_cr09)
us_array = np.array([0.24, 0.32, 0.40, 0.48, 0.56, 0.67])
Sh_array = np.array([0.012, 0.022, 0.035, 0.05, 0.068, 0.098])
Q_array = 1.0e-3 * np.array([5.25, 12.84, 21.09, 32, 48.85, 79.55])
Cr09_data = {
	"u_star": us_array,
	"u_star_err": 0.05 * us_array,
	"Sh": Sh_array,
	"Sh_err": 0.1 * Sh_array,
	"Q": Q_array,
	"Q_err": 0.05 * Q_array,
	"Q_star": Q_array/(rhop_cr09*d_cr09*np.sqrt(s_cr09*g_cr09*d_cr09)),
	"Q_star_err": 0.05 * Q_array
}

# Ho 2011 的数据
d_ho11 = 2.30e-4
rhop_ho11 = 2450
rhof_ho11 = 1.2
s_ho11 = rhop_ho11 / rhof_ho11
g_ho11 = 9.81
Sh_array = np.array([0.01776, 0.04970, 0.08710, 0.13516, 0.17804, 0.32099])
Q_array = np.array([0.35012, 1.12710, 2.34532, 3.64988, 5.59712, 9.11271])
Q_array = rhop_ho11 * np.sqrt(g_ho11*d_ho11) * d_ho11 * Q_array
Ho2011_data = {
	"Sh": Sh_array,
	"Q": Q_array,
	"Q_err": 0.05 * Q_array,
	"Q_star": Q_array/(rhop_ho11*d_ho11*np.sqrt(s_ho11*g_ho11*d_ho11)),
	"Q_star_err": 0.05 * Q_array
}

# 定义参照量字典
others_case_dict = {
	"Cr09": Cr09_data,
	"Ho2011": Ho2011_data
}

# 确定读取文件名及x,y轴标签
output_file_dict = {
	"Q_x": "mass_flux_x.dat",
	"Q_z": "mass_flux_z.dat",
	"M": "carrying_capacity.dat",
	"E": "total_energy.dat",
	"d": "d_in_air.dat"
}
output_list = ["Q_", "M", "E"]
direction_list = ["x", "z"]
output_file_key = output_list[y_axis]
dir_file_key = direction_list[direction]
if output_file_key == "Q_":
	output_file_key = "Q_" + dir_file_key
label_dict = {
	"u_star": r"$u^\ast \, (\mathrm{m/s})$",
	"Sh": r"$Sh$",
	"Q": r"$Q \, (\mathrm{kg/(ms)})$",
	"Q_star": r"$Q^\ast$",
	"M": r"$M \, (\mathrm{kg/m^2})$",
	"M_star": r"$M^\ast$",
	"E": r"$E \, (\mathrm{J})$",
}
x_entry, y_entry = get_entrys(y_axis, nondim)
x_label = label_dict[x_entry]
y_label = label_dict[y_entry]

tail_list = ["dm", "d50", "d90", "dair", "dair50", "dair90"]
rslts_dict = {}
for dictkey, case_dict in cases_dict.items():
	x_points = []
	y_points = []
	for case_name in case_dict.values():
		read_file_path = os.path.join(working_dir, case_name, output_file_dict[output_file_key])
		tuple_val = read_par_file(read_file_path)
		parameters = tuple_val[0]
		lines = tuple_val[2]
		data = process_data(lines, y_axis)
		x_point = parameters[x_entry]
		y_list = [case[y_entry] for case in data if case["t"] >= start and case["t"] <= end]
		y_array = np.array(y_list)

		if nondim_d >= 3:
			if nondim_d == 3:
				case_key = "d"
			elif nondim_d == 4:
				case_key = "d50"
			elif nondim_d == 5:
				case_key = "d90"
			read_file_path = os.path.join(working_dir, case_name, output_file_dict["d"])
			tuple_val = read_par_file(read_file_path)
			lines = tuple_val[2]
			data = process_data(lines, 3)
			d_list = [case[case_key] for case in data if case["t"] >= start and case["t"] <= end]
			d_array = np.array(d_list)

		if nondim and nondim_d == 1:
			d_old = parameters["d"]
			d = d_dict[dictkey]["d50"]
			x_point = x_point * d_old / d
			if y_axis == 0:
				nond_old = rhop * d_old * np.sqrt(s * g * d_old)
				nond = rhop * d * np.sqrt(s * g * d)
			else:
				nond_old = rhop * d_old
				nond = rhop * d
			y_array = y_array * nond_old / nond
		elif nondim and nondim_d == 2:
			d_old = parameters["d"]
			d = d_dict[dictkey]["d90"]
			x_point = x_point * d_old / d
			if y_axis == 0:
				nond_old = rhop * d_old * np.sqrt(s * g * d_old)
				nond = rhop * d * np.sqrt(s * g * d)
			else:
				nond_old = rhop * d_old
				nond = rhop * d
			y_array = y_array * nond_old / nond
		elif nondim and nondim_d >= 3:
			d_old = parameters["d"]
			d = d_array.mean()
			x_point = x_point * d_old / d
			if y_axis == 0:
				nond_old = rhop * d_old * np.sqrt(s * g * d_old)
				nond = rhop * d * np.sqrt(s * g * d)
			else:
				nond_old = rhop * d_old
				nond = rhop * d
			y_array = y_array * nond_old / nond
		y_point = np.mean(y_array)
		x_points.append(x_point)
		y_points.append(y_point)
	# 计算拟合曲线
	coeffs = np.polyfit(x_points, y_points, 1)
	rslts_dict[dictkey] = {"x": np.array(x_points), "y": np.array(y_points), "fit": coeffs}
	tail = tail_list[nondim_d]
	if y_axis == 0:
		file_name = f"Q_monobed_{dictkey}_{tail}.npz"
	elif y_axis == 1:
		file_name = f"M_monobed_{dictkey}_{tail}.npz"
	np.savez(file_name, **rslts_dict[dictkey])

# 绘图
fig = plt.figure(1, figsize=(8, 6), constrained_layout=True)
ax = fig.gca()

if comparison and nondim and y_entry == "Q_star":
	ax.errorbar(others_case_dict["Cr09"]["Sh"],
			 	others_case_dict["Cr09"]["Q_star"],
				xerr=others_case_dict["Cr09"]["Sh_err"],
				yerr=others_case_dict["Cr09"]["Q_star_err"],
				fmt='o',
				color='orange',
				label='Cr09',
				markersize=marker_size,
				markerfacecolor='none',
				markeredgewidth=marker_width
				)
	ax.errorbar(others_case_dict["Ho2011"]["Sh"],
				others_case_dict["Ho2011"]["Q_star"],
				yerr=others_case_dict["Ho2011"]["Q_star_err"],
				fmt='o',
				color='purple',
				label='Ho2011',
				markersize=marker_size,
				markerfacecolor='none',
				markeredgewidth=marker_width
				)

#legend_str = "Narrow"
#dictkey = "d300stdd50"
#style = "ro"
#arrow = True
#fitline = True
#rslt_plot(ax, rslts_dict, dictkey, legend_str, fitline, arrow, lims, sizes, style)
#
#legend_str = "Narrow 240"
#dictkey = "d240stdd50"
#style = "rv"
#arrow = False
#fitline = False
#rslt_plot(ax, rslts_dict, dictkey, legend_str, fitline, arrow, lims, sizes, style)
#
#legend_str = "Narrow 400"
#dictkey = "d400stdd50"
#style = "rv"
#arrow = False
#fitline = False
#rslt_plot(ax, rslts_dict, dictkey, legend_str, fitline, arrow, lims, sizes, style)
#
#legend_str = "Narrow 0.1"
#dictkey = "d250stdd25"
#style = "rv"
#arrow = False
#fitline = False
#rslt_plot(ax, rslts_dict, dictkey, legend_str, fitline, arrow, lims, sizes, style)
#
#legend_str = "Medium"
#dictkey = "d300stdd100"
#style = "g^"
#arrow = False
#fitline = False
#rslt_plot(ax, rslts_dict, dictkey, legend_str, fitline, arrow, lims, sizes, style)
#
#legend_str = "Medium 0.4"
#dictkey = "d271stdd121"
#style = "g^"
#arrow = False
#fitline = False
#rslt_plot(ax, rslts_dict, dictkey, legend_str, fitline, arrow, lims, sizes, style)
#
#legend_str = "Wide"
#dictkey = "d300stdd200"
#style = "bs"
#arrow = False
#fitline = False
#rslt_plot(ax, rslts_dict, dictkey, legend_str, fitline, arrow, lims, sizes, style)
#
#legend_str = "Wide 0.7"
#dictkey = "d317stdd252"
#style = "bs"
#arrow = False
#fitline = False
#rslt_plot(ax, rslts_dict, dictkey, legend_str, fitline, arrow, lims, sizes, style)
#
#legend_str = "Very wide"
#dictkey = "d300stdd300"
#style = "ch"
#arrow = True
#fitline = True
#rslt_plot(ax, rslts_dict, dictkey, legend_str, fitline, arrow, lims, sizes, style)
#
#legend_str = "Very wide 1.0"
#dictkey = "d347stdd537"
#style = "ch"
#arrow = False
#fitline = False
#rslt_plot(ax, rslts_dict, dictkey, legend_str, fitline, arrow, lims, sizes, style)
#
#legend_str = "Right Truncation"
#dictkey = "d269stdd100"
#style = "y*"
#arrow = False
#fitline = False
#rslt_plot(ax, rslts_dict, dictkey, legend_str, fitline, arrow, lims, sizes, style)
#
#legend_str = "Left Truncation"
#dictkey = "d321stdd100"
#style = "k*"
#arrow = False
#fitline = False
#rslt_plot(ax, rslts_dict, dictkey, legend_str, fitline, arrow, lims, sizes, style)
#
#legend_str = "Large Mu"
#dictkey = "d430stdd100"
#style = "cx"
#arrow = False
#fitline = False
#rslt_plot(ax, rslts_dict, dictkey, legend_str, fitline, arrow, lims, sizes, style)
#
#legend_str = "Small Mu"
#dictkey = "d167stdd100"
#style = "kx"
#arrow = False
#fitline = False
#rslt_plot(ax, rslts_dict, dictkey, legend_str, fitline, arrow, lims, sizes, style)


if use_lims:
    ax.set_xlim(lims["x_min"], lims["x_max"])
    ax.set_ylim(lims["y_min"], lims["y_max"])
ax.set_xlabel(x_label, fontsize=label_size)
ax.set_ylabel(y_label, fontsize=label_size)
ax.tick_params(axis='both', labelsize=ticks_size)
ax.legend(fontsize=ticks_size)

plt.show()