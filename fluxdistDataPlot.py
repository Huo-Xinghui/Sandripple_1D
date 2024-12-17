# -*-*-*- Author: EkalHxH -*-*-*-
# version: 1.0 (2024-12-17)
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

"""读取数据并画床面通量分布图"""

# ********************************************************************************************************

# 导入必要的库
import os # 用于文件操作
import re # 用于正则表达式
import matplotlib.pyplot as plt # 用于绘图
import numpy as np # 用于数值计算
from dataclasses import dataclass # 用于定义数据类
from typing import List # 用于定义列表
from typing import Dict # 用于定义字典

# 定义一个数据类，用于存储文件中的数据
@dataclass
class file_data:
	data: List[float]

# 定义一个函数来读取任意文件
def read_file(file_path) -> List[str]:
	if not os.path.exists(file_path):
		print(f"文件 {file_path} 不存在。")
		exit()

	with open(file_path, 'r', encoding='utf-8') as file:
		line = file.readlines()
		return line

# 读取廓线文件内容, 并将其存储到一个字典中
def read_data_file1(folder_path, file_name) -> Dict[int, List[file_data]]:
	file_path = os.path.join(folder_path, file_name) # 拼接文件路径
	all_lines = read_file(file_path) # 读取文件内容

	time_step_data_dict: Dict[int, List[file_data]] = {} # 定义一个字典，用于存储数据
	for line in all_lines:
		if "vs" in line:
			nums_in_first_line = [int(num) for num in re.findall(r'\d+\.?\d*', line)]
			current_time = nums_in_first_line[0]
			point_num = nums_in_first_line[1]
			data_list = []
			point_count = 0
		else:
			point_count += 1
			columns = line.split()
			columns_float = [float(column) for column in columns]
			current_data = file_data(
				data=columns_float
			)
			data_list.append(current_data)
			if point_count == point_num:
				time_step_data_dict[current_time] = data_list
	return time_step_data_dict

# 主程序
if __name__ == "__main__":
	# 判断操作系统
	sys_OS = "w" # "w" for windows, "l" for linux
	if sys_OS == "l":
		linux_flag = True
	elif sys_OS == "w":
		linux_flag = False
	else:
		print("Invalid input!")
		exit()

# ----------------------------------------------------------------------------------------
	# 固定参数
	nu = 1.51e-5 # 运动粘度
	interval = 30 # 源文件时间间隔
	# 控制参数
	output_num = 2 # 出图类型：0为床面廓线，1为相关性，2为波长, 3为波高, 4为波速, 5为床面粒径分布
	x_type = 0 # 算例对比图中x轴的类型：0为u*, 1为d_min, 2为d_max, 3为d, 4为stddev, 5为Sh, 6为Ga
	start = 30 # 时空图的起始时间或者变量随时间变化的起始时间
	end = 3000 # 时空图的结束时间或者变量随时间变化的结束时间
	average_start = 300 # 开始计算平均值的时间
	average_end = 600 # 结束计算平均值的时间
	profile_offset = 2e-4 # 廓线图的纵向偏移
	corr_offset = 1e-8 # 相关性图的纵向偏移
	diameter_offset = 2e-5 # 廓线图的纵向偏移
	case_num = 36 # 算例号，用于画时空图
# ----------------------------------------------------------------------------------------

	# 定义文件路径
	if linux_flag:
		working_dir = "/home/ekalhxh/ripple/coll1"
	else:
		working_dir = "E:/Data/Sandripples1DFluid/ripple/coll1"

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
		30: "uStar035_300_0_2650_3600",
		31: "uStar040_300_0_2650_3600",
		32: "uStar045_300_0_2650_3600",
		33: "uStar050_300_0_2650_3600",
		34: "uStar055_300_0_2650_3600",
		35: "uStar060_300_0_2650_3600",
		36: "uStar065_300_0_2650_3600",
		37: "uStar035_300stdd100_0_2650_3600",
		38: "uStar040_300stdd100_0_2650_3600",
		39: "uStar045_300stdd100_0_2650_3600",
		40: "uStar050_300stdd100_0_2650_3600",
		41: "uStar055_300stdd100_0_2650_3600",
		42: "uStar060_300stdd100_0_2650_3600",
		43: "uStar065_300stdd100_0_2650_3600",
	}