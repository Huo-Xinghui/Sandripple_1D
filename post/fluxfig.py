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
import matplotlib.pyplot as plt # 导入绘图库
import numpy as np
import matplotlib as mpl
import matplotlib.patches as patches
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import ScalarFormatter
from matplotlib.legend_handler import HandlerTuple
from scipy.optimize import curve_fit

# Set up matplotlib parameters
mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = r'\usepackage{amssymb}'

def rslt_plot(ax, rslts_dict, dictkey, legend_str, fitline, arrow, lims, sizes, style, axis_type):
	x = rslts_dict[dictkey]["x"]
	y = rslts_dict[dictkey]["y"]
	x_max = lims["x_max"]
	if axis_type == "loglog":
		ax.set_xscale('log')
		ax.set_yscale('log')
	if style["fill"] == "full":
		ax.scatter(x, y,
			color=style["color"],
			marker=style["marker"],
			s=sizes["marker"]**2,
			label=legend_str)
	else:
		ax.scatter(x, y,
			color=style["color"],
			marker=style["marker"],
			s=sizes["marker"]**2,
			linewidths=sizes["markerwidth"],
			label=legend_str,
			facecolors='none')
	if fitline:
		coeffs = rslts_dict[dictkey]["fit"]
		fit_line = np.poly1d(coeffs)
		# 计算y=0时的x值（x轴截距）
		x_intercept = -coeffs[1] / coeffs[0]
		# 扩展x轴范围到x轴截距
		x_fit = np.linspace(min(min(x), x_intercept), max(x_max, x_intercept), 10000)
		# 在x轴截距处添加标注
		x_loc = x_intercept * 1.1
		y_loc = fit_line(x_loc)
		if arrow:
			plt.annotate(f'$S_d$ = {x_intercept:.4f}',
				xy=(x_loc, y_loc),
				ha='center',
				va='top',
				xytext=(lims["x_loc"], lims["y_loc"]),
				fontsize=sizes["ticks"],
				arrowprops=dict(arrowstyle='->'))
		ax.plot(x_fit, fit_line(x_fit), '--', color=style["color"], linewidth=sizes["line"])

def Pahtz_Q_func(x, Sd, cM, mub):
	"""S, Sd"""
	kapa = 0.4
	cM = 1.7
	Q = 2.0*np.sqrt(Sd)/(kapa*mub)*(x - Sd)*(1.0 + cM/mub*(x - Sd))
	return Q

# 控制参数
ampf = 1/0.6 # amplification factor
label_size = 15*ampf
ticks_size = 10*ampf
ticks_size_in = 6*ampf
marker_size = 8*ampf
marker_size_in = 6*ampf
line_width = 1.5
marker_width = 2
sizes = {
	"marker": marker_size,
	"line": line_width,
	"ticks": ticks_size,
	"markerwidth": marker_width,
}
sizes_in = {
	"marker": marker_size_in,
	"line": line_width,
	"markerwidth": marker_width
}
use_lims = False # 是否使用坐标轴限制
lims = {
	"x_min": 1e-4,
	"x_max": 0.4,
	"y_min": 3e-3,
	"y_max": 0.2,
	"x_loc": 0.1,
	"y_loc": 0.1
}
lims_in = {
	"x_min": 1e-2,
	"x_max": 0.2,
	"y_min": 3e-3,
	"y_max": 0.2,
}
tailin = "dair"
tailout = "d90"
working_dir = "E:/Data/Q_on_flat_bed"
dmin = 1.0e-4
dmax = 10.0e-4
bin_num = 1000
dx = dmax / bin_num
ddx = dx * 0.5
mu_list0 = [-8.1254, -8.1644, -8.3221, -8.5521]
sigma_list0 = [0.1655, 0.3246, 0.6167, 0.87]
mu_list1 = [-8.3, -8.3, -8.3, -8.3]
sigma_list1 = [0.1, 0.4, 0.7, 1.0]

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

# Ho 2012 的数据
d_ho12 = 2.30e-4*1.63
rhop_ho12 = 2450
rhof_ho12 = 1.204
s_ho12 = rhop_ho12 / rhof_ho12
g_ho12 = 9.81
ghat_ho12 = (1.0 - 1.0/s_ho12) * g_ho12
Sh_array = np.array([0.01776, 0.04970, 0.08710, 0.13516, 0.17804, 0.32099])
Q_array = np.array([0.35012, 1.12710, 2.34532, 3.64988, 5.59712, 9.11271])
#Sh_array = np.array([0.01776, 0.04970, 0.08710, 0.13516, 0.17804])
#Q_array = np.array([0.35012, 1.12710, 2.34532, 3.64988, 5.59712])
Q_array = rhop_ho12 * np.sqrt(g_ho12*d_ho12) * d_ho12 * Q_array
ud_ho12 = 0.17
Shd_ho12 = rhof_ho12*ud_ho12**2 / (rhop_ho12*ghat_ho12*d_ho12)
Ho2012_data = {
	"Sh": Sh_array,
	"Q": Q_array,
	"Q_err": 0.05 * Q_array,
	"Q_star": Q_array/(rhop_ho12*d_ho12*np.sqrt(s_ho12*ghat_ho12*d_ho12)),
	"Q_star_err": 0.05 * Q_array,
	"Shd": Shd_ho12,
	"ud_star": ud_ho12
}

# Ho 2012 coarse 的数据
d_ho12 = 6.30e-4*2.01
rhop_ho12 = 2450
rhof_ho12 = 1.204
s_ho12 = rhop_ho12 / rhof_ho12
g_ho12 = 9.81
ghat_ho12 = (1.0 - 1.0/s_ho12) * g_ho12
Sh_array = np.array([0.02562, 0.02958, 0.03992, 0.04560, 0.05019, 0.05115, 0.06251, 0.07904, 0.10017])
Q_array = np.array([0.5800, 0.6744, 1.2739, 1.7121, 1.9586, 2.1716, 2.4790, 2.8351, 4.0949])
Q_array = rhop_ho12 * np.sqrt(g_ho12*d_ho12) * d_ho12 * Q_array
ud_ho12c = 0.28
Shd_ho12c = rhof_ho12*ud_ho12c**2 / (rhop_ho12*ghat_ho12*d_ho12)
Ho2012_c_data = {
	"Sh": Sh_array,
	"Q": Q_array,
	"Q_err": 0.05 * Q_array,
	"Q_star": Q_array/(rhop_ho12*d_ho12*np.sqrt(s_ho12*ghat_ho12*d_ho12)),
	"Q_star_err": 0.05 * Q_array,
	"Shd": Shd_ho12c,
	"ud_star": ud_ho12c
}

# Zhu 2019 的数
d50_zh19 = np.array([187, 251, 351, 438, 190, 529])*1e-6
d90_zh19 = np.array([240, 350, 489, 535, 424, 1300])*1e-6
#d50_zh19 = np.array([187, 251, 351, 438])*1e-6
#d90_zh19 = np.array([240, 350, 489, 535])*1e-6
d9050 = d90_zh19 / d50_zh19
ud_zh19 = np.array([16.2, 17.4, 17.6, 18.7, 27.4, 65.2])*1e-2
#ud_zh19 = np.array([16.2, 17.4, 17.6, 18.7])*1e-2
rhop_zh19 = 2650
rhof_zh19 = 1.2
s_zh19 = rhop_zh19 / rhof_zh19
g_zh19 = 9.81
ghat_zh19 = (1.0 - 1.0/s_zh19) * g_zh19
Shd_zh19 = rhof_zh19*ud_zh19**2 / (rhop_zh19*ghat_zh19*d90_zh19)
Zh19_data = {
	"d50": d50_zh19,
	"d90": d90_zh19,
	"d9050": d9050,
	"ud_star": ud_zh19,
	"Shd": Shd_zh19
}

# Martin 2018的数据
d50_ma18 = np.array([526, 533, 398])*1e-6
d90_ma18 = np.array([847, 839, 650])*1e-6
d9050_ma18 = d90_ma18 / d50_ma18
tau_ma18 = np.array([0.111, 0.110, 0.088])
ud_ma18 = np.array([0.309, 0.300, 0.269])
rhof_ma18 = tau_ma18 / ud_ma18**2
rhop_ma18 = 2650
s_ma18 = rhop_ma18 / rhof_ma18
g_ma18 = 9.81
ghat_ma18 = (1.0 - 1.0/s_ma18) * g_ma18
Shd_ma18 = rhof_ma18*ud_ma18**2 / (rhop_ma18*ghat_ma18*d90_ma18)
Ma18_data = {
	"d50": d50_ma18,
	"d90": d90_ma18,
	"d9050": d9050_ma18,
	"tau": tau_ma18,
	"ud": ud_ma18,
	"rhof": rhof_ma18,
	"rhop": rhop_ma18,
	"Shd": Shd_ma18
}

# Campmans & Wijnberg 2022 DEM data
d50_ca22 = 2.5e-4
rhop_ca22 = 2700
rhof_ca22 = 1.225
s_ca22 = rhop_ca22 / rhof_ca22
g_ca22 = 9.81
ghat_ca22 = (1.0 - 1.0/s_ca22) * g_ca22
mu_ca22 = np.log(d50_ca22)
sigma_phi_ca22 = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
sigma_ca22 = sigma_phi_ca22 * np.log(2)
sd_ca22 = np.exp(mu_ca22 + 0.5*sigma_ca22**2) * np.sqrt(np.exp(sigma_ca22**2) - 1)
dm_ca22 = d50_ca22 * np.exp(0.5*sigma_ca22**2)
Sh_array1 = np.array([0.02, 0.03, 0.04, 0.05, 0.06, 0.07])
Q_array1 = np.array([0.01000, 0.01318, 0.02206, 0.02008, 0.03000, 0.04106])
Q_array1 = Q_array1*np.sqrt((s_ca22 - 1)*g_ca22*d50_ca22**3)*rhop_ca22
Q_array1 = Q_array1/(rhop_ca22*d50_ca22*np.sqrt(s_ca22*ghat_ca22*d50_ca22))
Sh_array2 = np.array([0.03, 0.04, 0.05, 0.06, 0.07])
Q_array2 = np.array([0.01053, 0.01781, 0.02240, 0.02765, 0.03479])
Q_array2 = Q_array2*np.sqrt((s_ca22 - 1)*g_ca22*d50_ca22**3)*rhop_ca22
Q_array2 = Q_array2/(rhop_ca22*d50_ca22*np.sqrt(s_ca22*ghat_ca22*d50_ca22))
Sh_array3 = np.array([0.02, 0.04, 0.05, 0.06, 0.07])
Q_array3 = np.array([0.00928, 0.01593, 0.02211, 0.02963, 0.03132])
Q_array3 = Q_array3*np.sqrt((s_ca22 - 1)*g_ca22*d50_ca22**3)*rhop_ca22
Q_array3 = Q_array3/(rhop_ca22*d50_ca22*np.sqrt(s_ca22*ghat_ca22*d50_ca22))
Sh_array4 = np.array([0.02, 0.03, 0.04, 0.05, 0.06, 0.07])
Q_array4 = np.array([0.01082, 0.01318, 0.01820, 0.02051, 0.02317, 0.03108])
Q_array4 = Q_array4*np.sqrt((s_ca22 - 1)*g_ca22*d50_ca22**3)*rhop_ca22
Q_array4 = Q_array4/(rhop_ca22*d50_ca22*np.sqrt(s_ca22*ghat_ca22*d50_ca22))
Sh_array5 = np.array([0.02, 0.03, 0.04, 0.05, 0.06, 0.07])
Q_array5 = np.array([0.00778, 0.01222, 0.01434, 0.01984, 0.02572, 0.03064])
Q_array5 = Q_array5*np.sqrt((s_ca22 - 1)*g_ca22*d50_ca22**3)*rhop_ca22
Q_array5 = Q_array5/(rhop_ca22*d50_ca22*np.sqrt(s_ca22*ghat_ca22*d50_ca22))
Sh_array6 = np.array([0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07])
Q_array6 = np.array([0.00441, 0.00860, 0.01164, 0.01608, 0.01931, 0.02408, 0.03088])
Q_array6 = Q_array6*np.sqrt((s_ca22 - 1)*g_ca22*d50_ca22**3)*rhop_ca22
Q_array6 = Q_array6/(rhop_ca22*d50_ca22*np.sqrt(s_ca22*ghat_ca22*d50_ca22))
Ca22_data = {
	"Sh0": Sh_array1,
	"Q_star0": Q_array1,
	"Sh1": Sh_array2,
	"Q_star1": Q_array2,
	"Sh2": Sh_array3,
	"Q_star2": Q_array3,
	"Sh3": Sh_array4,
	"Q_star3": Q_array4,
	"Sh4": Sh_array5,
	"Q_star4": Q_array5,
	"Sh5": Sh_array6,
	"Q_star5": Q_array6
}

# 定义参照量字典
others_case_dict = {
	"Cr09": Cr09_data,
	"Ho2012": Ho2012_data,
	"Ho2012_c": Ho2012_c_data
}

Q_nr1_in = np.load(f"Q_d300stdd50_{tailin}.npz")
M_nr1_in = np.load(f"M_d300stdd50_{tailin}.npz")
Q_md1_in = np.load(f"Q_d300stdd100_{tailin}.npz")
M_md1_in = np.load(f"M_d300stdd100_{tailin}.npz")
Q_wd1_in = np.load(f"Q_d300stdd200_{tailin}.npz")
M_wd1_in = np.load(f"M_d300stdd200_{tailin}.npz")
Q_vwd1_in = np.load(f"Q_d300stdd300_{tailin}.npz")
M_vwd1_in = np.load(f"M_d300stdd300_{tailin}.npz")
Q_nr1 = np.load(f"Q_d300stdd50_{tailout}.npz")
M_nr1 = np.load(f"M_d300stdd50_{tailout}.npz")
Q_md1 = np.load(f"Q_d300stdd100_{tailout}.npz")
M_md1 = np.load(f"M_d300stdd100_{tailout}.npz")
Q_wd1 = np.load(f"Q_d300stdd200_{tailout}.npz")
M_wd1 = np.load(f"M_d300stdd200_{tailout}.npz")
Q_vwd1 = np.load(f"Q_d300stdd300_{tailout}.npz")
M_vwd1 = np.load(f"M_d300stdd300_{tailout}.npz")
Q_nr2 = np.load(f"Q_d250stdd25_{tailin}.npz")
M_nr2 = np.load(f"M_d250stdd25_{tailin}.npz")
Q_md2 = np.load(f"Q_d271stdd121_{tailin}.npz")
M_md2 = np.load(f"M_d271stdd121_{tailin}.npz")
Q_wd2 = np.load(f"Q_d317stdd252_{tailin}.npz")
M_wd2 = np.load(f"M_d317stdd252_{tailin}.npz")
Q_vwd2 = np.load(f"Q_d347stdd537_{tailin}.npz")
M_vwd2 = np.load(f"M_d347stdd537_{tailin}.npz")
Q_nr3 = np.load(f"Q_d240stdd50_{tailin}.npz")
M_nr3 = np.load(f"M_d240stdd50_{tailin}.npz")
Q_nr4 = np.load(f"Q_d400stdd50_{tailin}.npz")
M_nr4 = np.load(f"M_d400stdd50_{tailin}.npz")
Q_nr2_out = np.load(f"Q_d250stdd25_{tailout}.npz")
M_nr2_out = np.load(f"M_d250stdd25_{tailout}.npz")
Q_md2_out = np.load(f"Q_d271stdd121_{tailout}.npz")
M_md2_out = np.load(f"M_d271stdd121_{tailout}.npz")
Q_wd2_out = np.load(f"Q_d317stdd252_{tailout}.npz")
M_wd2_out = np.load(f"M_d317stdd252_{tailout}.npz")
Q_vwd2_out = np.load(f"Q_d347stdd537_{tailout}.npz")
M_vwd2_out = np.load(f"M_d347stdd537_{tailout}.npz")
Q_nr3_in = np.load(f"Q_mono_d250stdd25_{tailin}.npz")
M_nr3_in = np.load(f"M_mono_d250stdd25_{tailin}.npz")
Q_md3_in = np.load(f"Q_mono_d271stdd121_{tailin}.npz")
M_md3_in = np.load(f"M_mono_d271stdd121_{tailin}.npz")
Q_wd3_in = np.load(f"Q_mono_d317stdd252_{tailin}.npz")
M_wd3_in = np.load(f"M_mono_d317stdd252_{tailin}.npz")
Q_vwd3_in = np.load(f"Q_mono_d347stdd537_{tailin}.npz")
M_vwd3_in = np.load(f"M_mono_d347stdd537_{tailin}.npz")
Q_nr4_in = np.load(f"Q_mono_d300stdd50_{tailin}.npz")
M_nr4_in = np.load(f"M_mono_d300stdd50_{tailin}.npz")
Q_md4_in = np.load(f"Q_mono_d300stdd100_{tailin}.npz")
M_md4_in = np.load(f"M_mono_d300stdd100_{tailin}.npz")
Q_wd4_in = np.load(f"Q_mono_d300stdd200_{tailin}.npz")
M_wd4_in = np.load(f"M_mono_d300stdd200_{tailin}.npz")
Q_vwd4_in = np.load(f"Q_mono_d300stdd300_{tailin}.npz")
M_vwd4_in = np.load(f"M_mono_d300stdd300_{tailin}.npz")
Q_nr5_in = np.load(f"Q_Eeff_d300stdd50_{tailin}.npz")
M_nr5_in = np.load(f"M_Eeff_d300stdd50_{tailin}.npz")
Q_md5_in = np.load(f"Q_Eeff_d300stdd100_{tailin}.npz")
M_md5_in = np.load(f"M_Eeff_d300stdd100_{tailin}.npz")
Q_wd5_in = np.load(f"Q_Eeff_d300stdd200_{tailin}.npz")
M_wd5_in = np.load(f"M_Eeff_d300stdd200_{tailin}.npz")
Q_vwd5_in = np.load(f"Q_Eeff_d300stdd300_{tailin}.npz")
M_vwd5_in = np.load(f"M_Eeff_d300stdd300_{tailin}.npz")
Q_nr5_out = np.load(f"Q_Eeff_d300stdd50_{tailout}.npz")
M_nr5_out = np.load(f"M_Eeff_d300stdd50_{tailout}.npz")
Q_md5_out = np.load(f"Q_Eeff_d300stdd100_{tailout}.npz")
M_md5_out = np.load(f"M_Eeff_d300stdd100_{tailout}.npz")
Q_wd5_out = np.load(f"Q_Eeff_d300stdd200_{tailout}.npz")
M_wd5_out = np.load(f"M_Eeff_d300stdd200_{tailout}.npz")
Q_vwd5_out = np.load(f"Q_Eeff_d300stdd300_{tailout}.npz")
M_vwd5_out = np.load(f"M_Eeff_d300stdd300_{tailout}.npz")
Q_nr6_in = np.load(f"Q_mono_d150stdd50_{tailin}.npz")
M_nr6_in = np.load(f"M_mono_d150stdd50_{tailin}.npz")
Q_md6_in = np.load(f"Q_mono_d150stdd100_{tailin}.npz")
M_md6_in = np.load(f"M_mono_d150stdd100_{tailin}.npz")
Q_wd6_in = np.load(f"Q_mono_d150stdd200_{tailin}.npz")
M_wd6_in = np.load(f"M_mono_d150stdd200_{tailin}.npz")
Q_vwd6_in = np.load(f"Q_mono_d150stdd300_{tailin}.npz")
M_vwd6_in = np.load(f"M_mono_d150stdd300_{tailin}.npz")
Q_mu1_out = np.load(f"Q_d167stdd100_{tailout}.npz")
M_mu1_out = np.load(f"M_d167stdd100_{tailout}.npz")
Q_mu2_out = np.load(f"Q_d197stdd65_{tailout}.npz")
M_mu2_out = np.load(f"M_d197stdd65_{tailout}.npz")
Q_mu3_out = np.load(f"Q_d290stdd97_{tailout}.npz")
M_mu3_out = np.load(f"M_d290stdd97_{tailout}.npz")
Q_mu4_out = np.load(f"Q_d430stdd100_{tailout}.npz")
M_mu4_out = np.load(f"M_d430stdd100_{tailout}.npz")
Q_rt_out = np.load(f"Q_d269stdd100_{tailout}.npz")
M_rt_out = np.load(f"M_d269stdd100_{tailout}.npz")
Q_lt_out = np.load(f"Q_d321stdd100_{tailout}.npz")
M_lt_out = np.load(f"M_d321stdd100_{tailout}.npz")
Q_nr7_out = np.load(f"Q_d240stdd50_{tailout}.npz")
M_nr7_out = np.load(f"M_d240stdd50_{tailout}.npz")
Q_nr8_out = np.load(f"Q_d400stdd50_{tailout}.npz")
M_nr8_out = np.load(f"M_d400stdd50_{tailout}.npz")

dir_list_nr = [
	"uStar030_300log50_0_2650_300",
	"uStar035_300log50_0_2650_300",
	"uStar040_300log50_0_2650_300",
	"uStar045_300log50_0_2650_300",
	"uStar050_300log50_0_2650_300",
	"uStar055_300log50_0_2650_300",
	"uStar060_300log50_0_2650_300",
	"uStar065_300log50_0_2650_300"
]


dir_list_vwd = [
	"uStar030_300log300_0_2650_300",
	"uStar035_300log300_0_2650_300",
	"uStar040_300log300_0_2650_300",
	"uStar045_300log300_0_2650_300",
	"uStar050_300log300_0_2650_300",
	"uStar055_300log300_0_2650_300",
	"uStar060_300log300_0_2650_300",
	"uStar065_300log300_0_2650_300",
]


dir_list_md = [
	"uStar030_300log100_0_2650_300",
	"uStar035_300log100_0_2650_300",
	"uStar040_300log100_0_2650_300",
	"uStar045_300log100_0_2650_300",
	"uStar050_300log100_0_2650_300",
	"uStar055_300log100_0_2650_300",
	"uStar060_300log100_0_2650_300",
	"uStar065_300log100_0_2650_300",
]

dir_list_wd = [
	"uStar030_300log200_0_2650_300",
	"uStar035_300log200_0_2650_300",
	"uStar040_300log200_0_2650_300",
	"uStar045_300log200_0_2650_300",
	"uStar050_300log200_0_2650_300",
	"uStar055_300log200_0_2650_300",
	"uStar060_300log200_0_2650_300",
	"uStar065_300log200_0_2650_300",
]

dir_list_nr1 = [
	"uStar030_250log25_0_2650_300",
	"uStar040_250log25_0_2650_300",
	"uStar050_250log25_0_2650_300",
	"uStar060_250log25_0_2650_300",
	"uStar035_250log25_0_2650_300",
	"uStar045_250log25_0_2650_300",
	"uStar055_250log25_0_2650_300",
	"uStar065_250log25_0_2650_300"
]

dir_list_md1 = [
	"uStar030_271log121_0_2650_300",
	"uStar040_271log121_0_2650_300",
	"uStar050_271log121_0_2650_300",
	"uStar060_271log121_0_2650_300",
	"uStar035_271log121_0_2650_300",
	"uStar045_271log121_0_2650_300",
	"uStar055_271log121_0_2650_300",
	"uStar065_271log121_0_2650_300"
]

dir_list_wd1 = [
	"uStar030_317log252_0_2650_300",
	"uStar040_317log252_0_2650_300",
	"uStar050_317log252_0_2650_300",
	"uStar060_317log252_0_2650_300",
	"uStar035_317log252_0_2650_300",
	"uStar045_317log252_0_2650_300",
	"uStar055_317log252_0_2650_300",
	"uStar065_317log252_0_2650_300"
]

dir_list_vwd1 = [
	"uStar030_347log537_0_2650_300",
	"uStar040_347log537_0_2650_300",
	"uStar050_347log537_0_2650_300",
	"uStar060_347log537_0_2650_300",
	"uStar035_347log537_0_2650_300",
	"uStar045_347log537_0_2650_300",
	"uStar055_347log537_0_2650_300",
	"uStar065_347log537_0_2650_300",
]

dir_list_mu4 = [
	"uStar040_430log100_0_2650_300",
	"uStar050_430log100_0_2650_300",
	"uStar060_430log100_0_2650_300",
]

dir_list_mu1 = [
	"uStar030_167log100_0_2650_300",
	"uStar040_167log100_0_2650_300",
	"uStar050_167log100_0_2650_300",
]

dir_list_rt = [
	"uStar030_269log100_0_2650_300",
	"uStar040_269log100_0_2650_300",
	"uStar050_269log100_0_2650_300",
	"uStar035_269log100_0_2650_300",
	"uStar045_269log100_0_2650_300",
	"uStar055_269log100_0_2650_300",
]

dir_list_lt = [
	"uStar030_321log100_0_2650_300",
	"uStar040_321log100_0_2650_300",
	"uStar050_321log100_0_2650_300",
]

dir_list_nr3 = [
	"uStar030_240log50_0_2650_300",
	"uStar040_240log50_0_2650_300",
	"uStar050_240log50_0_2650_300",
	"uStar035_240log50_0_2650_300",
	"uStar045_240log50_0_2650_300",
	"uStar055_240log50_0_2650_300",
]

dir_list_nr4 = [
	"uStar040_400log50_0_2650_300",
	"uStar050_400log50_0_2650_300",
	"uStar060_400log50_0_2650_300"
]

dir_list_mu3 = [
	"uStar030_290log97_0_2650_300",
	"uStar040_290log97_0_2650_300",
	"uStar050_290log97_0_2650_300",
	"uStar060_290log97_0_2650_300"
]

dir_list_mu2 = [
	"uStar030_197log65_0_2650_300",
	"uStar040_197log65_0_2650_300",
	"uStar050_197log65_0_2650_300",
	"uStar060_197log65_0_2650_300"
]

cases_dict = {
	"NR1": dir_list_nr,
	"MD1": dir_list_md,
	"WD1": dir_list_wd,
	"VWD1": dir_list_vwd,
	"NR2": dir_list_nr1,
	"MD2": dir_list_md1,
	"WD2": dir_list_wd1,
	"VWD2": dir_list_vwd1,
	"MU4": dir_list_mu4,
	"MU1": dir_list_mu1,
	"RT": dir_list_rt,
	"LT": dir_list_lt,
	"NR3": dir_list_nr3,
	"NR4": dir_list_nr4,
	"MU3": dir_list_mu3,
	"MU2": dir_list_mu2
}

d_hist_addr_nr = f"{working_dir}/{dir_list_nr[-1]}/d_in_air_hist.npz"
d_hist_addr_md = f"{working_dir}/{dir_list_md[-1]}/d_in_air_hist.npz"
d_hist_addr_wd = f"{working_dir}/{dir_list_wd[-1]}/d_in_air_hist.npz"
d_hist_addr_vwd = f"{working_dir}/{dir_list_vwd[-1]}/d_in_air_hist.npz"
d_hist_addr_nr1 = f"{working_dir}/{dir_list_nr1[-1]}/d_in_air_hist.npz"
d_hist_addr_md1 = f"{working_dir}/{dir_list_md1[-1]}/d_in_air_hist.npz"
d_hist_addr_wd1 = f"{working_dir}/{dir_list_wd1[-1]}/d_in_air_hist.npz"
d_hist_addr_vwd1 = f"{working_dir}/{dir_list_vwd1[-1]}/d_in_air_hist.npz"
d_hist_nr_dict = np.load(d_hist_addr_nr)
d_hist_md_dict = np.load(d_hist_addr_md)
d_hist_wd_dict = np.load(d_hist_addr_wd)
d_hist_vwd_dict = np.load(d_hist_addr_vwd)
d_hist_nr1_dict = np.load(d_hist_addr_nr1)
d_hist_md1_dict = np.load(d_hist_addr_md1)
d_hist_wd1_dict = np.load(d_hist_addr_wd1)
d_hist_vwd1_dict = np.load(d_hist_addr_vwd1)
d_hist_nr = d_hist_nr_dict["hist"]
d_hist_md = d_hist_md_dict["hist"]
d_hist_wd = d_hist_wd_dict["hist"]
d_hist_vwd = d_hist_vwd_dict["hist"]
d_hist_nr1 = d_hist_nr1_dict["hist"]
d_hist_md1 = d_hist_md1_dict["hist"]
d_hist_wd1 = d_hist_wd1_dict["hist"]
d_hist_vwd1 = d_hist_vwd1_dict["hist"]
d_bin_centers_nr = d_hist_nr_dict["bin_centers"]
d_bin_centers_md = d_hist_md_dict["bin_centers"]
d_bin_centers_wd = d_hist_wd_dict["bin_centers"]
d_bin_centers_vwd = d_hist_vwd_dict["bin_centers"]
d_bin_centers_nr1 = d_hist_nr1_dict["bin_centers"]
d_bin_centers_md1 = d_hist_md1_dict["bin_centers"]
d_bin_centers_wd1 = d_hist_wd1_dict["bin_centers"]
d_bin_centers_vwd1 = d_hist_vwd1_dict["bin_centers"]
d_bin_edges_nr = d_hist_nr_dict["bin_edges"]
d_bin_edges_md = d_hist_md_dict["bin_edges"]
d_bin_edges_wd = d_hist_wd_dict["bin_edges"]
d_bin_edges_vwd = d_hist_vwd_dict["bin_edges"]
d_bin_edges_nr1 = d_hist_nr1_dict["bin_edges"]
d_bin_edges_md1 = d_hist_md1_dict["bin_edges"]
d_bin_edges_wd1 = d_hist_wd1_dict["bin_edges"]
d_bin_edges_vwd1 = d_hist_vwd1_dict["bin_edges"]
dm_nr = d_hist_nr_dict["d_mean"]
dm_md = d_hist_md_dict["d_mean"]
dm_wd = d_hist_wd_dict["d_mean"]
dm_vwd = d_hist_vwd_dict["d_mean"]
dm_nr1 = d_hist_nr1_dict["d_mean"]
dm_md1 = d_hist_md1_dict["d_mean"]
dm_wd1 = d_hist_wd1_dict["d_mean"]
dm_vwd1 = d_hist_vwd1_dict["d_mean"]
x = np.linspace(dmin+dx, dmax, bin_num)
x = x - ddx

center_list0 = []
hist_list0 = []
for mu, sigma in zip(mu_list0, sigma_list0):
	# 对数正态分布概率密度函数
	y = (1 / (x * sigma * np.sqrt(2 * np.pi))) * np.exp(-(np.log(x) - mu) ** 2 / (2 * sigma ** 2))
	# 归一化
	y = y / y.sum()
	center_list0.append(x)
	hist_list0.append(y)

center_list1 = []
hist_list1 = []
for mu, sigma in zip(mu_list1, sigma_list1):
	# 对数正态分布概率密度函数
	y = (1 / (x * sigma * np.sqrt(2 * np.pi))) * np.exp(-(np.log(x) - mu) ** 2 / (2 * sigma ** 2))
	# 归一化
	y = y / y.sum()
	center_list1.append(x)
	hist_list1.append(y)

d_in_air_dict = {}
for case_name, case_dir in cases_dict.items():
	dm_list = []
	d50_list = []
	d90_list = []
	for i in range(len(case_dir)):
		folder = case_dir[i]
		d_hist_nr_dict = np.load(f"{working_dir}/{folder}/d_in_air_hist.npz")
		dm_list.append(d_hist_nr_dict["d_mean"])
		d50_list.append(d_hist_nr_dict["d50"])
		d90_list.append(d_hist_nr_dict["d90"])
	d_in_air_dict[case_name] = {
		"dm": dm_list,
		"d50": d50_list,
		"d90": d90_list
	}

Q_rslt = {
	f"d300stdd50_{tailin}": Q_nr1_in,
	f"d300stdd100_{tailin}": Q_md1_in,
	f"d300stdd200_{tailin}": Q_wd1_in,
	f"d300stdd300_{tailin}": Q_vwd1_in,
	f"d300stdd50_{tailout}": Q_nr1,
	f"d300stdd100_{tailout}": Q_md1,
	f"d300stdd200_{tailout}": Q_wd1,
	f"d300stdd300_{tailout}": Q_vwd1,
	f"d250stdd25_{tailin}": Q_nr2,
	f"d271stdd121_{tailin}": Q_md2,
	f"d317stdd252_{tailin}": Q_wd2,
	f"d347stdd537_{tailin}": Q_vwd2,
	f"d250stdd25_{tailout}": Q_nr2_out,
	f"d271stdd121_{tailout}": Q_md2_out,
	f"d317stdd252_{tailout}": Q_wd2_out,
	f"d347stdd537_{tailout}": Q_vwd2_out,
	f"d250stdd25_m_{tailin}": Q_nr3_in,
	f"d271stdd121_m_{tailin}": Q_md3_in,
	f"d317stdd252_m_{tailin}": Q_wd3_in,
	f"d347stdd537_m_{tailin}": Q_vwd3_in,
	f"d300stdd50_m_{tailin}": Q_nr4_in,
	f"d300stdd100_m_{tailin}": Q_md4_in,
	f"d300stdd200_m_{tailin}": Q_wd4_in,
	f"d300stdd300_m_{tailin}": Q_vwd4_in,
	f"d300stdd50_e_{tailin}": Q_nr5_in,
	f"d300stdd100_e_{tailin}": Q_md5_in,
	f"d300stdd200_e_{tailin}": Q_wd5_in,
	f"d300stdd300_e_{tailin}": Q_vwd5_in,
	f"d300stdd50_e_{tailout}": Q_nr5_out,
	f"d300stdd100_e_{tailout}": Q_md5_out,
	f"d300stdd200_e_{tailout}": Q_wd5_out,
	f"d300stdd300_e_{tailout}": Q_vwd5_out,
	f"d150stdd50_m_{tailin}": Q_nr6_in,
	f"d150stdd100_m_{tailin}": Q_md6_in,
	f"d150stdd200_m_{tailin}": Q_wd6_in,
	f"d150stdd300_m_{tailin}": Q_vwd6_in,
	f"d167stdd100_{tailout}": Q_mu1_out,
	f"d197stdd65_{tailout}": Q_mu2_out,
	f"d290stdd97_{tailout}": Q_mu3_out,
	f"d430stdd100_{tailout}": Q_mu4_out,
	f"d269stdd100_{tailout}": Q_rt_out,
	f"d321stdd100_{tailout}": Q_lt_out,
	f"d240stdd50_{tailout}": Q_nr7_out,
	f"d400stdd50_{tailout}": Q_nr8_out,
}

# 绘图
fig = plt.figure(1, figsize=(8, 6), constrained_layout=True)
ax = fig.gca()

ax.plot(others_case_dict["Cr09"]["Sh"],
		others_case_dict["Cr09"]["Q_star"],
		'*',
		color='orange',
		label='Cr09',
		markersize=marker_size,
		markerfacecolor='none',
		markeredgewidth=marker_width
		)
ax.plot(others_case_dict["Ho2012"]["Sh"],
		others_case_dict["Ho2012"]["Q_star"],
		'x',
		color='purple',
		label='Ho12 Fine',
		markersize=marker_size,
		markerfacecolor='none',
		markeredgewidth=marker_width
		)
# 线性拟合
x_array = others_case_dict["Ho2012"]["Sh"]
y_array = others_case_dict["Ho2012"]["Q_star"]
fit = np.polyfit(x_array, y_array, 1)
x_array = np.linspace(lims["x_min"], lims["x_max"], 1000)
#ax.plot(x_array, fit[0]*x_array + fit[1], '--', color='purple', linewidth=line_width, label='Ho12 Fine Fit')
ksho12 = fit[0]

ax.plot(others_case_dict["Ho2012_c"]["Sh"],
		others_case_dict["Ho2012_c"]["Q_star"],
		'+',
		color='brown',
		label='Ho12 Coarse',
		markersize=marker_size,
		markerfacecolor='none',
		markeredgewidth=marker_width
		)
# 线性拟合
x_array = others_case_dict["Ho2012_c"]["Sh"]
y_array = others_case_dict["Ho2012_c"]["Q_star"]
fit = np.polyfit(x_array[3:], y_array[3:], 1)
x_array = np.linspace(lims["x_min"], lims["x_max"], 1000)
#ax.plot(x_array, fit[0]*x_array + fit[1], '--', color='brown', linewidth=line_width, label='Ho12 Coarse Fit')
ksho12c = fit[0]

axis_type = "loglog"

legend_str = "NR1"
dictkey = f"d300stdd50_{tailout}"
style = {"color": "r", "marker": "o", "fill": "full"}
arrow = False
fitline = False
rslt_plot(ax, Q_rslt, dictkey, legend_str, fitline, arrow, lims, sizes, style, axis_type)

legend_str = "MD1"
dictkey = f"d300stdd100_{tailout}"
style = {"color": "g", "marker": "^", "fill": "full"}
arrow = False
fitline = False
rslt_plot(ax, Q_rslt, dictkey, legend_str, fitline, arrow, lims, sizes, style, axis_type)

legend_str = "WD1"
dictkey = f"d300stdd200_{tailout}"
style = {"color": "b", "marker": "s", "fill": "full"}
arrow = False
fitline = False
rslt_plot(ax, Q_rslt, dictkey, legend_str, fitline, arrow, lims, sizes, style, axis_type)

legend_str = "VWD1"
dictkey = f"d300stdd300_{tailout}"
style = {"color": "c", "marker": "h", "fill": "full"}
arrow = False
fitline = False
rslt_plot(ax, Q_rslt, dictkey, legend_str, fitline, arrow, lims, sizes, style, axis_type)

legend_str = "NR2"
dictkey = f"d250stdd25_{tailout}"
style = {"color": "r", "marker": "o", "fill": "none"}
arrow = False
fitline = True
rslt_plot(ax, Q_rslt, dictkey, legend_str, fitline, arrow, lims, sizes, style, axis_type)

legend_str = "MD2"
dictkey = f"d271stdd121_{tailout}"
style = {"color": "g", "marker": "^", "fill": "none"}
arrow = False
fitline = False
rslt_plot(ax, Q_rslt, dictkey, legend_str, fitline, arrow, lims, sizes, style, axis_type)

legend_str = "WD2"
dictkey = f"d317stdd252_{tailout}"
style = {"color": "b", "marker": "s", "fill": "none"}
arrow = False
fitline = False
rslt_plot(ax, Q_rslt, dictkey, legend_str, fitline, arrow, lims, sizes, style, axis_type)

legend_str = "VWD2"
dictkey = f"d347stdd537_{tailout}"
style = {"color": "c", "marker": "h", "fill": "none"}
arrow = False
fitline = True
rslt_plot(ax, Q_rslt, dictkey, legend_str, fitline, arrow, lims, sizes, style, axis_type)

#rect = patches.Rectangle(
#	(0.315, 0.5), # (x, y)
#	0.1, # width
#	0.32, # height
#	transform=ax.transAxes,
#	angle = -45,
#	edgecolor='black',
#	facecolor='none',
#	linestyle='--',
#	linewidth=line_width*1.5
#)
#ax.add_patch(rect)

inset_ax = inset_axes(ax, width='42%', height='42%', loc='lower right')

sigma0 = np.array([0.1655, 0.3246, 0.6167, 0.87])
std0 = np.array([50, 100, 200, 300])
cv0 = std0 / 300
d500 = np.array([296, 285, 261, 249])
d900 = np.array([366, 432, 537, 589])
d90500 = d900 / d500
sigma1 = np.array([0.1, 0.4, 0.7, 1])
std1 = np.array([25, 121, 252, 537])
dpa = np.array([250, 271, 317, 347])
cv1 = std1 / dpa
d501 = np.array([248, 250, 265, 281])
d901 = np.array([282, 415, 584, 685])
d90501 = d901 / d501
d502 = np.array([157, 186, 275, 409])
d902 = np.array([233, 280, 416, 618])
d90502 = d902 / d502
d503 = np.array([261, 294])
d903 = np.array([356, 430])
d90503 = d903 / d503
d504 = np.array([235, 397])
d904 = np.array([306, 466])
d90504 = d904 / d504

x0, x1, x2, x3, x4 = d90500, d90501, d90502, d90503, d90504

ks0 = Q_rslt[f"d300stdd50_{tailout}"]["fit"][0]
ks1 = Q_rslt[f"d300stdd100_{tailout}"]["fit"][0]
ks2 = Q_rslt[f"d300stdd200_{tailout}"]["fit"][0]
ks3 = Q_rslt[f"d300stdd300_{tailout}"]["fit"][0]
ks10 = Q_rslt[f"d250stdd25_{tailout}"]["fit"][0]
ks11 = Q_rslt[f"d271stdd121_{tailout}"]["fit"][0]
ks12 = Q_rslt[f"d317stdd252_{tailout}"]["fit"][0]
ks13 = Q_rslt[f"d347stdd537_{tailout}"]["fit"][0]
ks20 = Q_rslt[f"d167stdd100_{tailout}"]["fit"][0]
ks21 = Q_rslt[f"d197stdd65_{tailout}"]["fit"][0]
ks22 = Q_rslt[f"d290stdd97_{tailout}"]["fit"][0]
ks23 = Q_rslt[f"d430stdd100_{tailout}"]["fit"][0]
ks30 = Q_rslt[f"d269stdd100_{tailout}"]["fit"][0]
ks31 = Q_rslt[f"d321stdd100_{tailout}"]["fit"][0]
ks40 = Q_rslt[f"d240stdd50_{tailout}"]["fit"][0]
ks41 = Q_rslt[f"d400stdd50_{tailout}"]["fit"][0]

y0 = np.array([ks0, ks1, ks2, ks3])
y1 = np.array([ks10, ks11, ks12, ks13])
y2 = np.array([ks20, ks21, ks22, ks23])
y3 = np.array([ks30, ks31])
y4 = np.array([ks40, ks41])
# merge x0~x4
x_array = np.concatenate([x0, x1, x2, x3, x4])
y_array = np.concatenate([y0, y1, y2, y3, y4])
fit = np.polyfit(x_array, y_array, 1)
x_array = np.linspace(0, 4, 100)
inset_ax.plot(x_array, fit[0] * x_array + fit[1], color='grey', linestyle='--', linewidth=line_width)


c_list0 = ['r', 'g', 'b', 'c']
c_list1 = ['r', 'g', 'b', 'c']
s_list0 = ['o', '^', 's', 'h']
s_list1 = ['o', '^', 's', 'h']
for x, y, s, c in zip(x0, y0, s_list0, c_list0):
    inset_ax.plot(x, y, s, color=c, markersize=marker_size)
for x, y, s, c in zip(x1, y1, s_list1, c_list1):
    inset_ax.plot(x, y, s, color=c, markersize=marker_size, markerfacecolor='none', markeredgewidth=marker_width)
inset_ax.plot(x2, y2, 'd', color='grey', markersize=marker_size, label='Other Groups')
inset_ax.plot(x3, y3, 'd', color='grey', markersize=marker_size)
inset_ax.plot(x4, y4, 'd', color='grey', markersize=marker_size)
inset_ax.plot(1.63, ksho12, 'x', color='purple', markersize=marker_size, markerfacecolor='none', markeredgewidth=marker_width)
inset_ax.plot(2.01, ksho12c, '+', color='brown', markersize=marker_size, markerfacecolor='none', markeredgewidth=marker_width)

ax.xaxis.set_major_formatter(ScalarFormatter())
ax.yaxis.set_major_formatter(ScalarFormatter())
ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
ax.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
ax.set_xlim(0.004, 0.2)
ax.set_ylim(0.0005, 0.2)
ax.set_xlabel('$S$', fontsize=label_size)
ax.set_ylabel(r'$\widetilde{Q}$', fontsize=label_size)
ax.tick_params(axis='both', labelsize=ticks_size)
ax.set_xticks([0.01, 0.05, 0.1])
ax.set_yticks([0.005, 0.01, 0.05, 0.1, 0.2])

inset_ax.xaxis.set_major_formatter(ScalarFormatter())
inset_ax.yaxis.set_major_formatter(ScalarFormatter())
inset_ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
inset_ax.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
# 把y的tick移到右边
#inset_ax.yaxis.tick_right()
#inset_ax.yaxis.set_label_position("right")
inset_ax.xaxis.tick_top()
inset_ax.xaxis.set_label_position("top")
inset_ax.set_xlabel('$\\eta_d$', fontsize=label_size)
inset_ax.set_ylabel('$k_s$', fontsize=label_size)
inset_ax.set_xlim(0.9, 2.6)
inset_ax.set_ylim(0.2, 1.1)
inset_ax.tick_params(axis='both', labelsize=ticks_size)

Cr_09 = ax.scatter([], [], marker='*', color='orange', s=marker_size**2, facecolors='none', linewidths=marker_width)
Ho_12 = ax.scatter([], [], marker='x', color='purple', s=marker_size**2, linewidths=marker_width)
Ho_12_c = ax.scatter([], [], marker='+', color='brown', s=marker_size**2, linewidths=marker_width)
proxy_ro = plt.Line2D([0], [0],
					 color='r',
					 marker='o',
					 linestyle='',
					 markersize=marker_size
					)
proxy_ro_h = plt.Line2D([0], [0],
					 color='r',
					 marker='o',
					 linestyle='',
					 markerfacecolor='none',
					 markersize=marker_size,
					 markeredgewidth=marker_width
					)
proxy_gt = plt.Line2D([0], [0],
					 color='g',
					 marker='^',
					 linestyle='',
					 markersize=marker_size
					)
proxy_gt_h = plt.Line2D([0], [0],
					 color='g',
					 marker='^',
					 linestyle='',
					 markerfacecolor='none',
					 markersize=marker_size,
					 markeredgewidth=marker_width
					)
proxy_bs = plt.Line2D([0], [0],
					 color='b',
					 marker='s',
					 linestyle='',
					 markersize=marker_size
					)
proxy_bs_h = plt.Line2D([0], [0],
					 color='b',
					 marker='s',
					 linestyle='',
					 markerfacecolor='none',
					 markersize=marker_size,
					 markeredgewidth=marker_width
					)
proxy_ch = plt.Line2D([0], [0],
					 color='c',
					 marker='h',
					 linestyle='',
					 markersize=marker_size
					)
proxy_ch_h = plt.Line2D([0], [0],
					 color='c',
					 marker='h',
					 linestyle='',
					 markerfacecolor='none',
					 markersize=marker_size,
					 markeredgewidth=marker_width
					)
ax.legend([(proxy_ro, proxy_ro_h),
		   (proxy_gt, proxy_gt_h),
		   (proxy_bs, proxy_bs_h),
		   (proxy_ch, proxy_ch_h),
		   Cr_09, Ho_12, Ho_12_c],
          ['NR1, NR2', 'MD1, MD2', 'WD1, WD2', 'VWD1, VWD2', 'Cr09 Unknown', 'Ho12 $\\eta_d=1.63$', 'Ho12 $\\eta_d=2.01$'],
		  handler_map={tuple: HandlerTuple(ndivide=None)},
          fontsize=ticks_size,
          loc='upper left',
          bbox_to_anchor=(-0.01, 1.02),
          frameon=True,
		  framealpha=1,
		  ncol=1
		  )
inset_ax.legend(fontsize=ticks_size,
				loc='upper left',
				bbox_to_anchor=(-0.1, 1.05),
				frameon=False,
				framealpha=1)

# 绘制床面和空中颗粒的CDF
fig = plt.figure(2, figsize=(8, 6), constrained_layout=True)
ax = fig.gca()
# 转为CDF
d_cdf_nr_bed = np.cumsum(hist_list0[0]) / np.sum(hist_list0[0])
d_cdf_md_bed = np.cumsum(hist_list0[1]) / np.sum(hist_list0[1])
d_cdf_wd_bed = np.cumsum(hist_list0[2]) / np.sum(hist_list0[2])
d_cdf_vwd_bed = np.cumsum(hist_list0[3]) / np.sum(hist_list0[3])
d_cdf_nr_air = np.cumsum(d_hist_nr) / np.sum(d_hist_nr)
d_cdf_md_air = np.cumsum(d_hist_md) / np.sum(d_hist_md)
d_cdf_wd_air = np.cumsum(d_hist_wd) / np.sum(d_hist_wd)
d_cdf_vwd_air = np.cumsum(d_hist_vwd) / np.sum(d_hist_vwd)
ax.plot(center_list0[0]*1e6, d_cdf_nr_bed, '-', color='r', label='NR1 Bed CDF', markersize=marker_size*0.8, markerfacecolor='none', markeredgewidth=marker_width)
ax.plot(center_list0[1]*1e6, d_cdf_md_bed, '-', color='g', label='MD1 Bed CDF', markersize=marker_size*0.8, markerfacecolor='none', markeredgewidth=marker_width)
ax.plot(center_list0[2]*1e6, d_cdf_wd_bed, '-', color='b', label='WD1 Bed CDF', markersize=marker_size*0.8, markerfacecolor='none', markeredgewidth=marker_width)
ax.plot(center_list0[3]*1e6, d_cdf_vwd_bed, '-', color='c', label='VWD1 Bed CDF', markersize=marker_size*0.8, markerfacecolor='none', markeredgewidth=marker_width)
ax.plot(d_bin_centers_nr*1e6, d_cdf_nr_air, 'o', color='r', label='NR1 Airborn CDF', markersize=marker_size*0.8, markerfacecolor='none', markeredgewidth=marker_width)
ax.plot(d_bin_centers_md*1e6, d_cdf_md_air, '^', color='g', label='MD1 Airborn CDF', markersize=marker_size*0.8, markerfacecolor='none', markeredgewidth=marker_width)
ax.plot(d_bin_centers_wd*1e6, d_cdf_wd_air, 's', color='b', label='WD1 Airborn CDF', markersize=marker_size*0.8, markerfacecolor='none', markeredgewidth=marker_width)
ax.semilogx(d_bin_centers_vwd*1e6, d_cdf_vwd_air, 'h', color='c', label='VWD1 Airborn CDF', markersize=marker_size*0.8, markerfacecolor='none', markeredgewidth=marker_width)
#ax.plot([], [], 'k-', label='Bed $\\mathbb{E}[d]$')
#ax.plot([], [], 'k--', label='Airborn $\\mathbb{E}[d]$')
#ax.axvline(dm_nr*1e6, color='r', linestyle='--')
#ax.axvline(dm_md*1e6, color='g', linestyle='--')
#ax.axvline(dm_wd*1e6, color='b', linestyle='--')
#ax.axvline(dm_vwd*1e6, color='c', linestyle='--')
#ax.axvline(3e-4*1e6, color='k', linestyle='-')
ax.axhline(0.5, color='k', linestyle='-')
ax.axhline(0.9, color='k', linestyle='-')
ax.text(146, 0.52, '50\\%', fontsize=ticks_size)
ax.text(146, 0.92, '90\\%', fontsize=ticks_size)

ax.xaxis.set_major_formatter(ScalarFormatter())
ax.yaxis.set_major_formatter(ScalarFormatter())
ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
ax.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
ax.tick_params(axis='both', labelsize=ticks_size)
ax.legend(fontsize=ticks_size,
		  loc='lower right',
		  #bbox_to_anchor=(1.0, 0.9),
		  framealpha=1,
		  )
ax.set_xlabel('$d$ ($\\mu$m)', fontsize=label_size)
ax.set_ylabel('CDF', fontsize=label_size)
ax.tick_params(axis='both', labelsize=ticks_size)
ax.set_xlim(100, 1000)

fig = plt.figure(6, figsize=(8, 6), constrained_layout=True)
ax = fig.gca()

eta_a0 = np.array(d_in_air_dict["NR1"]["d90"])/np.array(d_in_air_dict["NR1"]["d50"])
eta_a1 = np.array(d_in_air_dict["MD1"]["d90"])/np.array(d_in_air_dict["MD1"]["d50"])
eta_a2 = np.array(d_in_air_dict["WD1"]["d90"])/np.array(d_in_air_dict["WD1"]["d50"])
eta_a3 = np.array(d_in_air_dict["VWD1"]["d90"])/np.array(d_in_air_dict["VWD1"]["d50"])
eta_a10 = np.array(d_in_air_dict["NR2"]["d90"])/np.array(d_in_air_dict["NR2"]["d50"])
eta_a11 = np.array(d_in_air_dict["MD2"]["d90"])/np.array(d_in_air_dict["MD2"]["d50"])
eta_a12 = np.array(d_in_air_dict["WD2"]["d90"])/np.array(d_in_air_dict["WD2"]["d50"])
eta_a13 = np.array(d_in_air_dict["VWD2"]["d90"])/np.array(d_in_air_dict["VWD2"]["d50"])
eta_a20 = np.array(d_in_air_dict["MU1"]["d90"])/np.array(d_in_air_dict["MU1"]["d50"])
eta_a21 = np.array(d_in_air_dict["MU2"]["d90"])/np.array(d_in_air_dict["MU2"]["d50"])
eta_a22 = np.array(d_in_air_dict["MU3"]["d90"])/np.array(d_in_air_dict["MU3"]["d50"])
eta_a23 = np.array(d_in_air_dict["MU4"]["d90"])/np.array(d_in_air_dict["MU4"]["d50"])
eta_a30 = np.array(d_in_air_dict["RT"]["d90"])/np.array(d_in_air_dict["RT"]["d50"])
eta_a31 = np.array(d_in_air_dict["LT"]["d90"])/np.array(d_in_air_dict["LT"]["d50"])
eta_a40 = np.array(d_in_air_dict["NR3"]["d90"])/np.array(d_in_air_dict["NR3"]["d50"])
eta_a41 = np.array(d_in_air_dict["NR4"]["d90"])/np.array(d_in_air_dict["NR4"]["d50"])
meta_a0 = np.mean(eta_a0)
meta_a1 = np.mean(eta_a1)
meta_a2 = np.mean(eta_a2)
meta_a3 = np.mean(eta_a3)
meta_a10 = np.mean(eta_a10)
meta_a11 = np.mean(eta_a11)
meta_a12 = np.mean(eta_a12)
meta_a13 = np.mean(eta_a13)
meta_a20 = np.mean(eta_a20)
meta_a21 = np.mean(eta_a21)
meta_a22 = np.mean(eta_a22)
meta_a23 = np.mean(eta_a23)
meta_a30 = np.mean(eta_a30)
meta_a31 = np.mean(eta_a31)
meta_a40 = np.mean(eta_a40)
meta_a41 = np.mean(eta_a41)
meta_a_array0 = np.array([meta_a0, meta_a1, meta_a2, meta_a3])
meta_a_array1 = np.array([meta_a10, meta_a11, meta_a12, meta_a13])
meta_a_others = np.array([meta_a20, meta_a21, meta_a22, meta_a23, meta_a40, meta_a41])
meta_a_array3 = np.array([meta_a30, meta_a31])
std_a0 = np.std(eta_a0)
std_a1 = np.std(eta_a1)
std_a2 = np.std(eta_a2)
std_a3 = np.std(eta_a3)
std_a10 = np.std(eta_a10)
std_a11 = np.std(eta_a11)
std_a12 = np.std(eta_a12)
std_a13 = np.std(eta_a13)
std_a20 = np.std(eta_a20)
std_a21 = np.std(eta_a21)
std_a22 = np.std(eta_a22)
std_a23 = np.std(eta_a23)
std_a30 = np.std(eta_a30)
std_a31 = np.std(eta_a31)
std_a40 = np.std(eta_a40)
std_a41 = np.std(eta_a41)
std_a_array0 = np.array([std_a0, std_a1, std_a2, std_a3])
std_a_array1 = np.array([std_a10, std_a11, std_a12, std_a13])
std_a_others = np.array([std_a20, std_a21, std_a22, std_a23, std_a40, std_a41])
std_a_array3 = np.array([std_a30, std_a31])
# 计算床面颗粒的eta
d500 = np.array([296, 285, 261, 249])
d900 = np.array([366, 432, 537, 589])
d90500 = d900 / d500
d501 = np.array([248, 250, 265, 281])
d901 = np.array([282, 415, 584, 685])
d90501 = d901 / d501
d502 = np.array([157, 186, 275, 409])
d902 = np.array([233, 280, 416, 618])
d90502 = d902 / d502
d503 = np.array([261, 294])
d903 = np.array([356, 430])
d90503 = d903 / d503
d504 = np.array([235, 397])
d904 = np.array([306, 466])
d90504 = d904 / d504
eta_b0 = d90500
eta_b1 = d90501
eta_b3 = d90503
meta_b_others = np.concatenate([d90502, d90504])

c_list = ['r', 'g', 'b', 'c']
s_list = ['o', '^', 's', 'h']
for i in range(len(meta_a_others)):
	ax.errorbar(meta_b_others[i], meta_a_others[i], yerr=std_a_others[i],
				fmt='d',
				color='grey',
				markersize=marker_size,
				markeredgewidth=marker_width,
				capsize=marker_size*0.6,
				label='Other Groups' if i == 0 else "")
l_list0 = ['NR1', 'MD1', 'WD1', 'VWD1']
for xi, yi, ei, si, ci, li in zip(d90500, meta_a_array0, std_a_array0, s_list, c_list, l_list0):
	ax.errorbar(xi, yi, yerr=ei,
				fmt=si,
				color=ci,
				label=li,
				markersize=marker_size,
				markeredgewidth=marker_width,
				capsize=marker_size*0.6)
l_list1 = ['NR2', 'MD2', 'WD2', 'VWD2']
for xi, yi, ei, si, ci, li in zip(d90501, meta_a_array1, std_a_array1, s_list, c_list, l_list1):
	ax.errorbar(xi, yi, yerr=ei,
				fmt=si,
				color=ci,
				label=li,
				markersize=marker_size,
				markerfacecolor='none',
				markeredgewidth=marker_width,
				capsize=marker_size*0.6)
ax.errorbar(d90503[0], meta_a_array3[0], yerr=std_a_array3[0],
			fmt='P',
			color='y',
			label='RT',
			markersize=marker_size,
			markeredgewidth=marker_width,
			capsize=marker_size*0.6)
ax.errorbar(d90503[1], meta_a_array3[1], yerr=std_a_array3[1],
			fmt='P',
			color='y',
			label='LT',
			markersize=marker_size,
			markerfacecolor='none',
			markeredgewidth=marker_width,
			capsize=marker_size*0.6)
# 45度参考线
x_array = np.linspace(0, 3, 100)
y_array = x_array
ax.plot(x_array, y_array, color='k', linestyle='-', linewidth=line_width)
ax.annotate('1:1 Line', xy=(1.4, 1.4), xytext=(1.4, 1.2),
			arrowprops=dict(arrowstyle='->', lw=1.5),
			fontsize=ticks_size
			)
ax.set_xlabel('$\\eta_d^{\\mathrm{bed}}$', fontsize=label_size)
ax.set_ylabel('$\\eta_d^{\\mathrm{air}}$', fontsize=label_size)
ax.tick_params(axis='both', labelsize=ticks_size)
ax.set_xlim(1.0, 2.5)
ax.set_ylim(1.0, 1.8)
others = ax.scatter([], [], marker='d', color='gray', s=marker_size**2, linewidths=marker_width)
proxy_ro = plt.Line2D([0], [0],
					 color='r',
					 marker='o',
					 linestyle='',
					 markersize=marker_size
					)
proxy_ro_h = plt.Line2D([0], [0],
					 color='r',
					 marker='o',
					 linestyle='',
					 markerfacecolor='none',
					 markersize=marker_size,
					 markeredgewidth=marker_width
					)
proxy_gt = plt.Line2D([0], [0],
					 color='g',
					 marker='^',
					 linestyle='',
					 markersize=marker_size
					)
proxy_gt_h = plt.Line2D([0], [0],
					 color='g',
					 marker='^',
					 linestyle='',
					 markerfacecolor='none',
					 markersize=marker_size,
					 markeredgewidth=marker_width
					)
proxy_bs = plt.Line2D([0], [0],
					 color='b',
					 marker='s',
					 linestyle='',
					 markersize=marker_size
					)
proxy_bs_h = plt.Line2D([0], [0],
					 color='b',
					 marker='s',
					 linestyle='',
					 markerfacecolor='none',
					 markersize=marker_size,
					 markeredgewidth=marker_width
					)
proxy_ch = plt.Line2D([0], [0],
					 color='c',
					 marker='h',
					 linestyle='',
					 markersize=marker_size
					)
proxy_ch_h = plt.Line2D([0], [0],
					 color='c',
					 marker='h',
					 linestyle='',
					 markerfacecolor='none',
					 markersize=marker_size,
					 markeredgewidth=marker_width
					)
proxy_yD = plt.Line2D([0], [0],
					 color='y',
					 marker='P',
					 linestyle='',
					 markersize=marker_size
					)
proxy_yD_h = plt.Line2D([0], [0],
					 color='y',
					 marker='P',
					 linestyle='',
					 markerfacecolor='none',
					 markersize=marker_size,
					 markeredgewidth=marker_width
					)
ax.legend([(proxy_ro, proxy_ro_h),
		   (proxy_gt, proxy_gt_h),
		   (proxy_bs, proxy_bs_h),
		   (proxy_ch, proxy_ch_h),
		   (proxy_yD, proxy_yD_h),
		   others],
          ['NR1, NR2', 'MD1, MD2', 'WD1, WD2', 'VWD1, VWD2', 'RT, LT', 'Other Groups'],
		  handler_map={tuple: HandlerTuple(ndivide=None)},
          fontsize=ticks_size,
          loc='upper left',
          bbox_to_anchor=(0.0, 1.02),
          frameon=True,
		  framealpha=1,
		  ncol=2
		  )

inset_ax = inset_axes(ax, width='42%', height='42%', loc='lower right')

d90_a0 = np.array(d_in_air_dict["NR1"]["d90"])
d90_a1 = np.array(d_in_air_dict["MD1"]["d90"])
d90_a2 = np.array(d_in_air_dict["WD1"]["d90"])
d90_a3 = np.array(d_in_air_dict["VWD1"]["d90"])
d90_a10 = np.array(d_in_air_dict["NR2"]["d90"])
d90_a11 = np.array(d_in_air_dict["MD2"]["d90"])
d90_a12 = np.array(d_in_air_dict["WD2"]["d90"])
d90_a13 = np.array(d_in_air_dict["VWD2"]["d90"])
d90_a20 = np.array(d_in_air_dict["MU1"]["d90"])
d90_a21 = np.array(d_in_air_dict["MU2"]["d90"])
d90_a22 = np.array(d_in_air_dict["MU3"]["d90"])
d90_a23 = np.array(d_in_air_dict["MU4"]["d90"])
d90_a30 = np.array(d_in_air_dict["RT"]["d90"])
d90_a31 = np.array(d_in_air_dict["LT"]["d90"])
d90_a40 = np.array(d_in_air_dict["NR3"]["d90"])
d90_a41 = np.array(d_in_air_dict["NR4"]["d90"])
dm_a0 = np.array(d_in_air_dict["NR1"]["dm"])
dm_a1 = np.array(d_in_air_dict["MD1"]["dm"])
dm_a2 = np.array(d_in_air_dict["WD1"]["dm"])
dm_a3 = np.array(d_in_air_dict["VWD1"]["dm"])
dm_a10 = np.array(d_in_air_dict["NR2"]["dm"])
dm_a11 = np.array(d_in_air_dict["MD2"]["dm"])
dm_a12 = np.array(d_in_air_dict["WD2"]["dm"])
dm_a13 = np.array(d_in_air_dict["VWD2"]["dm"])
dm_a20 = np.array(d_in_air_dict["MU1"]["dm"])
dm_a21 = np.array(d_in_air_dict["MU2"]["dm"])
dm_a22 = np.array(d_in_air_dict["MU3"]["dm"])
dm_a23 = np.array(d_in_air_dict["MU4"]["dm"])
dm_a30 = np.array(d_in_air_dict["RT"]["dm"])
dm_a31 = np.array(d_in_air_dict["LT"]["dm"])
dm_a40 = np.array(d_in_air_dict["NR3"]["dm"])
dm_a41 = np.array(d_in_air_dict["NR4"]["dm"])
md90_a0 = np.mean(d90_a0)
md90_a1 = np.mean(d90_a1)
md90_a2 = np.mean(d90_a2)
md90_a3 = np.mean(d90_a3)
md90_a10 = np.mean(d90_a10)
md90_a11 = np.mean(d90_a11)
md90_a12 = np.mean(d90_a12)
md90_a13 = np.mean(d90_a13)
md90_a20 = np.mean(d90_a20)
md90_a21 = np.mean(d90_a21)
md90_a22 = np.mean(d90_a22)
md90_a23 = np.mean(d90_a23)
md90_a30 = np.mean(d90_a30)
md90_a31 = np.mean(d90_a31)
md90_a40 = np.mean(d90_a40)
md90_a41 = np.mean(d90_a41)
mdm_a0 = np.mean(dm_a0)
mdm_a1 = np.mean(dm_a1)
mdm_a2 = np.mean(dm_a2)
mdm_a3 = np.mean(dm_a3)
mdm_a10 = np.mean(dm_a10)
mdm_a11 = np.mean(dm_a11)
mdm_a12 = np.mean(dm_a12)
mdm_a13 = np.mean(dm_a13)
mdm_a20 = np.mean(dm_a20)
mdm_a21 = np.mean(dm_a21)
mdm_a22 = np.mean(dm_a22)
mdm_a23 = np.mean(dm_a23)
mdm_a30 = np.mean(dm_a30)
mdm_a31 = np.mean(dm_a31)
mdm_a40 = np.mean(dm_a40)
mdm_a41 = np.mean(dm_a41)
md90_a_array0 = np.array([md90_a0, md90_a1, md90_a2, md90_a3])*1e6
md90_a_array1 = np.array([md90_a10, md90_a11, md90_a12, md90_a13])*1e6
md90_a_others = np.array([md90_a20, md90_a21, md90_a22, md90_a23, md90_a40, md90_a41])*1e6
md90_a_array3 = np.array([md90_a30, md90_a31])*1e6
mdm_a_array0 = np.array([mdm_a0, mdm_a1, mdm_a2, mdm_a3])*1e6
mdm_a_array1 = np.array([mdm_a10, mdm_a11, mdm_a12, mdm_a13])*1e6
mdm_a_others = np.array([mdm_a20, mdm_a21, mdm_a22, mdm_a23, mdm_a40, mdm_a41])*1e6
mdm_a_array3 = np.array([mdm_a30, mdm_a31])*1e6
std90_a0 = np.std(d90_a0)
std90_a1 = np.std(d90_a1)
std90_a2 = np.std(d90_a2)
std90_a3 = np.std(d90_a3)
std90_a10 = np.std(d90_a10)
std90_a11 = np.std(d90_a11)
std90_a12 = np.std(d90_a12)
std90_a13 = np.std(d90_a13)
std90_a20 = np.std(d90_a20)
std90_a21 = np.std(d90_a21)
std90_a22 = np.std(d90_a22)
std90_a23 = np.std(d90_a23)
std90_a30 = np.std(d90_a30)
std90_a31 = np.std(d90_a31)
std90_a40 = np.std(d90_a40)
std90_a41 = np.std(d90_a41)
stdm_a0 = np.std(dm_a0)
stdm_a1 = np.std(dm_a1)
stdm_a2 = np.std(dm_a2)
stdm_a3 = np.std(dm_a3)
stdm_a10 = np.std(dm_a10)
stdm_a11 = np.std(dm_a11)
stdm_a12 = np.std(dm_a12)
stdm_a13 = np.std(dm_a13)
stdm_a20 = np.std(dm_a20)
stdm_a21 = np.std(dm_a21)
stdm_a22 = np.std(dm_a22)
stdm_a23 = np.std(dm_a23)
stdm_a30 = np.std(dm_a30)
stdm_a31 = np.std(dm_a31)
stdm_a40 = np.std(dm_a40)
stdm_a41 = np.std(dm_a41)
std90_a_array0 = np.array([std90_a0, std90_a1, std90_a2, std90_a3])*1e6
std90_a_array1 = np.array([std90_a10, std90_a11, std90_a12, std90_a13])*1e6
std90_a_others = np.array([std90_a20, std90_a21, std90_a22, std90_a23, std90_a40, std90_a41])*1e6
std90_a_array3 = np.array([std90_a30, std90_a31])*1e6
stdm_a_array0 = np.array([stdm_a0, stdm_a1, stdm_a2, stdm_a3])*1e6
stdm_a_array1 = np.array([stdm_a10, stdm_a11, stdm_a12, stdm_a13])*1e6
stdm_a_others = np.array([stdm_a20, stdm_a21, stdm_a22, stdm_a23, stdm_a40, stdm_a41])*1e6
stdm_a_array3 = np.array([stdm_a30, stdm_a31])*1e6

dm0 = np.array([300, 300, 300, 300])
dm1 = np.array([250, 271, 317, 347])
dm2 = np.array([167, 197, 290, 430])
dm3 = np.array([269, 313])
dm4 = np.array([240, 400])
d90_others = np.concatenate([d902, d904])
dm_others = np.concatenate([dm2, dm4])

d90900 = md90_a_array0/d900
d90901 = md90_a_array1/d901
d90903 = md90_a_array3/d903
d9090_others = md90_a_others/d90_others
d90900 = 1.0/d90900
d90901 = 1.0/d90901
d90903 = 1.0/d90903
d9090_others = 1.0/d9090_others
dmm0 = mdm_a_array0/dm0
dmm1 = mdm_a_array1/dm1
dmm3 = mdm_a_array3/dm3
dmm_others = mdm_a_others/dm_others
dmm0 = 1.0/dmm0
dmm1 = 1.0/dmm1
dmm3 = 1.0/dmm3
dmm_others = 1.0/dmm_others
st90900 = d900*std90_a_array0/md90_a_array0**2
st90901 = d901*std90_a_array1/md90_a_array1**2
st90903 = d903*std90_a_array3/md90_a_array3**2
st9090_others = d90_others*std90_a_others/md90_a_others**2
stmm0 = dm0*stdm_a_array0/mdm_a_array0**2
stmm1 = dm1*stdm_a_array1/mdm_a_array1**2
stmm3 = dm3*stdm_a_array3/mdm_a_array3**2
stmm_others = dm_others*stdm_a_others/mdm_a_others**2
c_list = ['m', 'm', 'm', 'm']
s_list = ['o', 'o', 'o', 'o']
scale = 1.0
for i in range(len(d9090_others)):
	inset_ax.errorbar(meta_b_others[i], d9090_others[i], yerr=st9090_others[i],
				fmt='o',
				color='m',
				markersize=marker_size_in*scale,
				markerfacecolor='none',
				markeredgewidth=marker_width,
				capsize=marker_size_in*0.6*scale)
l_list0 = ['NR1', 'MD1', 'WD1', 'VWD1']
for xi, yi, ei, si, ci, li in zip(d90500, d90900, st90900, s_list, c_list, l_list0):
	inset_ax.errorbar(xi, yi, yerr=ei,
				fmt=si,
				color=ci,
				markersize=marker_size_in*scale,
				markerfacecolor='none',
				markeredgewidth=marker_width,
				capsize=marker_size_in*0.6*scale)
l_list1 = ['NR2', 'MD2', 'WD2', 'VWD2']
for xi, yi, ei, si, ci, li in zip(d90501, d90901, st90901, s_list, c_list, l_list1):
	inset_ax.errorbar(xi, yi, yerr=ei,
				fmt=si,
				color=ci,
				markersize=marker_size_in*scale,
				markerfacecolor='none',
				markeredgewidth=marker_width,
				capsize=marker_size_in*0.6*scale)
inset_ax.errorbar(d90503[0], d90903[0], yerr=st90903[0],
			fmt='o',
			color='m',
			label='Ratios of $d_{90}$',
			markersize=marker_size_in*scale,
			markerfacecolor='none',
			markeredgewidth=marker_width,
			capsize=marker_size_in*0.6*scale)
inset_ax.errorbar(d90503[1], d90903[1], yerr=st90903[1],
			fmt='o',
			color='m',
			markersize=marker_size_in*scale,
			markerfacecolor='none',
			markeredgewidth=marker_width,
			capsize=marker_size_in*0.6*scale)

c_list = ['k', 'k', 'k', 'k']
s_list = ['s', 's', 's', 's']
for i in range(len(dmm_others)):
	inset_ax.errorbar(meta_b_others[i], dmm_others[i], yerr=stmm_others[i],
				fmt='s',
				color='k',
				markersize=marker_size_in*scale,
				markerfacecolor='none',
				markeredgewidth=marker_width,
				capsize=marker_size_in*0.6*scale)
l_list0 = ['NR1', 'MD1', 'WD1', 'VWD1']
for xi, yi, ei, si, ci, li in zip(d90500, dmm0, stmm0, s_list, c_list, l_list0):
	inset_ax.errorbar(xi, yi, yerr=ei,
				fmt=si,
				color=ci,
				markersize=marker_size_in*scale,
				markerfacecolor='none',
				markeredgewidth=marker_width,
				capsize=marker_size_in*0.6*scale)
l_list1 = ['NR2', 'MD2', 'WD2', 'VWD2']
for xi, yi, ei, si, ci, li in zip(d90501, dmm1, stmm1, s_list, c_list, l_list1):
	inset_ax.errorbar(xi, yi, yerr=ei,
				fmt=si,
				color=ci,
				markersize=marker_size_in*scale,
				markerfacecolor='none',
				markeredgewidth=marker_width,
				capsize=marker_size_in*0.6*scale)
inset_ax.errorbar(d90503[0], dmm3[0], yerr=stmm3[0],
			fmt='s',
			color='k',
			label='Ratios of $\\mathbb{E}[d]$',
			markersize=marker_size_in*scale,
			markerfacecolor='none',
			markeredgewidth=marker_width,
			capsize=marker_size_in*0.6*scale)
inset_ax.errorbar(d90503[1], dmm3[1], yerr=stmm3[1],
			fmt='s',
			color='k',
			markersize=marker_size_in*scale,
			markerfacecolor='none',
			markeredgewidth=marker_width,
			capsize=marker_size_in*0.6*scale)

x_array = np.concatenate([d90501, d90501, d90502, d90504, d90503])
y_array = np.concatenate([dmm0, dmm1, dmm_others, dmm3])
# 限制x=1时y=1
k, _, _, _ = np.linalg.lstsq((x_array-1).reshape(-1, 1), (y_array-1), rcond=None)
x_array = np.linspace(0, 4, 100)
y_array = k * (x_array - 1) + 1
inset_ax.plot(x_array, y_array, color='grey', linestyle='--', linewidth=line_width)

inset_ax.set_xlim(0.8, 3.2)
inset_ax.set_ylim(0.5, 3.8)
inset_ax.xaxis.tick_top()
inset_ax.xaxis.set_label_position("top")
inset_ax.set_xlabel('$\\eta_d^\\mathrm{bed}$', fontsize=label_size, labelpad=10)
inset_ax.set_ylabel('$d^\\mathrm{bed}/d^\\mathrm{air}$', fontsize=label_size)
inset_ax.tick_params(axis='both', labelsize=ticks_size)
inset_ax.legend(fontsize=ticks_size,
				loc='lower left',
				bbox_to_anchor=(-1.05, 0),
				frameon=True,
				framealpha=1,
				)


fig = plt.figure(3, figsize=(8, 6), constrained_layout=True)
ax = fig.gca()

M_rslt = {
	f"d300stdd50_{tailin}": M_nr1_in,
	f"d300stdd100_{tailin}": M_md1_in,
	f"d300stdd200_{tailin}": M_wd1_in,
	f"d300stdd300_{tailin}": M_vwd1_in,
	f"d300stdd50_{tailout}": M_nr1,
	f"d300stdd100_{tailout}": M_md1,
	f"d300stdd200_{tailout}": M_wd1,
	f"d300stdd300_{tailout}": M_vwd1,
	f"d250stdd25_{tailin}": M_nr2,
	f"d271stdd121_{tailin}": M_md2,
	f"d317stdd252_{tailin}": M_wd2,
	f"d347stdd537_{tailin}": M_vwd2,
	f"d240stdd50_{tailin}": M_nr3,
	f"d400stdd50_{tailin}": M_nr4,
	f"d430stdd100_{tailin}": M_mu4_out,
	f"d167stdd100_{tailin}": M_mu1_out,
	f"d269stdd100_{tailin}": M_rt_out,
	f"d321stdd100_{tailin}": M_lt_out,
	f"d250stdd25_m_{tailin}": M_nr3_in,
	f"d271stdd121_m_{tailin}": M_md3_in,
	f"d317stdd252_m_{tailin}": M_wd3_in,
	f"d347stdd537_m_{tailin}": M_vwd3_in,
	f"d300stdd50_m_{tailin}": M_nr4_in,
	f"d300stdd100_m_{tailin}": M_md4_in,
	f"d300stdd200_m_{tailin}": M_wd4_in,
	f"d300stdd300_m_{tailin}": M_vwd4_in,
	f"d300stdd50_e_{tailin}": M_nr5_in,
	f"d300stdd100_e_{tailin}": M_md5_in,
	f"d300stdd200_e_{tailin}": M_wd5_in,
	f"d300stdd300_e_{tailin}": M_vwd5_in,
	f"d300stdd50_e_{tailout}": M_nr5_out,
	f"d300stdd100_e_{tailout}": M_md5_out,
	f"d300stdd200_e_{tailout}": M_wd5_out,
	f"d300stdd300_e_{tailout}": M_vwd5_out,
	f"d150stdd50_m_{tailin}": M_nr6_in,
	f"d150stdd100_m_{tailin}": M_md6_in,
	f"d150stdd200_m_{tailin}": M_wd6_in,
	f"d150stdd300_m_{tailin}": M_vwd6_in,
	f"d250stdd25_{tailout}": M_nr2_out,
	f"d271stdd121_{tailout}": M_md2_out,
	f"d317stdd252_{tailout}": M_wd2_out,
	f"d347stdd537_{tailout}": M_vwd2_out,
	f"d167stdd100_{tailout}": M_mu1_out,
	f"d197stdd65_{tailout}": M_mu2_out,
	f"d290stdd97_{tailout}": M_mu3_out,
	f"d430stdd100_{tailout}": M_mu4_out,
	f"d269stdd100_{tailout}": M_rt_out,
	f"d321stdd100_{tailout}": M_lt_out,
	f"d240stdd50_{tailout}": M_nr7_out,
	f"d400stdd50_{tailout}": M_nr8_out,
}

legend_str = "NR300"
dictkey = f"d300stdd50_m_{tailin}"
style = {"color": "r", "marker": "o", "fill": "full"}
arrow = False
fitline = True
lims["x_loc"] = 0.006
lims["y_loc"] = 0.01
rslt_plot(ax, M_rslt, dictkey, legend_str, fitline, arrow, lims, sizes, style, axis_type)

legend_str = "MD300"
dictkey = f"d300stdd100_m_{tailin}"
style = {"color": "g", "marker": "^", "fill": "full"}
arrow = False
fitline = False
rslt_plot(ax, M_rslt, dictkey, legend_str, fitline, arrow, lims, sizes, style, axis_type)

legend_str = "WD300"
dictkey = f"d300stdd200_m_{tailin}"
style = {"color": "b", "marker": "s", "fill": "full"}
arrow = False
fitline = False
rslt_plot(ax, M_rslt, dictkey, legend_str, fitline, arrow, lims, sizes, style, axis_type)

legend_str = "VWD300"
dictkey = f"d300stdd300_m_{tailin}"
style = {"color": "c", "marker": "h", "fill": "full"}
arrow = False
fitline = True
lims["x_loc"] = 0.012
lims["y_loc"] = 0.04
rslt_plot(ax, M_rslt, dictkey, legend_str, fitline, arrow, lims, sizes, style, axis_type)

legend_str = "NR150"
dictkey = f"d150stdd50_m_{tailin}"
style = {"color": "r", "marker": "o", "fill": "full"}
arrow = False
fitline = True
S = np.array(M_rslt[dictkey]["x"])*1.5/3
M = np.array(M_rslt[dictkey]["y"])*1.5/3
#rslt_plot(ax, M_rslt, dictkey, legend_str, fitline, arrow, lims, sizes, style, axis_type)
ax.scatter(S, M, color=style["color"], marker=style["marker"], s=sizes["marker"]**2, label=legend_str, edgecolors='k', linewidths=marker_width)

legend_str = "MD150"
dictkey = f"d150stdd100_m_{tailin}"
style = {"color": "g", "marker": "^", "fill": "full"}
arrow = False
fitline = False
S = np.array(M_rslt[dictkey]["x"])*1.5/3
M = np.array(M_rslt[dictkey]["y"])*1.5/3
#rslt_plot(ax, M_rslt, dictkey, legend_str, fitline, arrow, lims, sizes, style, axis_type)
ax.scatter(S, M, color=style["color"], marker=style["marker"], s=sizes["marker"]**2, label=legend_str, edgecolors='k', linewidths=marker_width)

legend_str = "WD150"
dictkey = f"d150stdd200_m_{tailin}"
style = {"color": "b", "marker": "s", "fill": "full"}
arrow = False
fitline = False
S = np.array(M_rslt[dictkey]["x"])*1.5/3
M = np.array(M_rslt[dictkey]["y"])*1.5/3
#rslt_plot(ax, M_rslt, dictkey, legend_str, fitline, arrow, lims, sizes, style, axis_type)
ax.scatter(S, M, color=style["color"], marker=style["marker"], s=sizes["marker"]**2, label=legend_str, edgecolors='k', linewidths=marker_width)

legend_str = "VWD150"
dictkey = f"d150stdd300_m_{tailin}"
style = {"color": "c", "marker": "h", "fill": "full"}
arrow = False
fitline = True
S = np.array(M_rslt[dictkey]["x"])*1.5/3
M = np.array(M_rslt[dictkey]["y"])*1.5/3
#rslt_plot(ax, M_rslt, dictkey, legend_str, fitline, arrow, lims, sizes, style, axis_type)
ax.scatter(S, M, color=style["color"], marker=style["marker"], s=sizes["marker"]**2, label=legend_str, edgecolors='k', linewidths=marker_width)

ax.set_xlabel('$S$', fontsize=label_size)
ax.set_ylabel('$M^*$', fontsize=label_size)
ax.tick_params(axis='both', labelsize=ticks_size)
ax.set_xlim(0.004, 0.3)
ax.set_ylim(0.005, 0.3)
ax.legend(fontsize=ticks_size)

inset_ax = inset_axes(ax, width='45%', height='45%', loc='lower right')

ks0 = M_rslt[f"d300stdd50_m_{tailin}"]["fit"][0]
ks1 = M_rslt[f"d300stdd100_m_{tailin}"]["fit"][0]
ks2 = M_rslt[f"d300stdd200_m_{tailin}"]["fit"][0]
ks3 = M_rslt[f"d300stdd300_m_{tailin}"]["fit"][0]
k0 = -M_rslt[f"d300stdd50_m_{tailin}"]["fit"][1]/ks0
k1 = -M_rslt[f"d300stdd100_m_{tailin}"]["fit"][1]/ks1
k2 = -M_rslt[f"d300stdd200_m_{tailin}"]["fit"][1]/ks2
k3 = -M_rslt[f"d300stdd300_m_{tailin}"]["fit"][1]/ks3
ks10 = M_rslt[f"d250stdd25_m_{tailin}"]["fit"][0]
ks11 = M_rslt[f"d271stdd121_m_{tailin}"]["fit"][0]
ks12 = M_rslt[f"d317stdd252_m_{tailin}"]["fit"][0]
ks13 = M_rslt[f"d347stdd537_m_{tailin}"]["fit"][0]
k10 = -M_rslt[f"d250stdd25_m_{tailin}"]["fit"][1]/ks10
k11 = -M_rslt[f"d271stdd121_m_{tailin}"]["fit"][1]/ks11
k12 = -M_rslt[f"d317stdd252_m_{tailin}"]["fit"][1]/ks12
k13 = -M_rslt[f"d347stdd537_m_{tailin}"]["fit"][1]/ks13
ks21 = M_rslt[f"d150stdd50_m_{tailin}"]["fit"][0]
ks22 = M_rslt[f"d150stdd100_m_{tailin}"]["fit"][0]
ks23 = M_rslt[f"d150stdd200_m_{tailin}"]["fit"][0]
ks24 = M_rslt[f"d150stdd300_m_{tailin}"]["fit"][0]
k21 = -M_rslt[f"d150stdd50_m_{tailin}"]["fit"][1]/ks21
k22 = -M_rslt[f"d150stdd100_m_{tailin}"]["fit"][1]/ks22
k23 = -M_rslt[f"d150stdd200_m_{tailin}"]["fit"][1]/ks23
k24 = -M_rslt[f"d150stdd300_m_{tailin}"]["fit"][1]/ks24

ks0_array = np.array([ks0, ks1, ks2, ks3])
ks1_array = np.array([ks10, ks11, ks12, ks13])
ks2_array = np.array([ks21, ks22, ks23, ks24])
k0_array = np.array([k0, k1, k2, k3])
k1_array = np.array([k10, k11, k12, k13])
k2_array = np.array([k21, k22, k23, k24])

dm1 = np.array([250, 271, 317, 347])
std0 = np.array([50, 100, 200, 300])
std1 = np.array([25, 121, 252, 537])
cv0 = std0 / 300
cv1 = std1 / dm1
#k_list = [k0, k1, k2, k3]
#k_array = np.array(k_list)
#k1_list = [k10, k11, k12, k13]
#k1_array = np.array(k1_list)
#sigma = np.array([0.1655, 0.3246, 0.6167, 0.87])
#sd_from_tau = np.array([0.005245, 0.006080*1.2, 0.007008*1.4, 0.007238*1.6])
#sd_from_tau = sd_from_tau
#sigma1 = np.array([0.1, 0.4, 0.7, 1])
#sd_from_tau1 = np.array([0.005152, 0.005030, 0.005092, 0.004603])
#sd_from_tau1 = sd_from_tau1
#sigma2 = np.array([0.1245, 0.2061])
##k2_array = np.array([k20, k21])
##sigma3 = np.array([0.3246, 0.3246])
##k3_array = np.array([k30, k31])
##sigma4 = np.array([0.3246, 0.3246])
##k4_array = np.array([k40, k41])
#k5_array = np.array([k51, k52, k53, k54])
#k6_array = np.array([k61, k62, k63, k64])
#k7_array = np.array([k71, k72, k73, k74])

#x, x1, x2, x3, x4, x5, x6 = sigma, sigma1, sigma2, sigma3, sigma4, sigma5, sigma6
x0, x1 = cv0, cv1

inset_ax.plot(x0, ks0_array, 'r*-', markersize=marker_size_in, markeredgewidth=marker_width, label='Relative $k_s$')
inset_ax.plot(x1, ks1_array, 'g*-', markersize=marker_size_in, markeredgewidth=marker_width, label='Relative $k_s$')
inset_ax.plot(x0, ks2_array, 'k*', markersize=marker_size_in, markeredgewidth=marker_width, label='Relative $k_s$')
#inset_ax.plot(x0, k0_array, 'rx-', markersize=marker_size_in, markeredgewidth=marker_width, label='Group 1')
#inset_ax.plot(x1, k1_array, 'gx-', markersize=marker_size_in, markeredgewidth=marker_width, label='Group 2')
#inset_ax.plot(x0, k2_array, 'kx', markersize=marker_size_in, markeredgewidth=marker_width, label='Group 8')
##inset_ax.plot(x2, k2_array, 'ko', markersize=marker_size_in, markeredgewidth=marker_width, label='Group 3')
##inset_ax.plot(x3, k3_array, 'k^', markersize=marker_size_in, markeredgewidth=marker_width, label='Group 4')
##inset_ax.plot(x4, k4_array, 'ks', markersize=marker_size_in, markeredgewidth=marker_width, label='Group 5')
#inset_ax.plot(x1, k5_array, 'kD', markersize=marker_size_in, markeredgewidth=marker_width, label='Group 6')
#inset_ax.plot(x, k6_array, 'kX', markersize=marker_size_in, markeredgewidth=marker_width, label='Group 7')
##inset_ax.plot(x, k7_array, 'rP', markersize=marker_size_in, markeredgewidth=marker_width, label='Group 8')

inset_ax.xaxis.tick_top()
inset_ax.xaxis.set_label_position("top")
inset_ax.set_xlabel('$\\sigma_d$', fontsize=label_size)
inset_ax.set_ylabel('$S_d$', fontsize=label_size)
inset_ax.tick_params(axis='both', labelsize=ticks_size)

inset_ax.legend(fontsize=ticks_size,
		  		loc='upper right',
		  		bbox_to_anchor=(1.0, 2),
		  		)

fig = plt.figure(4, figsize=(8, 6), constrained_layout=True)
ax = fig.gca()

axis_type = "loglog"

dictkey = f"d300stdd50_{tailout}"
x = np.array(M_rslt[dictkey]["y"])
y = np.array(Q_rslt[dictkey]["y"])
ax.plot(x, y, 'o', color='r', label='NR1', markersize=marker_size, markerfacecolor='none', markeredgewidth=marker_width)

dictkey = f"d300stdd100_{tailout}"
x = np.array(M_rslt[dictkey]["y"])
y = np.array(Q_rslt[dictkey]["y"])
ax.plot(x, y, 'o', color='g', label='NR1', markersize=marker_size, markerfacecolor='none', markeredgewidth=marker_width)

dictkey = f"d300stdd200_{tailout}"
x = np.array(M_rslt[dictkey]["y"])
y = np.array(Q_rslt[dictkey]["y"])
ax.plot(x, y, 'o', color='b', label='NR1', markersize=marker_size, markerfacecolor='none', markeredgewidth=marker_width)

dictkey = f"d300stdd300_{tailout}"
x = np.array(M_rslt[dictkey]["y"])
y = np.array(Q_rslt[dictkey]["y"])
ax.plot(x, y, 'o', color='c', label='NR1', markersize=marker_size, markerfacecolor='none', markeredgewidth=marker_width)


ax.plot(Ca22_data["Sh0"],
		Ca22_data["Q_star0"],
		'X',
		color='red',
		label='Ca22 0',
		markersize=marker_size,
		markerfacecolor='none',
		markeredgewidth=marker_width
		)
ax.plot(Ca22_data["Sh1"],
		Ca22_data["Q_star1"],
		'X',
		color='green',
		label='Ca22 1',
		markersize=marker_size,
		markerfacecolor='none',
		markeredgewidth=marker_width
		)
ax.plot(Ca22_data["Sh3"],
		Ca22_data["Q_star3"],
		'X',
		color='blue',
		label='Ca22 3',
		markersize=marker_size,
		markerfacecolor='none',
		markeredgewidth=marker_width
		)
ax.plot(Ca22_data["Sh5"],
		Ca22_data["Q_star5"],
		'X',
		color='cyan',
		label='Ca22 5',
		markersize=marker_size,
		markerfacecolor='none',
		markeredgewidth=marker_width
		)
ax.legend(fontsize=ticks_size,
		  		loc='upper right',
		  		#bbox_to_anchor=(-1.3, -0.05),
		  		)

#rslt_dict_smonod3_ex = np.load('rb_vs_sigma_monod2noiter_ex.npz')
#rslt_dict_smonod3_ex_sd1 = np.load('rb_vs_sigma_monod2noiter_ex_smalld1.npz')
#rslt_dict_3D_ex = np.load('rb_vs_sigma_3D_ex.npz')
#ex_array_3D_ex_10 = rslt_dict_3D_ex['ex_10'][::4]
#ez_array_3D_ex_10 = rslt_dict_3D_ex['ez_10'][::4]
#sigma_array = rslt_dict_3D_ex['sigma'][::4]
#th = 10/180*np.pi
#mu_b_10 = (1.0 - ex_array_3D_ex_10)/(ez_array_3D_ex_10 - 1.0)/np.tan(th)
#mu_b_10_sp = (1.0 - rslt_dict_smonod3_ex["ex_10"])/(rslt_dict_smonod3_ex["ez_10"] - 1.0)/np.tan(th)
#mu_b_10_sp_sd1 = (1.0 - rslt_dict_smonod3_ex_sd1["ex_10"])/(rslt_dict_smonod3_ex_sd1["ez_10"] - 1.0)/np.tan(th)
#k = 1/mu_b_10
#k_sp = 1/mu_b_10_sp
#k_sp_sd1 = 1/mu_b_10_sp_sd1
#
#k0 = M_rslt[f"d300stdd50_{tailin}"]["fit"][0]
#k1 = M_rslt[f"d300stdd100_{tailin}"]["fit"][0]
#k2 = M_rslt[f"d300stdd200_{tailin}"]["fit"][0]
#k3 = M_rslt[f"d300stdd300_{tailin}"]["fit"][0]
#k10 = M_rslt[f"d250stdd25_{tailin}"]["fit"][0]
#k11 = M_rslt[f"d271stdd121_{tailin}"]["fit"][0]
#k12 = M_rslt[f"d317stdd252_{tailin}"]["fit"][0]
#k13 = M_rslt[f"d347stdd537_{tailin}"]["fit"][0]
#k20 = M_rslt[f"d300stdd50_m_{tailin}"]["fit"][0]
#k21 = M_rslt[f"d300stdd100_m_{tailin}"]["fit"][0]
#k22 = M_rslt[f"d300stdd200_m_{tailin}"]["fit"][0]
#k23 = M_rslt[f"d300stdd300_m_{tailin}"]["fit"][0]
#k30 = M_rslt[f"d250stdd25_m_{tailin}"]["fit"][0]
#k31 = M_rslt[f"d271stdd121_m_{tailin}"]["fit"][0]
#k32 = M_rslt[f"d317stdd252_m_{tailin}"]["fit"][0]
#k33 = M_rslt[f"d347stdd537_m_{tailin}"]["fit"][0]
#k_list = [k0, k1, k2, k3]
#k1_list = [k10, k11, k12, k13]
#k2_list = [k20, k21, k22, k23]
#k3_list = [k30, k31, k32, k33]
#sigma_array = np.array([0.1655, 0.3246, 0.6167, 0.87])
#sigma_array1 = np.array([0.1, 0.4, 0.7, 1])
#sd_array = np.array([50, 100, 200, 300])
#sd_array1 = np.array([25, 121, 252, 537])
#dm_array1 = np.array([250, 271, 317, 347])
#cv_array = sd_array/300
#cv_array1 = sd_array1/dm_array1
#x, x1 = cv_array, cv_array1
#
##ax.plot(sigma_array, k, 'r-')
##ax.plot(rslt_dict_smonod3_ex["sigma"], k_sp, 'r--')
##ax.plot(rslt_dict_smonod3_ex_sd1["sigma"], k_sp_sd1, 'r:')
#ax.plot(x, k_list, 'ko', markersize=marker_size, markerfacecolor='none', markeredgewidth=marker_width)
#ax.plot(x1, k1_list, 'yo', markersize=marker_size, markerfacecolor='none', markeredgewidth=marker_width)
#ax.plot(x, k2_list, 'go', markersize=marker_size, markerfacecolor='none', markeredgewidth=marker_width)
#ax.plot(x1, k3_list, 'bo', markersize=marker_size, markerfacecolor='none', markeredgewidth=marker_width)

fig = plt.figure(5, figsize=(8, 6), constrained_layout=True)
ax = fig.gca()


legend_str = "NR1"
dictkey = f"d300stdd50_{tailout}"
#style = {"color": "r", "marker": "o", "fill": "full"}
#arrow = True
#fitline = True
#lims["x_loc"] = 0.015
#lims["y_loc"] = 0.001
#rslt_plot(ax, M_rslt, dictkey, legend_str, fitline, arrow, lims, sizes, style, axis_type)
vxt = 1 #2.0/0.4*np.sqrt(ks0)
Qvxt0 = Q_rslt[dictkey]["y"]/vxt
#Qvxt = Q_rslt[dictkey]["y"]/M_rslt[dictkey]["y"]
ax.plot(M_rslt[dictkey]["y"], Qvxt0, 'o', color='r', label=legend_str, markersize=marker_size)

legend_str = "MD1"
dictkey = f"d300stdd100_{tailout}"
#style = {"color": "g", "marker": "^", "fill": "full"}
#arrow = False
#fitline = False
#rslt_plot(ax, M_rslt, dictkey, legend_str, fitline, arrow, lims, sizes, style, axis_type)
vxt = 1 #2.0/0.4*np.sqrt(ks1)
Qvxt1 = Q_rslt[dictkey]["y"]/vxt
ax.plot(M_rslt[dictkey]["y"], Qvxt1, '^', color='g', label=legend_str, markersize=marker_size)

legend_str = "WD1"
dictkey = f"d300stdd200_{tailout}"
#style = {"color": "b", "marker": "s", "fill": "full"}
#rslt_plot(ax, M_rslt, dictkey, legend_str, fitline, arrow, lims, sizes, style, axis_type)
vxt = 1 #2.0/0.4*np.sqrt(ks2)
Qvxt2 = Q_rslt[dictkey]["y"]/vxt
ax.plot(M_rslt[dictkey]["y"], Qvxt2, 's', color='b', label=legend_str, markersize=marker_size)

legend_str = "VWD1"
dictkey = f"d300stdd300_{tailout}"
#style = {"color": "c", "marker": "h", "fill": "full"}
#arrow = True
#fitline = True
#lims["x_loc"] = 0.007
#lims["y_loc"] = 0.007
#rslt_plot(ax, M_rslt, dictkey, legend_str, fitline, arrow, lims, sizes, style, axis_type)
vxt = 1 #2.0/0.4*np.sqrt(ks3)
Qvxt3 = Q_rslt[dictkey]["y"]/vxt
ax.plot(M_rslt[dictkey]["y"], Qvxt3, 'h', color='c', label=legend_str, markersize=marker_size)

legend_str = "NR2"
dictkey = f"d250stdd25_{tailout}"
#style = {"color": "r", "marker": "o", "fill": "none"}
#arrow = False
#fitline = False
#rslt_plot(ax, M_rslt, dictkey, legend_str, fitline, arrow, lims, sizes, style, axis_type)
vxt = 1 #2.0/0.4*np.sqrt(ks10)
Qvxt10 = Q_rslt[dictkey]["y"]/vxt
ax.plot(M_rslt[dictkey]["y"], Qvxt10, 'o', color='r', label=legend_str, markersize=marker_size, markerfacecolor='none', markeredgewidth=marker_width)

legend_str = "MD2"
dictkey = f"d271stdd121_{tailout}"
#style = {"color": "g", "marker": "^", "fill": "none"}
#rslt_plot(ax, M_rslt, dictkey, legend_str, fitline, arrow, lims, sizes, style, axis_type)
vxt = 1 #2.0/0.4*np.sqrt(ks11)
Qvxt11 = Q_rslt[dictkey]["y"]/vxt
ax.plot(M_rslt[dictkey]["y"], Qvxt11, '^', color='g', label=legend_str, markersize=marker_size, markerfacecolor='none', markeredgewidth=marker_width)

legend_str = "WD2"
dictkey = f"d317stdd252_{tailout}"
#style = {"color": "b", "marker": "s", "fill": "none"}
#rslt_plot(ax, M_rslt, dictkey, legend_str, fitline, arrow, lims, sizes, style, axis_type)
vxt = 1 #2.0/0.4*np.sqrt(ks12)
Qvxt12 = Q_rslt[dictkey]["y"]/vxt
ax.plot(M_rslt[dictkey]["y"], Qvxt12, 's', color='b', label=legend_str, markersize=marker_size, markerfacecolor='none', markeredgewidth=marker_width)

legend_str = "VWD2"
dictkey = f"d347stdd537_{tailout}"
#style = {"color": "c", "marker": "h", "fill": "none"}
#rslt_plot(ax, M_rslt, dictkey, legend_str, fitline, arrow, lims, sizes, style, axis_type)
vxt = 1 #2.0/0.4*np.sqrt(ks13)
Qvxt13 = Q_rslt[dictkey]["y"]/vxt
ax.loglog(M_rslt[dictkey]["y"], Qvxt13, 'h', color='c', label=legend_str, markersize=marker_size, markerfacecolor='none', markeredgewidth=marker_width)

legend_str = "MU1"
dictkey = f"d167stdd100_{tailout}"
Q_20 = Q_rslt[dictkey]["y"]
M_20 = M_rslt[dictkey]["y"]
ax.plot(M_20, Q_20, 'd', color='gray', label=legend_str, markersize=marker_size, markerfacecolor='none', markeredgewidth=marker_width)

legend_str = "MU2"
dictkey = f"d197stdd65_{tailout}"
Q_21 = Q_rslt[dictkey]["y"]
M_21 = M_rslt[dictkey]["y"]
ax.plot(M_21, Q_21, 'd', color='gray', label=legend_str, markersize=marker_size, markerfacecolor='none', markeredgewidth=marker_width)

legend_str = "MU3"
dictkey = f"d290stdd97_{tailout}"
Q_22 = Q_rslt[dictkey]["y"]
M_22 = M_rslt[dictkey]["y"]
ax.plot(M_22, Q_22, 'd', color='gray', label=legend_str, markersize=marker_size, markerfacecolor='none', markeredgewidth=marker_width)

legend_str = "MU4"
dictkey = f"d430stdd100_{tailout}"
Q_23 = Q_rslt[dictkey]["y"]
M_23 = M_rslt[dictkey]["y"]
ax.plot(M_23, Q_23, 'd', color='gray', label=legend_str, markersize=marker_size, markerfacecolor='none', markeredgewidth=marker_width)

legend_str = "RT"
dictkey = f"d269stdd100_{tailout}"
Q_30 = Q_rslt[dictkey]["y"]
M_30 = M_rslt[dictkey]["y"]
ax.plot(M_30, Q_30, 'd', color='gray', label=legend_str, markersize=marker_size, markerfacecolor='none', markeredgewidth=marker_width)

legend_str = "LT"
dictkey = f"d321stdd100_{tailout}"
Q_31 = Q_rslt[dictkey]["y"]
M_31 = M_rslt[dictkey]["y"]
ax.plot(M_31, Q_31, 'd', color='gray', label=legend_str, markersize=marker_size, markerfacecolor='none', markeredgewidth=marker_width)

legend_str = "NR3"
dictkey = f"d240stdd50_{tailout}"
Q_40 = Q_rslt[dictkey]["y"]
M_40 = M_rslt[dictkey]["y"]
ax.plot(M_40, Q_40, 'd', color='gray', label=legend_str, markersize=marker_size, markerfacecolor='none', markeredgewidth=marker_width)

legend_str = "NR4"
dictkey = f"d400stdd50_{tailout}"
Q_41 = Q_rslt[dictkey]["y"]
M_41 = M_rslt[dictkey]["y"]
ax.plot(M_41, Q_41, 'd', color='gray', label=legend_str, markersize=marker_size, markerfacecolor='none', markeredgewidth=marker_width)

x = np.concatenate([
					M_rslt[f"d300stdd50_{tailout}"]["y"],
                    M_rslt[f"d300stdd100_{tailout}"]["y"],
                    M_rslt[f"d300stdd200_{tailout}"]["y"],
                    M_rslt[f"d300stdd300_{tailout}"]["y"],
                    M_rslt[f"d250stdd25_{tailout}"]["y"],
                    M_rslt[f"d271stdd121_{tailout}"]["y"],
                    M_rslt[f"d317stdd252_{tailout}"]["y"],
                    M_rslt[f"d347stdd537_{tailout}"]["y"],
					M_rslt[f"d167stdd100_{tailout}"]["y"],
					M_rslt[f"d197stdd65_{tailout}"]["y"],
					M_rslt[f"d290stdd97_{tailout}"]["y"],
					M_rslt[f"d430stdd100_{tailout}"]["y"],
					M_rslt[f"d269stdd100_{tailout}"]["y"],
					M_rslt[f"d321stdd100_{tailout}"]["y"],
					M_rslt[f"d240stdd50_{tailout}"]["y"],
					M_rslt[f"d400stdd50_{tailout}"]["y"],
					])
y = np.concatenate([Qvxt0, Qvxt1, Qvxt2, Qvxt3, Qvxt10, Qvxt11, Qvxt12, Qvxt13, Q_20, Q_21, Q_22, Q_23, Q_30, Q_31, Q_40, Q_41])
fit = np.polyfit(x, y, 1)
x_array = np.linspace(lims["x_min"], lims["x_max"], 1000)
y_array = np.polyval(fit, x_array)
ax.plot(x_array, y_array, 'k--', label='fit', linewidth=line_width)

ax.set_xlabel(r'$\widetilde{M}$', fontsize=label_size)
ax.set_ylabel(r'$\widetilde{Q}$', fontsize=label_size)
ax.tick_params(axis='both', labelsize=ticks_size)

ax.set_xlim(0.01, 0.15)
ax.set_ylim(0.002, 0.05)
ax.xaxis.set_major_formatter(ScalarFormatter())
ax.yaxis.set_major_formatter(ScalarFormatter())
ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
ax.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
ax.set_xticks([0.01, 0.05, 0.1])
ax.set_yticks([0.005, 0.01, 0.05])


proxy_ro = plt.Line2D([0], [0],
					 color='r',
					 marker='o',
					 linestyle='',
					 markersize=marker_size
					)
proxy_ro_h = plt.Line2D([0], [0],
					 color='r',
					 marker='o',
					 linestyle='',
					 markerfacecolor='none',
					 markersize=marker_size,
					 markeredgewidth=marker_width
					)
proxy_gt = plt.Line2D([0], [0],
					 color='g',
					 marker='^',
					 linestyle='',
					 markersize=marker_size
					)
proxy_gt_h = plt.Line2D([0], [0],
					 color='g',
					 marker='^',
					 linestyle='',
					 markerfacecolor='none',
					 markersize=marker_size,
					 markeredgewidth=marker_width
					)
proxy_bs = plt.Line2D([0], [0],
					 color='b',
					 marker='s',
					 linestyle='',
					 markersize=marker_size
					)
proxy_bs_h = plt.Line2D([0], [0],
					 color='b',
					 marker='s',
					 linestyle='',
					 markerfacecolor='none',
					 markersize=marker_size,
					 markeredgewidth=marker_width
					)
proxy_ch = plt.Line2D([0], [0],
					 color='c',
					 marker='h',
					 linestyle='',
					 markersize=marker_size
					)
proxy_ch_h = plt.Line2D([0], [0],
					 color='c',
					 marker='h',
					 linestyle='',
					 markerfacecolor='none',
					 markersize=marker_size,
					 markeredgewidth=marker_width
					)
proxy_gd = plt.Line2D([0], [0],
					 color='gray',
					 marker='d',
					 linestyle='',
					 markersize=marker_size
					)
proxy_gd_h = plt.Line2D([0], [0],
					 color='gray',
					 marker='d',
					 linestyle='',
					 markerfacecolor='none',
					 markersize=marker_size,
					 markeredgewidth=marker_width
					)
Cr_09 = ax.scatter([], [], marker='x', color='purple', s=marker_size**2, linewidths=marker_width)
Ho_12 = ax.scatter([], [], marker='+', color='brown', s=marker_size**2, linewidths=marker_width)
Ho_12_c = ax.scatter([], [], marker='.', color='k', s=marker_size**2, linewidths=marker_width)
ax.legend([(proxy_ro, proxy_ro_h),
		   (proxy_gt, proxy_gt_h),
		   (proxy_bs, proxy_bs_h),
		   (proxy_ch, proxy_ch_h),
		   (proxy_gd, proxy_gd_h),
		   Cr_09,
		   Ho_12,
		   Ho_12_c],
          ['NR1, NR2',
		   'MD1, MD2',
		   'WD1, WD2',
		   'VWD1, VWD2',
		   'Other Groups',
		   'Ho12 Roughness',
		   'Ma18 TFEM',
		   'Zh19 Roughness'],
		  handler_map={tuple: HandlerTuple(ndivide=None)},
          fontsize=ticks_size,
          loc='lower right',
          bbox_to_anchor=(1.01, -0.01),
          frameon=True,
		  framealpha=1,
		  ncol=1
		  )

inset_ax = inset_axes(ax, width='45%', height='45%', loc='upper left')

ks0 = -M_rslt[f"d300stdd50_{tailout}"]["fit"][1]/M_rslt[f"d300stdd50_{tailout}"]["fit"][0]
ks1 = -M_rslt[f"d300stdd100_{tailout}"]["fit"][1]/M_rslt[f"d300stdd100_{tailout}"]["fit"][0]
ks2 = -M_rslt[f"d300stdd200_{tailout}"]["fit"][1]/M_rslt[f"d300stdd200_{tailout}"]["fit"][0]
ks3 = -M_rslt[f"d300stdd300_{tailout}"]["fit"][1]/M_rslt[f"d300stdd300_{tailout}"]["fit"][0]
ks10 = -M_rslt[f"d250stdd25_{tailout}"]["fit"][1]/M_rslt[f"d250stdd25_{tailout}"]["fit"][0]
ks11 = -M_rslt[f"d271stdd121_{tailout}"]["fit"][1]/M_rslt[f"d271stdd121_{tailout}"]["fit"][0]
ks12 = -M_rslt[f"d317stdd252_{tailout}"]["fit"][1]/M_rslt[f"d317stdd252_{tailout}"]["fit"][0]
ks13 = -M_rslt[f"d347stdd537_{tailout}"]["fit"][1]/M_rslt[f"d347stdd537_{tailout}"]["fit"][0]
ks20 = -M_rslt[f"d167stdd100_{tailout}"]["fit"][1]/M_rslt[f"d167stdd100_{tailout}"]["fit"][0]
ks21 = -M_rslt[f"d197stdd65_{tailout}"]["fit"][1]/M_rslt[f"d197stdd65_{tailout}"]["fit"][0]
ks22 = -M_rslt[f"d290stdd97_{tailout}"]["fit"][1]/M_rslt[f"d290stdd97_{tailout}"]["fit"][0]
ks23 = -M_rslt[f"d430stdd100_{tailout}"]["fit"][1]/M_rslt[f"d430stdd100_{tailout}"]["fit"][0]
ks30 = -M_rslt[f"d269stdd100_{tailout}"]["fit"][1]/M_rslt[f"d269stdd100_{tailout}"]["fit"][0]
ks31 = -M_rslt[f"d321stdd100_{tailout}"]["fit"][1]/M_rslt[f"d321stdd100_{tailout}"]["fit"][0]
ks40 = -M_rslt[f"d240stdd50_{tailout}"]["fit"][1]/M_rslt[f"d240stdd50_{tailout}"]["fit"][0]
ks41 = -M_rslt[f"d400stdd50_{tailout}"]["fit"][1]/M_rslt[f"d400stdd50_{tailout}"]["fit"][0]

x0, x1, x2, x3, x4 = d90500, d90501, d90502, d90503, d90504

y0 = np.array([ks0, ks1, ks2, ks3])
y1 = np.array([ks10, ks11, ks12, ks13])
y2 = np.array([ks20, ks21, ks22, ks23])
y3 = np.array([ks30, ks31])
y4 = np.array([ks40, ks41])

c_list0 = ['r', 'g', 'b', 'c']
c_list1 = ['r', 'g', 'b', 'c']
s_list0 = ['o', '^', 's', 'h']
s_list1 = ['o', '^', 's', 'h']
inset_ax.plot(x2, y2, 'd', color='grey', markersize=marker_size_in)
inset_ax.plot(x3, y3, 'd', color='grey', markersize=marker_size_in)
inset_ax.plot(x4, y4, 'd', color='grey', markersize=marker_size_in)
for x, y, s, c in zip(x0, y0, s_list0, c_list0):
    inset_ax.plot(x, y, s, color=c, markersize=marker_size_in)
for x, y, s, c in zip(x1, y1, s_list1, c_list1):
    inset_ax.plot(x, y, s, color=c, markersize=marker_size_in, markerfacecolor='none', markeredgewidth=marker_width)
inset_ax.plot([1.63, 2.01], [Shd_ho12, Shd_ho12c], 'x', color='purple', markersize=marker_size_in, markeredgewidth=marker_width, label='Ho12 Roughness')
inset_ax.plot(Zh19_data["d9050"], Zh19_data["Shd"], '.', color='k', markersize=marker_size_in, markeredgewidth=marker_width, label='Zh19 Roughness')
inset_ax.semilogy(Ma18_data["d9050"], Ma18_data["Shd"], '+', color='brown', markersize=marker_size_in, markeredgewidth=marker_width, label='Ma18 TFEM')

inset_ax.set_xlabel('$\\eta_d$', fontsize=label_size)
inset_ax.set_ylabel('$S_d$', fontsize=label_size)
inset_ax.tick_params(axis='both', labelsize=ticks_size)
# 把y的tick移到右边
inset_ax.yaxis.tick_right()
inset_ax.yaxis.set_label_position("right")
inset_ax.set_xlim(0.9, 2.6)
inset_ax.set_ylim(0.0008, 0.02)
#inset_ax.legend(
#          fontsize=ticks_size,
#          loc='upper right',
#          bbox_to_anchor=(2.22, 1.1),
#          frameon=True,
#		  framealpha=1,
#		  ncol=2
#		  )

plt.show()