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
				xytext=(x_intercept, max(y)/4),
				fontsize=sizes["ticks"],
				arrowprops=dict(arrowstyle='->'))
		ax.plot(x_fit, fit_line(x_fit), '--', color=style["color"], linewidth=sizes["line"])

# 控制参数
ampf = 1/0.6 # amplification factor
label_size = 15*ampf
ticks_size = 10*ampf
ticks_size_in = 6*ampf
marker_size = 10*ampf
marker_size_in = 6*ampf
line_width = 1.5
marker_width = 2
sizes = {
	"marker": marker_size,
	"line": line_width,
	"ticks": ticks_size,
	"markerwidth": marker_width
}
sizes_in = {
	"marker": marker_size_in,
	"line": line_width,
	"markerwidth": marker_width
}
use_lims = True # 是否使用坐标轴限制
lims = {
	"x_min": 1e-2,
	"x_max": 0.4,
	"y_min": 3e-3,
	"y_max": 0.2,
}
lims_in = {
	"x_min": 1e-2,
	"x_max": 0.2,
	"y_min": 3e-3,
	"y_max": 0.2,
}
tailin = "dair"
tailout = "d50"
working_dir = "E:/Data/Q_on_flat_bed"
dmin = 1.0e-4
dmax = 10.0e-4
bin_num = 1000
dx = dmax / bin_num
ddx = dx * 0.5
mu_nr = -8.1254
sigma_nr = 0.1655
mu_vwd = -8.5521
sigma_vwd = 0.87

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
d_ho12 = 2.30e-4
rhop_ho12 = 2450
rhof_ho12 = 1.2
s_ho12 = rhop_ho12 / rhof_ho12
g_ho12 = 9.81
Sh_array = np.array([0.01776, 0.04970, 0.08710, 0.13516, 0.17804, 0.32099])
Q_array = np.array([0.35012, 1.12710, 2.34532, 3.64988, 5.59712, 9.11271])
Q_array = rhop_ho12 * np.sqrt(g_ho12*d_ho12) * d_ho12 * Q_array
Ho2012_data = {
	"Sh": Sh_array,
	"Q": Q_array,
	"Q_err": 0.05 * Q_array,
	"Q_star": Q_array/(rhop_ho12*d_ho12*np.sqrt(s_ho12*g_ho12*d_ho12)),
	"Q_star_err": 0.05 * Q_array
}

# Ho 2012 coarse 的数据
d_ho12 = 6.30e-4
rhop_ho12 = 2450
rhof_ho12 = 1.2
s_ho12 = rhop_ho12 / rhof_ho12
g_ho12 = 9.81
Sh_array = np.array([0.02562, 0.02958, 0.03992, 0.04560, 0.05019, 0.05115, 0.06251, 0.07904])
Q_array = np.array([0.5800, 0.6744, 1.2739, 1.7121, 1.9586, 2.1716, 2.4790, 2.8351])
Q_array = rhop_ho12 * np.sqrt(g_ho12*d_ho12) * d_ho12 * Q_array
Ho2012_c_data = {
	"Sh": Sh_array,
	"Q": Q_array,
	"Q_err": 0.05 * Q_array,
	"Q_star": Q_array/(rhop_ho12*d_ho12*np.sqrt(s_ho12*g_ho12*d_ho12)),
	"Q_star_err": 0.05 * Q_array
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
Q_lm = np.load(f"Q_d430stdd100_{tailin}.npz")
M_lm = np.load(f"M_d430stdd100_{tailin}.npz")
Q_sm = np.load(f"Q_d167stdd100_{tailin}.npz")
M_sm = np.load(f"M_d167stdd100_{tailin}.npz")
Q_rt = np.load(f"Q_d269stdd100_{tailin}.npz")
M_rt = np.load(f"M_d269stdd100_{tailin}.npz")
Q_lt = np.load(f"Q_d321stdd100_{tailin}.npz")
M_lt = np.load(f"M_d321stdd100_{tailin}.npz")
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


dir_list_vwd1 = [
	"uStar030_347log537_0_2650_300",
	"uStar040_347log537_0_2650_300",
	"uStar050_347log537_0_2650_300",
	"uStar060_347log537_0_2650_300",
]

d_hist_addr_nr = f"{working_dir}/{dir_list_nr[-1]}/d_in_air_hist.npz"
d_hist_addr_vwd = f"{working_dir}/{dir_list_vwd[-1]}/d_in_air_hist.npz"
d_hist_nr_dict = np.load(d_hist_addr_nr)
d_hist_vwd_dict = np.load(d_hist_addr_vwd)
d_hist_nr = d_hist_nr_dict["hist"]
d_hist_vwd = d_hist_vwd_dict["hist"]
d_bin_centers_nr = d_hist_nr_dict["bin_centers"]
d_bin_centers_vwd = d_hist_vwd_dict["bin_centers"]
d_bin_edges_nr = d_hist_nr_dict["bin_edges"]
d_bin_edges_vwd = d_hist_vwd_dict["bin_edges"]
dm_nr = d_hist_nr_dict["d_mean"]
dm_vwd = d_hist_vwd_dict["d_mean"]
x = np.linspace(dmin+dx, dmax, bin_num)
x = x - ddx
# 对数正态分布概率密度函数
y = (1 / (x * sigma_nr * np.sqrt(2 * np.pi))) * np.exp(-(np.log(x) - mu_nr) ** 2 / (2 * sigma_nr ** 2))
# 归一化
y = y / y.sum()
d_bin_centers_nr_bed = x
d_hist_nr_bed = y
# 对数正态分布概率密度函数
y = (1 / (x * sigma_vwd * np.sqrt(2 * np.pi))) * np.exp(-(np.log(x) - mu_vwd) ** 2 / (2 * sigma_vwd ** 2))
# 归一化
y = y / y.sum()
d_bin_centers_vwd_bed = x
d_hist_vwd_bed = y

dm_nr_list = []
d50_nr_list = []
d90_nr_list = []
for i in range(len(dir_list_nr)):
	folder = dir_list_nr[i]
	d_hist_nr_dict = np.load(f"{working_dir}/{folder}/d_in_air_hist.npz")
	dm_nr_list.append(d_hist_nr_dict["d_mean"])
	d50_nr_list.append(d_hist_nr_dict["d50"])
	d90_nr_list.append(d_hist_nr_dict["d90"])

dm_vwd_list = []
d50_vwd_list = []
d90_vwd_list = []
for i in range(len(dir_list_vwd)):
	folder = dir_list_vwd[i]
	d_hist_vwd_dict = np.load(f"{working_dir}/{folder}/d_in_air_hist.npz")
	dm_vwd_list.append(d_hist_vwd_dict["d_mean"])
	d50_vwd_list.append(d_hist_vwd_dict["d50"])
	d90_vwd_list.append(d_hist_vwd_dict["d90"])

dm_vwd_list1 = []
d50_vwd_list1 = []
d90_vwd_list1 = []
for i in range(len(dir_list_vwd1)):
	folder = dir_list_vwd1[i]
	d_hist_vwd_dict = np.load(f"{working_dir}/{folder}/d_in_air_hist.npz")
	dm_vwd_list1.append(d_hist_vwd_dict["d_mean"])
	d50_vwd_list1.append(d_hist_vwd_dict["d50"])
	d90_vwd_list1.append(d_hist_vwd_dict["d90"])

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
}

# 绘图
fig = plt.figure(1, figsize=(8, 6), constrained_layout=True)
ax = fig.gca()

ax.plot(others_case_dict["Cr09"]["Sh"],
		others_case_dict["Cr09"]["Q_star"],
		'x',
		color='orange',
		label='Cr09',
		markersize=marker_size,
		markerfacecolor='none',
		markeredgewidth=marker_width
		)
ax.plot(others_case_dict["Ho2012"]["Sh"],
		others_case_dict["Ho2012"]["Q_star"],
		'*',
		color='purple',
		label='Ho12 Fine',
		markersize=marker_size,
		markerfacecolor='none',
		markeredgewidth=marker_width
		)
ax.plot(others_case_dict["Ho2012_c"]["Sh"],
		others_case_dict["Ho2012_c"]["Q_star"],
		'+',
		color='brown',
		label='Ho12 Coarse',
		markersize=marker_size,
		markerfacecolor='none',
		markeredgewidth=marker_width
		)

axis_type = "loglog"

legend_str = "NR1"
dictkey = f"d300stdd50_{tailout}"
style = {"color": "r", "marker": "o", "fill": "full"}
arrow = False
fitline = True
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
fitline = True
rslt_plot(ax, Q_rslt, dictkey, legend_str, fitline, arrow, lims, sizes, style, axis_type)

legend_str = "NR2"
dictkey = f"d250stdd25_{tailout}"
style = {"color": "r", "marker": "o", "fill": "none"}
arrow = False
fitline = False
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
fitline = False
rslt_plot(ax, Q_rslt, dictkey, legend_str, fitline, arrow, lims, sizes, style, axis_type)

rect = patches.Rectangle(
	(0.315, 0.5), # (x, y)
	0.1, # width
	0.32, # height
	transform=ax.transAxes,
	angle = -45,
	edgecolor='black',
	facecolor='none',
	linestyle='--',
	linewidth=line_width*1.5
)
ax.add_patch(rect)

inset_ax = inset_axes(ax, width='42%', height='42%', loc='lower right')

#inset_ax.plot(others_case_dict["Cr09"]["Sh"],
#			  others_case_dict["Cr09"]["Q_star"],
#			  'x',
#			  color='orange',
#			  label='Cr09',
#			  markersize=marker_size_in,
#			  markerfacecolor='none',
#			  markeredgewidth=marker_width
#			  )
#inset_ax.plot(others_case_dict["Ho2012"]["Sh"],
#			  others_case_dict["Ho2012"]["Q_star"],
#			  '*',
#			  color='purple',
#			  label='Ho12 Fine',
#			  markersize=marker_size_in,
#			  markerfacecolor='none',
#			  markeredgewidth=marker_width
#			  )
#inset_ax.plot(others_case_dict["Ho2012_c"]["Sh"],
#			  others_case_dict["Ho2012_c"]["Q_star"],
#			  '+',
#			  color='brown',
#			  label='Ho12 Coarse',
#			  markersize=marker_size_in,
#			  markerfacecolor='none',
#			  markeredgewidth=marker_width
#			  )

legend_str = "NR1"
dictkey = f"d300stdd50_{tailin}"
style = {"color": "r", "marker": "o", "fill": "full"}
arrow = False
fitline = True
rslt_plot(inset_ax, Q_rslt, dictkey, legend_str, fitline, arrow, lims_in, sizes_in, style, axis_type)

legend_str = "MD1"
dictkey = f"d300stdd100_{tailin}"
style = {"color": "g", "marker": "^", "fill": "full"}
arrow = False
fitline = False
rslt_plot(inset_ax, Q_rslt, dictkey, legend_str, fitline, arrow, lims_in, sizes_in, style, axis_type)

legend_str = "WD1"
dictkey = f"d300stdd200_{tailin}"
style = {"color": "b", "marker": "s", "fill": "full"}
arrow = False
fitline = False
rslt_plot(inset_ax, Q_rslt, dictkey, legend_str, fitline, arrow, lims_in, sizes_in, style, axis_type)

legend_str = "VWD1"
dictkey = f"d300stdd300_{tailin}"
style = {"color": "c", "marker": "h", "fill": "full"}
arrow = False
fitline = True
rslt_plot(inset_ax, Q_rslt, dictkey, legend_str, fitline, arrow, lims_in, sizes_in, style, axis_type)

legend_str = "NR2"
dictkey = f"d250stdd25_{tailin}"
style = {"color": "r", "marker": "o", "fill": "none"}
arrow = False
fitline = False
rslt_plot(inset_ax, Q_rslt, dictkey, legend_str, fitline, arrow, lims_in, sizes_in, style, axis_type)

legend_str = "MD2"
dictkey = f"d271stdd121_{tailin}"
style = {"color": "g", "marker": "^", "fill": "none"}
arrow = False
fitline = False
rslt_plot(inset_ax, Q_rslt, dictkey, legend_str, fitline, arrow, lims_in, sizes_in, style, axis_type)

legend_str = "WD2"
dictkey = f"d317stdd252_{tailin}"
style = {"color": "b", "marker": "s", "fill": "none"}
arrow = False
fitline = False
rslt_plot(inset_ax, Q_rslt, dictkey, legend_str, fitline, arrow, lims_in, sizes_in, style, axis_type)

legend_str = "VWD2"
dictkey = f"d347stdd537_{tailin}"
style = {"color": "c", "marker": "h", "fill": "none"}
arrow = False
fitline = False
rslt_plot(inset_ax, Q_rslt, dictkey, legend_str, fitline, arrow, lims_in, sizes_in, style, axis_type)


ax.xaxis.set_major_formatter(ScalarFormatter())
ax.yaxis.set_major_formatter(ScalarFormatter())
ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
ax.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
inset_ax.xaxis.set_major_formatter(ScalarFormatter())
inset_ax.yaxis.set_major_formatter(ScalarFormatter())
inset_ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
inset_ax.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
# 把y的tick移到右边
#inset_ax.yaxis.tick_right()
#inset_ax.yaxis.set_label_position("right")
inset_ax.xaxis.tick_top()
inset_ax.xaxis.set_label_position("top")
inset_ax.set_xlabel('$S^{\\mathrm{Ea}}$', fontsize=label_size)
inset_ax.set_ylabel(r'$\widetilde{Q}^{\mathrm{Ea}}$', fontsize=label_size)

if use_lims:
    ax.set_xlim(lims["x_min"], lims["x_max"])
    ax.set_ylim(lims["y_min"], lims["y_max"])
    inset_ax.set_xlim(lims_in["x_min"], lims_in["x_max"])
    inset_ax.set_ylim(lims_in["y_min"], lims_in["y_max"])
ax.set_xlabel('$S^{d50}$', fontsize=label_size)
ax.set_ylabel(r'$\widetilde{Q}^{d50}$', fontsize=label_size)
ax.tick_params(axis='both', labelsize=ticks_size)
ax.set_xticks([0.01, 0.05, 0.1, 0.4])
ax.set_yticks([0.005, 0.01, 0.05, 0.1, 0.3])
inset_ax.tick_params(axis='both', labelsize=ticks_size)

fill_ro = ax.plot([], [], 'ro', markersize=marker_size)[0]
fill_gt = ax.plot([], [], 'g^', markersize=marker_size)[0]
fill_bs = ax.plot([], [], 'bs', markersize=marker_size)[0]
fill_ch = ax.plot([], [], 'ch', markersize=marker_size)[0]
Cr_09 = ax.scatter([], [], marker='x', color='orange', s=marker_size**2, linewidths=marker_width)
Ho_12 = ax.scatter([], [], marker='*', color='purple', s=marker_size**2, facecolors='none', linewidths=marker_width)
Ho_12_c = ax.scatter([], [], marker='+', color='brown', s=marker_size**2, linewidths=marker_width)
ax.legend([fill_ro, fill_gt, fill_bs, fill_ch, Cr_09, Ho_12, Ho_12_c],
          ['NR', 'MD', 'WD', 'VWD', 'Cr09', 'Ho12 Fine', 'Ho12 Coarse'],
          fontsize=ticks_size,
          loc='upper left',
          bbox_to_anchor=(-0.01, 1.05),
          frameon=True,
		  framealpha=1,
		  ncol=1
		  )

# 绘制床面和空中颗粒的CDF
fig = plt.figure(2, figsize=(8, 6), constrained_layout=True)
ax = fig.gca()
# 转为CDF
d_cdf_nr_bed = np.cumsum(d_hist_nr_bed) / np.sum(d_hist_nr_bed)
d_cdf_vwd_bed = np.cumsum(d_hist_vwd_bed) / np.sum(d_hist_vwd_bed)
d_cdf_nr_air = np.cumsum(d_hist_nr) / np.sum(d_hist_nr)
d_cdf_vwd_air = np.cumsum(d_hist_vwd) / np.sum(d_hist_vwd)
ax.plot(d_bin_centers_nr_bed*1e6, d_cdf_nr_bed, '-', color='r', label='NR1 Bed CDF', markersize=marker_size, markerfacecolor='none', markeredgewidth=marker_width)
ax.plot(d_bin_centers_nr*1e6, d_cdf_nr_air, 'o', color='r', label='NR1 Air CDF', markersize=marker_size, markerfacecolor='none', markeredgewidth=marker_width)
ax.plot(d_bin_centers_vwd_bed*1e6, d_cdf_vwd_bed, '-', color='b', label='VWD1 Bed CDF', markersize=marker_size, markerfacecolor='none', markeredgewidth=marker_width)
ax.plot(d_bin_centers_vwd*1e6, d_cdf_vwd_air, 's', color='b', label='VWD1 Air CDF', markersize=marker_size, markerfacecolor='none', markeredgewidth=marker_width)
ax.axvline(dm_nr*1e6, color='r', linestyle='--')
ax.axvline(dm_vwd*1e6, color='b', linestyle='--')
ax.axvline(3e-4*1e6, color='k', linestyle='-')

ax.tick_params(axis='both', labelsize=ticks_size)
ax.legend(fontsize=ticks_size,
		  loc='upper right',
		  bbox_to_anchor=(1.0, 0.9),
		  )

inset_ax = inset_axes(ax, width='45%', height='45%', loc='lower right')
#dictkey = f"d300stdd50_{tailout}"
## x:Sh
#x = Q_rslt[dictkey]["x"]
#inset_ax.plot(x, dm_nr_list, 'o', color='r', label='NR1', markersize=marker_size_in, markerfacecolor='none', markeredgewidth=marker_width)
dictkey = f"d300stdd200_{tailin}"
x = np.array(Q_rslt[dictkey]["x"])
inset_ax.plot(x, np.array(dm_vwd_list)/3e-4, '-', color='b', label='VWD1', markersize=marker_size_in, markerfacecolor='none', markeredgewidth=marker_width)
#inset_ax.plot(x, np.array(dm_nr_list)/3e-4, '-', color='r', label='NR1', markersize=marker_size_in, markerfacecolor='none', markeredgewidth=marker_width)
x_arrow = x[4]
y_arrow = dm_vwd_list[4]/3e-4
inset_ax.annotate('Mean',
				  xy=(x_arrow, y_arrow),
				  xytext=(x_arrow, y_arrow-0.05),
				  ha='center',
				  va='top',
				  fontsize=ticks_size,
				  arrowprops=dict(arrowstyle='->'))
inset_ax.plot(x, np.array(d50_vwd_list)/2.49e-4, '--', color='b', label='VWD1', markersize=marker_size_in, markerfacecolor='none', markeredgewidth=marker_width)
#inset_ax.plot(x, np.array(d50_nr_list)/2.96e-4, '--', color='r', label='NR1', markersize=marker_size_in, markerfacecolor='none', markeredgewidth=marker_width)
x_arrow = x[6]
y_arrow = d50_vwd_list[6]/2.49e-4
inset_ax.annotate('d50',
				  xy=(x_arrow, y_arrow),
				  xytext=(x_arrow, y_arrow-0.1),
				  ha='center',
				  va='top',
				  fontsize=ticks_size,
				  arrowprops=dict(arrowstyle='->'))
inset_ax.plot(x, np.array(d90_vwd_list)/5.89e-4, ':', color='b', label='VWD1', markersize=marker_size_in, markerfacecolor='none', markeredgewidth=marker_width)
#inset_ax.plot(x, np.array(d90_nr_list)/3.66e-4, ':', color='r', label='NR1', markersize=marker_size_in, markerfacecolor='none', markeredgewidth=marker_width)
x_arrow = x[2]
y_arrow = d90_vwd_list[2]/5.89e-4
inset_ax.annotate('d90',
				  xy=(x_arrow, y_arrow),
				  xytext=(x_arrow, y_arrow+0.1),
				  ha = 'center',
				  va = 'bottom',
				  fontsize=ticks_size,
				  arrowprops=dict(arrowstyle='->'))

ax.set_xlabel('$d$ ($\\mu$m)', fontsize=label_size)
ax.set_ylabel('CDF', fontsize=label_size)
ax.tick_params(axis='both', labelsize=ticks_size)
inset_ax.xaxis.tick_top()
inset_ax.xaxis.set_label_position("top")
inset_ax.set_xlabel('$S$', fontsize=label_size)
inset_ax.set_ylabel('$d_{\\mathrm{air}}/d_{\\mathrm{bed}}$', fontsize=label_size)
inset_ax.tick_params(axis='both', labelsize=ticks_size)

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
	f"d430stdd100_{tailin}": M_lm,
	f"d167stdd100_{tailin}": M_sm,
	f"d269stdd100_{tailin}": M_rt,
	f"d321stdd100_{tailin}": M_lt,
	f"d250stdd25_m_{tailin}": M_nr3_in,
	f"d271stdd121_m_{tailin}": M_md3_in,
	f"d317stdd252_m_{tailin}": M_wd3_in,
	f"d347stdd537_m_{tailin}": M_vwd3_in,
}

legend_str = "NR1"
dictkey = f"d300stdd50_{tailin}"
style = {"color": "r", "marker": "o", "fill": "full"}
arrow = True
fitline = True
rslt_plot(ax, M_rslt, dictkey, legend_str, fitline, arrow, lims, sizes, style, axis_type)

legend_str = "MD1"
dictkey = f"d300stdd100_{tailin}"
style = {"color": "g", "marker": "^", "fill": "full"}
arrow = False
fitline = False
rslt_plot(ax, M_rslt, dictkey, legend_str, fitline, arrow, lims, sizes, style, axis_type)

legend_str = "WD1"
dictkey = f"d300stdd200_{tailin}"
style = {"color": "b", "marker": "s", "fill": "full"}
arrow = False
fitline = False
rslt_plot(ax, M_rslt, dictkey, legend_str, fitline, arrow, lims, sizes, style, axis_type)

legend_str = "VWD1"
dictkey = f"d300stdd300_{tailin}"
style = {"color": "c", "marker": "h", "fill": "full"}
arrow = True
fitline = True
rslt_plot(ax, M_rslt, dictkey, legend_str, fitline, arrow, lims, sizes, style, axis_type)

ax.set_xlabel('$S$', fontsize=label_size)
ax.set_ylabel('$M^*$', fontsize=label_size)
ax.tick_params(axis='both', labelsize=ticks_size)
ax.set_xlim(0.004, 0.3)
ax.legend(fontsize=ticks_size)

inset_ax = inset_axes(ax, width='45%', height='45%', loc='lower right')

ks0 = M_rslt[f"d300stdd50_{tailin}"]["fit"][0]
ks1 = M_rslt[f"d300stdd100_{tailin}"]["fit"][0]
ks2 = M_rslt[f"d300stdd200_{tailin}"]["fit"][0]
ks3 = M_rslt[f"d300stdd300_{tailin}"]["fit"][0]
k0 = -M_rslt[f"d300stdd50_{tailin}"]["fit"][1]/ks0
k1 = -M_rslt[f"d300stdd100_{tailin}"]["fit"][1]/ks1
k2 = -M_rslt[f"d300stdd200_{tailin}"]["fit"][1]/ks2
k3 = -M_rslt[f"d300stdd300_{tailin}"]["fit"][1]/ks3
k10 = -M_rslt[f"d250stdd25_{tailin}"]["fit"][1]/M_rslt[f"d250stdd25_{tailin}"]["fit"][0]
k11 = -M_rslt[f"d271stdd121_{tailin}"]["fit"][1]/M_rslt[f"d271stdd121_{tailin}"]["fit"][0]
k12 = -M_rslt[f"d317stdd252_{tailin}"]["fit"][1]/M_rslt[f"d317stdd252_{tailin}"]["fit"][0]
k13 = -M_rslt[f"d347stdd537_{tailin}"]["fit"][1]/M_rslt[f"d347stdd537_{tailin}"]["fit"][0]
k21 = -M_rslt[f"d240stdd50_{tailin}"]["fit"][1]/M_rslt[f"d240stdd50_{tailin}"]["fit"][0]
k20 = -M_rslt[f"d400stdd50_{tailin}"]["fit"][1]/M_rslt[f"d400stdd50_{tailin}"]["fit"][0]
k30 = -M_rslt[f"d430stdd100_{tailin}"]["fit"][1]/M_rslt[f"d430stdd100_{tailin}"]["fit"][0]
k31 = -M_rslt[f"d167stdd100_{tailin}"]["fit"][1]/M_rslt[f"d167stdd100_{tailin}"]["fit"][0]
k40 = -M_rslt[f"d269stdd100_{tailin}"]["fit"][1]/M_rslt[f"d269stdd100_{tailin}"]["fit"][0]
k41 = -M_rslt[f"d321stdd100_{tailin}"]["fit"][1]/M_rslt[f"d321stdd100_{tailin}"]["fit"][0]
k51 = -M_rslt[f"d250stdd25_m_{tailin}"]["fit"][1]/M_rslt[f"d250stdd25_m_{tailin}"]["fit"][0]
k52 = -M_rslt[f"d271stdd121_m_{tailin}"]["fit"][1]/M_rslt[f"d271stdd121_m_{tailin}"]["fit"][0]
k53 = -M_rslt[f"d317stdd252_m_{tailin}"]["fit"][1]/M_rslt[f"d317stdd252_m_{tailin}"]["fit"][0]
k54 = -M_rslt[f"d347stdd537_m_{tailin}"]["fit"][1]/M_rslt[f"d347stdd537_m_{tailin}"]["fit"][0]
k_list = [k0, k1, k2, k3]
ks_array = np.array([ks0, ks1, ks2, ks3])
ks_array = ks_array
k_array = np.array(k_list)
k1_list = [k10, k11, k12, k13]
k1_array = np.array(k1_list)
x = np.array([0.1655, 0.3246, 0.6167, 0.87])
sd_from_tau = np.array([0.005245, 0.006080*1.2, 0.007008*1.4, 0.007238*1.6])
sd_from_tau = sd_from_tau
x1 = np.array([0.1, 0.4, 0.7, 1])
sd_from_tau1 = np.array([0.005152, 0.005030, 0.005092, 0.004603])
sd_from_tau1 = sd_from_tau1
x2 = np.array([0.1245, 0.2061])
k2_array = np.array([k20, k21])
x3 = np.array([0.3246, 0.3246])
k3_array = np.array([k30, k31])
x4 = np.array([0.3246, 0.3246])
k4_array = np.array([k40, k41])
x5 = np.array([0.1, 0.4, 0.7, 1])
k5_array = np.array([k51, k52, k53, k54])
inset_ax.plot(x, k_array, 'kx-', markersize=marker_size_in, markeredgewidth=marker_width, label='Group 1')
#inset_ax.plot(x, ks_array, 'y*-', markersize=marker_size_in, markeredgewidth=marker_width, label='Relative $k_s$')
inset_ax.plot(x1, k1_array, 'y*-', markersize=marker_size_in, markeredgewidth=marker_width, label='Group 2')
#inset_ax.plot(x2, k2_array, 'ko', markersize=marker_size_in, markeredgewidth=marker_width, label='Group 3')
#inset_ax.plot(x3, k3_array, 'k^', markersize=marker_size_in, markeredgewidth=marker_width, label='Group 4')
#inset_ax.plot(x4, k4_array, 'ks', markersize=marker_size_in, markeredgewidth=marker_width, label='Group 5')
inset_ax.plot(x5, k5_array, 'kD', markersize=marker_size_in, markeredgewidth=marker_width, label='Group 6')

inset_ax.xaxis.tick_top()
inset_ax.xaxis.set_label_position("top")
inset_ax.set_xlabel('$\\sigma_d$', fontsize=label_size)
inset_ax.set_ylabel('$S_d$', fontsize=label_size)
inset_ax.tick_params(axis='both', labelsize=ticks_size)

inset_ax.legend(fontsize=ticks_size,
		  		loc='upper right',
		  		bbox_to_anchor=(1.0, 1.6),
		  		)

fig = plt.figure(4, figsize=(8, 6), constrained_layout=True)
ax = fig.gca()

rslt_dict_smonod3_ex = np.load('rb_vs_sigma_monod2noiter_ex.npz')
rslt_dict_smonod3_ex_sd1 = np.load('rb_vs_sigma_monod2noiter_ex_smalld1.npz')
rslt_dict_3D_ex = np.load('rb_vs_sigma_3D_ex.npz')
ex_array_3D_ex_10 = rslt_dict_3D_ex['ex_10'][::4]
ez_array_3D_ex_10 = rslt_dict_3D_ex['ez_10'][::4]
sigma_array = rslt_dict_3D_ex['sigma'][::4]
th = 10/180*np.pi
mu_b_10 = (1.0 - ex_array_3D_ex_10)/(ez_array_3D_ex_10 - 1.0)/np.tan(th)
mu_b_10_sp = (1.0 - rslt_dict_smonod3_ex["ex_10"])/(rslt_dict_smonod3_ex["ez_10"] - 1.0)/np.tan(th)
mu_b_10_sp_sd1 = (1.0 - rslt_dict_smonod3_ex_sd1["ex_10"])/(rslt_dict_smonod3_ex_sd1["ez_10"] - 1.0)/np.tan(th)
k = 1/mu_b_10
k_sp = 1/mu_b_10_sp
k_sp_sd1 = 1/mu_b_10_sp_sd1

k0 = M_rslt[f"d300stdd50_{tailin}"]["fit"][0]
k1 = M_rslt[f"d300stdd100_{tailin}"]["fit"][0]
k2 = M_rslt[f"d300stdd200_{tailin}"]["fit"][0]
k3 = M_rslt[f"d300stdd300_{tailin}"]["fit"][0]
k10 = M_rslt[f"d250stdd25_{tailin}"]["fit"][0]
k11 = M_rslt[f"d271stdd121_{tailin}"]["fit"][0]
k12 = M_rslt[f"d317stdd252_{tailin}"]["fit"][0]
k13 = M_rslt[f"d347stdd537_{tailin}"]["fit"][0]
k_list = [k0, k1, k2, k3]
k1_list = [k10, k11, k12, k13]
sigma_array = np.array([0.1655, 0.3246, 0.6167, 0.87])
sigma_array1 = np.array([0.1, 0.4, 0.7, 1])
sd_array = np.array([50, 100, 200, 300])
sd_array1 = np.array([25, 121, 252, 537])
dm_array1 = np.array([250, 271, 317, 347])
cv_array = sd_array/300
cv_array1 = sd_array1/dm_array1
x, x1 = cv_array, cv_array1

#ax.plot(sigma_array, k, 'r-')
#ax.plot(rslt_dict_smonod3_ex["sigma"], k_sp, 'r--')
#ax.plot(rslt_dict_smonod3_ex_sd1["sigma"], k_sp_sd1, 'r:')
ax.plot(x, k_list, 'ko', markersize=marker_size, markerfacecolor='none', markeredgewidth=marker_width)
ax.plot(x1, k1_list, 'yo', markersize=marker_size, markerfacecolor='none', markeredgewidth=marker_width)

fig = plt.figure(5, figsize=(8, 6), constrained_layout=True)
ax = fig.gca()

axis_type = "linear"

legend_str = "NR1"
dictkey = f"d300stdd50_{tailout}"
style = {"color": "r", "marker": "o", "fill": "full"}
arrow = True
fitline = True
rslt_plot(ax, M_rslt, dictkey, legend_str, fitline, arrow, lims, sizes, style, axis_type)

legend_str = "MD1"
dictkey = f"d300stdd100_{tailout}"
style = {"color": "g", "marker": "^", "fill": "full"}
arrow = True
fitline = True
rslt_plot(ax, M_rslt, dictkey, legend_str, fitline, arrow, lims, sizes, style, axis_type)

legend_str = "WD1"
dictkey = f"d300stdd200_{tailout}"
style = {"color": "b", "marker": "s", "fill": "full"}
arrow = True
fitline = True
rslt_plot(ax, M_rslt, dictkey, legend_str, fitline, arrow, lims, sizes, style, axis_type)

legend_str = "VWD1"
dictkey = f"d300stdd300_{tailout}"
style = {"color": "c", "marker": "h", "fill": "full"}
arrow = True
fitline = True
rslt_plot(ax, M_rslt, dictkey, legend_str, fitline, arrow, lims, sizes, style, axis_type)

ax.set_xlim([0, 0.175])
ax.set_ylim([0, 0.25])

plt.show()