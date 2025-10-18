import numpy as np

def get_exp_data(case):
    """
    获取相应case的实验数据
    """
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
    	"Ho2012_c": Ho2012_c_data,
        "Zh19": Zh19_data,
        "Ma18": Ma18_data,
        "Ca22": Ca22_data
    }
    if case in others_case_dict:
        return others_case_dict[case]
    else:
        raise ValueError(f"未知的实验数据案例: {case}")