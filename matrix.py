import numpy as np
from numpy.linalg import inv
from util import np_array_to_latex

matrix = np.array([
    [0.28591, 0.00001, 0.00000, 0.00000],
    [0.00000, 0.91794, 0.00577, 0.00000],
    [0.00303, 0.00913, 0.79342, 0.00405],
    [0.00102, 0.00000, 0.01198, 0.99471],


])

err_matrix = np.array([
    [0.0029,  0.0000,  0.0000,  0.0000],
    [0.0000,  0.0009,  0.0003,  0.0000],
    [0.0004,  0.0003,  0.0014,  0.0002],
    [0.0002,  0.0000,  0.0004,  0.0002],

])

print("Matrix")
print(matrix)
print(np_array_to_latex(matrix,2))
inv_matrix = inv(matrix)
print("inverse matrix")
print(inv_matrix)
print(np_array_to_latex(inv_matrix,2))


print("Error Matrix")
print(err_matrix)
inv_err_mat = - np.dot(inv_matrix, np.dot(err_matrix, inv_matrix))
print(inv_err_mat)
print(np_array_to_latex(inv_err_mat,2))


"""
crazy shit coming up here
"""
e = np.zeros((7, 4))
e[0] = [20,93,97,2500]
e[1] = [81,282,266,6650]
e[2] = [91,389,345,8770]
e[3] = [631,2879,2398,66863]
e[4] = [114,577,484,12995]
e[5] = [37,263,251,6257]
e[6] = [56,292,255,6982]

#electron from dataset2
#e[0] = [58,93,97,2500]
#e[1] = [120,282,266,6650]
#e[2] = [156,389,345,8770]
#e[3] = [1177,2879,2398,66863]
#e[4] = [214,577,484,12995]
#e[5] = [84,263,251,6257]
#e[6] = [84,292,255,6982]

#corrected
#e[0] = [70,93,97,2500]
#e[1] = [192,282,266,6650]
#e[2] = [444,389,345,8770]
#e[3] = [755,2879,2398,66863]
#e[4] = [789,577,484,12995]
#e[5] = [170,263,251,6257]
#e[6] = [93,292,255,6982]

print('data from experiment')
print(np_array_to_latex(e,0))

print("True events")
ntrue = inv_matrix.dot(e.T).T
print(inv_matrix.dot(e.T).T)
print(np_array_to_latex(inv_matrix.dot(e.T).T,3))
print("error")
dntrue = inv_err_mat.dot(e.T).T
print(dntrue)
print(np_array_to_latex(inv_err_mat.dot(e.T).T,3))

fixes = np.array([
    [0.09]*3+[2.0],
    [0.2]*3+[4.3],
    [0.36]*3+[7.7],
    [0.52]*3+[10.8],
    [0.22]*3+[4.7],
    [-0.01]*3+[-0.2],
    [-0.08]*3+[-1.6],
], dtype=np.float32)

print('fixes_latex')
print(np_array_to_latex(fixes,2))
print('L')
lumi_errs = [[4.249604],[5.691792],[4.454466],[16.43293],[4.848926],[4.276552],[6.104764]]
lumi = [[463.979], [667.5236], [486.7642], [2246.568], [535.908], [450.6], [709.698]]
print(np_array_to_latex(lumi,4))
print('ERROR')
print(np_array_to_latex(lumi_errs))
lumis = np.array([a*4 for a in lumi], dtype=np.float32)


monte_carlo_corrections = np.array([
    [1.1145911123,    1.060006996, 1.2627697592,    1.0152181196] for x in range(7)
])
#sigma = e/lumis+fixes

#monto carlo correction
sigma = (monte_carlo_corrections*ntrue)/lumis + fixes
print('sigma')
print(sigma)

print(np_array_to_latex(sigma,5))

dsigma = np.empty_like(sigma)
for i in range(0,7):
    for j in range(0,4):
        dsigma[i][j] = np.sqrt((1.0/lumi[i][0]*dntrue[i][j])**2+(-ntrue[i][j]/(lumi[i][0])**2*lumi_errs[i][0])**2)

print('dsigma')
print(dsigma)
print(np_array_to_latex(dsigma,4))



sigma_lept = np.zeros((7,2))
for i in range(7):
    sigma_lept[i][0] = sigma[i][0] + sigma[i][1] + sigma[i][2]
    print("y[%d] = %.5f;" % (i, sigma[i][0]))
    sigma_lept[i][1] = sigma[i][3]

print("sigma leptons sum")
print(np_array_to_latex(sigma_lept,4))


dsigma_lept = np.zeros((7,2))
for i in range(7):
    dsigma_lept[i][0] = np.sqrt(dsigma[i][0]**2 + dsigma[i][1]**2 + dsigma[i][2]**2)
    print("error[%d] = %.5f;" % (i, dsigma[i][0]))
    dsigma_lept[i][1] = dsigma[i][3]

print("dsigma leptons sum")
print(np_array_to_latex(dsigma_lept,4))



"""
determine partial width
"""
def p_width(peak, gz, ge, mz):
    return (peak*gz**2*mz**2)/(12*np.pi*ge)

def p_width_e(peak, gz, mz):
    return np.sqrt((peak*gz**2*mz**2)/(12*np.pi))



def gz_error_e(error, peak, mz):
    return np.sqrt(peak/(12*np.pi))*mz*error

def mz_error_e(error, peak, gz):
    return np.sqrt(peak/(12*np.pi))*gz*error

def peak_error_e(error, mz, gz, peak):
    return (gz*mz*np.sqrt(peak/(12*np.pi)))/(2*peak)*error



def gz_error(error, peak, gz, mz, ge):
    return (peak*gz*mz**2)/(6*np.pi*ge)*error

def mz_error(error, peak, gz, mz, ge):
    return (peak*gz**2*mz)/(6*np.pi*ge)*error

def peak_error(error, peak, gz, mz, ge):
    return (gz**2*mz**2)/(12*np.pi*ge)*error


Z_mass = [90.97, 91.18, 91.21, 91.13]
Z_mass_error = [0.0048, 0.004858, 0.003, 0.00669]

print("Z mass %f" % np.mean(Z_mass))
print("Z mass error %f" % (np.sqrt(np.sum(np.square(Z_mass_error)))))


mb_factor = 2.57
e_peak = 0.000001735
e_peak_error = 0.00000001165
e_gz = 2.098
e_gz_error = 0.012
e_mz = 90.97
e_mz_error = 0.00481

gamma_e_m = p_width_e(mb_factor*e_peak, e_gz, e_mz)
print("electron width: %f" % gamma_e_m)
e_width_error = np.sqrt(gz_error_e(e_gz_error, mb_factor*e_peak, e_mz)**2 + mz_error_e(e_mz_error, mb_factor*e_peak, e_gz)**2 + peak_error_e(mb_factor*e_peak_error, e_mz, e_gz, mb_factor*e_peak)**2)
print("electron width error: %f" % e_width_error)


##continue with litaerture value
gamma_e = 0.0838
#gamma_e = gamma_e_m

m_peak = 0.000001836 #(+/- 0.00903)
m_peak_error = 0.00000000903
m_gz = 2.398
m_gz_error = 0.01142
m_mz = 91.18
m_mz_error = 0.004858

muon_width = p_width(m_peak*mb_factor, m_gz, gamma_e, m_mz)
print("muon: %f" % muon_width)
m_width_error = np.sqrt(gz_error(m_gz_error, mb_factor*m_peak, m_gz, m_mz, gamma_e)**2 + mz_error(m_mz_error, mb_factor*m_peak, m_gz, m_mz, gamma_e)**2 + peak_error(m_peak_error, mb_factor*m_peak, m_gz, m_mz, gamma_e)**2)
print("muon width error: %f" % m_width_error)

t_peak = 0.000001612
t_peak_error = 0.00000001087
t_gz = 2.516
t_gz_error = 0.01656
t_mz = 91.13
t_mz_error = 0.006693

tau_width = p_width(m_peak*mb_factor, t_gz, gamma_e, t_mz)
print("tau: %f" % tau_width)
t_width_error = np.sqrt(gz_error(t_gz_error, mb_factor*t_peak, t_gz, t_mz, gamma_e)**2 + mz_error(t_mz_error, mb_factor*t_peak, t_gz, t_mz, gamma_e)**2 + peak_error(t_peak_error, mb_factor*t_peak, t_gz, t_mz, gamma_e)**2)
print("tau width error: %f" % t_width_error)


l_peak = 0.000004702
leptons_width = p_width(l_peak*mb_factor, 2.372, gamma_e, 91.11)
print("leptons: %f" % leptons_width)

leptons_width = gamma_e_m + tau_width + muon_width
leptons_width_error = np.sqrt(e_width_error**2 + t_width_error**2 + m_width_error**2)

q_peak = 0.00004022
q_peak_error = 0.0000001821
q_gz = 2.576
q_gz_error = 0.008535
q_mz = 91.21
q_mz_error = 0.003095

hadrons_width = p_width(q_peak*mb_factor, q_gz, gamma_e, q_mz)
print("quark: %f" % hadrons_width)
q_width_error = np.sqrt(gz_error(q_gz_error, mb_factor*q_peak, q_gz, q_mz, gamma_e)**2 + mz_error(q_mz_error, mb_factor*q_peak, q_gz, q_mz, gamma_e)**2 + peak_error(q_peak_error, mb_factor*q_peak, q_gz, q_mz, gamma_e)**2)
print("quark width error: %f" % q_width_error)



#forward backward
muon_monte = [44700,49639]

true_monte_muon_left = matrix[1][0]*muon_monte[0] + matrix[1][1]*muon_monte[0] + matrix[1][2]*muon_monte[0]
true_monte_muon_right = matrix[1][0]*muon_monte[1] + matrix[1][1]*muon_monte[1] + matrix[1][2]*muon_monte[1]
    
print("True muon monte events: backward: %.4f, forward: %.4f" % (true_monte_muon_left, true_monte_muon_right))
print("Monte carlo correction: backward: %.4f, forward: %.4f" % (true_monte_muon_left*monte_carlo_corrections[0][1], true_monte_muon_right*monte_carlo_corrections[0][1]))


opal = [
    [631, 0],
    [1363, 1516],
    [517,1881]
]

true_opal_muon_left = matrix[1][0]*opal[0][0] + matrix[1][1]*opal[1][0] + matrix[1][2]*opal[2][0]
true_opal_muon_right = matrix[1][0]*opal[0][1] + matrix[1][1]*opal[1][1] + matrix[1][2]*opal[2][1]

print("True muon opal events: backward: %.4f, forward: %.4f" % (true_opal_muon_left, true_opal_muon_right))
print("Monte carlo correction: backward: %.4f, forward: %.4f" % (true_opal_muon_left*monte_carlo_corrections[0][1], true_opal_muon_right*monte_carlo_corrections[0][1]))


# forward backward via A_fb = (n_right - n_left) / (n_right + n_left)
a_fb_monte = (true_monte_muon_right -true_monte_muon_left)/ (true_monte_muon_right + true_monte_muon_left)
print("muon monte fb asy: %f" % a_fb_monte)
a_fb_opal = (true_opal_muon_right -true_opal_muon_left)/ (true_opal_muon_right + true_opal_muon_left)
print("muon opal fb asy: %f" % a_fb_opal)

def sin_2_weinberg(a_fb):
    return (1 - np.sqrt(a_fb)/3.0)/4.0


print(sin_2_weinberg(a_fb_monte))
print(sin_2_weinberg(a_fb_opal))


# cross section ratios
total_cross_section_measured = e_peak + m_peak + t_peak + q_peak

print("Leptonic cross section %.4f" % ((e_peak + m_peak + t_peak) / total_cross_section_measured*100))
print("Hadronic cross section %.4f" % (q_peak / total_cross_section_measured*100))


# if we compare the width of leptons and hadrons we enc up with
z_width_measured = (hadrons_width + leptons_width)
print("Total width of hardrons and leptons %.4f" % z_width_measured)
z_width_measured_error = np.sqrt(q_width_error**2 + leptons_width_error**2)
print("With an error of %.4f" % z_width_measured_error)

g_z_fit = [m_gz, t_gz, q_gz] # excluded electron

z_width_theoretical = 2.5027 # calculated by us
z_width_theoretical = (3*83.8+299*2+378*3+3*167.6) / 1000.0 # from the mappe
z_width_theoretical = np.mean(g_z_fit) # from our measuerment
z_width_theoretical_error = np.sqrt(t_gz_error**2 + m_gz_error**2 + q_gz_error**2)

neutrino_width_deficiet = z_width_theoretical - z_width_measured
print("Missing %.4f +/- %f to get to theoretical value of %f +/- %f" % (neutrino_width_deficiet, np.sqrt(z_width_theoretical_error**2 + z_width_measured_error**2), z_width_theoretical, z_width_theoretical_error))
print("Split by 3 we end up with %f" % (neutrino_width_deficiet/3.0))
neutrino_width_theo = 0.1659
print("Or we split by %f and get %f families" % (neutrino_width_theo, neutrino_width_deficiet/neutrino_width_theo))



gz_std = np.std(g_z_fit)
print("Errors: %f" % (np.sqrt(np.sqrt(z_width_theoretical_error**2 + z_width_measured_error**2)**2 + gz_std**2)/neutrino_width_theo))

lept_univ = [65.6, 71.405, 78.519]

print("lept univ %f +/- %f" % (np.mean(lept_univ), np.std(lept_univ)))



