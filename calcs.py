#-*- coding: utf-8 -*-
import numpy as np


G_F = 1.1663*10**(-5) # GeV^-2

M_Z = 91.187 #GeV
weinberg = 28.741


def partial_width(color_factor, Q_f, I_3):
    g_A = I_3
    g_V = I_3 - 2*Q_f*(np.sin(weinberg)**2)

    return (color_factor*np.sqrt(2)/float(12*np.pi))*G_F*(M_Z**3)*(g_V**2 + g_A**2) # GeV


lepton = partial_width(1, -1, -0.5)
print("Leptons: %f GeV" % lepton)

neutrino = partial_width(1, 0, 0.5)
print("Neutrinos: %f GeV" % neutrino)

quark_u_c = partial_width(3, 2.0/3.0, 0.5)
print("Up and charm quark: %f GeV" % quark_u_c)

quark_d_s_b = partial_width(3, -1.0/3.0, -0.5)
print("Down, strange and bottom quark: %f GeV" % quark_d_s_b)


Z = 3*lepton + 2*quark_u_c + 3*quark_d_s_b +3*neutrino
print("Z_0: %f GeV" % Z)

print("Hadronisch: %f GeV, Rate: %.1f Prozent" % (2*quark_u_c + 3*quark_d_s_b, 100*(2*quark_u_c + 3*quark_d_s_b)/Z))
print("Leptonisch: %f GeV, Rate: %.1f Prozent" % (3*lepton, 100*(3*lepton)/Z))
print("Leptonisch (invisible (neutrinos)): %f GeV, Rate: %.1f Prozent" % (3*neutrino, 100*(3*neutrino)/Z))


# peak: s=M_Z^2
def peak_cross_section(partial_width):
    return (12*np.pi/(M_Z**2))*(lepton/Z)*(partial_width/Z)


print("Peak cross section for different fermion channels")
print("Lepton (For Bhaba-scatter e+e->e+e- additional Feynman diagram): %f MeV^-2" % (peak_cross_section(lepton)*1000))
print("Neutrino: %f MeV^-2" % (peak_cross_section(neutrino)*1000))
print("U and C: %f MeV^-2" % (peak_cross_section(quark_u_c)*1000))
print("D, S and B: %f MeV^-2" % (peak_cross_section(quark_d_s_b)*1000))


print("Change of partial width for Z in case of an additional channel")
print("Additional lepton %f GeV thats %.1f percent change" % (Z + lepton, lepton/Z*100))
print("Additional neutrino %f GeV thats %.1f percent change" % (Z + neutrino, neutrino/Z*100))
print("Additional u quark %f GeV thats %.1f percent change" % (Z + quark_u_c, quark_u_c/Z*100))
print("Additional d quark %f GeV thats %.1f percent change" % (Z + quark_d_s_b, quark_d_s_b/Z*100))




# forward-back asymmetry

def forward_backward(w_sin_2_value, s):
    weinberg = np.arcsin(np.sqrt(w_sin_2_value))*(360/(2*np.pi))

    def x(s):
        return s/((s-M_Z**2)+1j*s*Z/M_Z)


    g_e = -0.5 - 2*(-1)*(np.sin(weinberg)**2)
    v_e = g_e/(2*np.sin(weinberg)*np.cos(weinberg))
    a_e = (-0.5)/(2*np.sin(weinberg)*np.cos(weinberg))
    def F_1(s, Q_f, I_3):
        g_V = I_3 - 2*Q_f*(np.sin(weinberg)**2)
        v_f = g_V/(2*np.sin(weinberg)*np.cos(weinberg))
        a_f = I_3/(2*np.sin(weinberg)*np.cos(weinberg))
        xii = x(s)

        return Q_f**2 - 2*v_e*v_f*Q_f*np.real(xii) + (v_e**2 + a_e**2)*(v_f**2 + a_f**2)*(np.absolute(xii)**2)

    def F_2(s, Q_f, I_3):
        g_V = I_3 - 2*Q_f*(np.sin(weinberg)**2)
        v_f = g_V/(2*np.sin(weinberg)*np.cos(weinberg))
        a_f = I_3/(2*np.sin(weinberg)*np.cos(weinberg))
        xii = x(s)

        return -2*a_e*a_f*Q_f*np.real(xii) + 4*v_e*a_e*v_f*a_f*(np.absolute(xii)**2)


    print("Vorwärts rückwärs Anti %f (at s:%f, w:%f)" % (3.0/4.0*(F_2(s, -1, -0.5)/F_1(s, -1, -0.5)), np.sqrt(s), w_sin_2_value))

    return 3.0/4.0*(F_2(s, -1, -0.5)/F_1(s, -1, -0.5))



w_sin_2_value = 0.23122
s = 91.225**2
#s = 89.225**2
#s = 93.225**2

s_range = np.arange(88, 94, 0.1)
w_range = np.arange(0.20, 0.26, 0.001)


result = np.zeros((len(w_range), len(s_range)), dtype=np.float32)
for i, sqrt_s in enumerate(s_range):
    for l, w_sin_2_value in enumerate(w_range):

        result[l][i] = forward_backward(w_sin_2_value, (sqrt_s+0.225)**2)


import os,sys
import numpy as np
BASE = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(BASE, 'opengl_plot_prototype'))

from gllib.plot.app import plot2d 
from gllib.plot.domain import FieldDomain 
from gllib.plot.field import Field
from gllib.plot.color.schemes import ColorMap

COLOR_SCHEME = ColorMap('IDL_Blue-Pastel-Red', colorrange=[-.2,0.9])



def plot(plotter):
    plotter.graphs['foo'] = Field(FieldDomain.from_numpy(result),
        top_left=(s_range[0]+0.225, w_range[-1]),
        bottom_right=(s_range[-1]+0.225, w_range[0]),
        #data_kernel="fragment_color = texture(tex[0], vec2(x.y, x.x));",
        color_scheme=COLOR_SCHEME)

plot2d(plot, axis=[s_range[-1]-s_range[0], w_range[-1]-w_range[0]], origin=[-s_range[0]-0.225, w_range[0]],
    xlabel='sqrt(s)',
    colorlegend=COLOR_SCHEME,
    ylabel='sin^2(phi)',
    title='Forward Backward Asymmetry',
    height=300,
    width=500)


