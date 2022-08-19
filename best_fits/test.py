import ast
from asteroids import Asteroid
import matplotlib.pyplot as plt
import numpy as np


ast_1 = Asteroid("triaxial_05693")
ast_2 = Asteroid("triaxial_68950")


def comparison_chi_template(packed_name1, packed_name2, values_x_s, chis_s, values_x_t, chis_t,
                                    label_x, unit_x=None):
    
    values_x_s = [np.log10(value) for value in values_x_s]
    values_x_t = [np.log10(value) for value in values_x_t]
    
    absolute_x_minima = min(values_x_s)
    if min(values_x_t) < min(values_x_s):
        absolute_x_minima = min(values_x_t)
    
    absolute_x_maxima = max(values_x_s)
    if max(values_x_t) > max(values_x_s):
        absolute_x_maxima = max(values_x_t)
    
    absolute_y_minima = min(chis_s)
    if min(chis_t) < min(chis_s):
        absolute_y_minima = min(chis_t)
    
    absolute_y_maxima = max(chis_s)
    if max(chis_t) > max(chis_s):
        absolute_y_maxima = max(chis_t)
    
    fig, (ax1, ax2) = plt.subplots(1, 2)
    scatter_plot = ax1.scatter(values_x_s, chis_s, s=1, color="black", marker='o',linewidths=1)

    log_values_copy = values_x_s[:]
    values_count = len(values_x_s)
    log_values_copy.sort()
    sigma_1_low = log_values_copy[int(values_count * 0.16)]
    sigma_1_high = log_values_copy[int(values_count * 0.84)]
    median = log_values_copy[int(values_count * 0.5)]
    sigma_2_low = log_values_copy[int(values_count * 0.025)]
    sigma_2_high = log_values_copy[int(values_count * 0.975)]
    ax1.axvline(median, color='#cc0000')
    ax1.axvline(sigma_1_low, color='#cc0000', ls='dashed')
    ax1.axvline(sigma_1_high, color='#cc0000', ls='dashed')
    ax1.axvline(sigma_2_low, color='#cc0000', ls='dotted')
    ax1.axvline(sigma_2_high, color='#cc0000', ls='dotted')

    ax1.set_xlim(absolute_x_minima, absolute_x_maxima)
    ax1.set_ylim(absolute_y_minima, absolute_y_maxima)
    label_string_x = f"Log {label_x}"
    if unit_x != None:
        label_string_x = f"Log {label_x} ({unit_x})"
    #ax1.set_xlabel(label_string_x)
    label_string_y = r"fit $\chi^2$"
    #ax1.set_ylabel(label_string_y)
    ax1.set_title(f"{packed_name1} (Triaxial Ellipsoid)", loc="left")

    scatter_plot = ax2.scatter(values_x_t, chis_t, s=1, color="black", marker='o',linewidths=1)

    log_values_copy = values_x_t[:]
    values_count = len(values_x_t)
    log_values_copy.sort()
    sigma_1_low = log_values_copy[int(values_count * 0.16)]
    sigma_1_high = log_values_copy[int(values_count * 0.84)]
    median = log_values_copy[int(values_count * 0.5)]
    sigma_2_low = log_values_copy[int(values_count * 0.025)]
    sigma_2_high = log_values_copy[int(values_count * 0.975)]
    ax2.axvline(median, color='#cc0000')
    ax2.axvline(sigma_1_low, color='#cc0000', ls='dashed')
    ax2.axvline(sigma_1_high, color='#cc0000', ls='dashed')
    ax2.axvline(sigma_2_low, color='#cc0000', ls='dotted')
    ax2.axvline(sigma_2_high, color='#cc0000', ls='dotted')

    ax2.set_yticklabels([])
    ax2.set_xlim(absolute_x_minima, absolute_x_maxima)
    ax2.set_ylim(absolute_y_minima, absolute_y_maxima)
    label_string_x = f"Log {label_x}"
    if unit_x != None:
        label_string_x = f"Log {label_x} ({unit_x})"
    label_string_y = r"fit $\chi^2$"
    ax2.set_title(f"{packed_name2} (Triaxial Ellipsoid)", loc="left")

    fig.text(0.5, 0.04, label_string_x, ha='center', va='center')
    fig.text(0.04, 0.5, label_string_y, ha='center', va='center', rotation='vertical')
    
    title_file_name = f"{label_x.replace(' ', '_').lower()}_vs_chis"
    fig.savefig(f"./comparison_plots/{packed_name1}vs{packed_name2}.png", dpi=1000)
    fig.savefig(f"./comparison_plots/{packed_name1}vs{packed_name2}.pdf", dpi=1000)
    plt.close(fig)

comparison_chi_template(ast_1.packed_name, ast_2.packed_name, ast_1.diameters, ast_1.chis, ast_2.diameters, ast_2.chis, "Diameter", "km")