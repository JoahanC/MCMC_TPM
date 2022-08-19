"""
This file takes the inputs for Ned's MCMC model and converts them into
an astropy and LaTeX table format for use in publications.
Must be run in the same directory with the run*.csh file 
"""
import os
from astropy.table import Table
from asteroids import Asteroid
from helpers import *


def generate_objects():
    """
    Generates all Asteroid objects for legally set up MCMC folders.
    """
    # Generate list of all valid asteroid folders
    files = os.listdir("../")
    dirs = []
    for file in files:
        if os.path.isdir(f"../{file}"):
            dirs.append(file)
    dirs.remove(".git")
    dirs.remove("best_fits")
    dirs.remove("asteroid_ex")

    # Generate clusters to have asteroids processed in order
    clusters = set()
    for dir in dirs:
        name = dir.replace("triaxial_", '')
        clusters.add(name)
    objects = []
    print("Generating Asteroid objects.")
    for name in clusters:
        object = Asteroid(name)
        objects.append(object)
        object = Asteroid("triaxial_" + name)
        objects.append(object)
    print("Finished generating Asteroid objects.\n")
    return objects


def generate_total_table(objects):
    """
    Generates an astropy table of all thermophysical property outputs from a set of
    Asteroid objects.

    Parameters
    ----------
    objects : list
        A list of Asteroid objects which contain information about a given asteroid's
        modeling outputs.
    """
    # Initializing all collector lists
    object_names = []
    is_triaxial = []
    diameter_mean = []
    diameter_mean_sig = []
    diameter_median = []
    diameter_median_pos_sig = []
    diameter_median_neg_sig = []
    p_V_mean = []
    p_V_mean_sig = []
    p_V_median = []
    p_V_median_sig = []
    theta_mean = []
    theta_mean_sig = []
    theta_median = []
    theta_median_pos_sig = []
    theta_median_neg_sig = []
    period = []
    period_pos_sig = []
    period_neg_sig = []
    gamma = []
    gamma_pos_sig = []
    gamma_neg_sig = []
    crater_frac = []
    crater_frac_pos_sig = []
    crater_frac_neg_sig = []
    alb_ratio = []
    alb_ratio_pos_sig = []
    alb_ratio_neg_sig = []
    pole_peak_ra = []
    pole_peak_dec = []
    mean_pole_bar = []
    mean_pole_ra = []
    mean_pole_dec = []
    eigenvalue_1 = []
    eigenvalue_1_ra = []
    eigenvalue_1_dec = []
    eigenvalue_2 = []
    eigenvalue_2_ra = []
    eigenvalue_2_dec = []
    eigenvalue_3 = []
    eigenvalue_3_ra = []
    eigenvalue_3_dec = []
    h_magnitude = []
    h_magnitude_error = []
    epoch_count = []
    model_type = []
    
    # Aggregating all parameter data from each object
    for object in objects:
        object_names.append(object.packed_name)
        is_triaxial.append(object.is_triaxial)
        diameter_mean.append(object.diameter_mean)
        diameter_mean_sig.append(object.diameter_mean_16th_percentile)
        diameter_median.append(object.diameter_median)
        diameter_median_pos_sig.append(object.diameter_median_16th_percentile)
        diameter_median_neg_sig.append(object.diameter_median_84th_percentile)
        p_V_mean.append(object.albedo_mean)
        p_V_mean_sig.append(object.albedo_mean_16th_percentile)
        p_V_median.append(object.albedo_median)
        p_V_median_sig.append(object.albedo_median_16th_percentile)
        theta_mean.append(object.theta_mean)
        theta_mean_sig.append(object.theta_mean_16th_percentile)
        theta_median.append(object.theta_median)
        theta_median_pos_sig.append(object.theta_median_16th_percentile)
        theta_median_neg_sig.append(object.theta_median_84th_percentile)
        period.append(object.period_median)
        period_pos_sig.append(object.period_16th_percentile)
        period_neg_sig.append(object.period_84th_percentile)
        gamma.append(object.gamma_median)
        gamma_pos_sig.append(object.gamma_16th_percentile)
        gamma_neg_sig.append(object.gamma_84th_percentile)
        crater_frac.append(object.crater_fraction_median)
        crater_frac_pos_sig.append(object.crater_fraction_16th_percentile)
        crater_frac_neg_sig.append(object.crater_fraction_84th_percentile)
        alb_ratio.append(object.ir_fraction_median)
        alb_ratio_pos_sig.append(object.ir_fraction_16th_percentile)
        alb_ratio_neg_sig.append(object.ir_fraction_84th_percentile)
        pole_peak_ra.append(object.pole_peak_ra)
        pole_peak_dec.append(object.pole_peak_dec)
        mean_pole_ra.append(object.pole_mean_ra)
        mean_pole_dec.append(object.pole_mean_dec)
        mean_pole_bar.append(object.pole_bar)
        eigenvalue_1.append(object.eigenvalues["eigenvalue_1"][0])
        eigenvalue_1_ra.append(object.eigenvalues["eigenvalue_1"][1])
        eigenvalue_1_dec.append(object.eigenvalues["eigenvalue_1"][2])
        eigenvalue_2.append(object.eigenvalues["eigenvalue_2"][0])
        eigenvalue_2_ra.append(object.eigenvalues["eigenvalue_2"][1])
        eigenvalue_2_dec.append(object.eigenvalues["eigenvalue_2"][2])
        eigenvalue_3.append(object.eigenvalues["eigenvalue_3"][0])
        eigenvalue_3_ra.append(object.eigenvalues["eigenvalue_3"][1])
        eigenvalue_3_dec.append(object.eigenvalues["eigenvalue_3"][2])
        h_magnitude.append(object.h_magnitude)
        h_magnitude_error.append(object.h_error)
        epoch_count.append(object.epoch_count)
        if object.is_triaxial:
            model_type.append("Triaxial")
        else:
            model_type.append("Spherical")

    # Writing to astropy table
    object_table = Table()
    object_table['name'] = object_names
    object_table['diameter_mean'] = diameter_mean
    object_table['diameter_mean_sig'] = diameter_mean_sig
    object_table['diameter_median'] = diameter_median
    object_table['diameter_median_pos_sig'] = diameter_median_pos_sig
    object_table['diameter_median_neg_sig'] = diameter_median_neg_sig
    object_table['p_V_mean'] = p_V_mean
    object_table['p_V_mean_sig'] = p_V_mean_sig
    object_table['p_V_median'] = p_V_median
    object_table['p_V_median_sig'] = p_V_median_sig
    object_table['theta_mean'] = theta_mean
    object_table['theta_mean_sig'] = theta_mean_sig
    object_table['theta_median'] = theta_median
    object_table['theta_median_pos_sig'] = theta_median_pos_sig
    object_table['theta_median_neg_sig'] = theta_median_neg_sig
    object_table['period'] = period
    object_table['period_pos_sig'] = period_pos_sig
    object_table['period_neg_sig'] = period_neg_sig
    object_table['kappa_rho_c'] = gamma
    object_table['kappa_rho_c_pos_sig'] = gamma_pos_sig
    object_table['kappa_rho_c_neg_sig'] = gamma_neg_sig
    object_table['crater_fraction'] = crater_frac
    object_table['crater_fraction_pos_sig'] = crater_frac_pos_sig
    object_table['crater_fraction_neg_sig'] = crater_frac_neg_sig
    object_table['p_IRp_V'] = alb_ratio
    object_table['p_IRp_V_pos_sig'] = alb_ratio_pos_sig
    object_table['p_IRp_V_neg_sig'] = alb_ratio_neg_sig
    object_table['pole_peak_ra'] = pole_peak_ra
    object_table['pole_peak_dec'] = pole_peak_dec
    object_table['mean_pole_ra'] = mean_pole_ra
    object_table['mean_pole_dec'] = mean_pole_dec
    object_table['mean_pole_bar'] = mean_pole_bar
    object_table['eigenvector_1'] = eigenvalue_1
    object_table['eigenvector_1_ra'] = eigenvalue_1_ra
    object_table['eigenvector_1_dec'] = eigenvalue_1_dec
    object_table['eigenvector_2'] = eigenvalue_2
    object_table['eigenvector_2_ra'] = eigenvalue_2_ra
    object_table['eigenvector_2_dec'] = eigenvalue_2_dec
    object_table['eigenvector_3'] = eigenvalue_3
    object_table['eigenvector_3_ra'] = eigenvalue_3_ra
    object_table['eigenvector_3_dec'] = eigenvalue_3_dec
    object_table['h_magnitude'] = h_magnitude
    object_table['h_error'] = h_magnitude_error
    object_table['epoch_count'] = epoch_count
    object_table['model_type'] = model_type
    object_table.write(f"MCMC_outputs.tbl", format="ipac", 
                       overwrite=True)


def generate_MCMC_results():
    """
    Generates a short LaTeX file which contains table information about all of the
    final MCMC modeling parameters, as well as input information for each asteroid
    run.
    """
    # Begin document header
    data_object = Table.read("MCMC_outputs.tbl", format='ipac')
    tex_file = open("thermophysical_outputs.tex", 'w')
    tex_file.write("\makeatletter\n")
    tex_file.write("\declare@file@substitution{revtex4-1.cls}{revtex4-2.cls}\n")
    tex_file.write("\makeatother\n")
    tex_file.write("\\documentclass[linenumbers]{aastex631}\n\n")
    tex_file.write("\\usepackage{float}\n\n")
    tex_file.write("\\shorttitle{MCMC Thermophysical Modeling Results}\n")
    tex_file.write("\\begin{document}\n\n")
    tex_file.write("\\title{MCMC Thermophysical Modeling Results \\footnote{Summer 2022}}\n")
    tex_file.write("\\section{Introduction}\n")
    tex_file.write("Thermophysical modeling results for nine asteroids are documented" +
                   "in this file.\n")
    tex_file.write("Included is a series of tables with major" +
    " physical characteristics and a series of plots for each\n")
    tex_file.write("object modeled.\n")

    # Generate physical characteristics output table
    tex_file.write("\\begin{deluxetable*}{" + 'c'*7 + "}[h]\n")
    tex_file.write(' ' * 4 + "\\tablenum{1}\n")
    tex_file.write(' ' * 4 + "\\tablecaption{Physical characteristics from thermophysical modeling of various objects}\n")
    tex_file.write(' ' * 4 + "\\tablewidth{10pt}\n")
    tex_file.write(' ' * 4 + "\\tablehead{\n")
    col_names = ["Name", "Diameter", "Albedo", "Theta", "Period", "Crater Fraction", "Model Type"]
    table_header = ' ' * 8
    for name in col_names:
        table_header += "\\colhead{" + name + "} & "
    table_header = table_header[:len(table_header) - 2] + "\\\ \n" + ' ' * 8
    col_types = ['', 'km', '', 'deg', 'hr', '', '']
    for type in col_types:
        table_header += "\\colhead{" + type + "} & "
    tex_file.write(table_header[:len(table_header) - 2] + "\n")
    tex_file.write(' ' * 4 + '}\n')
    tex_file.write(' ' * 4 + "\startdata\n")

    for iter in range(len(data_object['name'])):
        line = ' ' * 8
        for idx, col in enumerate(data_object.colnames):
            if idx == 0:
                line += f"{data_object[col][iter]} & "
            elif idx in [3, 8, 12, 15, 21]:
                line += f"${data_object[col][iter]}^"
            elif idx in [4, 9, 13, 16]:
                line += "{+" + f"{data_object[col][iter]}" + "\%}_"
            elif idx in [22]:
                line += "{+" + f"{data_object[col][iter]}" + "}_"
            if idx in [5, 9, 14, 17]:
                line += "{-" + f"{data_object[col][iter]}" + "\%}$ & "
            elif idx == 23:
                line += "{-" + f"{data_object[col][iter]}" + "}$ & "
            elif idx == 44:
                line += f"{data_object[col][iter]} " + "\\\ \n"
        tex_file.write(line)
    tex_file.write(' ' * 4 + "\enddata\n")
    tex_file.write("\end{deluxetable*}\n\n")
    tex_file.write("\pagebreak\n")

    # Generate input parameter table
    tex_file.write("\\begin{deluxetable*}{" + 'c'*5 + "}[h]\n")
    tex_file.write(' ' * 4 + "\\tablenum{2}\n")
    tex_file.write(' ' * 4 + "\\tablecaption{Input paramaters for modeled objects}\n")
    tex_file.write(' ' * 4 + "\\tablewidth{10pt}\n")
    tex_file.write(' ' * 4 + "\\tablehead{\n")
    col_names = ["Name", "H Magnitude", "H Magnitude Uncertainty", "Number of Epochs", "Model Type"]
    table_header = ' ' * 8
    for name in col_names:
        table_header += "\\colhead{" + name + "} & "
    table_header = table_header[:len(table_header) - 2] + "\\\ \n" + ' ' * 8
    col_types = ['', 'mag', 'mag', '', '']
    for type in col_types:
        table_header += "\\colhead{" + type + "} & "
    tex_file.write(table_header[:len(table_header) - 2] + "\n")
    tex_file.write(' ' * 4 + '}\n')
    tex_file.write(' ' * 4 + "\startdata\n")
    
    for iter in range(len(data_object['name'])):
        line = ' ' * 8
        for idx, col in enumerate(data_object.colnames):
            if idx == 0:
                line += f"{data_object[col][iter]} & "
            if idx in [41, 42]:
                line += f"${data_object[col][iter]}$ & "
            if idx == 43:
                line += f"${data_object[col][iter]}$ & "
            elif idx == 44:
                line += f"{data_object[col][iter]} " + "\\\ \n"
        tex_file.write(line)
    tex_file.write(' ' * 4 + "\enddata\n")
    tex_file.write("\end{deluxetable*}\n\n")
    tex_file.write("\end{document}\n")
    tex_file.close()