"""
This file takes the inputs for Ned's MCMC model and converts them into
an astropy and LaTeX table format for use in publications.
Must be run in the same directory with the run*.csh file 
"""
import sys
import os
from matplotlib.pyplot import table
import numpy as np
from astropy.table import Table, Column


def conv_num2mpcformat(n):

    """

    TBD

    :param n:

    :return:

    """

    nn2c = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"

    if n < 100000:

        return "%05d" % n

    elif n < 620000:

        nn = int(n / 10000.0)

        c = nn2c[nn]

        return "{:1s}{:04d}".format(c, n - nn * 10000)

    else:

        # For numbers larger than 620,000 the MPC has defined a new packing

        # scheme. Code by J. Masiero

        nn = n - 620000

        dig4 = int(nn % 62)

        hold3 = (nn - dig4) / 62.0

        dig3 = int(hold3 % 62)

        hold2 = (hold3 - dig3) / 62.0

        dig2 = int(hold2 % 62)

        hold1 = (hold2 - dig2) / 62.0

        dig1 = int(hold1 % 62)

        return "~{:1s}{:1s}{:1s}{:1s}".format(

            nn2c[dig1], nn2c[dig2], nn2c[dig3], nn2c[dig4]

        )


def unpack_MPC_name(packed_name):
    """
    Generates the unpacked MPC designation for a given object.

    Arguments: packed_name (str) -- The packed name of the object.

    Returns: (str) -- The unpacked name of the object
    """

    char_map = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    char_map += "abcdefghijklmnopqrstuvwxyz"

    if len(packed_name) == 5:

        if packed_name[4] == 'P':
            # Periodic Comet
            unpacked_name = str(int(packed_name[0:4])) + 'P'

        elif packed_name[4] == 'S':
            # Satellite
            unpacked_name = packed_name[0] + '/'
            unpacked_name += str(int(packed_name[1:4])) + 'S'

        else:
            # Numbered Object
            input_1 = char_map.index(packed_name[0]) * 10000
            input_2 = int(packed_name[1:])
            unpacked_name = '(' + str(input_1 + input_2) + ')'

    elif len(packed_name) == 7:

        # PLS Object
        if packed_name[0:3] in ["PLS", "T1S", "T2S", "T3S"]:
            unpacked_name = packed_name[3:] + " " + packed_name[0]
            unpacked_name += '-' + packed_name[1]

        else:
            # Unnumbered Asteroid
            year = char_map.index(packed_name[0]) * 100
            year += int(packed_name[1:3])

            if packed_name[4:6] == "00":
                prov = packed_name[3] + packed_name[6]

            else:
                input_1 = char_map.index(packed_name[4]) * 10
                input_2 = int(packed_name[5])
                prov = packed_name[3] + packed_name[6]
                prov += str(input_1 + input_2)
    
            unpacked_name = str(year) + " " + prov
    
    elif len(packed_name) == 8:
        # Parabolic Comet
        year = char_map.index(packed_name[1]) * 100
        year += int(packed_name[2:4])
        
        if packed_name[7] not in ['0', '1', '2', '3', '4', '5', 
                                  '6', '7', '8', '9']:
            prov = packed_name[4] + packed_name[7]
            input_1 = char_map.index(packed_name[5]) * 10
            input_2 = int(packed_name[6])
            prov += str(input_1 + input_2)
            unpacked_name = packed_name[0] + '/'
            unpacked_name += str(year) + " " + prov
        
        else:
            unpacked_name = packed_name[0] + '/' + str(year)
            unpacked_name += " " + packed_name[4]
            unpacked_name += str(int(packed_name[5:7]))
    
    else:
        print("Did not recognize input format of:", packed_name)
        unpacked_name = packed_name
    
    return unpacked_name


def regex_loop(line, value, index, start_key, stop_key):
    while True:
        if value == start_key:
            value = ""
        if line[index] == stop_key:
            value = value.strip()
            break
        value += line[index]
        index += 1
    return value, index

def determine_mean_median_vals(line, median_sigma_type):
    """
    Returns the mean, median, and error values in a given format.

    Arguments: line (str) -- The line containing all values
    """
    
    index = 0
    mean = ""
    
    while True:
        if mean == "dia=" or mean == "p_V =" or mean == "theta1=" :
            mean = ""
        if line[index] == '+' or line[index] == '+/-':
            mean = mean.strip()
            break
        mean += line[index]
        index += 1
    
    mean_sigma, index = regex_loop(line, "", index, "+/-", 'm')
    median, index = regex_loop(line, "", index, "median", '+')
    
    if median_sigma_type == "multi":
        
        median_pos_sigma, index = regex_loop(line, "", index, '+', '-')
        median_neg_sigma, index = regex_loop(line, "", index, '-', '%')
        outputs = [mean, mean_sigma, median, median_pos_sigma, median_neg_sigma]
        
        return outputs
    
    if median_sigma_type == "single":
        
        median_sigma, index = regex_loop(line, "", index, "+/-", '%')
        outputs = [mean, mean_sigma, median, median_sigma]
        
        return mean, mean_sigma, median, median_sigma


def determine_period(line):
    period_pos_sigma = ""
    period_neg_sigma = ""
    period, index = regex_loop(line, "", 0, "Period [h] =", '+')
    while True:
        if period_pos_sigma.strip() == "0.0   0.0%":
            outputs = [period, "0.0", "0.0"]
            return outputs
        if period_pos_sigma == '+':
            period_pos_sigma = ""
        if line[index] == '-':
            period_pos_sigma = period_pos_sigma.strip()
            break
        period_pos_sigma += line[index]
        index += 1
    period_neg_sigma, index = regex_loop(line, "", index, '-', '%')

    outputs = [period, period_pos_sigma, period_neg_sigma]
    return outputs

    
def determine_square_vals(line):
    value, index = regex_loop(line, "", 0, "sqrt(kappa*rho*C)=", '+')
    value_pos_sigma, index = regex_loop(line, "", index, '+', '-')
    value_neg_sigma, index = regex_loop(line, "", index, '-', '%')
    outputs = [value, value_pos_sigma, value_neg_sigma]
    return outputs

def determine_crater_fraction(line):
    fraction = ""
    fraction_pos_sigma = ""
    fraction_neg_sigma = ""
    index = 0
    while True:
        if line[index] == '+':
            fraction = fraction.strip()[17:]
            break
        fraction += line[index]
        index += 1
    fraction_pos_sigma, index = regex_loop(line, "", index, '+', '-')
    while index < len(line):
        if fraction_neg_sigma == '-':
            fraction_neg_sigma = ""
        fraction_neg_sigma += line[index]
        index += 1
    outputs = [fraction, fraction_pos_sigma, fraction_neg_sigma.strip()]
    return outputs


def determine_p_V_ratio(line):
    ratio = ""
    ratio_pos_sigma = ""
    ratio_neg_sigma = ""
    index = 0
    while True:
        if line[index] == '+':
            ratio = ratio.strip()[9:]
            break
        ratio += line[index]
        index += 1
    ratio_pos_sigma, index = regex_loop(line, "", index, '+', '-')
    ratio_neg_sigma, index = regex_loop(line, "", index, '-', '%')
    outputs = [ratio, ratio_pos_sigma, ratio_neg_sigma]
    return outputs


def call_input_files():
    """
    Returns all the input files in this directory.

    Arguments: (None)

    Returns: (tup) A list of all the input irsa files[0] and the
            .csh file used to run the MCMC model
    """

    irsa_files = os.popen("ls irsa*.tbl").readlines()
    for i in range(len(irsa_files)):
        irsa_files[i] = irsa_files[i].replace('\n', '')
    csh_file = os.popen("ls run*.csh").readlines()[0].replace('\n', '')

    return irsa_files, csh_file


def read_csh_inputs(csh_file, triaxial=True):
    """
    Reads in all the MCMC inptu information from the .csh file

    Arguments: csh_file (str) -- the name of the file housing the 
               direct MCMC inputs

    Returns: (dict) -- contains all relevant information in the
             the firle
    """

    csh_inputs = {}

    with open(csh_file, 'r') as input_file:
        
        for i in range(3):
            input_file.readline()
        
        header_info = input_file.readline().split(',')
        csh_inputs["h_mag"] = header_info[0]
        csh_inputs["h_error"] = header_info[1]
        csh_inputs["period"] = header_info[2]
        csh_inputs["up_desig"] = header_info[3].replace('\n', '')
        
        if triaxial:
            input_file.readline()

        current_line = input_file.readline()
        index = 1
        while True:
            if "LAST" in current_line:
                break
            
            csh_inputs[f"epoch_{index}"] = [index]
            for datum in current_line.split(','):
                datum.replace('+', '').replace('\n', '')
                csh_inputs[f"epoch_{index}"].append(float(datum))
            
            current_line = input_file.readline()
            index += 1

        return csh_inputs


def generate_object_table():
    """
    Generates an astropy table of all thermophysical property
    outputs from existing objects in the directory.

    Arguments: csh_file (str) -- the name of the file housing the 
               direct MCMC inputs

    Returns: None -- Generates an .tbl file in the same directory
    """

    
    #files = list(os.popen('ls *.txt'))
    #files = ["01990.txt", "2002.txt", "2100.txt", "02212.txt", "5693.txt", "7735.txt", "23606.txt", "85713.txt", "G1819.txt"]
    files = []
    curr = os.listdir(".")
    files = [neo.replace("\n", '') for neo in curr if ".txt" in neo]
    
    object_files = []
    object_names = []
    for file in files:
        if "triaxial" in file:
            object_names.append(file.replace(".txt", '').replace('\n', '').replace("triaxial_", "Triaxial "))
            object_files.append(file[:len(file)])
        else:
            object_names.append(file.replace(".txt", '').replace('\n', ''))
            object_files.append(file[:len(file)])
    
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

    sqrt = []
    sqrt_pos_sig = []
    sqrt_neg_sig = []

    crater_frac = []
    crater_frac_pos_sig = []
    crater_frac_neg_sig = []

    ratio = []
    ratio_pos_sig = []
    ratio_neg_sig = []

    pole_peak_ra = []
    pole_peak_dec = []

    mean_pole_ra = []
    mean_pole_dec = []
    mean_pole_val = []

    eigenvector_1 = []
    eigenvector_1_ra = []
    eigenvector_1_dec = []

    eigenvector_2 = []
    eigenvector_2_ra = []
    eigenvector_2_dec = []

    eigenvector_3 = []
    eigenvector_3_ra = []
    eigenvector_3_dec = []

    for file in object_files:
        object_file = open(file, 'r')
        for i in range(5):
            object_file.readline()
        
        diameter_vals = determine_mean_median_vals(object_file.readline(), "multi")
        diameter_mean.append(diameter_vals[0])
        diameter_mean_sig.append(diameter_vals[1])
        diameter_median.append(diameter_vals[2])
        diameter_median_pos_sig.append(diameter_vals[3])
        diameter_median_neg_sig.append(diameter_vals[4])

        p_V_vals = determine_mean_median_vals(object_file.readline(), "single")
        p_V_mean.append(p_V_vals[0])
        p_V_mean_sig.append(p_V_vals[1])
        p_V_median.append(p_V_vals[2])
        p_V_median_sig.append(p_V_vals[3])

        theta_vals = determine_mean_median_vals(object_file.readline(), "multi")
        theta_mean.append(theta_vals[0])
        theta_mean_sig.append(theta_vals[1])
        theta_median.append(theta_vals[2])
        theta_median_pos_sig.append(theta_vals[3])
        theta_median_neg_sig.append(theta_vals[4])

        period_vals = determine_period(object_file.readline())
        period.append(period_vals[0])
        period_pos_sig.append(period_vals[1])
        period_neg_sig.append(period_vals[2])

        sqrt_vals = determine_square_vals(object_file.readline())
        sqrt.append(sqrt_vals[0])
        sqrt_pos_sig.append(sqrt_vals[1])
        sqrt_neg_sig.append(sqrt_vals[2])

        crater_vals = determine_crater_fraction(object_file.readline())
        crater_frac.append(crater_vals[0])
        crater_frac_pos_sig.append(crater_vals[1])
        crater_frac_neg_sig.append(crater_vals[2])

        ratio_vals = determine_p_V_ratio(object_file.readline())
        ratio.append(ratio_vals[0])
        ratio_pos_sig.append(ratio_vals[1])
        ratio_neg_sig.append(ratio_vals[2])

        if "triaxial" in file:
            object_file.readline()
            object_file.readline()

        pole_peak_line = object_file.readline().split()
        pole_peak_ra.append(pole_peak_line[4])
        pole_peak_dec.append(pole_peak_line[5])

        mean_pole_line = object_file.readline().split()
        mean_pole_ra.append(mean_pole_line[4])
        mean_pole_dec.append(mean_pole_line[5])
        mean_pole_val.append(mean_pole_line[7])

        eigenvector_1_line = object_file.readline().split()
        eigenvector_1.append(eigenvector_1_line[2])
        eigenvector_1_ra.append(eigenvector_1_line[5])
        eigenvector_1_dec.append(eigenvector_1_line[6])

        eigenvector_2_line = object_file.readline().split()
        eigenvector_2.append(eigenvector_2_line[2])
        eigenvector_2_ra.append(eigenvector_2_line[5])
        eigenvector_2_dec.append(eigenvector_2_line[6])

        eigenvector_3_line = object_file.readline().split()
        eigenvector_3.append(eigenvector_3_line[2])
        eigenvector_3_ra.append(eigenvector_3_line[5])
        eigenvector_3_dec.append(eigenvector_3_line[6])

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
    object_table['kappa_rho_c'] = sqrt
    object_table['kappa_rho_c_pos_sig'] = sqrt_pos_sig
    object_table['kappa_rho_c_neg_sig'] = sqrt_neg_sig
    object_table['crater_fraction'] = crater_frac
    object_table['crater_fraction_pos_sig'] = crater_frac_pos_sig
    object_table['crater_fraction_neg_sig'] = crater_frac_neg_sig
    object_table['p_IRp_V'] = ratio
    object_table['p_IRp_V_pos_sig'] = ratio_pos_sig
    object_table['p_IRp_V_neg_sig'] = ratio_neg_sig
    object_table['pole_peak_ra'] = pole_peak_ra
    object_table['pole_peak_dec'] = pole_peak_dec
    object_table['mean_pole'] = mean_pole_val
    object_table['mean_pole_ra'] = mean_pole_ra
    object_table['mean_pole_dec'] = mean_pole_dec
    object_table['eigenvector_1'] = eigenvector_1
    object_table['eigenvector_1_ra'] = eigenvector_1_ra
    object_table['eigenvector_1_dec'] = eigenvector_1_dec
    object_table['eigenvector_2'] = eigenvector_2
    object_table['eigenvector_2_ra'] = eigenvector_2_ra
    object_table['eigenvector_2_dec'] = eigenvector_2_dec
    object_table['eigenvector_3'] = eigenvector_3
    object_table['eigenvector_3_ra'] = eigenvector_3_ra
    object_table['eigenvector_3_dec'] = eigenvector_3_dec
    object_table.write(f"MCMC_outputs.tbl", format="ipac", 
                       overwrite=True)


def generate_LaTeX_table():

    data_object = Table.read("MCMC_outputs.tbl", format='ipac')

    tex_file = open("thermophysical_table.tex", 'w')    
    
    tex_file.write("\\begin{deluxetable*}{" + 'c'*6 + "}\n")
    tex_file.write(' ' * 4 + "\\tablenum{1}\n")
    tex_file.write(' ' * 4 + "\\tablecaption{Physical characteristics from thermophysical modeling of various objects}\n")
    tex_file.write(' ' * 4 + "\\tablewidth{0pt}\n")
    tex_file.write(' ' * 4 + "\\tablehead{\n")
    col_names = ["Name", "Diameter", "Albedo", "Theta", "Period", "Crater Fraction"]
    table_header = ' ' * 8
    for name in col_names:
        table_header += "\\colhead{" + name + "} & "
    table_header = table_header[:len(table_header) - 2] + "\\\ \n" + ' ' * 8
    col_types = ['', 'km', '', 'deg', 'hr', '']
    for type in col_types:
        table_header += "\\colhead{" + type + "} & "
    tex_file.write(table_header[:len(table_header) - 2] + "\n")
    tex_file.write(' ' * 4 + '}\n')
    tex_file.write(' ' * 4 + "\decimalcolnumbers\n" + ' ' * 4 + "\startdata\n")
    

    for iter in range(len(data_object['name'])):
        line = ' ' * 8
        for idx, col in enumerate(data_object.colnames):
            if idx == 0:
                line += f"{data_object[col][iter]} & "
            elif idx in [3, 8, 12, 15, 21]:
                line += f"${data_object[col][iter]}^"
            elif idx in [4, 9, 13, 16]:
                line += "{+" + data_object[col][iter] + "\%}_"
            elif idx in [22]:
                line += "{+" + data_object[col][iter] + "}_"
            if idx in [5, 9, 14, 17]:
                line += "{-" + data_object[col][iter] + "\%}$ & "
            elif idx == 23:
                line += "{-" + data_object[col][iter] + "}$ \\\ \n"
        tex_file.write(line)
    tex_file.write(' ' * 4 + "\enddata\n")
    tex_file.write("\end{deluxetable*}\n")
    tex_file.close()


def generate_MCMC_results():
    data_object = Table.read("MCMC_outputs.tbl", format='ipac')
    tex_file = open("thermophysical_outputs.tex", 'w')
    tex_file.write("\makeatletter\n")
    tex_file.write("\declare@file@substitution{revtex4-1.cls}{revtex4-2.cls}\n")
    tex_file.write("\makeatother\n")
    tex_file.write("\\documentclass[linenumbers]{aastex631}\n\n")
    tex_file.write("\\usepackage{float}\n\n")

    tex_file.write("\\shorttitle{MCMC Thermophysical Modeling Results}\n")
    tex_file.write("\\shortauthors{Castaneda Jaimes, Macias}\n")
    tex_file.write("\\begin{document}\n\n")
    tex_file.write("\\title{MCMC Thermophysical Modeling Results \\footnote{Summer 2022}}\n")
    tex_file.write("\\author{Joahan Castaneda Jaimes}\n")

    tex_file.write("\\begin{abstract}\n")
    tex_file.write("The Wide-field Infrared Survey Explorer (WISE)" +
    " telescope has been scanning the sky for over a decade," +
    " uncovering a vast population of Near Earth Asteroids (NEAs)" + 
    " as part of the NEOWISE mission. NEAs are of interest to the" + 
    " scientific community due to their proximity and chaotic" + 
    " nature, making them a potential threat to Earth while also" + 
    " promising targets for future missions. In this study, we" + 
    " assess asteroids by first recovering observational epochs" + 
    " missed by NEOWISE’s detection software and then using all" + 
    " valid observational epochs to run a triaxial ellipsoidal" + 
    " thermophysical model utilizing Monte Carlo Markov Chain (MCMC)" +
    " techniques. We present predictions of the diameter, albedo," + 
    " thermal inertia, and other physical characteristics for these" + 
    " asteroids. Additionally, we report newly discovered epochs for" +
    " these objects to the International Astronomical Union’s Minor" + 
    " Planet Center and develop software tools for discovering more" +
    " missed epochs by NEOWISE across the entire WISE database in an" +
    " automated fashion.\n")
    tex_file.write("\\end{abstract}\n")
    
    tex_file.write("\\section{Introduction}\n")

    tex_file.write("All the results from the spherical MCMC" + 
    " thermophysical application can be found in this file.\n")
    tex_file.write("Included is a series of tables with major" +
    " physical characteristics and a series of plots for each\n")
    tex_file.write("object modeled.\n")

    tex_file.write("\\begin{deluxetable*}{" + 'c'*6 + "}[H]\n")
    tex_file.write(' ' * 4 + "\\tablenum{1}\n")
    tex_file.write(' ' * 4 + "\\tablecaption{Physical characteristics from thermophysical modeling of various objects}\n")
    tex_file.write(' ' * 4 + "\\tablewidth{10pt}\n")
    tex_file.write(' ' * 4 + "\\tablehead{\n")
    col_names = ["Name", "Diameter", "Albedo", "Theta", "Period", "Crater Fraction"]
    table_header = ' ' * 8
    for name in col_names:
        table_header += "\\colhead{" + name + "} & "
    table_header = table_header[:len(table_header) - 2] + "\\\ \n" + ' ' * 8
    col_types = ['', 'km', '', 'deg', 'hr', '']
    for type in col_types:
        table_header += "\\colhead{" + type + "} & "
    tex_file.write(table_header[:len(table_header) - 2] + "\n")
    tex_file.write(' ' * 4 + '}\n')
    tex_file.write(' ' * 4 + "\decimalcolnumbers\n" + ' ' * 4 + "\startdata\n")
    

    for iter in range(len(data_object['name'])):
        line = ' ' * 8
        for idx, col in enumerate(data_object.colnames):
            if idx == 0:
                line += f"{data_object[col][iter]} & "
            elif idx in [3, 8, 12, 15, 21]:
                line += f"${data_object[col][iter]}^"
            elif idx in [4, 9, 13, 16]:
                line += "{+" + data_object[col][iter] + "\%}_"
            elif idx in [22]:
                line += "{+" + data_object[col][iter] + "}_"
            if idx in [5, 9, 14, 17]:
                line += "{-" + data_object[col][iter] + "\%}$ & "
            elif idx == 23:
                line += "{-" + data_object[col][iter] + "}$ \\\ \n"
        tex_file.write(line)
    tex_file.write(' ' * 4 + "\enddata\n")
    tex_file.write("\end{deluxetable*}\n\n")
    """plot_paths = return_all_image_files()

    for neo in plot_paths:
        tex_file.write("\\section*{\LARGE " + neo + "}\n")
        for path in plot_paths[neo]:
            plot_image(tex_file, path)"""
    tex_file.write("\end{document}\n")
    tex_file.close()


def plot_image(tex_file, image_pdf):
    tex_file.write("\\begin{figure}[H]\n")
    tex_file.write(' ' * 4 + "\\plotone{" + image_pdf + "}\n")
    tex_file.write("\\end{figure}\n\n")


def return_all_image_files():
    """
    Returns all of the existing image files located in the MCMC TPM folder.
    """
    curr = os.listdir(".")
    current_neos = [neo.replace(".txt", '') for neo in curr if ".txt" in neo]
    neo_paths = {}
    for neo in current_neos:
        neo_paths[neo] = f"../{neo}/general_plots/"

    plot_names = ["bestfit_SED.pdf", "diameter_histogram.pdf", 
    "diameter_vs_albedo.pdf", "diameter_vs_chi.pdf", 
    "diameter_vs_gamma.pdf", "diameter_vs_period.pdf",
    "gamma_vs_chi.pdf"]
    plot_paths = {}

    for neo in neo_paths:
        neo_dir = os.listdir(neo_paths[neo])
        for plot in plot_names:
            if plot in neo_dir:
                if neo in plot_paths:
                    plot_paths[neo].append(neo_paths[neo] + plot)
                else:
                    plot_paths[neo] = [neo_paths[neo] + plot]
    return plot_paths


#generate_object_table()
#generate_LaTeX_table()
generate_MCMC_results()
#os.popen("pdflatex thermophysical_outputs.tex")
#return_all_image_files()