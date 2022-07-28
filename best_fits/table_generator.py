"""
This file takes the inputs for Ned's MCMC model and converts them into
an astropy and LaTeX table format for use in publications.
Must be run in the same directory with the run*.csh file 
"""
import os
from anyio import current_effective_deadline
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

    files = list(os.popen('ls *.txt'))
    object_files = []
    object_names = []
    for file in files:
        object_names.append(file.replace(".txt", '').replace('\n', ''))
        object_files.append(file[:len(file) - 1])
    
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
    
    initializer = " " * 4 + "\\begin{tabular}{"
    initializer_cols = ""
    

    for i in range(6):
        initializer_cols += "|c"
    initializer = initializer + initializer_cols[1:] + "}\n"

    tex_file = open("thermophysical_outputs.tex", 'w')
    tex_file.write("\\documentclass[linenumbers]{aastex631}\n\n")
    tex_file.write("\\begin{document}\n\n")
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
    tex_file.write("\end{document}\n")
    tex_file.close()

generate_object_table()
generate_LaTeX_table()