"""
This file takes the inputs for Ned's MCMC model and converts them into
an astropy and LaTeX table format for use in publications.
Must be run in the same directory with the run*.csh file 
"""
import os
from anyio import current_effective_deadline
import numpy as np
from astropy.table import Table
from sklearn.metrics import mean_absolute_error


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


def generate_MCMC_input_astropy_table(csh_file, triaxial=True):
    """
    Generates an astropy table of the MCMC input file with the .csh
    file input as a reference.

    Arguments: csh_file (str) -- the name of the file housing the 
               direct MCMC inputs

    Returns: None -- Generates an .tbl file in the same directory
    """

    csh_inputs = read_csh_inputs(csh_file, triaxial)
    epochs = []
    for key in csh_inputs:
        if "epoch" in key:
            epochs.append(csh_inputs[key])

    header_names = ["Epoch", "RA", "DEC", "MJD", "Dist"]
    for i in range(4):
        header_names.append(f"W{i + 1}_Flux")
        header_names.append(f"W{i + 1}_Sigma")
        header_names.append(f"W{i + 1}_Amp")

    object_table = Table(rows=epochs, names=header_names)
    unpacked_designation = csh_inputs["up_desig"]
    object_table.write(f"{unpacked_designation}_MCMC_inputs.tbl", format="ipac", 
                       overwrite=True)


def generate_LaTeX_Obs_input_table():

    input_files = os.popen("ls *MCMC_inputs.tbl").readlines()
    if len(input_files) != 1:
        print("")
    table_file = input_files[0].replace('\n', '')
    epochs = Table.read(table_file, format="ipac")

    header = "    Epoch & RA & DEC & MJD & Source ID & "
    subheader = "     & & & & & "

    for i in range(4):
        subheader += "Mag & Mag $\sigma$ & "
        subheader += "Flux & Flux $\sigma$ & $\chi$ & "
        header += f" & & & W{i + 1} & & "

    print(header)
    print(subheader)

    with open("table.tex", 'w') as tex_file:
        tex_file.write("\\begin{table}[H]\n")
        tex_file.write(' ' * 4 + "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|}\n")
        tex_file.write(' ' * 4 + "\\centering\n")
        tex_file.write(' ' * 4 + header + "\\\ \n")
        tex_file.write(' ' * 4 + subheader[:len(subheader) - 3] + "\\\ \n")
        tex_file.write(' ' * 4 + "\\end{tabular}\n")
        tex_file.write("\end{table}")
    

def generate_LaTeX_thermo_results_table():

    files = os.popen("ls *.txt").readlines()
    for i in range(len(files)):
        files[i] = files[i][:len(files[i]) - 1]
    values = {"diameters" : [], "p_V_vals" : [], "thetas" : [], "periods" : [], "sqrt_vals" : [],
              "crater_vals" : [], "ratio_vals" : [], "pole_peak" : [], "mean_pole" : [],
              "moment_eigenvector_1" : [], "moment_eigenvector_2" : [], "moment_eigenvector_3" : []}
    for fit_file in files:
        with open(fit_file, 'r') as file:
            
            file.readline()
            file.readline()
            file.readline()
            file.readline()
            file.readline()
            line_6 = file.readline()
            diameter_vals = determine_mean_median_vals(line_6, "multi")
            for value in diameter_vals:
                values["diameters"].append(value)
            line_8 = file.readline()
            p_V_vals = determine_mean_median_vals(line_8, "single")
            for value in p_V_vals:
                values["p_V_vals"].append(value)
            line_10 = file.readline()
            theta_vals = determine_mean_median_vals(line_10, "multi")
            for value in theta_vals:
                values["thetas"].append(value)
            line_12 = file.readline()
            period_vals = determine_period(line_12)
            for value in period_vals:
                values["periods"].append(value)
            line_13 = file.readline()
            sqrt_vals = determine_square_vals(line_13)
            for value in sqrt_vals:
                values["sqrt_vals"].append(value)
            line_14 = file.readline()
            crater_vals = determine_crater_fraction(line_14)
            for value in crater_vals:
                values["crater_vals"].append(value)
            line_15 = file.readline()
            ratio_vals = determine_p_V_ratio(line_15)
            for value in ratio_vals:
                values["ratio_vals"].append(value)
            line_16 = file.readline().split()
            values["pole_peak"].append(line_16[4])
            values["pole_peak"].append(line_16[5])
            line_17 = file.readline().split()
            values["mean_pole"].append(line_17[4])
            values["mean_pole"].append(line_17[5])
            values["mean_pole"].append(line_17[7])
            line_18 = file.readline().split()
            values["moment_eigenvector_1"].append(line_18[2])
            values["moment_eigenvector_1"].append(line_18[5])
            values["moment_eigenvector_1"].append(line_18[6])
            line_19 = file.readline().split()
            values["moment_eigenvector_2"].append(line_19[2])
            values["moment_eigenvector_2"].append(line_19[5])
            values["moment_eigenvector_2"].append(line_19[6])
            line_20 = file.readline().split()
            values["moment_eigenvector_3"].append(line_20[2])
            values["moment_eigenvector_3"].append(line_20[5])
            values["moment_eigenvector_3"].append(line_20[6])

    initializer = " " * 4 + "\\begin{tabular}{"
    initializer_cols = ""
    for i in range(len(files) + 1):
        initializer_cols += "|c"
    initializer = initializer + initializer_cols[1:] + "}\n"
    header = " " * 8 + "MPC Designation & "

    objects = []
    for file in files:
        name = conv_num2mpcformat(int(file.replace(".txt", '')))
        objects.append(name)
    for object in objects:
        header += f"{object} & "

    with open("thermophysical_outputs.tex", 'w') as tex_file:
        tex_file.write("\\begin{table}[H]\n")
        tex_file.write(' ' * 4 + "\\centering\n")
        tex_file.write(initializer)
        tex_file.write(" " * 8 + "\hline\n")
        tex_file.write(" " * 8 + "\hline\n")
        tex_file.write(header[:len(header) - 2] + '\\\ \n')
        tex_file.write(" " * 8 + "\hline\n")

        mean_diameters = " " * 8 + "Mean Diameter & "
        median_diameters = " " * 8 + "Median Diameter & "
        mean_p_V = " " * 8 + "Mean p$_V$ & "
        median_p_V = " " * 8 + "Median p$_V$ & "
        mean_thetas = " " * 8 + "Mean Theta & "
        median_thetas = " " * 8 + "Median Theta & "
        periods = " " * 8 + "Period & "
        sqrt_vals = " " * 8 + "$\sqrt{\kappa\cdot\\rho\cdot C}$ & "
        crater_fraction = " " * 8 + "Crater Fraction & "
        p_IR_p_V = " " * 8 + "p$_{IR}$/p$_V$ & "
        pole_peak_RA = " " * 8 + "Pole Peak RA & "
        pole_peak_DEC = " " * 8 + "Pole Peak DEC & "
        mean_pole_RA = " " * 8 + "Mean Pole RA & "
        mean_pole_DEC = " " * 8 + "Mean Pole DEC & "
        mean_vector = " " * 8 + "|<p>| & "
        eigenvector_1_val = " " * 8 + "Eigenvector 1 & "
        eigenvector_2_val = " " * 8 + "Eigenvector 2 & "
        eigenvector_3_val = " " * 8 + "Eigenvector 3 & "
        eigenvector_1_RA = " " * 8 + "Eigenvector 1 RA & "
        eigenvector_2_RA = " " * 8 + "Eigenvector 2 RA & "
        eigenvector_3_RA = " " * 8 + "Eigenvector 3 RA & "
        eigenvector_1_DEC = " " * 8 + "Eigenvector 1 DEC & "
        eigenvector_2_DEC = " " * 8 + "Eigenvector 2 DEC & "
        eigenvector_3_DEC = " " * 8 + "Eigenvector 3 DEC & "
        for i in range(0, len(objects)):
            mean_diameters += f"${values['diameters'][i * 5]}^" + '{' f"+{values['diameters'][i * 5 + 1]}" + "}_{" + f"-{values['diameters'][i * 5 + 1]}" + "}$ & "
            median_diameters += f"${values['diameters'][i * 5 + 2]}^" + '{' f"+{values['diameters'][i * 5 + 3]}\%" + "}_{" + f"-{values['diameters'][i * 5 + 4]}\%" + "}$ & "
            mean_p_V += f"${values['p_V_vals'][i * 4]}^" + '{' f"+{values['p_V_vals'][i * 4 + 1]}\%" + "}_{" + f"-{values['p_V_vals'][i * 4 + 1]}\%" + "}$ & "
            median_p_V += f"${values['p_V_vals'][i * 4 + 2]}^" + '{' f"+{values['p_V_vals'][i * 4 + 3]}\%" + "}_{" + f"-{values['p_V_vals'][i * 4 + 3]}\%" + "}$ & "
            mean_thetas += f"${values['thetas'][i * 5]}^" + '{' f"+{values['thetas'][i * 5 + 1]}" + "}_{" + f"-{values['thetas'][i * 5 + 1]}" + "}$ & "
            median_thetas += f"${values['thetas'][i * 5 + 2]}^" + '{' f"+{values['thetas'][i * 5 + 3]}\%" + "}_{" + f"-{values['thetas'][i * 5 + 4]}\%" + "}$ & "
            periods += f"${values['periods'][i * 3]}^" + '{' f"+{values['periods'][i * 3 + 1]}" + "}_{" + f"-{values['periods'][i * 3 + 2]}" + "}$ & "
            sqrt_vals += f"${values['sqrt_vals'][i * 3]}^" + '{' f"+{values['sqrt_vals'][i * 3 + 1]}" + "}_{" + f"-{values['sqrt_vals'][i * 3 + 2]}" + "}$ & "
            crater_fraction += f"${values['crater_vals'][i * 3]}^" + '{' f"+{values['crater_vals'][i * 3 + 1]}" + "}_{" + f"-{values['crater_vals'][i * 3 + 2]}" + "}$ & "
            p_IR_p_V += f"${values['ratio_vals'][i * 3]}^" + '{' f"+{values['ratio_vals'][i * 3 + 1]}" + "}_{" + f"-{values['ratio_vals'][i * 3 + 2]}" + "}$ & "
            pole_peak_RA += f"${values['pole_peak'][i * 2]}$ & " 
            pole_peak_DEC += f"${values['pole_peak'][i * 2 + 1]}$ & "
            mean_pole_RA += f"${values['mean_pole'][i * 3]}$ & " 
            mean_pole_DEC += f"${values['mean_pole'][i * 3 + 1]}$ & "
            mean_vector += f"${values['mean_pole'][i * 3 + 2]}$ & "
            eigenvector_1_val += f"${values['moment_eigenvector_1'][i * 3]}$ & "
            eigenvector_2_val += f"${values['moment_eigenvector_2'][i * 3]}$ & "
            eigenvector_3_val += f"${values['moment_eigenvector_3'][i * 3]}$ & "
            eigenvector_1_RA += f"${values['moment_eigenvector_1'][i * 3 + 1]}$ & "
            eigenvector_2_RA += f"${values['moment_eigenvector_2'][i * 3 + 1]}$ & "
            eigenvector_3_RA += f"${values['moment_eigenvector_3'][i * 3 + 1]}$ & "
            eigenvector_1_DEC += f"${values['moment_eigenvector_1'][i * 3 + 2]}$ & "
            eigenvector_2_DEC += f"${values['moment_eigenvector_2'][i * 3 + 2]}$ & "
            eigenvector_3_DEC += f"${values['moment_eigenvector_3'][i * 3 + 2]}$ & "
        tex_file.write(mean_diameters[:len(mean_diameters) - 2] + '\\\ \n')
        tex_file.write(median_diameters[:len(median_diameters) - 2] + '\\\ \n')
        tex_file.write(mean_p_V[:len(mean_p_V) - 2] + '\\\ \n')
        tex_file.write(median_p_V[:len(median_p_V) - 2] + '\\\ \n')
        tex_file.write(mean_thetas[:len(mean_thetas) - 2] + '\\\ \n')
        tex_file.write(median_thetas[:len(median_thetas) - 2] + '\\\ \n')
        tex_file.write(periods[:len(periods) - 2] + '\\\ \n')
        tex_file.write(sqrt_vals[:len(sqrt_vals) - 2] + '\\\ \n')
        tex_file.write(pole_peak_RA[:len(pole_peak_RA) - 2] + '\\\ \n')
        tex_file.write(pole_peak_DEC[:len(pole_peak_DEC) - 2] + '\\\ \n')
        tex_file.write(mean_pole_RA[:len(mean_pole_RA) - 2] + '\\\ \n')
        tex_file.write(mean_pole_DEC[:len(mean_pole_DEC) - 2] + '\\\ \n')
        tex_file.write(eigenvector_1_val[:len(eigenvector_1_val) - 2] + '\\\ \n')
        tex_file.write(eigenvector_1_RA[:len(eigenvector_1_RA) - 2] + '\\\ \n')
        tex_file.write(eigenvector_1_DEC[:len(eigenvector_1_DEC) - 2] + '\\\ \n')
        tex_file.write(eigenvector_2_val[:len(eigenvector_2_val) - 2] + '\\\ \n')
        tex_file.write(eigenvector_2_RA[:len(eigenvector_2_RA) - 2] + '\\\ \n')
        tex_file.write(eigenvector_2_DEC[:len(eigenvector_2_DEC) - 2] + '\\\ \n')
        tex_file.write(eigenvector_3_val[:len(eigenvector_3_val) - 2] + '\\\ \n')
        tex_file.write(eigenvector_3_RA[:len(eigenvector_3_RA) - 2] + '\\\ \n')
        tex_file.write(eigenvector_3_DEC[:len(eigenvector_3_DEC) - 2] + '\\\ \n')
        tex_file.write(' ' * 4 + "\\end{tabular}\n")
        caption_line = ' ' * 4 + "\\caption{Thermophysical modeling results"
        caption_line += " generated by the spherical MCMC model.}\n"
        tex_file.write(caption_line)
        tex_file.write("\end{table}")
        

#irsa_files, csh_file = call_input_files()
#generate_MCMC_input_astropy_table(csh_file)
generate_LaTeX_thermo_results_table()
