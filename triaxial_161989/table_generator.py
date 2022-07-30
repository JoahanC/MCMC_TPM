"""
This file takes the inputs for Ned's MCMC model and converts them into
an astropy and LaTeX table format for use in publications.
Must be run in the same directory with the run*.csh file 
"""
import os
from anyio import current_effective_deadline
import numpy as np
from astropy.table import Table


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
    
    with open("thermophysical_outputs.tex", 'w') as tex_file:
        tex_file.write("\\begin{table}[H]\n")
        tex_file.write(' ' * 4 + "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|}\n")
        tex_file.write(' ' * 4 + "\\centering\n")
        tex_file.write(' ' * 4 + header + "\\\ \n")
        tex_file.write(' ' * 4 + subheader[:len(subheader) - 3] + "\\\ \n")
        tex_file.write(' ' * 4 + "\\end{tabular}\n")
        tex_file.write("\end{table}")

irsa_files, csh_file = call_input_files()
generate_MCMC_input_astropy_table(csh_file)
generate_LaTeX_Obs_input_table()
