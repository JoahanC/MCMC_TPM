"""
This file generates a pdf with all the plots and findings from the MCMC
modeling runs.
"""
import os
from table_generator import generate_LaTeX_table



def create_pdf():

    os.popen("cp thermophysical_outputs.tex results.tex")

generate_LaTeX_table()
create_pdf()