"""
python3 parser for Ned's MCMC model for spherical objects.
"""
import itertools
import os
from re import L
import matplotlib.pyplot as plt
import numpy as np
from math import modf, floor


def julian_days_utc_converter(jd):
    """
    Returns a date in utc format corresponding to the given 
    Julian day number.

    Parameters
    ----------
    jd : float
        A julian date entry.

    Returns
    -------
    A tuple representation of a utc format entry with year[0], month[1], and
    day[2].
    """
    jd_adjusted = float(jd) + 0.5
    decimal_day = modf(jd_adjusted)[0]
    x_2 = floor(jd - 1721119.5)
    c_2 = floor((4 * x_2 + 3) / 146097)
    x_1 = x_2 - floor(146097 * c_2 / 4)
    c_1 = floor((100 * x_1 + 99) / 36525)
    x_0 = x_1 - floor(36525 * c_1 / 100)
    year = 100 * c_2 + c_1
    month = floor((5 * x_0 + 461) / 153)
    day = x_0 - floor((153 * month - 457) / 5) + 1

    if month > 12:
        month = month - 12
        year = year + 1    
    
    return year, month, day + decimal_day


def retrieve_output_data(object):
    """
    Generates the diameter and albedo values for an object based on the fort.21 file.

    Parameters
    ----------
    object : str
        The folder name of the object being analyzed.

    Returns
    -------
    A tuple containing the diameters[0], albedos[1], gammas[2], 
    chi^2 values[3], and period[4] values for each solution point. 
    """
    
    output_file = open(f"../{object}/PJDFC.out")
    output_file.readline()

    diameters = []
    chis = []
    gammas = []
    periods = []
    albedos = []
    index = 0
    for line in output_file.readlines():
        albedo = np.e ** float(line.strip()[17:25].strip())
        
        period = np.e ** float(line.strip()[25:34].strip())
        gamma = np.e ** float(line.strip()[34:43].strip())
        diameter = np.e ** float(line.strip()[43:52].strip())
        chi = float(line.strip()[70:87].strip())
        diameters.append(diameter)
        chis.append(chi)
        gammas.append(gamma)
        periods.append(period)
        index += 1

    with open("fort.4", 'r') as albedo_file:
        for line in albedo_file.readlines():
            albedos.append(float(line.split()[1]))

    return diameters, albedos, gammas, chis, periods


def retrieve_MCMC_data(object):
    """
    Retrieves all of the best fit parameters and SED plotting data.

    Parameters
    ----------
    object : str
        The folder name of the object being analyzed.

    Returns
    -------
    A set of plotting parameters and epoch conditions for the modeled object.
    """
    
    os.system(f"/bin/cp ../{object}/fort.21 ../{object}/PJDFC.out")
    os.system(f"/bin/cp ../{object}/fort.22 ../{object}/SED_data.out")

    SED_file = open(f"../{object}/SED_data.out")
    print("Reading in best fit solution used for SED.")
    datum_minor = SED_file.readline().rstrip().split()
    datum_major = SED_file.readline().rstrip().split()
    
    best_fit = {}
    best_fit["pole_ra"] = np.degrees(float(datum_minor[1]))
    best_fit["pole_dec"] = np.degrees(float(datum_minor[2]))
    best_fit["pv_median"] = np.e ** (float(datum_minor[3]))
    best_fit["period"] = np.e ** (float(datum_minor[4]))
    best_fit["gamma"] = np.e ** (float(datum_minor[5]))
    best_fit["pjdfc"] = (float(datum_minor[6])) #period * thermal inertia Gamma * diameter * crater frac * color
    x = float(datum_minor[7])
    best_fit["crater_frac"] = np.e ** (x) / (1 + np.e ** x)
    best_fit["ir_albedo_ratio"] = np.e ** (float(datum_minor[8]))
    best_fit["chisq"] = float(datum_minor[10])
    best_fit["penalties"] = float(datum_minor[9]) - best_fit['chisq']
    best_fit["pv"] = float(datum_major[1])
    best_fit["diameter"] = float(datum_major[2])
    best_fit["theta"] = float(datum_major[3])

    print("Reading in all epochs")
    SED_file.readline()
    epoch_condition = {}
    epoch_condition["delta"] = []
    epoch_condition["r_helio"] = []
    epoch_condition["phase"] = []
    epoch_condition["sub_sun_lat"] = []
    epoch_condition["sub_earth_lat"] = []

    epoch = SED_file.readline()
    while epoch[0:4] != "/wvl":
        epoch_data = epoch.split()
        epoch_condition["delta"].append(float(epoch_data[1]))
        epoch_condition["r_helio"].append(float(epoch_data[2]))
        epoch_condition["phase"].append(float(epoch_data[3]))
        epoch_condition["sub_sun_lat"].append(float(epoch_data[4]))
        epoch_condition["sub_earth_lat"].append(float(epoch_data[5]))
        epoch = SED_file.readline()

    # Reads in all wavelength data for corresponding epochs
    print("Reading all wavelengths")
    wavelengths = []
    wave_line = SED_file.readline()
    while "def" not in wave_line:
        wavelength_data = wave_line.strip().split()
        for wavelength in wavelength_data:
            wavelengths.append(float(wavelength))
        wave_line = SED_file.readline()
    
    color = []
    esed = []
    data_x = []
    data_y = []
    data_y_error = []
    first_line = True

    while True:

        line = SED_file.readline()
        if "SRGB" not in line and first_line:
            color.append("#444444")
            data_x.append([])
            data_y.append([])
            data_y_error.append([])
        # Ends cycle
        if line == "":
            break
        # Gather fluxes for entire epoch
        if "/ft" in line:
            esed.append([])
            line = SED_file.readline()
            while "def doit" not in line:
                flux = float(line.strip())
                esed[-1].append(flux)
                line = SED_file.readline()
        # Acquire x and y data for epoch
        elif "QQ" in line:
            x_datum, error_datum, y_datum = line.rstrip().split()[0:3]
            data_x[-1].append(float(x_datum))
            data_y[-1].append(float(y_datum))
            data_y_error[-1].append(float(error_datum))
        elif "SRGB" in line:
            rr, gg, bb = line.rstrip().split()[0:3]
            red = hex(int(float(rr) * 255))[2:]
            green = hex(int(float(gg) * 255))[2:]
            blue = hex(int(float(bb) * 255))[2:]
            if len(red) == 1:
                red = '0' + red
            if len(green) == 1:
                green = '0' + green
            if len(blue) == 1:
                blue ='0' + blue
            
            color.append('#%2s%2s%2s'%(red, green, blue))
            data_x.append([])
            data_y.append([])
            data_y_error.append([])
        elif "setgray" in line:
            color.append("#444444")
            data_x.append([])
            data_y.append([])
            data_y_error.append([])

        first_line = False

    datelabels = []
    cshfile = os.popen(f"ls ../{object}/run*.csh").readline().rstrip()
    print("\n*** Using "+cshfile+" to get MJDs. Make sure this is right*** \n")
    cshlines = open(cshfile)
    for line in cshlines.readlines():
        epoch_data = line.rstrip().split(',')
        if len(epoch_data) != 12:
            continue
        mjd = float(epoch_data[0])
        year, month, day = julian_days_utc_converter(2400000.5 + mjd)
        date = "%4i-%02i-%02i"%(year, month, day)
        datelabels.append(date)

    best_fit_plotters = [esed, wavelengths, datelabels, data_x, data_y, data_y_error]
    return best_fit, best_fit_plotters, epoch_condition, wavelengths


def generate_SED_plot(object, esed, wavelengths, datelabels, data_x, data_y, data_y_error):
    """
    Generates a plot of best fit for wavelengths across different epochs.
    The outputted graph is in log 10 scale.

    Parameters
    ----------
    esed : list 
        A series of flux lists for each wavelengths.
    
    wavelengths : list  
        A series of wavelengths to plot for each flux list entry.
    
    datelabels : list
        A set of utc time entries for each epoch
    
    data_x : list 
        A set of positional points on the x axis
    
    data_y : list
        A set of positional points on the y axis
    
    data_y_error : list 
        Error values for the y axis data.

    Returns
    -------
    None, generates two plots titled best_fit_SED.pdf/png
    """
    print("Generating best fit plot for SED.")
    colors = ["#3498db", "#229954", "#c0392b", "#8e44ad", "#f1c40f", "#ec7063", "#34495e", "#6e2c00"]

    plt.figure(figsize=[8,6])
    plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15)
    ax = plt.gca()
    y_low = 1e99
    y_high = 0

    marker = itertools.cycle(('.', 'v', '^', 's', 'D')) 
    lines = itertools.cycle(('-', ':', '--', '-.'))
    for i in range(len(esed)):
        # Self corrects x and y limits
        if min(esed[i]) < y_low:
            y_low = min(esed[i])
        if max(esed[i]) > y_high:
            y_high = max(esed[i])

        # Plot flux for each wavelength
        plt.plot(wavelengths, esed[i], label=datelabels[i], color=colors[i], lw=0.75, ls=next(lines))    

        # Adjust y values and set up error ranges
        lower_errors = []
        higher_errors = []
        data_y_adjusted=[]
        for j in range(len(data_y[i])):
            lower_errors.append(data_y[i][j] * (1 - 1 / (10 ** (data_y_error[i][j] / 2.5))))
            higher_errors.append(data_y[i][j] * (10 ** (data_y_error[i][j] / 2.5) - 1))

            if data_y_error[i][j] > 0:    
                data_y_adjusted.append(data_y[i][j])
            else:
                # Measured flux is negative, so plot at an artifically
                # low point, and make error bars correct
                data_y_adjusted.append(1e-99)

        # Set scale and include error bars
        ax.set_xscale("log")
        ax.set_yscale("log")
        limtest = ax.get_ylim()
        
        if limtest[0] < 1e-99:
            plt.ylim(y_low / 10., y_high * 10)

        plt.errorbar(data_x[i], data_y_adjusted, 
                    yerr=[lower_errors, higher_errors],
                    ecolor=colors[i],
                    fmt=next(marker),
                    elinewidth=0.5,
                    color=colors[i],
                    ms=4,
                    capsize=2)
    ax = plt.gca()
    ax.set_facecolor("white")
    plt.legend(loc=2)
    plt.xlabel("Wavelength (microns)", fontsize=12)
    plt.ylabel(r"$\nu$F$_\nu$", fontsize=12)
    plt.title(f"{OBJECT_MAP[object]} ({object})", loc="left")
    plt.savefig(f"../{object}/general_plots/bestfit_SED.png")
    plt.savefig(f"../{object}/general_plots/bestfit_SED.pdf")
    plt.close()


def generate_diameter_histogram(object, diameters):
    """
    Generates a histogram of diameter solutions for the MCMC output.

    Parameters
    ----------
    object : str
        The folder name of the object being analyzed.

    diameters : list
        A list of diameter values for each MCMC solution

    Returns: None, generates two plots titled diameter_histogram.png/pdf
    """

    print("Generating diameter histogram.")
    # Set up log scale plot arguments and standard deviation lines
    log_diameters = [np.log10(diameter) for diameter in diameters]
    log_diameters_copy = log_diameters[:]
    diameters_count = len(log_diameters)
    log_diameters_copy.sort()

    sigma_1_low = log_diameters_copy[int(diameters_count * 0.16)]
    sigma_1_high = log_diameters_copy[int(diameters_count * 0.84)]
    sigma_2_low = log_diameters_copy[int(diameters_count * 0.025)]
    sigma_2_high = log_diameters_copy[int(diameters_count * 0.975)]

    # Adjust step size based on spacing
    if sigma_2_high - sigma_2_low > 0.2:
        hist_step = 0.01
    elif sigma_2_high - sigma_2_low > 0.02:
        hist_step = 0.001
    else:
        hist_step = 0.0001

    plt.figure(figsize=[8, 6])
    
    hist_low_limit = int(min(log_diameters) * 1000) / 1000.
    hist_high_limit = (int(max(log_diameters) * 1000) + 1) / 1000.
    plt.hist(log_diameters, bins=np.arange(hist_low_limit, hist_high_limit, hist_step),
            histtype="step",
            color="black")

    # Plotting standard deviation limits on plot
    plt.axvline(sigma_1_low, color='#cc0000')
    plt.axvline(sigma_1_high, color='#cc0000')
    plt.axvline(sigma_2_low, color='#cc0000',ls='dotted')
    plt.axvline(sigma_2_high, color='#cc0000',ls='dotted')

    # Labeling plot
    plt.xlabel("Log Diameter (km)")
    plt.ylabel("Number of Monte Carlo Results")
    plt.title(f"{OBJECT_MAP[object]} ({object})", loc="left")
    plt.savefig(f"../{object}/general_plots/diameter_histogram.png")
    plt.savefig(f"../{object}/general_plots/diameter_histogram.pdf")


def generate_albedo_histogram(object, albedos):
    """
    Generates a histogram of diameter solutions for the MCMC output.

    Parameters
    ----------
    object : str
        The folder name of the object being analyzed.

    albedos : list
        A list of albedo values for each MCMC solution

    Returns: None, generates two plots titled albedo_histogram.png/pdf
    """

    print("Generating albedo histogram.")
    # Set up log scale plot arguments and standard deviation lines
    log_albedos = [np.log10(albedo) for albedo in albedos]
    log_albedos_copy = log_albedos[:]
    albedos_count = len(log_albedos)
    log_albedos_copy.sort()

    sigma_1_low = log_albedos_copy[int(albedos_count * 0.16)]
    sigma_1_high = log_albedos_copy[int(albedos_count * 0.84)]
    sigma_2_low = log_albedos_copy[int(albedos_count * 0.025)]
    sigma_2_high = log_albedos_copy[int(albedos_count * 0.975)]

    # Adjust step size based on spacing
    if sigma_2_high - sigma_2_low > 0.2:
        hist_step = 0.01
    elif sigma_2_high - sigma_2_low > 0.02:
        hist_step = 0.001
    else:
        hist_step = 0.0001

    plt.figure(figsize=[8, 6])
    
    hist_low_limit = int(min(log_albedos) * 1000) / 1000.
    hist_high_limit = (int(max(log_albedos) * 1000) + 1) / 1000.
    plt.hist(log_albedos, bins=np.arange(hist_low_limit, hist_high_limit, hist_step),
            histtype="step",
            color="black")

    # Plotting standard deviation limits on plot
    plt.axvline(sigma_1_low, color='#cc0000')
    plt.axvline(sigma_1_high, color='#cc0000')
    plt.axvline(sigma_2_low, color='#cc0000',ls='dotted')
    plt.axvline(sigma_2_high, color='#cc0000',ls='dotted')

    # Labeling plot
    plt.xlabel("Log Albedo")
    plt.ylabel("Number of Monte Carlo Results")
    plt.title(f"{OBJECT_MAP[object]} ({object})", loc="left")
    plt.savefig(f"../{object}/general_plots/albedo_histogram.png")
    plt.savefig(f"../{object}/general_plots/albedo_histogram.pdf")


def generate_gamma_histogram(object, gammas):
    """
    Generates a histogram of diameter solutions for the MCMC output.

    Parameters
    ----------
    object : str
        The folder name of the object being analyzed.

    gammas : list
        A list of gamma values for each MCMC solution

    Returns: None, generates two plots titled gamma_histogram.png/pdf
    """

    print("Generating gamma histogram.")
    # Set up log scale plot arguments and standard deviation lines
    log_gammas = [np.log10(gamma) for gamma in gammas]
    log_gammas_copy = log_gammas[:]
    gammas_count = len(log_gammas)
    log_gammas_copy.sort()

    sigma_1_low = log_gammas_copy[int(gammas_count * 0.16)]
    sigma_1_high = log_gammas_copy[int(gammas_count * 0.84)]
    sigma_2_low = log_gammas_copy[int(gammas_count * 0.025)]
    sigma_2_high = log_gammas_copy[int(gammas_count * 0.975)]

    # Adjust step size based on spacing
    if sigma_2_high - sigma_2_low > 0.2:
        hist_step = 0.01
    elif sigma_2_high - sigma_2_low > 0.02:
        hist_step = 0.001
    else:
        hist_step = 0.0001

    plt.figure(figsize=[8, 6])
    
    hist_low_limit = int(min(log_gammas) * 1000) / 1000.
    hist_high_limit = (int(max(log_gammas) * 1000) + 1) / 1000.
    plt.hist(log_gammas, bins=np.arange(hist_low_limit, hist_high_limit, hist_step),
            histtype="step",
            color="black")

    # Plotting standard deviation limits on plot
    plt.axvline(sigma_1_low, color='#cc0000')
    plt.axvline(sigma_1_high, color='#cc0000')
    plt.axvline(sigma_2_low, color='#cc0000',ls='dotted')
    plt.axvline(sigma_2_high, color='#cc0000',ls='dotted')

    # Labeling plot
    plt.xlabel(r"Log Thermal inertia ($J~m^{-2}~s^{-0.5}~K^{-1})$)", fontsize=13)
    plt.ylabel("Number of Monte Carlo Results")
    plt.title(f"{OBJECT_MAP[object]} ({object})", loc="left")
    plt.savefig(f"../{object}/general_plots/gamma_histogram.png")
    plt.savefig(f"../{object}/general_plots/gamma_histogram.pdf")


def generate_period_histogram(object, periods):
    """
    Generates a histogram of diameter solutions for the MCMC output.

    Parameters
    ----------
    object : str
        The folder name of the object being analyzed.

    periods : list
        A list of period values for each MCMC solution

    Returns: None, generates two plots titled period_histogram.png/pdf
    """

    print("Generating period histogram.")
    if max(periods) - min(periods) == 0:
        print("Fixed period used, skipping period histogram.")
    # Set up log scale plot arguments and standard deviation lines
    log_periods = [np.log10(period) for period in periods]
    log_periods_copy = log_periods[:]
    periods_count = len(log_periods)
    log_periods_copy.sort()

    sigma_1_low = log_periods_copy[int(periods_count * 0.16)]
    sigma_1_high = log_periods_copy[int(periods_count * 0.84)]
    sigma_2_low = log_periods_copy[int(periods_count * 0.025)]
    sigma_2_high = log_periods_copy[int(periods_count * 0.975)]

    # Adjust step size based on spacing
    if sigma_2_high - sigma_2_low > 0.2:
        hist_step = 0.01
    elif sigma_2_high - sigma_2_low > 0.02:
        hist_step = 0.001
    else:
        hist_step = 0.0001

    plt.figure(figsize=[8, 6])
    
    hist_low_limit = int(min(log_periods) * 1000) / 1000.
    hist_high_limit = (int(max(log_periods) * 1000) + 1) / 1000.
    plt.hist(log_periods, bins=np.arange(hist_low_limit, hist_high_limit, hist_step),
            histtype="step",
            color="black")

    # Plotting standard deviation limits on plot
    plt.axvline(sigma_1_low, color='#cc0000')
    plt.axvline(sigma_1_high, color='#cc0000')
    plt.axvline(sigma_2_low, color='#cc0000',ls='dotted')
    plt.axvline(sigma_2_high, color='#cc0000',ls='dotted')

    # Labeling plot
    plt.xlabel(r"Log Period (hr)", fontsize=13)
    plt.ylabel("Number of Monte Carlo Results")
    plt.title(f"{OBJECT_MAP[object]} ({object})", loc="left")
    plt.savefig(f"../{object}/general_plots/period_histogram.png")
    plt.savefig(f"../{object}/general_plots/period_histogram.pdf")


def generate_diameter_vs_albedo_plot(object, diameters, albedos, chis):
    """
    Generates a plot of diameter and albedo solutions for the MCMC output.

    Parameters
    ----------
    object : str
        The folder name of the object being analyzed.

    diameters : list
        A list of diameter values for each MCMC solution.

    albedos : list
        A list of albedo values for each MCMC solution.

    chis : list
        A list of chi^2 values for each MCMC solution.

    Returns
    -------
    None, generates two plots titled diameter_vs_albedo.png/pdf
    """

    print("Generating diameter vs albedo plots.")
    for i in range(48600):
        if albedos[i] > 1:
            print(albedos[i])
    pairings = []
    diameter_break = 0.1
    print(len(diameters))

    for i in range(len(diameters)):
        pairings.append([diameters[i], albedos[i], chis[i]])
    pairings.sort()
    sorted_diameters = []
    sorted_albedos = []
    sorted_chis = []
    for pairing in pairings:
        sorted_diameters.append(pairing[0])
        sorted_albedos.append(pairing[1])
        sorted_chis.append(pairing[2])

    diameter_epochs = []
    diameter_start = sorted_diameters[0]
    for i in range(len(sorted_diameters) - 1):
        if (sorted_diameters[i + 1] - sorted_diameters[i]) > diameter_break:
            diameter_epochs.append((diameter_start, sorted_diameters[i]))
            diameter_start = sorted_diameters[i + 1]
        if i == len(sorted_diameters) - 2:
            diameter_epochs.append((diameter_start, sorted_diameters[i]))

    
    if "diameter_albedo_segments" not in os.listdir(f'../{object}/'):
        os.popen(f"mkdir ../{object}/diameter_albedo_segments")
        print("Creating diameter_vs_albedo subplot directory")
    
    for idx, epoch in enumerate(diameter_epochs):
        segment_diameters = []
        segment_albedos = []
        segment_chis = []
        for i in range(len(pairings)):
            if pairings[i][0] >= epoch[0] and pairings[i][0] <= epoch[1]:
                segment_diameters.append(pairings[i][0])
                segment_albedos.append(pairings[i][1])
                segment_chis.append(pairings[i][2])
        plt.figure(figsize=[8,6])
        ax = plt.gca()
        ax.set_xscale("log")
        ax.set_yscale("log")
        plt.scatter(segment_diameters, segment_albedos, s=1,marker='o',c=segment_chis,linewidths=1, cmap=plt.cm.get_cmap('plasma'))
        plt.colorbar(label=r"fit $\chi^2$")
        plt.xlabel("Log Diameter (km)")
        plt.ylabel("Log Albedo")
        plt.title(f"{OBJECT_MAP[object]} ({object})", loc="left")
        plt.title(f"{len(segment_diameters)}/{len(diameters)} points displayed", loc="right")
        plt.savefig(f"../{object}/diameter_albedo_segments/{idx}_diameter_vs_albedo.png")
        plt.savefig(f"../{object}/diameter_albedo_segments/{idx}_diameter_vs_albedo.pdf")
        plt.close()

    # Generate base figure
    plt.figure(figsize=[8,6])
    ax = plt.gca()
    ax.set_xscale("log")
    ax.set_yscale("log")
    plt.scatter(diameters, albedos, s=1,marker='o',c=chis,linewidths=1, cmap=plt.cm.get_cmap('plasma'))
    plt.colorbar(label=r"fit $\chi^2$")
    plt.xlabel("Log Diameter (km)")
    plt.ylabel("Log Albedo")
    plt.title(f"{OBJECT_MAP[object]} ({object})", loc="left")
    plt.title(f"{len(diameters)} points displayed", loc="right")
    plt.savefig(f"../{object}/general_plots/diameter_vs_albedo.png")
    plt.savefig(f"../{object}/general_plots/diameter_vs_albedo.pdf")
    plt.close()


def generate_diameter_vs_albedo_hex(object, diameters, albedos):
    """
    Generates a hexplot of the diameter and albedo solutions for an MCMC run.

    Parameters
    ----------
    object : str
        The folder specific name of the object being modeled.

    diameters : list
        

    """
    fig, ax = plt.subplots()
    binplot = ax.hexbin(diameters, albedos, gridsize=100, bins="log", cmap='plasma')#, cmap=plt.cm.get_cmap('plasma'))
    cbar = fig.colorbar(binplot, ax=ax, label="Number of Monte Carlo Results", pad=0.08)
    cbar.set_label("Number of Monte Carlo Results", rotation=270, labelpad=12)
    cbar.ax.yaxis.set_ticks_position("left")
    ax.set(xlim=min(diameters), ylim=min(albedos))
    ax.set_xlabel("Diameter (km)") 
    ax.set_ylabel("Albedo")
    ax.set_title(f"{OBJECT_MAP[object]} ({object})", loc="left")
    fig.savefig(f"../{object}/general_plots/diameter_vs_albedo_hex.png")
    fig.savefig(f"../{object}/general_plots/diameter_vs_albedo_hex.pdf")
    plt.close()

def generate_diameter_vs_period_plot(object, diameters, periods, chis):
    """
    Generates a plot of diameter and period solutions for the MCMC output.

    Parameters
    ----------
    object : str
        The folder name of the object being analyzed.

    diameters : list
        A list of diameter values for each MCMC solution.

    periods : list
        A lsit of period values for each MCMC solution.

    Returns
    -------
    None, generates two plots titled diameter_vs_period.png/pdf.
    """
    print("Generating diameter vs period plots.")

    if max(periods) - min(periods) == 0:
        print("Fixed period used, skipping diameter vs period plot")
    else:
        pairings = []
        diameter_break = 0.1
        for i in range(len(diameters)):
            pairings.append([diameters[i], periods[i], chis[i]])
        pairings.sort()
        sorted_diameters = []
        sorted_periods = []
        sorted_chis = []
        for pairing in pairings:
            sorted_diameters.append(pairing[0])
            sorted_periods.append(pairing[1])
            sorted_chis.append(pairing[2])
        
        # Segments diameter clusters to create solution cluster subplots.
        diameter_epochs = []
        diameter_start = sorted_diameters[0]
        for i in range(len(sorted_diameters) - 1):
            if (sorted_diameters[i + 1] - sorted_diameters[i]) > diameter_break:
                diameter_epochs.append((diameter_start, sorted_diameters[i]))
                diameter_start = sorted_diameters[i + 1]
            if i == len(sorted_diameters) - 2:
                diameter_epochs.append((diameter_start, sorted_diameters[i]))

        if "diameter_period_segments" not in os.listdir(f'../{object}/'):
            os.popen(f"mkdir ../{object}/diameter_period_segments")
            print("Creating diameter_vs_period subplot directory")
        
        for idx, epoch in enumerate(diameter_epochs):
            segment_diameters = []
            segment_periods = []
            segment_chis = []
            for i in range(len(pairings)):
                if pairings[i][0] >= epoch[0] and pairings[i][0] <= epoch[1]:
                    segment_diameters.append(pairings[i][0])
                    segment_periods.append(pairings[i][1])
                    segment_chis.append(pairings[i][2])
            plt.figure(figsize=[8,6])
            ax = plt.gca()
            ax.set_xscale("log")
            ax.set_yscale("log")
            plt.scatter(segment_diameters, segment_periods, s=1,marker='o',c=segment_chis,linewidths=1, cmap=plt.cm.get_cmap('plasma'))
            plt.colorbar(label=r"fit $\chi^2$")
            plt.text(0.5, .97,f"{len(segment_diameters)}/{len(diameters)} points displayed", ha='center', va='center', transform=ax.transAxes) 
            plt.xlabel("Diameter (km)")
            plt.ylabel("Rotation period")
            plt.title(f"{OBJECT_MAP[object]} ({object})", loc="left")
            plt.title(f"{len(segment_diameters)}/{len(diameters)} points displayed", loc="right")
            plt.savefig(f"../{object}/diameter_period_segments/{idx}_diameter_vs_period.png")
            plt.savefig(f"../{object}/diameter_period_segments/{idx}_diameter_vs_period.pdf")
            plt.close()
            
        # Plot the figure with all points
        plt.figure(figsize=[8,6])
        ax = plt.gca()
        ax.set_xscale("log")
        ax.set_yscale("log")
        plt.scatter(diameters, periods, s=1,marker='o',c=chis,linewidths=1, cmap=plt.cm.get_cmap('plasma'))
        plt.colorbar(label=r"fit $\chi^2$")
        plt.text(0.5, .97,f"{len(diameters)} points displayed", ha='center', va='center', transform=ax.transAxes) 
        plt.xlabel("Log Diameter (km)")
        plt.ylabel("Log Rotation period (hr)")
        plt.title(f"{OBJECT_MAP[object]} ({object})", loc="left")
        plt.title(f"{len(diameters)} points displayed", loc="right")
        plt.savefig(f"../{object}/general_plots/diameter_vs_period.png")
        plt.savefig(f"../{object}/general_plots/diameter_vs_period.pdf")
        plt.close()


def generate_diameter_vs_gamma_plot(object, diameters, gammas, chis):
    """
    Generates a plot of diameter and gamma solutions for the MCMC output.

    Parameters
    ----------
    object : str
        The folder name of the object being analyzed.

    diameters : list
        A list of diameter values for each MCMC solution.

    gammas : list
        A list of gamma values for each MCMC solution.

    chis : list
        A list of chi^2 values for each MCMC solution.

    Returns
    -------
    None, generates two plots titled diameter_vs_gamma.png/pdf
    """
    print("Generating diameter vs gamma plots.")
    pairings = []
    diameter_break = 0.1

    for i in range(len(diameters)):
        pairings.append([diameters[i], gammas[i], chis[i]])
    
    pairings.sort()
    sorted_diameters = []
    sorted_gammas = []
    sorted_chis = []
    for pairing in pairings:
        sorted_diameters.append(pairing[0])
        sorted_gammas.append(pairing[1])
        sorted_chis.append(pairing[2])

    # Segments diameter clusters to create solution cluster subplots.
    diameter_epochs = []
    diameter_start = sorted_diameters[0]
    for i in range(len(sorted_diameters) - 1):
        if (sorted_diameters[i + 1] - sorted_diameters[i]) > diameter_break:
            diameter_epochs.append((diameter_start, sorted_diameters[i]))
            diameter_start = sorted_diameters[i + 1]
        if i == len(sorted_diameters) - 2:
            diameter_epochs.append((diameter_start, sorted_diameters[i]))

    if "diameter_gamma_segments" not in os.listdir(f'../{object}/'):
        os.popen(f"mkdir ../{object}/diameter_gamma_segments")
        print("Creating diameter_vs_gamma subplot directory") 
    for idx, epoch in enumerate(diameter_epochs):
        segment_diameters = []
        segment_gammas = []
        segment_chis = []
        for i in range(len(pairings)):
            if pairings[i][0] >= epoch[0] and pairings[i][0] <= epoch[1]:
                segment_diameters.append(pairings[i][0])
                segment_gammas.append(pairings[i][1])
                segment_chis.append(pairings[i][2])
        plt.figure(figsize=[8,6])
        ax = plt.gca()
        ax.set_xscale("log")
        ax.set_yscale("log")
        plt.scatter(segment_diameters, segment_gammas, s=1,marker='o',c=segment_chis,linewidths=1, cmap=plt.cm.get_cmap('plasma'))
        plt.colorbar(label=r"fit $\chi^2$")
        plt.xlabel("Diameter (km)")
        plt.ylabel("Rotation period")
        plt.title(f"{OBJECT_MAP[object]} ({object})", loc="left")
        plt.title(f"{len(segment_diameters)}/{len(diameters)} points displayed", loc="right")
        plt.savefig(f"../{object}/diameter_gamma_segments/{idx}_diameter_vs_gamma.png")
        plt.savefig(f"../{object}/diameter_gamma_segments/{idx}_diameter_vs_gamma.pdf")
        plt.close()

    # Plot the figure with all points
    plt.figure(figsize=[8,6])
    ax = plt.gca()
    plt.scatter(diameters, gammas, s=1, marker='o', c=chis, linewidths=1, cmap=plt.cm.get_cmap('plasma'))
    plt.colorbar(label=r"fit $\chi^2$")
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.xlabel("Log Diameter (km)")
    plt.ylabel(r"Log Thermal inertia ($J~m^{-2}~s^{-0.5}~K^{-1})$)", fontsize=13)
    plt.title(f"{OBJECT_MAP[object]} ({object})", loc="left")
    plt.title(f"{len(diameters)} points displayed", loc="right")
    plt.savefig(f"../{object}/general_plots/diameter_vs_gamma.png")
    plt.savefig(f"../{object}/general_plots/diameter_vs_gamma.pdf")
    plt.close()


def generate_diameter_vs_chi_plots(object, diameters, chis):
    """
    Generates a plot of diameter and chi solutions for the MCMC output.

    Parameters
    ----------
    object : str
        The folder name of the object being analyzed.

    diameters : list
        A list of diameter values for each MCMC solution.

    chis : list
        A list of chis^2 values for each MCMC solution.

    Returns
    -------
    None, generates two plots titled diameter_vs_chi.png/pdf and two zoomed plots 
    titled diameter_vs_chi_zoom.png/pdf.
    """

    print("Generating diameter vs chi plots.")
    plt.figure(figsize=[8, 6])
    plt.subplots_adjust(left=0.15, right=0.95, top=0.95)
    plt.loglog(diameters, chis, color='k', marker='o', ms=0.5, ls="None")
    plt.xlabel("Log Diameter (km)", fontsize=13)
    plt.ylabel(r"Log fit $\chi^2$", fontsize=13)
    plt.title(f"{OBJECT_MAP[object]} ({object})", loc="left")
    plt.savefig(f"../{object}/general_plots/diameter_vs_chi.png")
    plt.savefig(f"../{object}/general_plots/diameter_vs_chi.pdf")
    plt.close()

    plt.figure(figsize=[8, 6])
    plt.subplots_adjust(left=0.15, right=0.95, top=0.95)
    plt.loglog(diameters, chis, color='k', marker='o', ms=0.5, ls='None')
    plt.xlabel("Log Diameter (km)", fontsize=13)
    plt.ylabel(r"Log fit $\chi^2$", fontsize=13)
    plt.title(f"{OBJECT_MAP[object]} ({object})", loc="left")
    plt.ylim(0.9 * min(chis), 2 * min(chis))  
    plt.savefig(f"../{object}/general_plots/diameter_vs_chi_zoom.png")
    plt.savefig(f"../{object}/general_plots/diameter_vs_chi_zoom.pdf")
    plt.close()

def generate_gamma_vs_chi_plots(object, gammas, chis):
    """
    Generates a plot of gamma and chi solutions for the MCMC output.

    Arguments
    ---------
    object : str
        The folder name of the object being analyzed.
    
    gammas : list
        A list of gamma values for each MCMC solution.
    
    chis : list
        A list of chi^2 values for each MCMC solution. 

    Returns
    -------
    None, generates two plots titled gamma_vs_chi.png/pdf
    """
    print("Generating gamma vs chi plots.")
    plt.figure(figsize=[8, 6])
    plt.subplots_adjust(left=0.15, right=0.95, top=0.95)
    plt.loglog(gammas, chis, color='k', marker='o', ms=0.5, ls='None')
    plt.xlabel(r"Log Thermal inertia ($J~m^{-2}~s^{-0.5}~K^{-1})$)", fontsize=13)
    plt.ylabel(r"Log fit $\chi^2$", fontsize=13)
    plt.title(f"{OBJECT_MAP[object]} ({object})", loc="left")
    plt.savefig(f"../{object}/general_plots/gamma_vs_chi.png")
    plt.savefig(f"../{object}/general_plots/gamma_vs_chi.pdf")
    plt.close()


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
    mean_sigma = ""
    while True:
        if mean_sigma == "+/-":
            mean_sigma = ""
        if line[index] == "m":
            mean_sigma = mean_sigma.strip()
            break
        mean_sigma += line[index]
        index += 1
    median = ""
    while True:
        if median == "median":
            median = ""
        if line[index] == "+":
            median = median.strip()
            break
        median += line[index]
        index += 1
    
    if median_sigma_type == "multi":
        median_pos_sigma = ""
        median_neg_sigma = ""
        while True:
            if median_pos_sigma == "+":
                median_pos_sigma = ""
            if line[index] == '-':
                median_pos_sigma = median_pos_sigma.strip()
                break
            median_pos_sigma += line[index]
            index += 1
        while True:
            if median_neg_sigma == "-":
                median_neg_sigma = ""
            if line[index] == '%':
                median_neg_sigma = median_neg_sigma.strip()
                break
            median_neg_sigma += line[index]
            index += 1
        outputs = [mean, mean_sigma, median, median_pos_sigma, median_neg_sigma]
        return outputs
    
    if median_sigma_type == "single":
        median_sigma = ""
        while True:
            if median_sigma == "+/-":
                median_sigma = ""
            if line[index] == '%':
                median_sigma = median_sigma.strip()
                break
            median_sigma += line[index]
            index += 1
        outputs = [mean, mean_sigma, median, median_sigma]
        return mean, mean_sigma, median, median_sigma


def determine_period(line):
    period = ""
    period_pos_sigma = ""
    period_neg_sigma = ""
    index = 0
    while True:
        if period == "Period [h] =":
            period = ""
        if line[index] == '+':
            period = period.strip()
            break
        period += line[index]
        index += 1
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
    while True:
        if period_neg_sigma == "-":
            period_neg_sigma = ""
        if line[index] == '%':
            period_neg_sigma = period_neg_sigma.strip()
            outputs = [period, period_pos_sigma, period_neg_sigma]
            return outputs
        period_neg_sigma += line[index]
        index += 1

    
def determine_square_vals(line):
    value = ""
    value_pos_sigma = ""
    value_neg_sigma = ""
    index = 0
    while True:
        if value == "sqrt(kappa*rho*C)=":
            value = ""
        if line[index] == '+':
            value = value.strip()
            break
        value += line[index]
        index += 1
    while True:
        if value_pos_sigma == '+':
            value_pos_sigma = ""
        if line[index] == '-':
            value_pos_sigma = value_pos_sigma.strip()
            break
        value_pos_sigma += line[index]
        index += 1
    while True:
        if value_neg_sigma == "-":
            value_neg_sigma = ""
        if line[index] == '%':
            value_neg_sigma = value_neg_sigma.strip()
            break
        value_neg_sigma += line[index]
        index += 1
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
    while True:
        if fraction_pos_sigma == '+':
            fraction_pos_sigma = ""
        if line[index] == '-':
            fraction_pos_sigma = fraction_pos_sigma.strip()
            break
        fraction_pos_sigma += line[index]
        index += 1
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
    while True:
        if ratio_pos_sigma == '+':
            ratio_pos_sigma = ""
        if line[index] == '-':
            ratio_pos_sigma = ratio_pos_sigma.strip()
            break
        ratio_pos_sigma += line[index]
        index += 1
    while True:
        if ratio_neg_sigma == '-':
            ratio_neg_sigma = ""
        if line[index] == '%':
            ratio_neg_sigma = ratio_neg_sigma.strip()
            break
        ratio_neg_sigma += line[index]
        index += 1
    outputs = [ratio, ratio_pos_sigma, ratio_neg_sigma]
    return outputs


def display_MCMC_results(mpc_name):
    mpc_data = []
    with open(f"../best_fits/{mpc_name}.txt", 'r') as file:
        for i in range(5):
            file.readline()
        line_6 = file.readline()
        diameter_vals = determine_mean_median_vals(line_6, "multi")
        output_6 = f"Diameter (Mean): {diameter_vals[0]} +/- {diameter_vals[1]}"
        print(output_6)
        output_7 = f"Diameter (Median): {diameter_vals[2]} +{diameter_vals[3]}%/-{diameter_vals[4]}%"
        print(output_7)
        line_8 = file.readline()
        p_V_vals = determine_mean_median_vals(line_8, "single")
        output_8 = f"p_V (Mean): {p_V_vals[0]} +/- {p_V_vals[1]}"
        print(output_8)
        output_9 = f"p_V (Median): {p_V_vals[2]} +/- {p_V_vals[3]}%"
        print(output_9)
        line_10 = file.readline()
        theta_vals = determine_mean_median_vals(line_10, "multi")
        output_10 = f"Theta_1 (Mean): {theta_vals[0]} +/- {theta_vals[1]}"
        print(output_10)
        output_11 = f"Theta_1 (Median): {theta_vals[2]} +{theta_vals[3]}%/-{theta_vals[4]}%"
        print(output_11)
        
        line_12 = file.readline()
        period_vals = determine_period(line_12)
        output_12 = f"Period (hours): {period_vals[0]} +{period_vals[1]}%/-{period_vals[2]}%"
        print(output_12)

        line_13 = file.readline()
        sqrt_vals = determine_square_vals(line_13)
        output_13 = f"Sqrt(Kappa*Rho*C): {sqrt_vals[0]} +{sqrt_vals[1]}%/{sqrt_vals[2]}%"
        print(output_13)

        line_14 = file.readline()
        crater_vals = determine_crater_fraction(line_14)
        output_14 = f"Crater Fraction: {crater_vals[0]} +{crater_vals[1]}/-{crater_vals[2]}"
        print(output_14)

        line_15 = file.readline()
        ratio_vals = determine_p_V_ratio(line_15)
        output_15 = f"p_IR/p_V: {ratio_vals[0]} +{ratio_vals[1]}%/-{ratio_vals[2]}%"
        print(output_15)

        line_16 = file.readline().split()
        output_16 = f"Pole peak at: {line_16[4]} {line_16[5]}"
        print(output_16)
        line_17 = file.readline().split()
        output_17 = f"Mean pole at: {line_17[4]} {line_17[5]} {line_17[6][0:5]} = {line_17[7]}"
        print(output_17)
        line_18 = file.readline().split()
        output_18 = f"Moment eigenvector: {line_18[2]} at RA,DEC: "
        output_18 += f"{line_18[5]} {line_18[6]}"
        print(output_18)
        line_19 = file.readline().split()
        output_19 = f"Moment eigenvector: {line_19[2]} at RA,DEC: "
        output_19 += f"{line_19[5]} {line_19[6]}"
        print(output_19)
        line_20 = file.readline().split()
        output_20 = f"Moment eigenvector: {line_20[2]} at RA,DEC: "
        output_20 += f"{line_20[5]} {line_20[6]}\n"
        print(output_20)


OBJECT_MAP = {"02100": "Ra-Shalom", 
              "02212": "Hephaistos", 
              "05189": "1990 UQ",
              "05693": "1993 EA", 
              "07335": "1989 JA", 
              "23606": "1996 AS1",
              "68950": "2002 QF15", 
              "85713": "1998 SS49",
              "G1989": "Cacus"}

objects = ["02100", "02212", "05189", "05693", "07335", "23606", "68950", "85713", "G1989"]

for object in objects[0:1]:
    print(f"{OBJECT_MAP[object]} ({object})")
    print("Echoing relevant files\n")
    os.system(f"echo '../{object}/PJDFC.out' | ../read-WISE-rc-MCMC-PJDFC") 
    os.system(f"/bin/cp ../{object}/fort.2 ../{object}/Dhist.dat")
    os.system(f"/bin/cp ../{object}/fort.32 ../{object}/Dhist_fine.dat")
    os.system(f"/bin/cp ../{object}/fort.3 ../{object}/DvsPeriod.dat")
    os.system(f"/bin/cp ../{object}/fort.4 ../{object}/DvsAlb.dat")

    best_fit, best_fit_plotters, epoch_condition, wavelengths = retrieve_MCMC_data(object)
    diameters, albedos, gammas, chis, periods = retrieve_output_data(object)
    if "general_plots" not in os.listdir(f'../{object}/.'):
        os.popen(f"mkdir ../{object}/general_plots")
        print("Creating plotting directory")
    generate_SED_plot(object, best_fit_plotters[0], best_fit_plotters[1], 
                    best_fit_plotters[2], best_fit_plotters[3], 
                    best_fit_plotters[4], best_fit_plotters[5])
    #generate_diameter_histogram(object, diameters)
    #generate_albedo_histogram(object, albedos)
    #generate_gamma_histogram(object, gammas)
    #generate_period_histogram(object, periods)
    #generate_diameter_vs_albedo_plot(object, diameters, albedos, chis)
    generate_diameter_vs_albedo_hex(object, diameters, albedos)
    #generate_diameter_vs_period_plot(object, diameters, periods, chis)
    #generate_diameter_vs_gamma_plot(object, diameters, gammas, chis)
    #generate_diameter_vs_chi_plots(object, diameters, chis)
    #generate_gamma_vs_chi_plots(object, gammas, chis)