from math import modf, floor
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, NullFormatter, LogFormatter


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

def determine_mean_median_vals(line, median_sigma_type):
    """
    Returns the mean, median, and error values in a given format.

    Parameters
    ----------
    
    line : str
        The line containing information about mean and median values.

    median_sigma_type : str
        Whether an line has concatenated '+/-' values. Accepted values are
        'single' or 'multi'.

    Returns
    -------
    A four or five length tuple with mean, mean uncertainty, median, and 
    median uncertainty values.
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
    """
    Returns all of the information regarding the best period solution.

    Parameters
    ----------

    line : str
        The line containing all the period information.
    
    """

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

    
def determine_gamma_vals(line):
    """
    Returns all of the information regarding the best gamma solution.

    Parameters
    ----------

    line : str
        The line containing all the gamma information.
    
    """
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
    """
    Returns all of the information regarding the best crater fraction solution.

    Parameters
    ----------

    line : str
        The line containing all the crater fraction information.
    
    """
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
    """
    Returns all of the information regarding the best albedo solution.

    Parameters
    ----------

    line : str
        The line containing all the albedo information.
    
    """
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


def histogram_template(directory, packed_name, values, label, unit=None):
    """
    A template function for generating histograms.

    Parameters
    ----------

    directory : str
        The name of the folder this object's results are located in.

    packed_name : str
        The packed MPC designated name.

    values : list
        The lost of plottable values.

    label : str
        The values being plotted.

    unit : str
        A unit if the quantity being plotted has units.
    
    """
    # Set up log scale plot arguments and standard deviation lines
    log_values = [np.log10(value) for value in values]
    log_values_copy = log_values[:]
    values_count = len(log_values)
    log_values_copy.sort()

    sigma_1_low = log_values_copy[int(values_count * 0.16)]
    sigma_1_high = log_values_copy[int(values_count * 0.84)]
    median = log_values_copy[int(values_count * 0.5)]
    sigma_2_low = log_values_copy[int(values_count * 0.025)]
    sigma_2_high = log_values_copy[int(values_count * 0.975)]

    # Adjust step size based on spacing
    if sigma_2_high - sigma_2_low > 0.2:
        hist_step = 0.01
    elif sigma_2_high - sigma_2_low > 0.02:
        hist_step = 0.001
    else:
        hist_step = 0.0001

    fig, ax = plt.subplots()

    hist_low_limit = int(min(log_values) * 1000) / 1000.
    hist_high_limit = (int(max(log_values) * 1000) + 1) / 1000.
    histo = ax.hist(log_values, bins=np.arange(hist_low_limit, hist_high_limit, hist_step),
            histtype="step",
            color="black")

    # Plotting standard deviation limits on plot
    ax.axvline(median, color='#cc0000')
    ax.axvline(sigma_1_low, color='#cc0000', ls='dashed')
    ax.axvline(sigma_1_high, color='#cc0000', ls='dashed')
    ax.axvline(sigma_2_low, color='#cc0000', ls='dotted')
    ax.axvline(sigma_2_high, color='#cc0000', ls='dotted')

    sigma_text = "+84% = " + f"{str(round(sigma_2_high, 3)).ljust(6, '0')}\n"
    sigma_text += "+16% = " + f"{str(round(sigma_1_high, 3)).ljust(6, '0')}\n"
    sigma_text += r"med = " + f"{str(round(median, 3)).ljust(6, '0')}\n"
    sigma_text += "-16% = " + f"{str(round(sigma_1_low, 3)).ljust(6, '0')}\n"
    sigma_text += "-84% = " + f"{str(round(sigma_2_low, 3)).ljust(6, '0')}"
    fig.text(0.82, 0.72, sigma_text, ha="center")

    n_count = len(values)
    if label == "albedo":
        ax.set_xlim(right=-0.0001)

    # Labeling plot
    label_string = f"Log {label}"
    if unit != None:
        label_string = f"Log {label} ({unit})"
    ax.set_xlabel(label_string)
    ax.set_ylabel("Number of Monte Carlo Results")
    ax.set_title(f"{packed_name}", loc="left")
    ax.set_title(f"n = {n_count}", loc="right")
    fig.savefig(f"../{directory}/general_plots/{label.replace(' ', '_').lower()}_histogram.png")
    fig.savefig(f"../{directory}/general_plots/{label.replace(' ', '_').lower()}_histogram.pdf")
    plt.close(fig)


def chi_scatterplot_template(directory, packed_name, values_x, values_y, chis, label_x, label_y, unit_x=None, unit_y=None, value_break=0.1):
    """
    A template function for generating chi valued scatterplots.

    Parameters
    ----------

    directory : str
        The name of the folder this object's results are located in.

    packed_name : str
        The packed MPC designated name.

    values_x : list
        The list of plottable values to be on the x-axis.

    values_y : list
        The list of plottable values to be on the y-axis.

    chis : list
        The list of chi^2 fit values associated with each point

    label_x : str
        The name of the values plotted on the x-axis.

    label_y : str
        The name of the values plotted on the y-axis.

    unit_x : str
        A unit of the quantity being plotted on the x-axis. OPTIONAL.
    
    unit_y : str
        A unit of the quantity being plotted on the y-axis. OPTIONAL.

    value_break : float
        The separation value needed to create a new 'segment'.
    
    """
    pairings = []

    for i in range(len(values_x)):
        pairings.append([values_x[i], values_y[i], chis[i]])
    pairings.sort()
    sorted_values_x = []
    sorted_values_y = []
    sorted_chis = []
    for pairing in pairings:
        sorted_values_x.append(pairing[0])
        sorted_values_y.append(pairing[1])
        sorted_chis.append(pairing[2])

    value_epochs = []
    value_start = sorted_values_x[0]
    for i in range(len(sorted_values_x) - 1):
        if (sorted_values_x[i + 1] - sorted_values_x[i]) > value_break:
            value_epochs.append((value_start, sorted_values_x[i]))
            value_start = sorted_values_x[i + 1]
        if i == len(sorted_values_x) - 2:
            value_epochs.append((value_start, sorted_values_x[i]))
    
    for idx, epoch in enumerate(value_epochs):
        segment_values_x = []
        segment_values_y = []
        segment_chis = []
        for i in range(len(pairings)):
            if pairings[i][0] >= epoch[0] and pairings[i][0] <= epoch[1]:
                segment_values_x.append(pairings[i][0])
                segment_values_y.append(pairings[i][1])
                segment_chis.append(pairings[i][2])
        
        fig, ax = plt.subplots()
        ax.set_xscale("log")
        ax.set_yscale("log")
        label_string_x = f"Log {label_x}"
        if unit_x != None:
            label_string_x = f"Log {label_x} ({unit_x})"
        ax.set_xlabel(label_string_x)
        label_string_y = f"Log {label_y}"
        if unit_y != None:
            label_string_y = f"Log {label_y} ({unit_y})"
        
        scatter_plot = ax.scatter(segment_values_x, segment_values_y, s=1,marker='o',c=segment_chis,linewidths=1, cmap=plt.cm.get_cmap('plasma'))
        cbar = fig.colorbar(scatter_plot, ax=ax, label=r"fit $\chi^2$")
        ax.set_title(f"{packed_name}", loc="left")
        ax.set_title(f"n = {len(segment_values_x)}", loc="right")
        title_file_name = f"{label_x.replace(' ', '_').lower()}_vs_{label_y.replace(' ', '_').lower()}"
        title_directory_name = f"{label_x.replace(' ', '_').lower()}_{label_y.replace(' ', '_').lower()}"
        fig.savefig(f"../{directory}/{title_directory_name}_segments/{idx}_{title_file_name}.png")
        fig.savefig(f"../{directory}/{title_directory_name}_segments/{idx}_{title_file_name}.pdf")
        plt.close(fig)

    # Generate base figure
    fig, ax = plt.subplots()
    ax.set_xscale("log")
    ax.set_yscale("log")
    label_string_x = f"Log {label_x.capitalize()}"
    if unit_x != None:
        label_string_x = f"Log {label_x} ({unit_x})"
    ax.set_xlabel(label_string_x)
    label_string_y = f"Log {label_y}"
    if unit_y != None:
        label_string_y = f"Log {label_y} ({unit_y})"
    ax.set_ylabel(label_string_y)
    for axis in [ax.xaxis]:
        axis.set_major_formatter(LogFormatter())
        axis.set_minor_formatter(LogFormatter(minor_thresholds=(15,0.4)))
    for axis in [ax.yaxis]:
        axis.set_major_formatter(LogFormatter())
        axis.set_minor_formatter(LogFormatter(minor_thresholds=(15,0.4)))
    scatter_plot = ax.scatter(values_x, values_y, s=1,marker='o',c=chis,linewidths=1, cmap=plt.cm.get_cmap('plasma'))
    cbar = fig.colorbar(scatter_plot, ax=ax, label=r"fit $\chi^2$")
    ax.set_title(f"{packed_name}", loc="left")
    ax.set_title(f"n = {len(values_x)}", loc="right")
    title_file_name = f"{label_x.replace(' ', '_').lower()}_vs_{label_y.replace(' ', '_').lower()}"
    fig.savefig(f"../{directory}/general_plots/{title_file_name}.png")
    fig.savefig(f"../{directory}/general_plots/{title_file_name}.pdf")
    plt.close(fig)


def hexbin_template(directory, packed_name, values_x, values_y, label_x, label_y, unit_x=None, unit_y=None, log_scale=True):
    """
    Template function for generating hexbin diagrams of two physical properties.

    Parameters
    ----------

    directory : str
        The name of the folder this object's results are located in.

    packed_name : str
        The packed MPC designated name.

    values_x : list
        The list of plottable values to be on the x-axis.

    values_y : list
        The list of plottable values to be on the y-axis.

    label_x : str
        The name of the values plotted on the x-axis.

    label_y : str
        The name of the values plotted on the y-axis.

    unit_x : str
        A unit of the quantity being plotted on the x-axis. OPTIONAL.
    
    unit_y : str
        A unit of the quantity being plotted on the y-axis. OPTIONAL.

    log_scale : bool
        Whether the plot should be created in logspace 10.
    """
    fig, ax = plt.subplots()
    ax.set_facecolor("#0e0783")
    if log_scale:
        binplot = ax.hexbin(values_x, values_y, gridsize=100, bins="log", cmap="plasma")
    elif not log_scale:
        binplot = ax.hexbin(values_x, values_y, gridsize=100, cmap="plasma")
    cbar = fig.colorbar(binplot, ax=ax, label="Number of Monte Carlo Results", pad=0.08)
    cbar.set_label("Number of Monte Carlo Results", rotation=270, labelpad=12)
    cbar.ax.yaxis.set_ticks_position("left")
    ax.set(xlim=min(values_x), ylim=min(values_y))
    label_string_x = f"Log {label_x}"
    if unit_x != None:
        label_string_x = f"Log {label_x} ({unit_x})"
    ax.set_xlabel(label_string_x)
    label_string_y = f"Log {label_y}"
    if unit_y != None:
        label_string_y = f"Log {label_y} ({unit_y})"
    ax.set_ylabel(label_string_y)
    ax.set_title(f"{packed_name}", loc="left")
    
    title_file_name = f"{label_x.replace(' ', '_').lower()}_vs_{label_y.replace(' ', '_').lower()}"
    fig.savefig(f"../{directory}/general_plots/{title_file_name}_hex.png")
    fig.savefig(f"../{directory}/general_plots/{title_file_name}_hex.pdf")
    plt.close(fig)