from math import modf, floor
import numpy as np
import matplotlib.pyplot as plt


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

    sigma_text = r"+2$\sigma$ = " + f"{str(round(sigma_2_high, 3)).ljust(6, '0')}\n"
    sigma_text += r"+1$\sigma$ = " + f"{str(round(sigma_1_high, 3)).ljust(6, '0')}\n"
    sigma_text += r"med = " + f"{str(round(median, 3)).ljust(6, '0')}\n"
    sigma_text += r"-1$\sigma$ = " + f"{str(round(sigma_1_low, 3)).ljust(6, '0')}\n"
    sigma_text += r"-2$\sigma$ = " + f"{str(round(sigma_2_low, 3)).ljust(6, '0')}"
    fig.text(0.82, 0.72, sigma_text, ha="center")

    n_count = len(values)
    if label == "albedo":
        ax.set_xlim(right=-0.0001)

    # Labeling plot
    label_string = f"Log {label.capitalize()}"
    if unit != None:
        label_string = f"Log {label} ({unit})"
    ax.set_xlabel(label_string)
    ax.set_ylabel("Number of Monte Carlo Results")
    ax.set_title(f"{packed_name}", loc="left")
    ax.set_title(f"n = {n_count}", loc="right")
    fig.savefig(f"../{directory}/general_plots/{label.replace(' ', '_').lower()}_histogram.png")
    fig.savefig(f"../{directory}/general_plots/{label.replace(' ', '_').lower()}_histogram.pdf")
    plt.close(fig)