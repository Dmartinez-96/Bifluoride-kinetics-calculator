# -*- coding: utf-8 -*-
"""
Bifluoride salt thermal decomposition kinetics curve fitting program.

Author: Dakotah Martinez
"""

import csv
import warnings
import pandas as pd
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")


def zero_order(inputs, arrhen, e_activ):
    """
    Return zero order fit of data.

    Parameters
    ----------
    inputs : Input parameters for fit (time (seconds), Temp. (Kelvin))
    arrhen : Arrhenius prefactor.
    e_activ : Activation energy in kcal/mol.

    """
    time_in, temp = inputs
    my_k = arrhen * np.exp(-e_activ / (ID_GAS_CONST * temp))
    return my_k * time_in


def first_order(inputs, arrhen, e_activ):
    """
    Return first order fit of data.

    Parameters
    ----------
    inputs : Input parameters for fit (time (seconds), Temp. (Kelvin))
    arrhen : Arrhenius prefactor.
    e_activ : Activation energy in kcal/mol.

    """
    time_in, temp = inputs
    my_k = arrhen * np.exp(-e_activ / (ID_GAS_CONST * temp))
    return 1 - np.exp(-my_k * time_in)


def second_order(inputs, arrhen, e_activ):
    """
    Return second order fit of data.

    Parameters
    ----------
    inputs : Input parameters for fit (time (seconds), Temp. (Kelvin))
    arrhen : Arrhenius prefactor.
    e_activ : Activation energy in kcal/mol.

    """
    time_in, temp = inputs
    my_k = arrhen * np.exp(-e_activ / (ID_GAS_CONST * temp))
    return (my_k * time_in)/(1 + my_k * time_in)


def third_order(inputs, arrhen, e_activ):
    """
    Return third order fit of data.

    Parameters
    ----------
    inputs : Input parameters for fit (time (seconds), Temp. (Kelvin))
    arrhen : Arrhenius prefactor.
    e_activ : Activation energy in kcal/mol.

    """
    time_in, temp = inputs
    my_k = arrhen * np.exp(-1 * e_activ / (ID_GAS_CONST * temp))
    return 1 - (1 / np.sqrt(np.abs(1 + 2 * my_k * time_in)))


def avrami_erofeyev_1(inputs, arrhen, e_activ):
    """
    Return Avrami Erofeyev 1 fit of data, alpha = 1- exp(-(k * t) ^ 2).

    Parameters
    ----------
    inputs : Input parameters for fit (time (seconds), Temp. (Kelvin))
    arrhen : Arrhenius prefactor.
    e_activ : Activation energy in kcal/mol.

    """
    time_in, temp = inputs
    my_k = arrhen * np.exp(-1 * e_activ / (ID_GAS_CONST * temp))
    return 1 - np.exp(-(my_k * time_in) ** 2)


def avrami_erofeyev_2(inputs, arrhen, e_activ):
    """
    Return Avrami Erofeyev 2 fit of data, alpha = 1- exp(-(k * t) ^ 3).

    Parameters
    ----------
    inputs : Input parameters for fit (time (seconds), Temp. (Kelvin))
    arrhen : Arrhenius prefactor.
    e_activ : Activation energy in kcal/mol.

    """
    time_in, temp = inputs
    my_k = arrhen * np.exp(-1 * e_activ / (ID_GAS_CONST * temp))
    return 1 - np.exp(-(my_k * time_in) ** 3)


def avrami_erofeyev_3(inputs, arrhen, e_activ):
    """
    Return Avrami Erofeyev 3 fit of data, alpha = 1- exp(-(k * t) ^ 4).

    Parameters
    ----------
    inputs : Input parameters for fit (time (seconds), Temp. (Kelvin))
    arrhen : Arrhenius prefactor.
    e_activ : Activation energy in kcal/mol.

    """
    time_in, temp = inputs
    my_k = arrhen * np.exp(-1 * e_activ / (ID_GAS_CONST * temp))
    return 1 - np.exp(-(my_k * time_in) ** 4)


def avrami_erofeyev_4(inputs, arrhen, e_activ):
    """
    Return Avrami Erofeyev 4 fit of data, alpha = 1- exp(-(k * t) ^ (2 / 3)).

    Parameters
    ----------
    inputs : Input parameters for fit (time (seconds), Temp. (Kelvin))
    arrhen : Arrhenius prefactor.
    e_activ : Activation energy in kcal/mol.

    """
    time_in, temp = inputs
    my_k = arrhen * np.exp(-e_activ / (ID_GAS_CONST * temp))
    return 1 - np.exp(-(my_k * time_in) ** (2 / 3))


def twothird_power_law(inputs, arrhen, e_activ):
    """
    Return (2/3) power law fit of data.

    Parameters
    ----------
    inputs : Input parameters for fit (time (seconds), Temp. (Kelvin))
    arrhen : Arrhenius prefactor.
    e_activ : Activation energy in kcal/mol.

    """
    time_in, temp = inputs
    my_k = arrhen * np.exp(-1 * e_activ / (ID_GAS_CONST * temp))
    return (my_k * time_in) ** (2 / 3)


def quadratic_power_law(inputs, arrhen, e_activ):
    """
    Return quadratic power law fit of data.

    Parameters
    ----------
    inputs : Input parameters for fit (time (seconds), Temp. (Kelvin))
    arrhen : Arrhenius prefactor.
    e_activ : Activation energy in kcal/mol.

    """
    time_in, temp = inputs
    my_k = arrhen * np.exp(-1 * e_activ / (ID_GAS_CONST * temp))
    return (my_k * time_in) ** (2)


def cubic_power_law(inputs, arrhen, e_activ):
    """
    Return cubic power law fit of data.

    Parameters
    ----------
    inputs : Input parameters for fit (time (seconds), Temp. (Kelvin))
    arrhen : Arrhenius prefactor.
    e_activ : Activation energy in kcal/mol.

    """
    time_in, temp = inputs
    my_k = arrhen * np.exp(-1 * e_activ / (ID_GAS_CONST * temp))
    return (2 * my_k * time_in) ** (3 / 2)


def quartic_power_law(inputs, arrhen, e_activ):
    """
    Return quartic power law fit of data.

    Parameters
    ----------
    inputs : Input parameters for fit (time (seconds), Temp. (Kelvin))
    arrhen : Arrhenius prefactor.
    e_activ : Activation energy in kcal/mol.

    """
    time_in, temp = inputs
    my_k = arrhen * np.exp(-1 * e_activ / (ID_GAS_CONST * temp))
    return (3 * my_k * time_in) ** (4 / 3)


def contracting_area(inputs, arrhen, e_activ):
    """
    Return contracting area fit of data.

    Parameters
    ----------
    inputs : Input parameters for fit (time (seconds), Temp. (Kelvin))
    arrhen : Arrhenius prefactor.
    e_activ : Activation energy in kcal/mol.

    """
    time_in, temp = inputs
    my_k = arrhen * np.exp(-1 * e_activ / (ID_GAS_CONST * temp))
    return 2 * my_k * time_in - ((my_k * time_in) ** 2)


def contracting_volume(inputs, arrhen, e_activ):
    """
    Return contracting volume fit of data.

    Parameters
    ----------
    inputs : Input parameters for fit (time (seconds), Temp. (Kelvin))
    arrhen : Arrhenius prefactor.
    e_activ : Activation energy in kcal/mol.

    """
    time_in, temp = inputs
    my_k = arrhen * np.exp(-1 * e_activ / (ID_GAS_CONST * temp))
    return 1 - (np.abs(1 - 2 * my_k * time_in) ** (3 / 2))


def onedimension_diffusion(inputs, arrhen, e_activ):
    """
    Return 1D diffusion fit of data.

    Parameters
    ----------
    inputs : Input parameters for fit (time (seconds), Temp. (Kelvin))
    arrhen : Arrhenius prefactor.
    e_activ : Activation energy in kcal/mol.

    """
    time_in, temp = inputs
    my_k = arrhen * np.exp(-1 * e_activ / (ID_GAS_CONST * temp))
    return np.sqrt(my_k * time_in)


def plot_model(selected_model_name):
    """Produce plot of fit and print fit parameters."""
    plt.figure(figsize=(12, 9))
    plt.clf()
    plt.scatter(time_data, plotalpha)
    plt.plot(mytime, mysol, color='b', linewidth=1.5, linestyle='-')
    plt.axis([0, np.amax(mytime) * 1.1, 0, np.amax(mysol) + 0.025])
    plt.xlabel('t (seconds)')
    plt.ylabel('Reaction progress fraction')
    plt.grid(True)
    plt.title(selected_model_name + " Curve Fit")
    plt.savefig(selected_model_name + "_CurveFit.png", dpi=400)
    plt.draw()
    plt.show()
    print("A = ", params[0], "\nEa = ", params[1], " kcal/mol")
    return print("Var(A) = ", params_covariance[0, 0],
                 "\nVar(Ea) = ", params_covariance[1, 1],
                 "\nCov(A, Ea) = Cov (Ea, A) = ", params_covariance[0, 1])


ID_GAS_CONST = 1.985877534 * (10 ** (-3))
direc = input('Enter the full directory for your data CSV file: ')
with open(direc) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    NUM_POINTS = 0
    for row in csv_reader:
        if NUM_POINTS == 0:
            NUM_POINTS += 1
        else:
            NUM_POINTS += 1
my_csv_data = pd.read_csv(direc)
time_data = []
for i in range(0, NUM_POINTS - 1):
    time_data.append(my_csv_data['time_data'][i])
time_data = np.array(time_data)
mytime = np.linspace(0., np.amax(time_data), NUM_POINTS)
initmol_data = []
for i in range(0, NUM_POINTS - 1):
    initmol_data.append(my_csv_data['initmol_data'][i])
initmol_data = np.array(initmol_data)
temp_data = []
for i in range(0, NUM_POINTS - 1):
    temp_data.append(my_csv_data['temp_data'][i])
temp_data = np.array(temp_data)
avgtemp = np.average(temp_data)
plottemp_data = np.empty_like(temp_data)
plottemp_data[:] = avgtemp
finalmol_data = []
for i in range(0, NUM_POINTS - 1):
    finalmol_data.append(my_csv_data['finalmol_data'][i])
finalmol_data = np.array(finalmol_data)
alpha = finalmol_data / initmol_data
plotalpha = np.empty_like(alpha)
mysol = np.empty_like(alpha)
myy = np.empty_like(finalmol_data)
USERCONTINUE = True
print('')

while USERCONTINUE:
    cmd = int(input("Now input the number of the model you want to use with "
                    + "your data, e.g., enter 1 for zero-order: "
                    + "\n1 zero-order, 2 first-order, 3 second-order, "
                    + "4 third-order, 5 Avrami-Erofeyev 1,"
                    + "\n6 Avrami-Erofeyev 2, 7 Avrami-Erofeyev 3, "
                    + "8 Avrami-Erofeyev 4, 9 2/3 power law, "
                    + "10 quadratic power law, \n11 cubic power law, "
                    + "12 quartic power law, 13 contracting area, "
                    + "14 contracting volume, 15 1-D diffusion: "))
    if cmd == 1:
        MODEL_NAME = 'Zero-order'
        params, params_covariance = optimize.curve_fit(zero_order,
                                                       (time_data, temp_data),
                                                       alpha, p0=[100, 10],
                                                       maxfev=100000)
        for i in range(0, NUM_POINTS - 1):
            plotalpha[i] = alpha[i] * np.exp((params[1] / ID_GAS_CONST)
                                             * ((1 / temp_data[i])
                                                - (1 / avgtemp)))
        mysol = zero_order((mytime, avgtemp), params[0], params[1])
        plot_model(MODEL_NAME)
    elif cmd == 2:
        MODEL_NAME = 'First-order'
        params, params_covariance = optimize.curve_fit(first_order,
                                                       (time_data, temp_data),
                                                       alpha, p0=[100, 10],
                                                       maxfev=100000)
        for i in range(0, NUM_POINTS - 1):
            plotalpha[i] = alpha[i] * np.exp((params[1] / ID_GAS_CONST)
                                             * ((1 / temp_data[i])
                                                - (1 / avgtemp)))
        mysol = first_order((mytime, avgtemp), params[0], params[1])
        plot_model(MODEL_NAME)
    elif cmd == 3:
        MODEL_NAME = 'Second-order'
        params, params_covariance = optimize.curve_fit(second_order,
                                                       (time_data, temp_data),
                                                       alpha, p0=[100, 10],
                                                       maxfev=100000)
        for i in range(0, NUM_POINTS - 1):
            plotalpha[i] = alpha[i] * np.exp((params[1] / ID_GAS_CONST)
                                             * ((1 / temp_data[i])
                                                - (1 / avgtemp)))
        mysol = second_order((mytime, avgtemp), params[0], params[1])
        plot_model(MODEL_NAME)
    elif cmd == 4:
        MODEL_NAME = 'Third-order'
        params, params_covariance = optimize.curve_fit(third_order,
                                                       (time_data, temp_data),
                                                       alpha, p0=[100, 10],
                                                       maxfev=100000)
        for i in range(0, NUM_POINTS - 1):
            plotalpha[i] = alpha[i] * np.exp((params[1] / ID_GAS_CONST)
                                             * ((1 / temp_data[i])
                                                - (1 / avgtemp)))
        mysol = third_order((mytime, avgtemp), params[0], params[1])
        plot_model(MODEL_NAME)
    elif cmd == 5:
        MODEL_NAME = 'Avrami-Erofeyev 1'
        params, params_covariance = optimize.curve_fit(avrami_erofeyev_1,
                                                       (time_data, temp_data),
                                                       alpha, p0=[100, 10],
                                                       maxfev=100000)
        for i in range(0, NUM_POINTS - 1):
            plotalpha[i] = alpha[i] * np.exp((params[1] / ID_GAS_CONST)
                                             * ((1 / temp_data[i])
                                                - (1 / avgtemp)))
        mysol = avrami_erofeyev_1((mytime, avgtemp), params[0], params[1])
        plot_model(MODEL_NAME)
    elif cmd == 6:
        MODEL_NAME = 'Avrami-Erofeyev 2'
        params, params_covariance = optimize.curve_fit(avrami_erofeyev_2,
                                                       (time_data, temp_data),
                                                       alpha, p0=[100, 10],
                                                       maxfev=100000)
        for i in range(0, NUM_POINTS - 1):
            plotalpha[i] = alpha[i] * np.exp((params[1] / ID_GAS_CONST)
                                             * ((1 / temp_data[i])
                                                - (1 / avgtemp)))
        mysol = avrami_erofeyev_2((mytime, avgtemp), params[0], params[1])
        plot_model(MODEL_NAME)
    elif cmd == 7:
        MODEL_NAME = 'Avrami-Erofeyev 3'
        params, params_covariance = optimize.curve_fit(avrami_erofeyev_3,
                                                       (time_data, temp_data),
                                                       alpha, p0=[100, 10],
                                                       maxfev=100000)
        for i in range(0, NUM_POINTS - 1):
            plotalpha[i] = alpha[i] * np.exp((params[1] / ID_GAS_CONST)
                                             * ((1 / temp_data[i])
                                                - (1 / avgtemp)))
        mysol = avrami_erofeyev_3((mytime, avgtemp), params[0], params[1])
        plot_model(MODEL_NAME)
    elif cmd == 8:
        MODEL_NAME = 'Avrami-Erofeyev 4'
        params, params_covariance = optimize.curve_fit(avrami_erofeyev_4,
                                                       (time_data, temp_data),
                                                       alpha, p0=[100, 10],
                                                       maxfev=100000)
        for i in range(0, NUM_POINTS - 1):
            plotalpha[i] = alpha[i] * np.exp((params[1] / ID_GAS_CONST)
                                             * ((1 / temp_data[i])
                                                - (1 / avgtemp)))
        mysol = avrami_erofeyev_4((mytime, avgtemp), params[0], params[1])
        plot_model(MODEL_NAME)
    elif cmd == 9:
        MODEL_NAME = 'Two-third power law'
        params, params_covariance = optimize.curve_fit(twothird_power_law,
                                                       (time_data, temp_data),
                                                       alpha, p0=[100, 10],
                                                       maxfev=100000)
        for i in range(0, NUM_POINTS - 1):
            plotalpha[i] = alpha[i] * np.exp((params[1] / ID_GAS_CONST)
                                             * ((1 / temp_data[i])
                                                - (1 / avgtemp)))
        mysol = twothird_power_law((mytime, avgtemp), params[0], params[1])
        plot_model(MODEL_NAME)
    elif cmd == 10:
        MODEL_NAME = 'Quadratic power law'
        params, params_covariance = optimize.curve_fit(quadratic_power_law,
                                                       (time_data, temp_data),
                                                       alpha, p0=[100, 10],
                                                       maxfev=100000)
        for i in range(0, NUM_POINTS - 1):
            plotalpha[i] = alpha[i] * np.exp((params[1] / ID_GAS_CONST)
                                             * ((1 / temp_data[i])
                                                - (1 / avgtemp)))
        mysol = quadratic_power_law((mytime, avgtemp), params[0], params[1])
        plot_model(MODEL_NAME)
    elif cmd == 11:
        MODEL_NAME = 'Cubic power law'
        params, params_covariance = optimize.curve_fit(cubic_power_law,
                                                       (time_data, temp_data),
                                                       alpha, p0=[100, 10],
                                                       maxfev=100000)
        for i in range(0, NUM_POINTS - 1):
            plotalpha[i] = alpha[i] * np.exp((params[1] / ID_GAS_CONST)
                                             * ((1 / temp_data[i])
                                                - (1 / avgtemp)))
        mysol = cubic_power_law((mytime, avgtemp), params[0], params[1])
        plot_model(MODEL_NAME)
    elif cmd == 12:
        MODEL_NAME = 'Quartic power law'
        params, params_covariance = optimize.curve_fit(quartic_power_law,
                                                       (time_data, temp_data),
                                                       alpha, p0=[100, 10],
                                                       maxfev=100000)
        for i in range(0, NUM_POINTS - 1):
            plotalpha[i] = alpha[i] * np.exp((params[1] / ID_GAS_CONST)
                                             * ((1 / temp_data[i])
                                                - (1 / avgtemp)))
        mysol = quartic_power_law((mytime, avgtemp), params[0], params[1])
        plot_model(MODEL_NAME)
    elif cmd == 13:
        MODEL_NAME = 'Contracting area'
        params, params_covariance = optimize.curve_fit(contracting_area,
                                                       (time_data, temp_data),
                                                       alpha, p0=[100, 10],
                                                       maxfev=100000)
        for i in range(0, NUM_POINTS - 1):
            plotalpha[i] = alpha[i] * np.exp((params[1] / ID_GAS_CONST)
                                             * ((1 / temp_data[i])
                                                - (1 / avgtemp)))
        mysol = contracting_area((mytime, avgtemp), params[0], params[1])
        plot_model(MODEL_NAME)
    elif cmd == 14:
        MODEL_NAME = 'Contracting volume'
        params, params_covariance = optimize.curve_fit(contracting_volume,
                                                       (time_data, temp_data),
                                                       alpha, p0=[100, 10],
                                                       maxfev=100000)
        for i in range(0, NUM_POINTS - 1):
            plotalpha[i] = alpha[i] * np.exp((params[1] / ID_GAS_CONST)
                                             * ((1 / temp_data[i])
                                                - (1 / avgtemp)))
        mysol = contracting_volume((mytime, avgtemp), params[0], params[1])
        plot_model(MODEL_NAME)
    elif cmd == 15:
        MODEL_NAME = '1D diffusion'
        params, params_covariance = optimize.curve_fit(onedimension_diffusion,
                                                       (time_data, temp_data),
                                                       alpha, p0=[100, 10],
                                                       maxfev=100000)
        for i in range(0, NUM_POINTS - 1):
            plotalpha[i] = alpha[i] * np.exp((params[1] / ID_GAS_CONST)
                                             * ((1 / temp_data[i])
                                                - (1 / avgtemp)))
        mysol = onedimension_diffusion((mytime, avgtemp), params[0], params[1])
        plot_model(MODEL_NAME)
    else:
        print("\nInvalid entry.\n")
    print("\nYour plot should be saved in your current working directory.\n")
    checkcontinue = input("Would you like to try again with the same data "
                          + "or quit? Enter Y to continue with the same "
                          + "data or N to stop, or enter C to change the "
                          + "input file: ")
    if checkcontinue in ('Y', 'y'):
        USERCONTINUE = True
        print('')
    elif checkcontinue in ('N', 'n'):
        USERCONTINUE = False
    elif checkcontinue in ('C', 'c'):
        direc = input('Enter the full directory for your data CSV file: ')
        with open(direc) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            NUM_POINTS = 0
            for row in csv_reader:
                if NUM_POINTS == 0:
                    NUM_POINTS += 1
                else:
                    NUM_POINTS += 1
        my_csv_data = pd.read_csv(direc)
        time_data = []
        for i in range(0, NUM_POINTS - 1):
            time_data.append(my_csv_data['time_data'][i])
        time_data = np.array(time_data)
        mytime = np.linspace(0., np.amax(time_data), NUM_POINTS)
        initmol_data = []
        for i in range(0, NUM_POINTS - 1):
            initmol_data.append(my_csv_data['initmol_data'][i])
        initmol_data = np.array(initmol_data)
        temp_data = []
        for i in range(0, NUM_POINTS - 1):
            temp_data.append(my_csv_data['temp_data'][i])
        temp_data = np.array(temp_data)
        avgtemp = np.average(temp_data)
        plottemp_data = np.empty_like(temp_data)
        plottemp_data[:] = avgtemp
        finalmol_data = []
        for i in range(0, NUM_POINTS - 1):
            finalmol_data.append(my_csv_data['finalmol_data'][i])
        finalmol_data = np.array(finalmol_data)
        alpha = finalmol_data / initmol_data
        plotalpha = np.empty_like(alpha)
        mysol = np.empty_like(alpha)
        myy = np.empty_like(finalmol_data)
    else:
        USERCONTINUE = True
        print("Invalid user input. Returning to model selection screen with "
              + "same data.\n")
