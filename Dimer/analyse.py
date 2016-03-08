import matplotlib.pyplot as plt
import json, colorsys
import numpy
from scipy.ndimage import convolve1d


def moving_average(data, window_size):

    # window = []
    # means = []
    #
    # if window_size % 2 == 0:
    #     back = int((window_size / 2) - 1)
    #     forward = int(window_size / 2)
    # else:
    #     back = int((window_size - 1) / 2)
    #     forward = back
    #
    # for z in range(back):
    #     window.append(-(z + 1))
    #     means.append(0)
    # for z in range(forward + 1):
    #     window.append(z)
    #
    # for q in range(back, len(data) - forward):
    #     current = 0
    #     for i in window:
    #         current += data[q + i]
    #     means.append(current/window_size)
    #
    # for z in range(forward):
    #     means.append(0)

    # Or, do it in two lines. Ya dingus.

    kernel = (1/window_size)*numpy.ones(window_size)
    means = convolve1d(data, kernel, mode="nearest")

    return means


def find_minimums(data_x, data_y):
    v = []
    pos = []
    for i in range(1, len(data_y) - 1):
        current = data_y[i]
        prev = data_y[i - 1]
        next = data_y[i + 1]

        if current != 0 and next != 0 and prev !=0:
            if current < prev and current < next:
                v.append(data_x[i])
                pos.append(data_y[i])
    return v, pos


def find_maximums(data_x, data_y):
    v = []
    pos = []
    for i in range(1, len(data_y) - 1):
        current = data_y[i]
        prev = data_y[i - 1]
        next = data_y[i + 1]

        if current != 0 and next != 0 and prev !=0:
            if current > prev and current > next:
                v.append(data_x[i])
                pos.append(data_y[i])
    return v, pos


def calculate_diffs(min_vals, max_vals):
    diff = []
    for i in range(len(min_vals)):
        try:
            min = min_vals[i]
            max = max_vals[i]
            diff.append(max - min)
        except IndexError:
            pass

    return diff


def read_full_results(filename):

    v_list = []
    w_means = []
    w_errbar_size = []

    with open(filename, 'r') as q:
        results = json.load(q)

    nconfig = results["params"][7]
    f_stiffness = results["params"][6]
    f_aco = results["params"][5]
    f_opt = results["params"][4]
    mode_fraction = f_opt/(f_opt + f_aco)
    T_approx = round(max([f_aco, f_opt]) / 1.874e-5)
    v_vals = list(results.keys())
    v_vals.remove("params")
    print(mode_fraction)

    for v in sorted(v_vals):
        v_list.append(v)
        w_means.append(results[v][0])
        w_errbar_size.append(results[v][1]/(results[v][2]**0.5))

    m_av = moving_average(w_means, 5)
    v_mins, min_vals = find_minimums(v_list, m_av)
    v_maxes, max_vals = find_maximums(v_list, m_av)
    diffs = calculate_diffs(min_vals, max_vals)

    diff_markers = []
    diff_markers_v = []

    for i in range(len(diffs)):
        for q in range(1, 10):
            diff_markers_v.append(v_mins[i])
            diff_markers.append(min_vals[i] + diffs[i]*(q/10))

    plt.errorbar(v_list, w_means, yerr=w_errbar_size, ls="none", marker="+")
    plt.plot(v_list, m_av, ls="-", color=(0.5, 0, 0.5))
    plt.xlabel("$\\frac{v}{v_s}$", size=24)
    plt.ylabel("$\\bar W$ / $\epsilon$", size=18)
    param_string = "Approximate Temperature = " + str(T_approx) + "K" + "\nopt = " + str(format(f_opt, ".3g")) +\
                   "\naco = " + str(format(f_aco, ".3g")) + "\nStiffness = " + str(f_stiffness)
    plt.title(param_string, loc="left", fontsize=12)
    plt.show()


def combine_results(filename_list):

    colours = []
    n = len(filename_list)
    for q in range(n):
        z = 0.00001 + q/n
        c = colorsys.hls_to_rgb(z, 0.5, 1)
        colours.append(c)

    low_dif = []
    low_dif_f = []

    for name in filename_list:

        i = filename_list.index(name)

        v_list = []
        w_means = []
        w_errbar_size = []

        with open(name, 'r') as q:
            results = json.load(q)

        nconfig = results["params"][7]
        f_stiffness = results["params"][6]
        f_aco = results["params"][5]
        f_opt = results["params"][4]
        mode_fraction = f_opt/(f_opt + f_aco)
        T_approx = round(max([f_aco, f_opt]) / 1.874e-5)

        v_vals = list(results.keys())
        v_vals.remove("params")

        for v in sorted(v_vals):
            v_list.append(v)
            w_means.append(results[v][0])
            w_errbar_size.append(results[v][1]/(results[v][2]**0.5))

        m_av = moving_average(w_means, 5)
        v_mins, min_vals = find_minimums(v_list, m_av)
        v_maxes, max_vals = find_maximums(v_list, m_av)
        diffs = calculate_diffs(min_vals, max_vals)

        # for g in range(len(v_mins)):
        #     v = float(v_mins[g])
        #     if  v < 0.1:#and v > 0.2:
        #         low_dif.append(diffs[g])
        #         low_dif_f.append(round(mode_fraction - 0.02, 2))
        #
        # diff_markers = []
        # diff_markers_v = []
        #
        # for z in range(len(diffs)):
        #     for q in range(1, 10):
        #         diff_markers_v.append(v_mins[z])
        #         diff_markers.append(min_vals[z] + diffs[z]*(q/10))


        plt.plot(v_list, m_av, ls="-", color=colours[i])
        #
        # plt.plot(diff_markers_v, diff_markers, ls="none", marker="o", color=colours[i])

    #plt.bar(low_dif_f, low_dif, 0.04, color=(0.5, 0, 0.75))

    # plt.xlabel("$\\frac{v}{v_s}$", size=24)
    # plt.xlabel("$\\frac{f_{opt}}{f_{opt} + f_{aco}}$", size=24)
    # plt.ylabel("$\Delta\\bar W$ / $\epsilon$", size=18)
    # plt.suptitle("Depth of resonance at $\\frac{v}{v_s}\\approx$0.05", fontsize=12)
    plt.show()

combine_results(["1200_s1.2_aco.json", "600_s1.2_mix.json", "1200_s1.2_opt.json"])

#read_full_results("E:\\Ben Vosper\\My Documents\\Silicon\\Dimer_3\\300_0.1_test.json")
#read_full_results("1200_s1.2_opt.json")
