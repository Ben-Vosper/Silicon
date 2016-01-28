import matplotlib.pyplot as plt
import json, colorsys


def combine_results(filename_list):

    colours = []
    n = len(filename_list)
    for q in range(n):
        z = 0.00001 + q/n
        c = colorsys.hls_to_rgb(z, 0.5, 1)
        colours.append(c)

    for name in filename_list:

        i = filename_list.index(name)

        v_list = []
        w_means = []
        w_errbar_size = []

        with open(name, 'r') as q:
            results = json.load(q)

        if i == 0:
            nconfig = results["params"][7]
            f_stiffness = results["params"][6]
            f_aco = results["params"][5]
            f_opt = results["params"][4]
            T_approx = round(((f_opt + f_aco)/2) * 11600)

        v_vals = list(results.keys())
        v_vals.remove("params")

        for v in v_vals:
            v_list.append(v)
            w_means.append(results[v][0])
            w_errbar_size.append(results[v][1]/(nconfig**0.5))

        plt.errorbar(v_list, w_means, yerr=w_errbar_size, ls="none", marker="+", color=colours[i])

    plt.xlabel("$\\frac{v}{v_s}$", size=24)
    plt.ylabel("$\\bar W$ / $\epsilon$", size=18)
    param_string = "Approximate Temperature = " + str(T_approx) + "K" + "\nopt = " + str(format(f_opt, ".3g")) +\
                   "\naco = " + str(format(f_aco, ".3g")) + "\nStiffness = " + str(f_stiffness)
    plt.title(param_string, loc="left", fontsize=12)
    plt.show()

combine_results(["Results\\300_f0.1.json", "Results\\300_f0.2_new.json", "Results\\300_f0.3_new.json",
                 "Results\\300_f0.4_new.json", "Results\\300_f0.5_new.json", "Results\\300_f0.6_new.json",
                 "Results\\300_f0.7_new.json", "Results\\300_f0.8_new.json", "Results\\300_f0.9_new.json",
                 "Results\\300_f1_new.json"])

#"E:\\Ben Vosper\\My Documents\\Silicon\\Dimer_3\\300_0.1_test.json"
