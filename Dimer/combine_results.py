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
#
# combine_results(["Old\\300_triple_f0.json", "Old\\300_triple_f0.1.json", "Old\\300_triple_f0.2.json", "Old\\300_triple_f0.3.json",
#                  "Old\\300_triple_f0.4.json", "Old\\300_triple_f0.5.json", "Old\\300_triple_f0.6.json",
#                  "Old\\300_triple_f0.7.json", "Old\\300_triple_f0.8.json", "Old\\300_triple_f0.9.json",
#                  "Old\\300_triple_f1.json"])

combine_results(["Results\\cold_s1.1.json", "cold_s1.1_vel.json"])

#"E:\\Ben Vosper\\My Documents\\Silicon\\Dimer_3\\300_0.1_test.json"
