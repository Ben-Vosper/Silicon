import matplotlib.pyplot as plt
import json


def read_full_results_old(filename):

    v_list = []
    w_means = []
    w_errbar_size = []

    with open(filename, 'r') as q:
        results = json.load(q)

    v_vals = list(results.keys())
    nconfig = 100

    for v in v_vals:
        v_list.append(v)
        w_means.append(results[v][0])
        w_errbar_size.append(results[v][1]/(nconfig**0.5))

    plt.errorbar(v_list, w_means, yerr=w_errbar_size, ls="none", marker="+")
    plt.xlabel("$\\frac{v}{v_s}$", size=24)
    plt.ylabel("$\\bar W$ / $\epsilon$", size=18)
    plt.show()


def read_full_results(filename):

    v_list = []
    w_means = []
    w_errbar_size = []

    with open(filename, 'r') as q:
        results = json.load(q)

    nconfig = results["params"][7]
    v_vals = list(results.keys())
    v_vals.remove("params")

    for v in v_vals:
        v_list.append(v)
        w_means.append(results[v][0])
        w_errbar_size.append(results[v][1]/(nconfig**0.5))

    plt.errorbar(v_list, w_means, yerr=w_errbar_size, ls="none", marker="+")
    plt.xlabel("$\\frac{v}{v_s}$", size=24)
    plt.ylabel("$\\bar W$ / $\epsilon$", size=18)
    plt.show()

read_full_results("Results\\opt0_aco0.1.json")