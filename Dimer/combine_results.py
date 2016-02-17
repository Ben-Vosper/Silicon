import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import json, colorsys, numpy, statistics


def combine_results(filename_list):

    patches = []
    custom_names = ["1.2", "1.1"]

    colours = []
    colour_offset = 0.1
    n = len(filename_list)
    for q in range(n):
        z = colour_offset + q/n
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
            T_approx = round(((f_opt + f_aco)/2) * 26682.1)

        v_vals = list(results.keys())
        v_vals.remove("params")

        for v in v_vals:
            v_list.append(v)
            w_means.append(results[v][0])
            w_errbar_size.append(results[v][1]/(results[v][2]**0.5))

        patches.append(mpatches.Patch(color=colours[i], label=custom_names[i]))

        plt.errorbar(v_list, w_means, yerr=w_errbar_size, ls="none", marker="+", color=colours[i])

    plt.xlabel("$\\frac{v}{v_s}$", size=24)
    plt.ylabel("$\\bar W$ / $\epsilon$", size=18)
    # param_string = "Approximate Temperature = " + str(T_approx) + "K" + "\nopt = " + str(format(f_opt, ".3g")) +\
    #                "\naco = " + str(format(f_aco, ".3g")) + "\nStiffness = " + str(f_stiffness)
    param_string = "Approximate Temperature = " + str(T_approx) + "K" + "\nopt = " + str(format(f_opt, ".3g")) +\
               "\naco = " + str(format(f_aco, ".3g"))
    plt.title(param_string, loc="left", fontsize=12)
    plt.legend(handles=patches, title="Stiffness", loc="upper left", fontsize=12)
    plt.show()


def combine_runs(filename_list):

    unique_vs = []

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
            T_approx = round(((f_opt + f_aco)/2) * 26682.1)

        v_vals = list(results.keys())
        v_vals.remove("params")

        for v in v_vals:
            v_list.append(v)
            w_means.append(results[v][0])
            w_errbar_size.append(results[v][1]/(results[v][2]**0.5))

        unique_vs += (numpy.unique(v_list).tolist())

    unique_vs = (numpy.unique(unique_vs)).tolist()
    data = numpy.zeros((len(filename_list), len(unique_vs), 4), float)

    for name in filename_list:

        i = filename_list.index(name)

        with open(name, 'r') as q:
            results = json.load(q)

        v_vals = list(results.keys())
        v_vals.remove("params")

        for v in v_vals:
            z = unique_vs.index(v)
            data[i, z, 0] = v
            data[i, z, 1] = results[v][0]
            data[i, z, 2] = results[v][1]
            data[i, z, 3] = results[v][2]

    v_list = []
    w_means = []
    w_errbar_size = []

    for z in range(len(unique_vs)):
        v_list.append(unique_vs[z])

        mean_to_combine = []
        sd_to_combine = []
        nc = []
        for i in range(len(filename_list)):
            if data[i, z, 3] != 0:
                mean_to_combine.append(data[i, z, 1])
                sd_to_combine.append(data[i, z, 2])
                nc.append(data[i, z, 3])

        w_means.append(statistics.mean(mean_to_combine))

        while len(sd_to_combine) > 1:
            s1 = sd_to_combine[0]**2
            s2 = sd_to_combine[1]**2

            m1 = mean_to_combine[0]
            m2 = mean_to_combine[1]

            n1 = nc[0]
            n2 = nc[1]

            new_sd = (n1**2)*s1 + (n2**2)*s2 - n1*s2 - n2*s1 - n1*s1 - n2*s2 + (n1*n2*s1) + (n1*n2*s2) + (n1*n2)*((m1-m2)**2)
            new_sd /= (n1 + n2 - 1)*(n1 + n2)

            new_sd **= 0.5

            sd_to_combine[0] = new_sd
            nc[0] += nc[1]
            mean_to_combine[0] = (m1 + m2)/2

            sd_to_combine.pop(1)
            nc.pop(1)
            mean_to_combine.pop(1)

        w_errbar_size.append(sd_to_combine[0]/(nc[0]**0.5))

    plt.errorbar(v_list, w_means, yerr=w_errbar_size, ls="none", marker="+")
    plt.show()






#
# combine_results(["Old\\300_triple_f0.json", "Old\\300_triple_f0.1.json", "Old\\300_triple_f0.2.json", "Old\\300_triple_f0.3.json",
#                  "Old\\300_triple_f0.4.json", "Old\\300_triple_f0.5.json", "Old\\300_triple_f0.6.json",
#                  "Old\\300_triple_f0.7.json", "Old\\300_triple_f0.8.json", "Old\\300_triple_f0.9.json",
#                  "Old\\300_triple_f1.json"])

# combine_results(["cold_s1.2.json", "cold_s1.1.json"])

combine_runs(["t1.json", "t2.json"])

#"E:\\Ben Vosper\\My Documents\\Silicon\\Dimer_3\\300_0.1_test.json"
