import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import json, colorsys, numpy, statistics


def combine_results(filename_list):

    patches = []
    custom_names = ["1.05", "1.10", "1.20"]

    colours = []
    colour_offset = 0.07
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

    def combine_sds(mean_to_combine, sd_to_combine, nc):
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
        return sd_to_combine[0]/(nc[0]**0.5)

    all_vs = []
    data = []

    for name in filename_list:

        i = filename_list.index(name)
        with open(name, 'r') as q:
            results = json.load(q)

        if i == 0:
            nconfig = results["params"][7]
            f_stiffness = results["params"][6]
            f_aco = results["params"][5]
            f_opt = results["params"][4]
            T_approx = round(((f_opt + f_aco)/2) * 26682.1)

        del results["params"]

        data.append(results)

        for key in list(results.keys()):
            all_vs.append((i, key, round(float(key), 6)))

    v_list = []
    w_means = []
    w_errbar_size = []

    for item in all_vs:
        v_to_check = item[2]
        other_refs = []
        for check_item in all_vs:
            if check_item != item:
                if check_item[2] == v_to_check:
                    other_refs.append((check_item[0], check_item[1]))
                    all_vs.pop(all_vs.index(check_item))

        mean = data[item[0]][item[1]][0]
        sd = data[item[0]][item[1]][1]
        nc = data[item[0]][item[1]][2]
        if not other_refs:
            v_list.append(v_to_check)
            w_means.append(mean)
            w_errbar_size.append(sd/(nc**0.5))
        else:
            temp_means = [mean]
            temp_sds = [sd]
            temp_ncs = [nc]
            for other in other_refs:
                temp_means.append(data[other[0]][other[1]][0])
                temp_sds.append(data[other[0]][other[1]][1])
                temp_ncs.append(data[other[0]][other[1]][2])
            v_list.append(v_to_check)
            w_means.append(statistics.mean(temp_means))
            w_errbar_size.append(combine_sds(temp_means, temp_sds, temp_ncs))

    #print(len(v_list), len(w_means), len(w_errbar_size))

    plt.errorbar(v_list, w_means, yerr=w_errbar_size, ls="none", marker="+")
    plt.xlabel("$\\frac{v}{v_s}$", size=24)
    plt.ylabel("$\\bar W$ / $\epsilon$", size=18)
    param_string = "Approximate Temperature = " + str(T_approx) + "K" + "\nopt = " + str(format(f_opt, ".3g")) +\
               "\naco = " + str(format(f_aco, ".3g"))
    plt.title(param_string, loc="left", fontsize=12)
    plt.show()






#
# combine_results(["Old\\300_triple_f0.json", "Old\\300_triple_f0.1.json", "Old\\300_triple_f0.2.json", "Old\\300_triple_f0.3.json",
#                  "Old\\300_triple_f0.4.json", "Old\\300_triple_f0.5.json", "Old\\300_triple_f0.6.json",
#                  "Old\\300_triple_f0.7.json", "Old\\300_triple_f0.8.json", "Old\\300_triple_f0.9.json",
#                  "Old\\300_triple_f1.json"])

# combine_results(["cold_s1.2.json", "cold_s1.1.json"])


#combine_runs(["Results\\300_f1.json", "300_f1_ext.json"])
combine_results(["cold_s1.05.json", "Results\\cold_s1.1.json", "Results\\cold_s1.2.json"])

#"E:\\Ben Vosper\\My Documents\\Silicon\\Dimer_3\\300_0.1_test.json"
