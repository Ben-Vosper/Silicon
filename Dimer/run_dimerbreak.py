from Dimer.dimerbreak import DimerBreak
import matplotlib.pyplot as plt
import json, statistics

results_temp = "cold_s1.1.json"

v = 0.001
v_max = 0.3
v_inc = 0.001

# Setting temperature. Using E = kB/T, E = 2.3eV, kB = 8.62e-5eVK-1

p = 2.3 / 8.62e-5

T = 0.00001
sum_of_modes = T/p * 2
mode_fraction = 0.5                # f_opt/sum_of_modes

f_opt = sum_of_modes*mode_fraction
f_aco = sum_of_modes - f_opt

f_timestep = 0.002
f_stiffness = 1.1

run_info = (v, v_max, v_inc)

nconfig = 1

v_list = []
w_means = []
w_errbar_size = []

results = {"params": (v, v_max, v_inc, f_timestep, f_opt, f_aco, f_stiffness, nconfig)}

total_breaks = ((v_max - v)/v_inc)*nconfig
breaks = 0

eks = []
ek_ebs = []

eps = []
ep_ebs = []

while v < v_max:

    t = DimerBreak(v, f_timestep, nconfig, f_opt, f_aco, f_stiffness, results_temp.replace(".json", ""),  run_info)
    works = t.run()

    n_successful_runs = len(works)

    if n_successful_runs == 0:
        print("Stiff bond(s) broke for every config at v = " + str(v))
        v = v_max
    else:
        w_mean = statistics.mean(works)
        w_dev = statistics.pstdev(works)

        results[v] = (w_mean, w_dev, n_successful_runs)

        w_means.append(w_mean)
        w_errbar_size.append(w_dev/(n_successful_runs**0.5))

        v_list.append(v)
        v += v_inc

        breaks += 1

        with open(results_temp, 'w') as q:
            json.dump(results, q, indent=2)

T_approx = round(((f_opt + f_aco)/2) * 11600)

plt.errorbar(v_list, w_means, yerr=w_errbar_size, ls="none", marker="+")#, color=(0, 0.5, 0.5))
# plt.errorbar(v_list, eks, yerr=ek_ebs, ls="none", marker="+", color=(0.5, 0, 0.5))
# plt.errorbar(v_list, eps, yerr=ep_ebs, ls="none", marker="+", color=(0.5, 0.5, 0))
plt.xlabel("$\\frac{v}{v_s}$", size=24)
plt.ylabel("$\\bar W$ / $\epsilon$", size=18)
param_string = "Approximate Temperature = " + str(T_approx) + "K" + "\nopt = " + str(format(f_opt, ".3g")) +\
                   "\naco = " + str(format(f_aco, ".3g")) + "\nStiffness = " + str(f_stiffness)
plt.title(param_string, loc="left", fontsize=12)
plt.show()
