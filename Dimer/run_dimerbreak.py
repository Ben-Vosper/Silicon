from Dimer.dimerbreak import DimerBreak
import matplotlib.pyplot as plt
import json

results_temp = "cold_wierddrift.json"

v = 0.05
v_max = 1
v_inc = 0.005

# Setting temperature. Using 1eV/K = 11600kB

T = 0.0000001
sum_of_modes = T/11600 * 2
mode_fraction = 0.5                # f_opt/sum_of_modes

f_opt = sum_of_modes*mode_fraction
f_aco = sum_of_modes - f_opt

f_timestep = 0.002
f_stiffness = 1.2

run_info = (v, v_max, v_inc)

nconfig = 1

v_list = []
w_means = []
w_errbar_size = []

results = {"params": (v, v_max, v_inc, f_timestep, f_opt, f_aco, f_stiffness, nconfig)}

total_breaks = ((v_max - v)/v_inc)*nconfig
breaks = 0

while v < v_max:

    t = DimerBreak(v, f_timestep, nconfig, f_opt, f_aco, f_stiffness, results_temp.replace(".json", ""),  run_info)
    w_mean, w_dev = t.run()

    results[v] = (w_mean, w_dev)

    w_means.append(w_mean)
    w_errbar_size.append(w_dev/(nconfig**0.5))

    v_list.append(v)
    v += v_inc

    breaks += 1

    with open(results_temp, 'w') as q:
        json.dump(results, q, indent=2)

T_approx = round(((f_opt + f_aco)/2) * 11600)

plt.errorbar(v_list, w_means, yerr=w_errbar_size, ls="none", marker="+")
plt.xlabel("$\\frac{v}{v_s}$", size=24)
plt.ylabel("$\\bar W$ / $\epsilon$", size=18)
param_string = "Approximate Temperature = " + str(T_approx) + "K" + "\nopt = " + str(format(f_opt, ".3g")) +\
                   "\naco = " + str(format(f_aco, ".3g")) + "\nStiffness = " + str(f_stiffness)
plt.title(param_string, loc="left", fontsize=12)
plt.show()
