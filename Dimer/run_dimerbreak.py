from Dimer.dimerbreak import DimerBreak
import matplotlib.pyplot as plt
import json

results_temp = "300_f0.4_new.json"

v = 0.01
v_max = 1
v_inc = 0.005

# Setting temperature. Using 1eV = 11600kB

T = 300
sum_of_modes = T/11600 * 2
mode_fraction = 0.4                # f_opt/sum_of_modes

f_opt = sum_of_modes*mode_fraction
f_aco = sum_of_modes - f_opt

f_timestep = 0.002
f_stiffness = 1.1

run_info = (v, v_max, v_inc)

nconfig = 100

v_list = []
w_means = []
w_errbar_size = []

results = {}
results["params"] = (v, v_max, v_inc, f_timestep, f_opt, f_aco, f_stiffness, nconfig)

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
    #print(str(round((breaks/total_breaks)*100, 2)) + "%   Complete.")

    with open(results_temp, 'w') as q:
        json.dump(results, q, indent=2)

plt.errorbar(v_list, w_means, yerr=w_errbar_size, ls="none", marker="+")
plt.xlabel("$\\frac{v}{v_s}$", size=24)
plt.ylabel("$\\bar W$ / $\epsilon$", size=18)
plt.show()