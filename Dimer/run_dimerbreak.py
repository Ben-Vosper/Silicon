from Dimer.dimerbreak import DimerBreak
import matplotlib.pyplot as plt
import json, statistics

results_temp = "600_s1.1_aco.json"

v_start = 0.005
v = v_start
v_max = 0.2
v_inc = 0.001

# Setting temperature. Using E = kB/T, E = 2.3eV, kB = 8.62e-5eVK-1

e_si_si = 2.3
kB = 8.62e-5

T = 600

e_av = (kB * T)/e_si_si
optical_fraction = 0
acoustic_fraction = 1

f_opt = e_av * optical_fraction
f_aco = e_av * acoustic_fraction

f_timestep = 0.002
f_stiffness = 1.1

run_info = (v, v_max, v_inc)

nconfig = 100

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
    if v == v_start:
        print(" -!-!-!- Mode = ", t.i_potential, " -!-!-!-")
        input("Continue?")
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

T_approx = round(max([f_aco, f_opt]) / 3.747e-5)

plt.errorbar(v_list, w_means, yerr=w_errbar_size, ls="none", marker="+")#, color=(0, 0.5, 0.5))
# plt.errorbar(v_list, eks, yerr=ek_ebs, ls="none", marker="+", color=(0.5, 0, 0.5))
# plt.errorbar(v_list, eps, yerr=ep_ebs, ls="none", marker="+", color=(0.5, 0.5, 0))
plt.xlabel("$\\frac{v}{v_s}$", size=24)
plt.ylabel("$\\bar W$ / $\epsilon$", size=18)
param_string = "Approximate Temperature = " + str(T_approx) + "K" + "\nopt = " + str(format(f_opt, ".3g")) +\
                   "\naco = " + str(format(f_aco, ".3g")) + "\nStiffness = " + str(f_stiffness)
plt.title(param_string, loc="left", fontsize=12)
plt.show()
