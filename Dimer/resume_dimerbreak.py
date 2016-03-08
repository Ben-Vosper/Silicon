from Dimer.dimerbreak import DimerBreak
import matplotlib.pyplot as plt
import json, os
import statistics

main_results = "600_s1.1_opt.json"
inc_results = "600_s1.1_opt_0.018.json"

with open(main_results, 'r') as q:
    results = json.load(q)

v = results["params"][0]
v_max = results["params"][1]
v_inc = results["params"][2]

f_timestep = results["params"][3]
f_opt = results["params"][4]
f_aco = results["params"][5]
f_stiffness = results["params"][6]

run_info = (v, v_max, v_inc)

nconfig = results["params"][7]

v_list = []
w_means = []
w_errbar_size = []

total_breaks = ((v_max - v)/v_inc)*nconfig
breaks = 0

v_vals = list(results.keys())
v_vals.remove("params")
v_vals_f = []
for item in v_vals:
    v_vals_f.append(float(item))

current_v = max(v_vals_f) + v_inc
v = current_v

if inc_results:
    with open(inc_results, 'r') as q:
        inc = json.load(q)
        remaining_runs = nconfig - len(inc)

    t = DimerBreak(current_v, f_timestep, remaining_runs, f_opt, f_aco, f_stiffness, main_results.replace(".json", ""), run_info)
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

        with open(main_results, 'w') as q:
            json.dump(results, q, indent=2)

    v = current_v + v_inc

while v < v_max:

    t = DimerBreak(v, f_timestep, nconfig, f_opt, f_aco, f_stiffness, main_results.replace(".json", ""),  run_info)
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

        with open(main_results, 'w') as q:
            json.dump(results, q, indent=2)


plt.errorbar(v_list, w_means, yerr=w_errbar_size, ls="none", marker="+")
plt.show()

os.remove(inc_results)