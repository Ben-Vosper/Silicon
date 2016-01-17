from Dimer.dimerbreak import DimerBreak
import matplotlib.pyplot as plt
import json, os

main_results = "resume_test.json"
inc_results = "test_0.07.json"

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

        t = DimerBreak(current_v, f_timestep, remaining_runs, f_opt, f_aco, f_stiffness, inc_results, run_info)
        w_mean, w_dev = t.run()

        results[v] = (w_mean, w_dev)

        w_means.append(w_mean)
        w_errbar_size.append(w_dev/(nconfig**0.5))

        v_list.append(v)
        v += v_inc

        breaks += 1
        #print(str(round((breaks/total_breaks)*100, 2)) + "%   Complete.")

        with open(main_results, 'w') as q:
            json.dump(results, q, indent=2)

    v = current_v + v_inc

while v < v_max:

    t = DimerBreak(v, f_timestep, nconfig, f_opt, f_aco, f_stiffness, "test",  run_info)
    w_mean, w_dev = t.run()

    results[v] = (w_mean, w_dev)

    w_means.append(w_mean)
    w_errbar_size.append(w_dev/(nconfig**0.5))

    v_list.append(v)
    v += v_inc

    breaks += 1
    #print(str(round((breaks/total_breaks)*100, 2)) + "%   Complete.")

    with open(main_results, 'w') as q:
        json.dump(results, q, indent=2)


plt.errorbar(v_list, w_means, yerr=w_errbar_size, ls="none", marker="+")
plt.show()

os.remove(inc_results)