from Dimer.dimerbreak import DimerBreak
import matplotlib.pyplot as plt
import json

results_temp = "test.json"

v = 0.01
v_max = 1
v_inc = 0.01

nconfig = 100

v_list = []
w_means = []
w_errbar_size = []

results = {}

total_breaks = ((v_max - v)/v_inc)*nconfig
breaks = 0

while v < v_max:

    t = DimerBreak(v, 0.002, nconfig, 0.0002, 0.0002, 1.2, "test")
    w_mean, w_dev = t.run()

    results[v] = (w_mean, w_dev)

    w_means.append(w_mean)
    w_errbar_size.append(w_dev/(nconfig**0.5))

    v_list.append(v)
    v += v_inc

    breaks += 1
    print(str(round((breaks/total_breaks)*100, 2)) + "   Complete.")

    with open(results_temp, 'w') as q:
        json.dump(results, q, indent=2)


plt.errorbar(v_list, w_means, yerr=w_errbar_size, ls="none", marker="+")
plt.show()