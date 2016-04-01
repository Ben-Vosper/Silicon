import matplotlib.pyplot as plt
from pylab import *
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
    f_stiffness = results["params"][6]
    f_aco = results["params"][5]
    f_opt = results["params"][4]
    T_approx = round(max([f_aco, f_opt]) / 3.747e-5)
    v_vals = list(results.keys())
    v_vals.remove("params")

    rescaled_v_vals = rescale(v_vals, f_stiffness, "stiffer_lj")
    #rescaled_v_vals = v_vals

    for v in v_vals:
        v_list.append(rescaled_v_vals[v_vals.index(v)])
        w_means.append(results[v][0])
        w_errbar_size.append(results[v][1]/(results[v][2]**0.5))

    print(len(v_list))

    plt.errorbar(v_list, w_means, yerr=w_errbar_size, ls="none", marker="+")
    plt.xlabel("$\\frac{v}{v_R}$", size=24)
    plt.ylabel("$\\bar W$ / $\epsilon$", size=18)
    param_string = "Approximate Temperature = " + str(T_approx) + "K" + "\nopt = " + str(format(f_opt, ".3g")) +\
                   "\naco = " + str(format(f_aco, ".3g")) + "\nStiffness = " + str(f_stiffness)
    plt.title(param_string, loc="left", fontsize=12)
    plt.ylim(-0.05)
    plt.tight_layout()
    plt.show()


def rescale(v_vals, stiffness, mode):

    mode = mode
    sig = 1
    eps = zeros(3)
    eps[0] = 1*stiffness
    eps[1] = 1
    eps[2] = 1*stiffness
    mass = 2
    red_mass = 1
    if mode == "stiffer_lj":
        per_opt = (2.0*pi*sig/(6.0*2.**(1./3)))*(red_mass/(0.5*eps[0]+eps[1]))**0.5
    elif mode == "triple_back":
        per_opt = (2.0*pi*sig/(6.0*2.**(1./3)))*(red_mass/(0.5*eps[0]/3.+eps[1]))**0.5
    speed_of_sound = 6.0*(2.*1/mass)**0.5
    lattice_const = sig*2.0**(1./6.0) * 4 * sin(70.5)
    t_vel = (lattice_const/per_opt)

    true_per_opt = 1/15.56e12
    print(true_per_opt)
    true_lattice_const = 5.43e-10
    true_t_vel = true_lattice_const/true_per_opt

    print(speed_of_sound * (true_t_vel/t_vel))

    # true_speed_of_sound = 8433
    # raleigh_fraction = 0.5
    # raleigh_speed = raleigh_fraction * true_speed_of_sound
    raleigh_speed = 4680
    pull_velocity_ratio = 0.25

    v_list = []
    for v in v_vals:
        true_v = float(v) * speed_of_sound * (true_t_vel/t_vel)
        v_list.append(true_v/raleigh_speed)

    return v_list

#read_full_results("E:\\Ben Vosper\\My Documents\\Silicon\\Dimer_3\\cold_v3_vel_drift.json")
read_full_results("cold_s1.01.json")