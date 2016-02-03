import numpy as np
import matplotlib.image as mpimg
import os, json, statistics, time
from pylab import *

class DimerBreak:

    def __init__(self, f_velocity, f_timestep, nconfig, f_optical, f_acoustic, stiffness_ratio, results_temp, run_info):

        # (1) opening velocity will be a fraction frac of the speed of sound in LJ 1D Lattice
        self.frac_vel = f_velocity
        print("velocity:  ", self.frac_vel)
        #
        # (2) time step will be a fraction of the period of a LJ dimer
        self.frac_tim = f_timestep
        print("time step: ", self.frac_tim)
        #
        # (3) number of configurations to average upon
        self.nconfig = nconfig
        print("aver. on : ", self.nconfig)
        #
        # (4) energy of the optical and acoustic modes as a fraction of the LJ binding energy
        self.frac_opt = f_optical
        self.frac_aco = f_acoustic
        print("phon. ens. ", self.frac_opt, self.frac_aco)
        #
        # (5) stiffnes ratio: energy of the back_bonds relative to the central bond
        self.frac_stif = stiffness_ratio
        print("stiffnes: ", self.frac_stif)
        #
        # (6) number of LJ sigmas the pulling will be taken to (must be enough to fully break the bond)
        #     the number of MD steps will be:   nstep = int(n_sigmas*sig/(vel*dt))
        self.n_sigmas = 7.0
        print("brk dist: ", self.n_sigmas)

        self.i_potential = "stiffer_lj"
        #self.i_potential = "triple_back"

        self.results_temp = results_temp + "_" + str(round(self.frac_vel, 4)) + ".json"
        self.run_info = run_info
        self.total_runs = ((run_info[1] - run_info[0])/run_info[2])*nconfig
        self.current_runs = ((self.frac_vel/run_info[2])*nconfig) - nconfig

    def en_and_for(self, xdum, nat, sig, eps):

        if self.i_potential == "stiffer_lj":
            en, fc = self.vdw_stiffer(xdum, nat, sig, eps)

        elif self.i_potential == "triple_back":
            en, fc = self.vdw_triple(xdum, nat, sig, eps)

        return en, fc

    def vdw_stiffer(self, xdum, nat, sig, eps):
        # returns the energy and forces for the single stiffer backbond case

        en = 0.0
        fc = zeros(nat)

        ene=zeros(3)
        di1   = zeros(3)
        di6   = zeros(3)
        di7   = zeros(3)
        di13  = zeros(3)
        di12  = zeros(3)

        di1 = (roll(xdum,-1)-xdum)[:3]
        di1 =  sig/di1
        di6 =  di1**6
        di7 =  di6*di1

        di12 = di6*di6
        di13 = di12*di1

        ene = 4.*eps*(di12-di6)
        fc[:3] = 4.*eps*(-12.*di13+6.*di7)/sig

        fc = fc - roll(fc,1)
        en = sum(ene)

    #    (amazingly,) the above does the following (2% faster ah ah ):
    #    en = 0.0
    #    fc = zeros(nat)
    #    for i in range(1,nat):
    #      r = abs( xdum[i] - xdum[i-1] )
    #      en = en + 4.0*eps[i-1]*( (sig/r)**12 - (sig/r)**6)
    #      fc[i-1] = fc[i-1]+ 4.0*eps[i-1]*(-12.0*(sig/r)**13 +6.0*(sig/r)**7)/sig
    #      fc[i] = fc[i] - 4.0*eps[i-1]*(-12.0*(sig/r)**13 +6.0*(sig/r)**7)/sig

        return en, fc

    def vdw_triple(self, xdum, nat, sig, eps):

    # returns the energy and forces for the triple backbond case
    # d_eq is the LJ eq distance

          en = 0.0
          fc = zeros(nat)

          d_eq=sig*2.**(1./6.)
          d_lat = d_eq*2.*(2**.5)/3.0    # tetrahedron base radius, if d_eq is the tetrahedron radius

          r =  xdum[2] - xdum[1]
          en = en + 4.0*eps[1]*( (sig/r)**12 - (sig/r)**6)
          fc[1] = 4.0*eps[1]*(-12.0*(sig/r)**13 +6.0*(sig/r)**7)/sig
          fc[2] = -fc[1]

          r = xdum[1] - xdum[0]

          den = d_lat**2 + r**2
          en = en + 12.0*eps[0]*( (sig**2/den)**6 - (sig**2/den)**3)

          fc[0] = 72.0*eps[0]*r*(1.-2.0*sig**6/den**3)/den**4
          fc[1] = fc[1] - fc[0]

          r = xdum[3] - xdum[2]

          den = d_lat**2 + r**2
          en = en + 12.0*eps[2]*( (sig**2/den)**6 - (sig**2/den)**3)

          fc[3] = -72.0*eps[2]*r*(1.-2.0*sig**6/den**3)/den**4
          fc[2] = fc[2] -fc[3]

          return en, fc

    def verlet(self, x0, xm, dt, force, mass, nat, vel):

        dstep = vel*dt/2.0  # vel is the opening speed
        xnew=zeros(nat)
        dwork = 0.0
        xnew[:] = x0[:]

        xnew[0] = xnew[0] - dstep
        xnew[3] = xnew[3] + dstep
        dwork = -dstep*force[0] + dstep*force[3]
        xnew[1:3] = x0[1:3] + (x0[1:3]-xm[1:3]) + force[1:3]*(dt**2)/mass

        return xnew, dwork

    def run(self):

        nat = 4  # number of atoms first and last mimic walls.
        x0 = zeros(nat)
        y0 = zeros(nat)
        xm = zeros(nat)
        t = 0.0

        #--------- potential setting
        epsilon = 1.0
        eps = zeros(nat-1)
        eps[0] = epsilon*self.frac_stif
        eps[1] = epsilon
        eps[2] = epsilon*self.frac_stif
        sig = 1.0
        d_eq = sig*2.0**(1./6.0)     # vdW eq. distance
        mass = 2.0  # mass of each atom
        red_mass = mass/2.0  # reduced mass (mass of the optical oscillator)
        aco_mass = mass*2.0  # mass of the acoustic oscillator

        period = (2.0*pi*sig/(6.0*2.**(1./3)))*(red_mass/epsilon)**0.5                 #--------- period of a free LJ dimer

        if self.i_potential == "stiffer_lj":
            per_opt = (2.0*pi*sig/(6.0*2.**(1./3)))*(red_mass/(0.5*eps[0]+eps[1]))**0.5    #--------- period of the optical mode (stiffer single backbond)
            per_aco = (2.0*pi*sig/(6.0*2.**(1./3)))*(aco_mass/(2.0*eps[0]))**0.5           #--------- period of the acoustic mode

        elif self.i_potential == "triple_back":
            per_opt = (2.0*pi*sig/(6.0*2.**(1./3)))*(red_mass/(0.5*eps[0]/3.+eps[1]))**0.5    #--------- period of the optical mode
            per_aco = (2.0*pi*sig/(6.0*2.**(1./3)))*(aco_mass/(2.*eps[0]/3.))**0.5

        speed_of_sound = 6.0*(2.*epsilon/mass)**0.5

        w_lj = 2*pi/period        # --------- angular frequencies
        w_opt = 2*pi/per_opt
        w_aco = 2*pi/per_aco

        # --------- opening velocity & time step

        vel = self.frac_vel * speed_of_sound    # velocity  as fraction of the speed of sound of a LJ 1D lattice
        dt = self.frac_tim * period            # time step as fraction of the period of a LJ dimer

        works = []
        for i_sample in range(0, self.nconfig):

            if i_sample == 0:
                start_time = time.time()

            # ---------equilibrium atoms positions
            if self.i_potential == "stiffer_lj":
                for i in range(nat):
                    x0[i] = (-float(nat-1)/2.0 + float(i))*d_eq
            elif self.i_potential == "triple_back":
                    x0[1] = -0.5*d_eq
                    x0[2] = 0.5*d_eq
                    x0[0] = x0[1]-d_eq/3.
                    x0[3] = x0[2]+d_eq/3.

            # PHONON SETUPS

            en_opt = self.frac_opt*epsilon   # ---------energy of the opt phonon in units of epsilon
            en_aco = self.frac_aco*epsilon   # ---------energy of the opt phonon in units of epsilon

            # -------------

            dis_opt = (2.0*en_opt/(red_mass * w_opt**2))**0.5  # displacement for the opt phonon to have that energy
            dis_aco = (2.0*en_aco/(2.*mass * w_aco**2))**0.5   # displacement for the opt phonon to have that energy

            phs_opt = 2.0*pi*rand()    # Present positions in relative phon coord (with random phases)
            xopt = dis_opt*cos(w_opt*t + phs_opt)
            xoptm = dis_opt*cos(w_opt*(t-dt) + phs_opt)

            phs_aco = 2.0*pi*rand()   # past positions in relative phon coord (with random phases)
            xaco = dis_aco*cos(w_aco*t + phs_aco)
            xacom = dis_aco*cos(w_aco*(t-dt) + phs_aco)

            xdtot=zeros(2)  # total current displ vector in relative phon coord
            xdtot[0]=+xopt/2.0 + xaco
            xdtot[1]=-xopt/2.0 + xaco

            xdmtot=zeros(2) # total past displ vector in relative phon coords
            xdmtot[0]=+xoptm/2.0 + xacom
            xdmtot[1]=-xoptm/2.0 + xacom

            # setting up atomic positions in absolute coords

            xm[:] = x0[:]

            x0[1] = x0[1]+xdtot[0]   # present
            x0[2] = x0[2]+xdtot[1]
            xm[1] = xm[1]+xdmtot[0]   # past time step
            xm[2] = xm[2]+xdmtot[1]

            delxd = zeros(2)

            delxd[0] = xm[1]-x0[1]
            delxd[1] = xm[2]-x0[2]

            # correction section for anharmonicity essentially this will be microcanonic dynamics with kin + pot fixed at the desired target
            # while epot and ekin won't necessarily be equal, so that half the target does not immediately give temperature (just a very good estimate)
            #
            # estimated potential: -(eps[0]+eps[1]+eps[2]) + phonon energy

            if self.i_potential == "stiffer_lj":
                e_bottom = -(eps[0]+eps[1]+eps[2])

            elif self.i_potential == "triple_back":
                e_bottom = -(3.0*eps[0]+eps[1]+3.*eps[2])

            e_est = 0.5*red_mass*w_opt**2*xopt**2 + 0.5*aco_mass*w_aco**2*xaco**2
            epot, fc = self.en_and_for(x0, nat, sig, eps)

            #print("estimated and real epots: ",e_est,epot + (eps[0]+eps[1]+eps[2]))

            for i in range(1, 7):  # iterative trick.

                rat_corr = e_est/(epot - e_bottom)
                x0[1] = x0[1]-xdtot[0]  # present
                x0[2] = x0[2]-xdtot[1]

                xdtot[:] = xdtot[:]*rat_corr**0.5

                x0[1] = x0[1]+xdtot[0]
                x0[2] = x0[2]+xdtot[1]

                epot, fc = self.en_and_for(x0, nat, sig, eps)
                #print("estimated and real epots: ",e_est,epot + (eps[0]+eps[1]+eps[2]))

            xm[1] = x0[1] + delxd[0]
            xm[2] = x0[2] + delxd[1]

            nstep = int(self.n_sigmas*sig/(vel*dt))
            work = 0.0

            veldrift = zeros(4)
            veldrift[0] = -vel/2.
            veldrift[1] = -vel/2.
            veldrift[2] = vel/2.
            veldrift[3] = vel/2.

            clf()

            for i in range(1, nstep):

                epot, fc = self.en_and_for(x0, nat, sig, eps)
                xp, dwork = self.verlet(x0, xm, dt, fc, mass, nat, vel)

                work += dwork
                velo = (xp - xm)/(2.0*dt)
                ekin = 0.5*mass*(velo[1]**2+velo[2]**2)

                etot = epot + ekin
                #   ekinotto = ekinotto + epot + 0.25*mass*vel**2 + eps[0]+eps[2] -en_opt - en_aco
                #
                #  A DEFINITION OF WORK DONE BY THE TWO FORCES:

                ekinotto = 0.5*mass*((velo[1]-veldrift[1])**2 + (velo[2]-veldrift[2])**2)   # kinetic energy of the two oscillating fragments
                                                                                            # in the translating ref. frames
                ekinotto += (epot - e_bottom)                  # total oscillator energy in the translating ref. frame
                ekinotto -= (en_opt + en_aco)                  # ..minus the initial oscillator energy before the pull
                ekinotto -= eps[1]                             # ..minus the energy it took to break the central bond
                                                               #   (this correction can be excluded as as an offset)
                #   ekinotto =  ekinotto + 2.*(0.5*mass*(vel/2.)**2)  # ..plus translation kin. energy of the two atoms
                                                                      #  (this last correction is temperature independent and
                                                                      #  could be excluded from the definition)

                # so this is the energy to break the bond at the net of the static bond energy and the translational kinetic energy
                # it can be negative (in which case, initial oscillations helped the bond breaking)
                # or positive (in which case, initial oscillations -temperature- made the bond breaking more difficult)
                # or zero (typical of zero initial oscillation and fracture at very low velocity)

                #if i%1000==1: print(("%6d"+"%12.8f"*6) % (i,ekin,epot,etot, work, etot + work, ekinotto))

                xm[:] = x0[:]
                x0[:] = xp[:]
                t += dt

            if i_sample == 0:
                end_time = time.time()
                run_time = end_time - start_time
                # if (vel/speed_of_sound) > 0.2:
                #     plot(x0, y0, marker='o', markersize=50)
                #     axis([-3.0*float(nat)/2.0, 3.0*float(nat)/2.0, -1, 1])
                #     plt.show()
                #     clf()

            self.current_runs += 1

            progress = (self.current_runs/self.total_runs)*100
            est_time_remaining = round(((self.total_runs - self.current_runs)*run_time)/60, 2)

            print(str(round(progress, 2)), "% Complete.      Max time remaining = ", str(est_time_remaining), " Minutes.")

            try:
                with open(self.results_temp, 'rw') as q:
                    works = json.load(q)

                    works.append(ekinotto)

                    json.dump(works, q)
            except:
                works.append(ekinotto)
                with open(self.results_temp, 'w') as q:
                    json.dump(works, q)

        os.remove(self.results_temp)
        return mean(works), statistics.pstdev(works)
