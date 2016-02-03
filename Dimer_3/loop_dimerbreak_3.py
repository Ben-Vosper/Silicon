""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
(1) Force integration
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#
import numpy as np
import matplotlib.pyplot as plt
import sys
from pylab import *


#------- RUN PARAMETERS:
# (1) opening velocity will be a fraction frac_vel of the speed of sound in LJ 1D Lattice
# frac_vel will be in n_velocities steps up to dfrac_vel*n_velocities
dfrac_vel    = 0.01
n_velocities = 100
veloc_start  = 0.05

print('max velocity:', dfrac_vel*float(n_velocities) )
#

# (2) energy of the optical and acoustic modes as a fraction of the LJ binding energy
# frac_opt (frac_aco) will grow in n_energies steps of dfrac_opt (dfrac_aco)

n_energies   = 1
dfrac_opt = 0.000000001
dfrac_aco = 0.000000001
ener_start_opt = 0.00002
ener_start_aco = 0.00002

print('phonon max energies. ',float(n_energies)*dfrac_opt,float(n_energies)*dfrac_aco)
#

# (3) number of configurations to average upon (for which velocity & energy combination) 
nconfig = 10
print('average will be taken on : ', nconfig,'  runs with different initial configurations ')

#
# (4) time step will be a fixed fraction of the period of a LJ dimer
frac_tim = 0.002
print('time step (as a fraction of the period of a LJ dimer): ',frac_tim)
#
# (5) number of LJ sigmas the pulling will be taken to (must be enough to fully break the bond)
#     the number of MD steps will be:   nstep = int(n_sigmas*sig/(vel*dt))
n_sigmas = 7.0
print('brk distance (in units of LJ sigmas): ',n_sigmas) 
#
#---BACKBOND FORCE FIELD CHOICES, and stiffnes ratios [energy of the back_bond(s) relative to the central bond]
#   the backbond force field can be either a stiffer Lennard-Jones or a triple (each not stiffer) LJ backbond proxected along z
#
#i_potential = "triple_back"
#frac_stif = 1.0   
i_potential = "stiffer_lj"
frac_stif = 1.2
#
print('backbond potential linking the two atoms to the opening walls: ', i_potential)     
print('fractional stiffnes factor applied to the backbond potential: ',frac_stif)
#
#

# If True, plot atom positions every 1000 steps
graphics = False
works = zeros(nconfig)


def en_and_for(xdum,nat,sig,eps,i_potential):

    if(i_potential=="stiffer_lj"): 
        en,fc = vdw_stiffer(xdum,nat,sig,eps)   

    elif (i_potential=="triple_back"):
        en,fc =  vdw_triple(xdum,nat,sig,eps)

    return(en,fc);


def velocity(time,tau,vel):
# this builds up to a the velocity vell with gaussian speed, as time grows from zero and exceeds tau
     
    vell = vel*(1.- exp(-time**2.0/(2.0*tau**2.0)))
    return(vell)


def vdw_stiffer(xdum,nat,sig,eps):
# returns the energy and forces for the single stiffer backbond case
#
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

    return(en,fc);


def vdw_triple(xdum,nat,sig,eps):

# returns the energy and forces for the triple backbond case
# d_eq is the LJ eq distance 


      en = 0.0
      fc = zeros(nat)

      d_eq=sig*2.**(1./6.) 
      d_lat = d_eq*2.*(2**.5)/3.0    # tetrahedron base radius, if d_eq is the tetrahedron radius  
#
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

      return(en,fc);


def verlet(x0,xm,dt,force,mass,nat,vell):

    dstep = vell*dt/2.0  # vel is the opening speed
    xnew=zeros(nat)
    dwork = 0.0
    xnew[:] = x0[:]

    xnew[0] = xnew[0] - dstep
    xnew[3] = xnew[3] + dstep
    dwork = -dstep*force[0] + dstep*force[3]
    xnew[1:3] = x0[1:3] + (x0[1:3]-xm[1:3]) + force[1:3]*(dt**2)/mass

    return (xnew,dwork);


def energy_of_breaking(frac_vel,frac_opt,frac_aco,frac_stif,nconfig,n_sigmas,frac_tim,i_vel,i_ener,en_matrix,en_mask,i_potential):

    nat = 4  # number of atoms first and last mimic walls.
    x0 = zeros(nat)
    xm = zeros(nat)
    xp = zeros(nat)
    y0 = zeros(nat)
    fc = zeros(nat)
    work=0.0

    #--------- potential setting
    epsilon = 1.0
    eps = zeros(3)

    if(i_potential=="stiffer_lj"): 
        eps[0] = epsilon*frac_stif
        eps[1] = epsilon
        eps[2] = eps[0]

    elif (i_potential=="triple_back"):  #---same stiffness for now, here there are 3 backbonds on each side
        eps[0] = epsilon*frac_stif
        eps[1] = epsilon
        eps[2] = eps[0]


    sig = 1.0                  # vdW sigma
    d_eq=sig*2.0**(1./6.0)     # vdW eq. distance

    mass = 2.0  # mass of each atom
    red_mass = mass/2.0 # reduced mass (mass of the optical oscillator)
    aco_mass = mass*2.0 # mass of the acoustic oscillator

    period  = (2.0*pi*sig/(6.0*2.**(1./3)))*(red_mass/epsilon)**0.5                #--------- period of a free LJ dimer

    if(i_potential=="stiffer_lj"): 
        per_opt = (2.0*pi*sig/(6.0*2.**(1./3)))*(red_mass/(0.5*eps[0]+eps[1]))**0.5    #--------- period of the optical mode (stiffer single backbond)
        per_aco = (2.0*pi*sig/(6.0*2.**(1./3)))*(aco_mass/(2.0*eps[0]))**0.5           #--------- period of the acoustic mode 

    elif (i_potential=="triple_back"):

#        per_opt = 2.0*pi*sig*( 2**(1./6.)/(2.0* 21.**0.5 )*(red_mass/(eps[1]))**0.5    #--------- period of the optical mode (triple backbond)
        per_opt = (2.0*pi*sig/(6.0*2.**(1./3)))*(red_mass/(0.5*eps[0]/3.+eps[1]))**0.5    #--------- period of the optical mode 

#        per_aco = 2.0*pi*sig*(1./(2.**(4./3)*(3.**0.5)))*(aco_mass/(2.0*eps[0]))**0.5    #--------- period of the acooustic mode 
        per_aco = (2.0*pi*sig/(6.0*2.**(1./3)))*(aco_mass/(2.*eps[0]/3.))**0.5           #--------- period of the acoustic mode 


    speed_of_sound = 6.0*(2.*epsilon/mass)**0.5
    #print(“speed of sound in LJ 1D: ", speed_of_sound)

    w_lj  = 2*pi/period        #--------- angular frequencies
    w_opt = 2*pi/per_opt
    w_aco = 2*pi/per_aco

    # --------- opening velocity & time step

    vel = frac_vel * speed_of_sound    # velocity  as fraction of the speed of sound of a LJ 1D lattice
    dt  = frac_tim * period            # time step as fraction of the period of a LJ dimer

    for i_sample in range(0,nconfig):
        
        print('iteration: ', i_sample) 
        t = 0.0

        # ---------initial equilibrium atoms positions

        if(i_potential=="stiffer_lj"): 
            for i in range(0,nat):
                x0[i]=(-float(nat-1)/2.0 + float(i))*d_eq
        elif (i_potential=="triple_back"):
                x0[1]= -0.5*d_eq
                x0[2]=  0.5*d_eq
                x0[0]=  x0[1]-d_eq/3.
                x0[3]=  x0[2]+d_eq/3.

        # PHONON SETUPS

        en_opt = frac_opt*epsilon # ---------target energy of the opt phonon in units of epsilon
        en_aco = frac_aco*epsilon # ---------target energy of the opt phonon in units of epsilon

        # -------------

        dis_opt  = ( 2.0*en_opt/(red_mass * w_opt**2) )**0.5 # displacement for the opt phonon to have that energy
        dis_aco  = ( 2.0*en_aco/(aco_mass * w_aco**2) )**0.5 # displacement for the opt phonon to have that energy

        phs_opt = 2.0*pi*rand()    # Present and past positions in relative optical phon coord (with random phases)
        xopt  = dis_opt*cos(w_opt*t + phs_opt)
        xoptm = dis_opt*cos(w_opt*(t-dt) + phs_opt)

        phs_aco = 2.0*pi*rand() # Present and past positions in relative acoustic phon coord (with random phases)
        xaco = dis_aco*cos(w_aco*t + phs_aco)
        xacom = dis_aco*cos(w_aco*(t-dt) + phs_aco)

        # total current and past displ vector in relative phon coord
        # (only for the two central atoms; they shift the same acoustically, half in opposite ways optically) 

        # current
        xdtot=zeros(2)  
        xdtot[0]=+xopt/2.0 + xaco
        xdtot[1]=-xopt/2.0 + xaco
        # past 
        xdmtot=zeros(2) 
        xdmtot[0]=+xoptm/2.0 + xacom
        xdmtot[1]=-xoptm/2.0 + xacom

        # setting up initial atomic positions in absolute coords

        xm[:]=x0[:]

        x0[1] = x0[1]+xdtot[0]    # present
        x0[2] = x0[2]+xdtot[1]
        xm[1] = xm[1]+xdmtot[0]   # initial "past" time step
        xm[2] = xm[2]+xdmtot[1]

        # stores the opposite of the initial displacement in a vector for both atoms (delta x diplacement)

        delxd = zeros(2)

        delxd[0] = xm[1]-x0[1]
        delxd[1] = xm[2]-x0[2]

        # correction section for anharmonicity essentially this will be microcanonic dynamics with kin + pot fixed at the desired target
        # while epot and ekin won't necessarily be equal, so that half the target does not immediately give temperature (just a very good estimate)
        #
        # estimated potential: relaxed ground state energy  + phonon energy

        e_bottom=0.0 
        if(i_potential=="stiffer_lj"): 
            e_bottom =  -(eps[0]+eps[1]+eps[2])

        elif (i_potential=="triple_back"):
            e_bottom =  -(3.0*eps[0]+eps[1]+3.*eps[2])

        e_est = 0.5*red_mass*w_opt**2*xopt**2 +0.5*aco_mass*w_aco**2*xaco**2
        epot,fc = en_and_for(x0,nat,sig,eps,i_potential)
        print('estimated and real epots: ',e_est,epot -e_bottom) 

        # iterative trick: to converge to the displacement which will be associated 
        # with exactly the target elastic e_est in spite of the system being anharmonic. A few steps suffice. 
        
        for i in range(1,7):  
            rat_corr = e_est/(epot -e_bottom)
            x0[1] = x0[1]-xdtot[0] # takes away the presently enforced displacement 
            x0[2] = x0[2]-xdtot[1]

            xdtot[:]=xdtot[:]*rat_corr**0.5  # corrects it with the sqrt of the target ratio, so a beter displacement is now on 

            x0[1] = x0[1]+xdtot[0] # enforces the new displacement 
            x0[2] = x0[2]+xdtot[1]

            epot,fc = en_and_for(x0,nat,sig,eps,i_potential) # computes the new exact epot at the new displacement
            
            print('estimated and real epots: ',e_est,epot -e_bottom) 

        xm[1] = x0[1] + delxd[0]
        xm[2] = x0[2] + delxd[1]

        
        if(vel > 0.000000000001): 
            nstep = int(n_sigmas*sig/(vel*dt))
        else: 
            nstep = 10000       # some standard for 0 velocity runs

        ekkek = zeros(nstep)
        ekavg = zeros(nstep)
        timev = zeros(nstep)
        worko = zeros(nstep)
        worvg = zeros(nstep)
        avg_worko = zeros(nstep)
        avg2_worko = zeros(nstep)
        avg3_worko = zeros(nstep)
        tau = per_aco*10.0

        dwork= 0.0
        work = 0.0
        en = 0.0
        epot=0.0

        veldrift= zeros(4)

        # standard expected drift for the four atoms

        veldrift[0]=-vel/2.
        veldrift[1]=0.0
        veldrift[2]=0.0
        veldrift[3]= vel/2.
        clf()

        for i in range(1,nstep):

            epot,fc = en_and_for(x0,nat,sig,eps,i_potential)
            time= dt*float(i)
            tau = 0.2*period
            vell = velocity(time,tau,vel)
            xp, dwork = verlet(x0,xm,dt,fc,mass,nat,vell)

            work = work + dwork

        #   worko[i]=work

        #   avg_worko[i] = avg_worko[i-1]*exp(-dt/tau) + (1.0 - exp(-dt/tau))*work
            if graphics and i%1000==0:
                plt.figure(1)
                clf()
                plot(x0,y0, marker='o', markersize=50)
                axis([-3.0*float(nat)/2.0, 3.0*float(nat)/2.0, -1, 1])
                plt.show()
                plt.draw()
                
            #    plt.figure(2)
            #    plot(t,fc[0],'bo')
            #    draw()

            velo = zeros(nat)
            ekin=0.0
            velo = (xp -xm)/(2.0*dt)
            ekin = 0.5*mass*(velo[1]**2+velo[2]**2)

        #  ekkek[i]=ekin
        #  timev[i]=t/dt
        #
        #           if i%1==0:
        #             plt.figure(3)
        #             plot(t,ekin,'+')
        #             draw()
        # if i%10==1:
        #   plt.figure(4)
        #   plot(t,work,'go')
        #   draw()

            etot = epot + ekin
        #
        #  A DEFINITION OF WORK DONE BY THE TWO FORCES:
        #
        # for best convergence, we first define the drifting ref. frame

            ekdrift =0.0
            veldrift[0] = -vel/2.
            veldrift[1] = -vel/2.
            veldrift[2] = vel/2.
            veldrift[3] = vel/2.

            ekdrift = 0.5*mass*(veldrift[1]**2 + veldrift[2]**2)  # total drifting kinetic energy for the two atoms

            ekinotto = 0.0
            ekinotto = 0.5*mass*((velo[1]-veldrift[1])**2 + (velo[2]-veldrift[2])**2) # kinetic energy of the oscillating fragments
                                                                                      # in their correct drifting ref. frames
            ekinotto =  ekinotto + (epot - e_bottom)               # total oscillators energy in the drifting ref. frames 

            #ekinotto =  ekinotto + ekdrift                         # ..plus drifting kin. energy of the two atoms

            ekinotto =  ekinotto -en_opt -en_aco                   # ..minus the initial oscillator energy before the pull
            ekinotto =  ekinotto - eps[1]                          # and we choose as zero the energy it took to break the central bond

# so this is the energy to break the bond at the net of the static bond energy and the translational kinetic energy
# it can be negative (in which case, initial oscillations helped the bond breaking)
# or positive (in which case, initial oscillations -temperature- made the bond breaking more difficult)
# or zero (typical of zero initial oscillation and fracture at very low velocity)

# Note that in the static frame the oscillator bound to a drifting wall is an open system, with the wall exerting an oscillating force
# which makes on average zero work. However, the total energy in the static frame has a m*vdrift*(v-vdrift) component which we filter out here
# to achieve fast convergence 


            if i%5000==1: print(('%6d'+'%12.8f'*6) % (i,ekin,epot,etot, work, etot + work, ekinotto)) 

            xm[:]=x0[:]
            x0[:]=xp[:]
            t = t + dt


        if ((x0[2]-x0[1]) > sig*(float(n_sigmas)-3.0)) :
            print('LOST COHERENCE at frac_vel,frac_opt,frac_aco', frac_vel,frac_opt,frac_aco) 

            if(  ((x0[1] < 0) and (x0[2] < 0)) or ((x0[1] > 0) and (x0[2] > 0)) ):    # a stiff, not a weak bond was broken
                print('STIFF BOND BROKEN, POSITIONS: ',x0) 


            elif(  ((x0[1] < 0) and (x0[2] > 0)) or ((x0[1] > 0) and (x0[2] < 0)) ):  # TWO stiff bonds were broken
                if(  ((x0[1]-x0[0]) > 2*sig) and ((x0[3]-x0[2]) > 2*sig)):                # (making sure)
                    print('TWO STIFF BONDS BROKEN, POSITIONS: ',x0) 

            else:
                print('DISASTER, POSITIONS: ',x0)
                sys.exit(0)

        works[i_sample] = ekinotto
        en_matrix[i_vel,i_ener,i_sample]=ekinotto


    # averages ekin
    # nwhile = int(per_aco/dt) # averaged over the period of bouncing over once broken (same as acoustic mode)
    # print(nwhile)

    #for i in range(0,nstep-nwhile):
    # dummy = 0.0
    # for j in range(0,nwhile):
    #   dummy = dummy + ekkek[i+j]
    # dummy = dummy/float(nwhile)
    # ekavg[i+nwhile/2]=dummy

    # for i in range(nstep-5*nwhile,nstep-nwhile):
    #   dummy = 0.0
    #   for j in range(0,nwhile):
    #     dummy = dummy + worko[i+j]
    #   dummy = dummy/float(nwhile)
    #   worvg[i+nwhile/2]=dummy
    # avg3_worko = (worvg + roll(worvg, nwhile/2))/2

    # dummy = 0.0
    # for j in range(1,nwhile+1):
    #    dummy = dummy + worko[nstep-j]
    # dummy = dummy/float(nwhile)
    # print "config, work done: ", i_sample, dummy
    # print "heat: ", -(dummy+eps[1])
    # works[i_sample] = -(dummy+eps[1])

    # avg2_worko = (avg_worko + roll(avg_worko, nwhile/2))/2
    # print "heat 2", avg_worko[nstep-1]

    xav = average(works)

    #---------printing out
    #print "LJ period, optical and acoustic periods: ",period, per_opt, per_aco
    #print "atomic mass:",mass

    print('nsteps:      ',nstep)
    print('time step:   ',dt)
    #print('vdW sigma:   ',sig)
    #print('vdw eps:     ',epsilon)
    print('    ')
    print('velocity:  ',frac_vel)
    print('phon. ens. ',frac_opt,frac_aco)
    print('    ')
    print('values: ',works)

    var = average(works*works) - xav**2
    sigmamed = (var/float(nconfig))**0.5
    print('nconfig, average work, sigmamed: ',nconfig,xav,sigmamed)

    #plt.figure(5)
    #plot(timev,worvg,'+')
    #plot(timev,ekkek,'+')
    #plot(timev,ekavg,'+')
    #draw()

    return (xav,sigmamed)


# -----------------------------------------------------------------------------------
# MAIN LOOP HERE, over equispaced velocities of pull, equispaced oscillator energies, each for a
# number nconfig of initial random conditions (phases of the two oscillators) 
#

velocities = zeros(n_velocities)
energies   = zeros((n_velocities,n_energies))
en_error   = zeros((n_velocities,n_energies))

en_matrix = zeros((n_velocities,n_energies,nconfig))
en_mask   = zeros((n_velocities,n_energies,nconfig),dtype=np.int)
ener_heat = zeros(n_energies)  


for i_vel in range(0,n_velocities):
    frac_vel = veloc_start + dfrac_vel*float(i_vel)
    velocities[i_vel]= frac_vel
    for i_ener in range(0,n_energies):

        frac_opt =  ener_start_opt + dfrac_opt*float(i_ener)
        frac_aco =  ener_start_aco + dfrac_aco*float(i_ener)
        ener_heat[i_ener] = frac_opt + frac_aco

        x_average = 0.0
        sig_avera = 0.0
        x_average, sig_avera = energy_of_breaking(frac_vel,frac_opt,frac_aco,frac_stif,nconfig,n_sigmas,frac_tim,i_vel,i_ener,en_matrix,en_mask,i_potential)
        print(x_average, sig_avera) 
        energies[i_vel,i_ener] = x_average
        en_error[i_vel,i_ener] = sig_avera


savetxt('energies.dat', energies)
savetxt('en_error.dat', en_error)
#savetxt('en_matrix.dat', en_matrix)
#savetxt('en_mask.dat', en_mask)

#works[:]=0.0
jdum = 0
xav  = 0.0

energies_corr = zeros((n_velocities,n_energies))
en_error_corr = zeros((n_velocities,n_energies))

for i_vel in range(0,n_velocities):
    for i_ener in range(0,n_energies):
        jdum = 0
        for i_sample in range(0,nconfig):
            if(en_mask[i_vel,i_ener,i_sample]==0):
                works[jdum]=en_matrix[i_vel,i_ener,i_sample]
                jdum = jdum +1
            work_sh=works[:jdum]
            xav = average(work_sh)
            var = average(work_sh*work_sh) - xav**2
            sigmamed = (var/float(jdum))**0.5

        energies_corr[i_vel,i_ener] = xav
        en_error_corr[i_vel,i_ener] = sigmamed


plt.errorbar(velocities, energies_corr, yerr=en_error_corr, ls="none", marker="+")
plt.xlabel("$\\frac{v}{v_s}$", size=24)
plt.ylabel("$\\bar W$ / $\epsilon$", size=18)
plt.show()

print("energies, corrected:")
print(energies_corr)
print("  ")
print('en_error, corrected:') 
print(en_error_corr) 
