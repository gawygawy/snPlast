#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
python version of the snPlast model from Brzosko et al., 2017 
adapted for the circular open field maze experiment
"""

import os
import numpy as np
from multiprocessing import Pool
import argparse
import matplotlib.pyplot as plt
#import pdb

from help_func import make_lattice_pts, find_wall_cells, action_directions, \
weights_walls, weight_matrix_init, initialise_variables, convolution, \
neuron, update_weights, reflect_vec

def simulate_agent(ACh_flag, n_trials, switch_tr, eta_DA, eta_ACh, w_max, w_init, plot):
    
    np.random.seed()

#%% task parameters 

    rew_flag1 = 1 # first reward 
    rew_flag2 = 0
    
    t_max = 15*10**3 # maximum duration of the trial 
    t_end = t_max
    
    c = [-0.43, 0.43] # coordinates of first reward 
    r_goal = 0.3 # radius of the rewarded area
    
    c2 = [0.43, -0.43]  # coordinates of the second reward 
    r_goal2 = 0.3 
    
    #%% place cells 
    n_pc = 121
    
    # create a grid of place cells tiling the circular open field maze 
    rad_lattice = 6.1 
    lattice_pts = make_lattice_pts(rad_lattice) 
    boundary_indx = find_wall_cells(lattice_pts) # boundary cells 
    
    sigma_pc = 0.4 # spacing between place cells 
    pc = lattice_pts*sigma_pc # place cell coordinates in the circular arena 
    
    rho_pc_max = 400*10**(-3)
    
    #%% Action neurons 
    
    # Neuron model - epsp of action neurons 
    eps0 = 20 # epsp scaling constant 
    tau_m = 20 # membrane time constant 
    tau_s = 5 # synaptic rise time 
    
    n_action = 40 # number of action neurons 
    a0 = .08 #velocity per rate unit 
    
    # spike train is filtered with kernel with parameters: 
    tau_gamma = 50 # kernel rise time (ms)
    v_gamma = 20 # kernel decay (ms)
    
    actions = action_directions(n_action, a0) # directions coded by the ring of action neurons 
    
    #%% 
    store_pos = np.zeros((t_max*n_trials,2))
    firing_rate_store = np.zeros((n_action, t_max*n_trials))
    
    #%% initial weights 
    
    #w_max = 3 # bounds of feed-forward weights 
    w_min = 1 
    w_walls = weights_walls(pc, boundary_indx, n_action, n_action+n_pc)
    w_tot = weight_matrix_init(w_init, n_action, n_pc)
    
    #%% synaptic plasticity parameters 
    
    a_plus = 1 # amplitude of the STDP window, pre-post 
    a_minus = 1 # post-pre 
    tau_pre_post = 10 # time-constant for STDP window 
    tau_post_pre = 10
    tau_e = 2*10**3 # time-constant of the eligibility trace 
    #eta_DA = 0.01 # parameter controlling amount of dopamine 
    #eta_ACh = 10**-3*2 # parameter controlling amount of acetylcholine 
    
    #%%
    starting_pos = np.array([[-1.6, -1.2], [1.6, 1.2]]) # coordinates of the two rewards 
    i = 0 # time counter 
    tr = 0 # trial number
    
    w_tot_old = w_tot[:, :n_pc]
    
    #%% for plotting
    rad_out = rad_lattice*sigma_pc
    rad_in = rad_out*0.5
    
    if plot:
        plt.ion()
        fig, ax = plt.subplots(ncols=2, nrows=2)
        
        # trajectory of the agent in the maze 
        circ_step = 2*np.pi/100
        x_circ = np.cos(np.arange(-np.pi, np.pi+circ_step, circ_step))
        y_circ = np.sin(np.arange(-np.pi, np.pi+circ_step, circ_step))
        new_col = 'g'
        old_col = 'k'
        
        ring_inner, = ax[0, 0].plot(rad_in*x_circ, rad_in*y_circ, 'k')
        ring_outer, = ax[0, 0].plot(rad_out*x_circ, rad_out*y_circ, 'k--')
        rew_circ1, = ax[0, 0].plot(c[0] + r_goal*x_circ, c[1]+r_goal*y_circ, color=new_col)
        rew_circ2, = ax[0, 0].plot(c2[0] + r_goal2*x_circ, c2[1]+r_goal2*y_circ, color=old_col)
        start_point, = ax[0, 0].plot([], [], '.r', markersize = 10)
        ax[0, 0].axhline(0, color='black', lw=2)
        ax[0, 0].axvline(0, color='black', lw=2)
        line, = ax[0, 0].plot([], [], lw=1)
        
        # activity of action neurons on a single trial. Firing rate of neurons 
        # shown as color map against time 
        firing_rate_plot = firing_rate_store[:,int((np.floor((i)/t_max))*t_max):int((np.floor((i)/t_max))*t_max+t_end)]
        plt1 = ax[0, 1].imshow(firing_rate_plot, \
                 aspect='auto', cmap=plt.cm.jet, vmin=0, vmax=1)
        plt.colorbar(plt1, ax=ax[0, 1], ticks=range(6))
        
        # weights over the maze, averaged over action neurons 
        w_plot = np.zeros((13, 13))
        grid_indices = lattice_pts+6
        plt2 = ax[1, 0].imshow(w_plot, aspect='auto', cmap=plt.cm.jet, vmin=w_min, vmax=w_max, origin='lower')
        plt.colorbar(plt2, ax=ax[1, 0], ticks=range(6))
        
        # vector field (averaging synaptic weights from each place cell to the action neurons) - agent's policy preference map 
        plt3 = ax[1, 1].quiver(pc[:,0], pc[:,1], np.ones(n_pc), np.ones(n_pc), scale_units='xy', pivot='mid', scale=20)
        
        plt.draw()
        plt.pause(0.001)
    
    #%% begin trial 
    
    # store for trial outcome 
    # 1 - 1st reward found; -1 - wrong visit to 1st rewarded site 
    # 2 - 2nd reward found; -2 - wrong visit to 2nd rewarded site 
    rew_visits = np.zeros(n_trials) 
    
    while i < t_max*n_trials:
               
        i +=1
        t = i % t_max
        
        if t == 1: # beginning a new trial
            rew_found = 0 # set reward flag to 0
            well_found = 0
            pos = starting_pos[np.random.randint(2), :] # start positions randomized
            first_pos = pos
        
            tr += 1 
            t_well = t_max # time of reward - initialized at t_max at the start  
            t_extreme = t_well + 300 
            if plot:
                firing_rate_plot[:, :] = 0
        
            y_action_neurons, canc,\
            epsp_rise, epsp_decay, epsp,\
            action_firing_rate, action_firing_decay, action_firing_rise,\
            last_post_spike, pre_post_trace, post_pre_trace, trace,\
            eligibility_trace = initialise_variables(n_action, n_pc)
            
            # switch the location of the reward 
            if tr in switch_tr: 
                if plot:
                    new_col, old_col = old_col, new_col
                    plt.setp(rew_circ1, color=new_col)
                    plt.setp(rew_circ2, color=old_col)
                rew_flag1, rew_flag2 = rew_flag2, rew_flag1  
        
    #%% place cell activity 
        
        # rate - inhomogeneous poisson process 
        rho_pc = rho_pc_max*np.exp(-np.sum((pos-pc)**2, axis=1)/sigma_pc**2)
        
        if t > t_well:
            rho_pc = 0 # place cells turned of if a well is encountered 
        
        x = np.random.rand(n_pc) <= rho_pc # realization of  place cell spikes 
        
    #%% reward 
        
        if np.sum((pos-c)**2) <= r_goal**2 and well_found == 0: # agent enters 1st rew site 
            well_found = 1
            t_well = t
            if rew_flag1 == 1:
                rew_found = 1
                t_extreme = t_well + 300 # trial terminates 300s after reward 
                rew_visits[tr-1] = 1
            elif rew_flag2 == 1: # wrong well 
                t_extreme = t_well
                rew_visits[tr-1] = -1
        elif np.sum((pos-c2)**2) <= r_goal2**2 and well_found == 0: # agent enters 2nd rew site
            well_found = 1
            t_well = t
            if rew_flag2 == 1:
                rew_found = 1
                t_extreme = t_well + 300 # trial terminates 300s after reward 
                rew_visits[tr-1] = 2
            elif rew_flag1 == 1: # wrong well 
                t_extreme = t_well
                rew_visits[tr-1] = -2
    
    #%% input to action neurons 
        
        input_spikes = np.vstack((x[:,np.newaxis], y_action_neurons))
        input_spikes = input_spikes*canc # reset if postsynaptic neuron spiked 
    
    #%% activity of the action neuron   
        epsp_rise = epsp_rise*canc
        epsp_decay = epsp_decay*canc 
        
        epsp, epsp_decay, epsp_rise = convolution(epsp_decay, epsp_rise, tau_m,\
                                                      tau_s, eps0, input_spikes, w_tot*w_walls)
        
        # spike train of the action neuron 
        y_action_neurons, last_post_spike, canc = neuron(epsp, last_post_spike, tau_m, i, n_action)
        
        blank_w = np.ones((1, n_action))
        
        # filter spike train with kernel gamma 
        action_firing_rate, \
        action_firing_decay, \
        action_firing_rise = convolution(action_firing_decay, action_firing_rise,\
                                         tau_gamma, v_gamma, 1, y_action_neurons, blank_w)
        
        firing_rate_store[:, int(i-1)] = action_firing_rate.flatten()
        action_now = np.dot(action_firing_rate.T, actions)/n_action
     
    #%% STDP
        
        pre_post_trace, post_pre_trace, eligibility_trace, trace, W = \
        update_weights(a_plus, a_minus, tau_pre_post, tau_post_pre, x, y_action_neurons,\
                       pre_post_trace, post_pre_trace, trace, tau_e)
        
        w_tot[:,:n_pc] = w_tot[:,:n_pc] - eta_ACh*W*(ACh_flag) # depress weights only if ACh is present 
        w_tot[:,:n_pc] = np.clip(w_tot[:,:n_pc], w_min, w_max) # limit weights 
    
    #%% update agent position 
    
        store_pos[int(i-1),:] = pos # for plotting
        pos = pos + action_now 
        
        origin_distance = np.sqrt(pos[0][0]**2 + pos[0][1]**2)
        
        if origin_distance >= rad_out: # check if agent is out of the boundary
            
            bounce = reflect_vec(pos, action_now)
            pos = pos + bounce
    
        if t > t_extreme and t < t_max:
            i = tr*t_max-1
            t_end = t_extreme
    
    #%%
        if t == 0: # end of the trial, t_max reached
            
            t = t_max
        
            if rew_found == 1:
                # retroactively potentiate weights through an eligibility trace 
                w_tot[:,:n_pc] = w_tot_old + eta_DA*eligibility_trace 
                w_tot[:,:n_pc] = np.clip(w_tot[:,:n_pc], w_min, w_max)
            
            #store weights for the beginning of the next trial 
            w_tot_old = w_tot[:, :n_pc]
            
            # to plot the policy (vector field)
            ac = np.dot(actions.T, (w_tot_old*w_walls[:, :n_pc])/a0) # preferred action for each pc
            ac[:, np.sort(boundary_indx)] = 0
            #%%
            
            ind_start = int((tr-1)*t_max)
            ind_end = ind_start+t_end
            path_x = store_pos[ind_start:ind_end,0]
            path_y = store_pos[ind_start:ind_end,1]
            
            if plot:
                line.set_data(path_x, path_y)
                start_point.set_data(first_pos[0], first_pos[1])
                    
                firing_rate_trial = firing_rate_store[:,ind_start:ind_end]
                firing_rate_plot[:, :t_end] = firing_rate_trial/firing_rate_trial.max() # normalize to 0 to 1
                plt1.set_data(firing_rate_plot)
                    
                w_vals = np.mean(w_tot[:,0:n_pc], axis=0)
                w_plot[grid_indices[:,0], grid_indices[:,1]] = w_vals
                plt2.set_data(w_plot.T)
                    
                plt3.set_UVC(ac[0,:],ac[1,:])
                    
                ax[0, 0].set_title('Trial No: %i' %tr)
                    
                fig.canvas.update()
                fig.canvas.flush_events()
        
            t_end = t_max
   
    return rew_visits
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-cr", "--cores", help="number of CPU cores to use", type=int, default=1)
    parser.add_argument("-it", "--niter", help="number of iterations to run model", type=int, default=100)
    parser.add_argument("-ACh", "--ACh_flag", help="ACh_flag", type=float, default=1)
    parser.add_argument("-wm", "--wmax", help ="weight ceiling", type=float, default=3)
    parser.add_argument("-wi", "--w_init", help ="initialized weights", type=float, default=2)
    parser.add_argument("-tr", "--trials", help="number of trials per experiment", type=int, default=40)
    parser.add_argument("-sw", "--switch_tr", help="switch trial", type=int, nargs='+', default=[20])
    parser.add_argument("-eD", "--eta_DA", help="learning rate with dopamine", type=float, default=0.01)
    parser.add_argument("-eA", "--eta_ACh", help="learning rate with ACh", type=float, default=10**-3*2)
    parser.add_argument("-out", "--output", help="out folder")
    parser.add_argument("-p", "--plot", help="flag to plot figure", default=False, action="store_true")
    
    args = parser.parse_args()
    
    p = Pool(args.cores)
    if args.plot and args.niter == 1:
        simulate_agent(args.ACh_flag, args.trials, args.switch_tr, args.eta_DA, args.eta_ACh, args.wmax, args.w_init, args.plot)
    else:
        result_ACh = p.starmap(simulate_agent, [(args.ACh_flag, args.trials, args.switch_tr, args.eta_DA, args.eta_ACh, args.wmax, args.w_init, args.plot)]*args.niter)
        rew_visits = np.stack(result_ACh)
        filename = os.path.join(args.output, "DA" + str(args.eta_DA) + "_ACh" + str(args.eta_ACh) + "_wmax" + str(args.wmax) + "_winit" + str(args.w_init) + ".csv")
        np.savetxt(filename, rew_visits.astype(int), fmt="%i", delimiter=",")