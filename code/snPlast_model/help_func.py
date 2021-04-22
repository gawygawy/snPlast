#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
helper functions to run snPlast model 
"""
import numpy as np
import collections 
import itertools

#%% functions to create place cells 

def points_in_circle(radius): 
    a = np.arange(radius + 1)
    for x, y in zip(*np.where(a[:,np.newaxis]**2 + a**2 <= radius**2)):
        yield from set(((x, y), (x, -y), (-x, y), (-x, -y),))
        
def make_lattice_pts(radius_lat):
    pts = list(points_in_circle(radius_lat))
    x, y = zip (*pts)
    ind = np.lexsort((x, y)) #sort by y, then by x
    X = np.array(x)[ind]
    Y = np.array(y)[ind]
    lattice_pts = np.vstack((X, Y)).T
    return lattice_pts

def find_wall_cells(lattice_pts):
    dist_to_o = np.sqrt(lattice_pts[:,0]**2 + lattice_pts[:,1]**2)
    ind = np.argpartition(dist_to_o, -24)[-24:] # because there are 24 pts that lie near the boundary
    return ind

#%% actions 

def theta_action(n_action):
    """
    angles of actions 
    """
    thetas = 2*np.pi*np.arange(1,n_action+1)/n_action
    
    return thetas
    
def action_directions(n_action, a0):
    """
    directions of action neurons 
    """
    thetas = theta_action(n_action)
    actions = a0*np.stack((np.sin(thetas), np.cos(thetas)))
    
    return actions.T

def reflect_vec(pos, action):
    """
    if agent encounters the boundary of the arena, it is "bounced"
    back
    """    
    normal_hat = pos/np.linalg.norm(pos)
    bounce = action - 2*np.dot(action, normal_hat.T)*normal_hat
    
    return bounce

#%% initialize at the start of each trial 

def lateral_connectivity(thetas, n_action):
    """ 
    winner-take-all mechabnism on action neurons
    by imposing lateral connectivity between action neurons 
    """
    psi = 20
    w_minus = -300
    w_plus = 100
    
    diff_theta = np.tile(thetas, (n_action,1)) - np.tile(thetas, (n_action,1)).T
    f = np.exp(psi*np.cos(diff_theta)) # lateral connectivity function 
    f = f - f*np.eye(n_action)
    
    normalised = sum(f)
    normalised = normalised[0]
    w_lateral = w_minus/n_action + w_plus*f/normalised
    
    return w_lateral 

def weight_matrix_init(w0, n_action, n_pc): 
    """
    initialise weight matrix
    """
    w_lateral = lateral_connectivity(theta_action(n_action), n_action)
    w_feedforward_init = np.ones((n_action, n_pc))*w0 # between place and action cells 
    w_tot = np.concatenate((w_feedforward_init, w_lateral), axis=1)
    
    return w_tot

def weights_walls(pc, boundary_indx, n_action, n_neurons): 
    """
    only allow 144 degrees of movement at each boundary point (=14 segments)
    """
    # find index on a 40 point circle 
    boundary_cells = pc[boundary_indx, :]
    angles = (np.arctan2(boundary_cells[:, 0], boundary_cells[:, 1]) + 2*np.pi) % (2*np.pi)
    
    # find index of boundary cell (in a ring of 40 action neurons)
    thetas = theta_action(n_action)
    diff = np.abs(angles[:, np.newaxis] - thetas[np.newaxis, :])
    indx = np.argmin(diff, axis=1) # the closest match

    w_walls = np.ones((n_action, n_neurons))
    
    for ii in np.arange(len(boundary_indx)):
        cell_ind = boundary_indx[ii]
        
        d = collections.deque(np.arange(n_action))
        w_walls[:, cell_ind] = 0 # first set all weights to 0
        d.rotate(n_action-(indx[ii]+ 13)) # 1st permitted action is 13 segments away
        actions_ind = list(itertools.islice(d, 0, 14))
        w_walls[actions_ind, cell_ind] = 1
    
    return w_walls

def initialise_variables(n_action, n_pc):
    """
    reset on each trial 
    """
    y_action_neurons = np.zeros((n_action,1))
    
    n_neurons = n_pc+n_action
    
    epsp_rise = np.zeros((n_neurons, n_action))
    epsp_decay = np.zeros((n_neurons, n_action))
    epsp_tot = np.zeros((n_neurons, n_action))
    
    action_firing_rate = np.zeros((n_action,1))
    action_firing_rise = np.zeros((n_action,1))
    action_firing_decay = np.zeros((n_action,1))
    
    canc = np.ones((1, n_action))
    
    last_post_spike = np.zeros((n_action,1)) - 1000 # last spike time of postsynaptic neuron 
    trace_pre_post = np.zeros((n_action, n_pc))
    trace_post_pre = np.zeros((n_action, n_pc))
    trace_tot = np.zeros((n_action, n_pc))
    eligibility_trace = np.zeros((n_action, n_pc))
    
    return y_action_neurons, canc, \
            epsp_rise, epsp_decay, epsp_tot, \
            action_firing_rate, action_firing_decay, action_firing_rise, \
            last_post_spike, trace_pre_post, trace_post_pre, trace_tot, eligibility_trace

#%% Spike Response Model 
            
def neuron(epsp, last_post_spike_time, tau_m, i, n_action):
    """
    for firing of the action neuron 
    """ 
    rho0 = 60*10**(-3) # maximum firing rate
    chi = -5 # scales the refractory effect 
    theta = 16 # 
    delta_u = 2 # randomness of spiking behaviour 
    
    refractory = chi*np.exp((-i + last_post_spike_time)/tau_m) 
    u = np.sum(epsp, axis=0)[:,np.newaxis] + refractory # membrane potential 
    
    # action cell activity 
    rho_action = rho0*np.exp((u - theta)/delta_u) # probabilty of emitting a spike 
    y = np.random.rand(n_action, 1) <= rho_action # realization of spike train 
    
    last_post_spike_time[y] = i 
    canc = 1 - y # 0 if postsynaptic neuron spiked, 1 if not 
    
    return y, last_post_spike_time, canc.T

#%% STDP 

def convolution(conv_decay, conv_rise, tau_m, tau_s, eps0, x, w):

    conv_decay = conv_decay + (-conv_decay)/tau_m + x*w.T
    conv_rise = conv_rise + (-conv_rise)/tau_s + x*w.T
    conv = eps0 * (conv_decay - conv_rise)/(tau_m - tau_s)
    
    return conv, conv_decay, conv_rise

def convolution2(conv1, tau, A, spikes):
    
    conv1 = conv1 + (-conv1)/tau + spikes # each spike leaves a trace 
    conv_scaled = A * conv1 
    
    return conv_scaled, conv1 

def update_weights(a_plus, a_minus, tau_plus, tau_minus, X, Y, pre_post_trace, post_pre_trace, trace, tau_e):
    """
    weight change due to STDP 
    """ 
    # pre trace without spikes - for coincident spikes 
    conv_pre_old, _ = convolution2(pre_post_trace, tau_plus, a_plus, 0) 
    # post trace without spikes - for coincident spikes 
    conv_post_old, _ = convolution2(post_pre_trace, tau_minus, a_minus, 0)
    
    # presynaptic neuron trace 
    conv_pre_scaled, pre_post_trace = convolution2(pre_post_trace, tau_plus, a_plus, X)
    # postynaptic neuron trace 
    conv_post_scaled, post_pre_trace = convolution2(post_pre_trace, tau_minus, a_minus, Y)
    
    # total synaptic change due to STDP 
    W = (conv_pre_scaled*Y + conv_post_scaled*X)* ~(X&Y) + \
    ((conv_pre_old*Y + conv_post_old*X)+(a_plus + a_minus)/2)*(X&Y)
    
    ## weight change is convoluted with eligibility trace 
    eligibility_trace, trace = convolution2(trace, tau_e, 1, W)
    
    return pre_post_trace, post_pre_trace, eligibility_trace, trace, W