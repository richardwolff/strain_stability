import pandas as pd
import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import time
from scipy.stats import gamma
from numba import jit
#from numba.typed import List
import config
from sklearn.metrics import auc

@jit(nopython=True)
def run_sim_jit_traj(num_iters,num_reps,sigma, K, tau,delta_t,sqrt_delta_t,sig_tau,init_val):

    psi_list=np.random.normal(loc=0,scale=sqrt_delta_t,size=(num_iters,num_reps))       

    trajectory = []
    #trajectory = List()
    trajectory.append(init_val*np.ones(num_reps))

    for i in range(num_iters):

        x_i = trajectory[-1]
        x_i = x_i + delta_t*(x_i/tau)*(1 - x_i/K) + sig_tau*x_i*psi_list[i]
        trajectory.append(x_i)

    return trajectory

@jit(nopython=True)
def run_sim_jit_final(num_iters,num_reps,sigma, K, tau,delta_t,sqrt_delta_t,sig_tau,init_val):

    psi_list=np.random.normal(loc=0,scale=sqrt_delta_t,size=(num_iters,num_reps))       

    x_i = init_val*np.ones(num_reps)
    
    for i in range(num_iters):

        x_i = x_i + delta_t*(x_i/tau)*(1 - x_i/K) + sig_tau*x_i*psi_list[i]

    return x_i


class slm:

    def __init__(self,sigma, K, tau,delta_t=1/10000,init_val=None):
        
        self.sigma = sigma
        self.K = K
        self.tau = tau
        self.delta_t = delta_t
        self.sqrt_delta_t = np.sqrt(delta_t)
        self.sig_tau = np.sqrt(1.0*sigma/tau)
        self.x_i = init_val
        self.alpha = (2/(1.0*sigma)) - 1
        self.beta = (2/(1.0*sigma*K))
        self.afd = gamma(a=self.alpha, scale=1/self.beta)
    
    def set_init_val(self,init_val):
        self.x_i = init_val
    
    def set_stationary_val(self,num_reps):
        
        self.x_i = self.afd.rvs(size=num_reps)
    
    def update(self,psi):

        self.x_i = self.x_i + self.delta_t*(self.x_i/self.tau)*(1 - self.x_i/self.K) + self.sig_tau*self.x_i*psi
    
    def run_sim(self,num_iters,num_reps,record_steps=False):
        
        if isinstance(self.x_i, type(None)):
            self.set_stationary_val(num_reps)

        if record_steps:
            
            self.trajectory = run_sim_jit_traj(num_iters,num_reps,self.sigma, self.K, self.tau,self.delta_t,self.sqrt_delta_t,self.sig_tau,self.x_i)
        
        else:
            
            self.x_i = run_sim_jit_final(num_iters,num_reps,self.sigma, self.K, self.tau,self.delta_t,self.sqrt_delta_t,self.sig_tau,self.x_i)

def fit_SLM_params(obs_data,n=None,strain_dates=None):
    
    if n == None:
        xbar = obs_data.values.mean()
        std_dev = obs_data.values.std()
        beta = (xbar/std_dev)**2
        
        sigma = 2/(beta+1)
        
        K = xbar/(1 - sigma/2)
           
        return {"K":K,"beta":beta,"xbar":xbar,"sigma":sigma}
    
    else:
        
        obs_data = obs_data[:n]
        
        if strain_dates is None:
            xbar = obs_data.values.mean()
            
        else:
            xbar = auc(strain_dates[:n],obs_data)/(strain_dates[n] - strain_dates[0])
        
        std_dev = obs_data.values.std()
        beta = (xbar/std_dev)**2
        
        sigma = 2/(beta+1)
        
        K = xbar/(1- sigma/2)
           
        return {"K":K,"beta":beta,"xbar":xbar,"sigma":sigma} 
    
    
def plot(S):
    
    fig,axs = plt.subplots(1,2,figsize = (30, 13),gridspec_kw={'width_ratios': [6, 1]})
    axs = axs.ravel()
    ax = axs[0]
    ax.plot(np.linspace(0,S.delta_t*T,int(T)+1),S.trajectory)
    ax.axhline(S.K,zorder=10,color="red");
    ax.grid(True)
    ax.set_xlabel("Time",size=20)
    ax.set_ylabel("Abundance",size=20)
    ax.set_ylim([S.afd.ppf(1e-4),S.afd.ppf(1-1e-4)])

    ax.plot(S.trajectory[::int(1/S.delta_t)],"o",color="red")

    axhist = axs[1]

    x = np.linspace(0,1,T)

    axhist.plot(S.afd.pdf(x),x,lw=3)
    axhist.hist(S.trajectory,bins=30,alpha=.3,density=True,orientation="horizontal");
    axhist.set_ylim([S.afd.ppf(1e-4),S.afd.ppf(1-1e-4)])
    axhist.set_yticks([])
    axhist.set_xticks([])

    fig.tight_layout()    
    