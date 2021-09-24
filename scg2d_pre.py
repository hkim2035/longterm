import math
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import itasca as it
it.command("python-reset-state false")

from itasca import zonearray as za
from itasca import gridpointarray as gpa


def strainv(bot,top,length):
    global gpa
    zdis_top = gpa.disp()[:,2][top-1].sum()/len(top)
    zdis_bot = gpa.disp()[:,2][bot-1].sum()/len(bot)
    axial_disp = abs(zdis_top - zdis_bot)
    axial_strain = axial_disp / float(length)
    return axial_disp, axial_strain



#=========== User setup ===========
gr = {'x': 50, 'y': 1, 'z': 110}
rd = {'x': 0.05, 'y': 0.001, 'z': 0.11}   
rp = {'den': 2700.0, 'bulk': 32.69e9, 'shear': 20.56e9, 'coh': 20.0e6, 'fric': 40.0, 'ten': 10.0e6, 'dil': 10.0}

# property random seed on/off
seed = [8, 5, 1974]

#====== initialize variables ======
#axial load: 64.0 MPa
aload = -64.0e6
eq_ratio = 0.001

#initial damage = 0
D0 = 0.0 

f_res_II = 36.0
f_res_I = 36.0
c_res = 1.0e6

#activation stess or stress threshold
s_act = 0.15 

#Charles equation
I_a0 = 0.1322    # 1/s
I_n = 29.39
II_a0 = 0.01424    # 1/s
II_n = 46.56
  
I_A = (1.0 - D0)/(math.exp(I_a0)**I_n)
II_A = (1.0 - D0)/(math.exp(II_a0)**II_n)

first_damage_on_I = 0
first_damage_on_II = 0
first_damage_time_I = 10.0e100
first_damage_time_II = 10.0e100
first_damage_i_I = 0
first_damage_j_I = 0
first_damage_i_II = 0
first_damage_j_II = 0
 
first_damage_type = 'none'

#scg
n_mode2 = 0
n_mode1 = 0

#plastic state
n_pl = 0
  
tot_time = 0.0
time_current = 0.0
time_pre = 0.0



title = "Long-term stability analysis based on damage mechanics: {} MPa,random properties".format(-aload)

it.command("""
model new
model title \'{}\'
""".format(title))

it.command("""
zone create brick size {} {} {} point 0 ({},{},{}) point 1 ({},{},{}) point 2 ({},{},{}) point 3 ({},{},{})
""".format(gr['x'],gr['y'],gr['z'], 0,0,0, rd['x'],0,0, 0,rd['y'],0, 0,0,rd['z']))

it.command("""                             
zone cmodel assign mohr-coulomb
zone property den {} bulk {} shear {} cohesion {} friction {} tension {} dilation {}
""".format(rp['den'], rp['bulk'],rp['shear'],rp['coh'],rp['fric'],rp['ten'],rp['dil']))


it.command("""
zone gridpoint initialize velocity-z 0 range position-z -0.01,0.01 
zone gridpoint fix velocity-z range position-z -0.01,0.01
""")

noz = len(za.live_mechanical())

# top, bottom ids()
gp_top = gpa.ids()[gpa.pos()[:,2] == rd['z']]
gp_bot = gpa.ids()[gpa.pos()[:,2] == 0]


c_mean = rp['coh']
c_std  = c_mean*0.15
f_mean = rp['fric']
f_std = f_mean*0.0125
t_mean = rp['ten']
t_std = t_mean*0.015

np.random.seed(seed[0])
za.set_prop_scalar('cohesion', c_mean + np.random.uniform(-1,1,noz)*c_std)
np.random.seed(seed[1])
za.set_prop_scalar('friction', f_mean + np.random.uniform(-1,1,noz)*f_std)
np.random.seed(seed[2])
za.set_prop_scalar('tension', t_mean + np.random.uniform(-1,1,noz)*t_std)

za.set_extra(1, np.zeros(noz))
za.set_extra(2, np.ones(noz)*D0)
za.set_extra(3, np.ones(noz)*D0)
za.set_extra(4, np.ones(noz)*D0)

c_min_id = za.prop_scalar('cohesion').argmin()
c_min = za.prop_scalar('cohesion')[c_min_id]
f_min = za.prop_scalar('friction')[c_min_id]
t_min = za.prop_scalar('tension')[c_min_id]


# dt scaling factor for random properties
s_factor = 0.5 
Nphi  = (1.0+math.sin(math.radians(f_mean)))/(1.0-math.sin(math.radians(f_mean)))
m_UCS = 2.0*c_mean*math.sqrt(Nphi)

# adaptitve values
if abs(aload/m_UCS) > 0.9:
    dt_min = 0.1
    n_step = 50
else:
    if abs(aload/m_UCS) > 0.8:
        dt_min = 1.0
        n_step = 20
    else:
        dt_min = 10.0
        n_step = 10

load_time = 60.0 #60 sec
ttf_I = (abs(aload)/m_UCS/math.exp(I_a0))**(-I_n) 
ttf_II = (abs(aload)/m_UCS/math.exp(II_a0))**(-II_n) 
ttf = s_factor * min(ttf_I, ttf_II)
stop_time = (1.5 * max(ttf_I, ttf_II))
dt_ini = 10**int(math.log(ttf/float(n_step)))
dt = min(dt_ini,dt_min)

axial_disp, ev = strainv(gp_bot,gp_top,rd['z'])
# s_zz 
sigmav = za.stress()[:,2,2].sum()/noz

# array: tot_time, unbal, ev, n_mode1, n_mode2, n_pl, AE_count, dt, axial_disp, sigmav






















