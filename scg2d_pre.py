import math
import os

import matplotlib.pyplot as plt
import numpy as np

import itasca as it

it.command("python-reset-state false")

from itasca import zonearray as za
from itasca import gridpointarray as gpa




def strainv(gpa,bot,top,length):
    
    zdis_top = gpa.disp()[:,2][top-1].sum()/len(top)
    zdis_bot = gpa.disp()[:,2][bot-1].sum()/len(bot)
    axial_disp = abs(zdis_top - zdis_bot)
    axial_strain = axial_disp / float(length)
    
    return axial_disp, axial_strain


def common_set(tmpzz, coh, fric, ten, group):
    
    tmpzz.set_extra(1,1)
    tmpzz.set_prop('cohesion', coh)
    tmpzz.set_prop('friction', fric)
    tmpzz.set_prop('tension', ten)
    tmpzz.set_group(group)

def hist_writing(tot_time, unbal, dt):
    
    global his_time, hist_unbal, hist_dt

    hist_time.append(tot_time)
    hist_unbal.append(unbal)
    hist_dt.append(dt)

def damage_model():

    global sI, sIc, sII, sIIc, d_I_sum, d_II_sum
    global tot_time, unbal, ev, n_mode1, n_mode2, n_pl, AE_count, dt, axial_disp, sigmav
    global his_time, hist_unbal, hist_ev, hist_n_mode1, hist_n_mode2
    global his_n_pl, hist_AE_count, hist_dt, hist_axial_disp, hist_sigmav

    eachZone = [dZone for id, dZone in enumerate(it.zone.list())]

    npstate = np.array([t.state(True) for t in eachZone])

    nps1 = np.array([np.sort(t.stress_prin())[0] for t in eachZone])
    nps2 = np.array([np.sort(t.stress_prin())[1] for t in eachZone])
    nps3 = np.array([np.sort(t.stress_prin())[2] for t in eachZone])

    npcoh = za.prop_scalar('cohesion')
    npfric = za.prop_scalar('friction')
    npten = za.prop_scalar('tension')

    old_ex1 = za.extra(1)
    old_DI = za.extra(2)
    old_DII = za.extra(3)

    sII = nps3 - nps1
    tmp1 = lambda x: math.tan(math.radians(x))
    vfunc1 = np.vectorize(tmp1)
    tfric = vfunc1(npfric)
    tmp2 = lambda x: math.sin(math.radians(x))
    vfunc2 = np.vectorize(tmp2)
    tsphi = vfunc2(npfric)
    tmp3 = lambda x: (1.0+x)/(1.0-x)
    vfunc3 = np.vectorize(tmp3)
    tanphi = vfunc3(npfric)
    tmp4 = lambda a,b: 2.0*a*math.sqrt(b)
    vfunc4 = np.vectorize(tmp4)
    tcn2 = vfunc4(npcoh,tanphi)

    ttcut = npten if npfric == 0 else min(npten,npcoh/npfric)
    ts1b = tanphi*nps3-tcn2
    sIIc = nps3-ts1b
    sI = nps3
    sIc = min(npten,npcoh/npfric)


   
#==================================================================== 
#=========== User setup =============================================
sav_folder = "D:\\Itasca\\Flac3d600_application_data\\My Projects\\longterm\\"

#grid & properties
gr = {'x': 50, 'y': 1, 'z': 110}
rd = {'x': 0.05, 'y': 0.001, 'z': 0.11}   
rp = {'den': 2700.0, 'bulk': 32.69e9, 'shear': 20.56e9, 'coh': 20.0e6, 'fric': 40.0, 'ten': 10.0e6, 'dil': 10.0}

# property random seed on/off (3 numbers)
seed = [8, 5, 1974]

#====== initialize variables ======
#axial load: 64.0 MPa
aload = -64.0e6

title = "Long-term stability analysis based on damage mechanics: {} MPa,random properties".format(-aload)

#unbal < eq_ratio -> OK
eq_ratio = 0.001

#activation stess or stress threshold
s_act = 0.15 

#initial damage = 0
D0 = 0.0 

#properties_residual
f_res_I = 36.0
f_res_II = 36.0
c_res = 1.0e6

#Charles equation
I_a0 = 0.1322    # 1/s
I_n = 29.39
II_a0 = 0.01424    # 1/s
II_n = 46.56
  
I_A = (1.0 - D0)/(math.exp(I_a0)**I_n)
II_A = (1.0 - D0)/(math.exp(II_a0)**II_n)

# [check, time, id, type]
first_dI_info = [False, 10.0e100, None, ""]
first_dII_info = [False, 10.0e100, None, ""]
 
#scg
n_mode1 = 0
n_mode2 = 0

#plastic state
n_pl = 0

d_I_sum = 0
d_II_sum = 0


tot_time = 0.0
time_current = 0.0
time_pre = 0.0


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

#number of zones
nozone = len(za.live_mechanical())

# top, bottom ids()
gp_top = gpa.ids()[gpa.pos()[:,2] == rd['z']]
gp_bot = gpa.ids()[gpa.pos()[:,2] == 0]

#mean and std
c_mean = rp['coh']
c_std  = c_mean*0.15
f_mean = rp['fric']
f_std = f_mean*0.0125
t_mean = rp['ten']
t_std = t_mean*0.015

np.random.seed(seed[0])
za.set_prop_scalar('cohesion', c_mean + np.random.uniform(-1,1,nozone)*c_std)
np.random.seed(seed[1])
za.set_prop_scalar('friction', f_mean + np.random.uniform(-1,1,nozone)*f_std)
np.random.seed(seed[2])
za.set_prop_scalar('tension', t_mean + np.random.uniform(-1,1,nozone)*t_std)

za.set_extra(1,np.zeros(len(za.ids())))  # 0-elastic 1-plastic
za.set_extra(2,np.zeros(len(za.ids())))  # DI
za.set_extra(3,np.zeros(len(za.ids())))  # DII
    
#properties of cohesion-minimum zone. why?
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
#ttf: time to failure
ttf_I = (abs(aload)/m_UCS/math.exp(I_a0))**(-I_n) 
ttf_II = (abs(aload)/m_UCS/math.exp(II_a0))**(-II_n) 
ttf = s_factor * min(ttf_I, ttf_II)

stop_time = (1.5 * max(ttf_I, ttf_II))

dt_ini = 10**int(math.log(ttf/float(n_step)))
dt = min(dt_ini,dt_min)

axial_disp, ev = strainv(gpa, gp_bot,gp_top,rd['z'])
# s_zz 
sigmav = za.stress()[:,2,2].sum()/nozone

# pre
plot_time_dt = 2.0*dt_ini
pl_id = 1
pl_str = str(int(abs(aload)/1.0e6))+"MPa"
ini_converged = False
second_converged = False

hist_tubdt = list()
append1 = hist_tubdt.append
hist_time = list()
hist_unbal = list()
hist_ev = list()
hist_n_mode1 = list()
hist_n_mode2 = list()
hist_n_pl = list()
hist_AE_count = list()
hist_dt = list()
hist_axial_disp = list()
hist_sigmav = list()


# ###### 1st stage: initial 60 sec #########
it.command("""
plot create 'ini'
plot 'ini' item create zone-boundary 
plot 'ini' item create zone label state average
plot 'ini' item create axes
""")



increment = 5   #sec
for ii in range (increment, 60+increment, increment):
    app_load = float(ii)/60.0*aload
    tot_time += increment
    
    it.command("zone face apply stress-normal {} range position-z {},{}".format(app_load,rd['z']-0.01,rd['z']+0.01))
    ini_converged = False
    for jj in range (1,501,1):
        it.command("model step 100")
        unbal = it.zone.unbal()
        print("{} sec unbal: {}".format(tot_time, unbal))
        if unbal < eq_ratio:
            ini_converged = True
            append1([tot_time, unbal, dt])
            break 
    
plt_fname = sav_folder + pl_str + "_initial_loading.svg"
sav_fname = sav_folder + pl_str + "_initial_loading.f3sav"
it.command("""
model save '{}'
plot 'ini' export svg filename '{}'
""".format(sav_fname, plt_fname))

#exit check
if ini_converged == False:
    print("===== Result =====")
    print("Failure during initial stage")
    print(" - Failure time (sec): {}".format(tot_time))
    print(" - Failure stress (Pa): {}".format(app_load))
    exit()

if abs(d_I_sum+d_II_sum) <= 1.0e-30:
    print("===== Result =====")
    print("Stresses of all zones are below the activation stress")
    exit()


#***** 2nd stage: until stop_time *****
high_unbal = 1000.0
low_unbal = 1.0

it.command("""
plot create 'second'
plot 'second' item create zone-boundary 
plot 'second' item create zone label state average
plot 'second' item create axes
""")

while tot_time < stop_time:
    if second_converged == False:
        break
    tot_time += dt
    damage_model()
    it.command("model step 10")
    if it.zone.unbal() > high_unbal:
        dt = max(dt/10.0,dt_min)
    if it.zone.unbal() < low_unbal:
        dt = min(dt*2.0,dt_ini)        
    
    second_converged = False
    for jj in range (1,501,1):
        it.command("model step 100")
        unbal = it.zone.unbal()
        print("{} sec unbal: {}".format(tot_time, unbal))
        if it.zone.unbal() < eq_ratio:
            second_converged = True
            hist_writing(tot_time, unbal, dt)
            break 

plt_fname = sav_folder + pl_str + "_second_loading.svg"
sav_fname = sav_folder + pl_str + "_second_loading.f3sav"
it.command("""
model save '{}'
plot 'second' export svg filename '{}'
""".format(sav_fname, plt_fname))

#===== making graph and csv using history data =====

# workflow automation -> simple report pdf maker ?