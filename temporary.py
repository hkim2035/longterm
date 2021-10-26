zmf = np.logical_and(za.live_mechanical(), za.extra(1)==1)
za.set_prop_scalar('cohesion', c_res*zmf)
za.set_prop_scalar('friction', f_res_II*zmf)
za.set_prop_scalar('tension', 0.0*zmf)

zme = np.logical_and(za.live_mechanical(), za.extra(1)!=1)
zsig = za.ids()[zme]
    
zsig = zz.stress_prin()
ts1 = -min(zsig[0], zsig[1], zsig[2])
ts3 = -max(zsig[0], zsig[1], zsig[2])
tc_coh = zz.prop('cohesion')
tc_fric = zz.prop('friction')
tc_ten = zz.prop('tension')
told_DI = zz.extra(2)
told_DII = zz.extra(3)

sII= ts3-ts1
tfric = math.tan(math.radians(tc_fric))
tsphi = math.sin(math.radians(tc_fric))
tanphi = (1.0+tsphi)/(1.0-tsphi)
tcn2 = 2.0*tc_coh*math.sqrt(tanphi)
ttcut = tc_ten if tfric == 0 else min(tc_ten,tc_coh/tfric)
ts1b = tanphi*ts3-tcn2
sIIc = ts3-ts1b
sI = ts3
sIc = min(tc_ten,tc_coh/tfric)

if sI > 0.:
    tdDI = I_A*((sI/sIc)**I_n) if sI>(s_act*sIc) else 0.
else:
    tdDI = 0.
tdDI = max(tdDI, 0.)
DI = told_DI + tdDI*(time_current - time_pre)

if sII > 0:
    tdDII = II_A*((sII/sIIc)**II_n) if sII>(s_act*sIIc) else 0.
else:
    tdDII = 0.
tdDII = max(tdDII, 0.)
DII = told_DII + tdDII*(time_current - time_pre)

zz.set_extra(2,DI)
zz.set_extra(3,DII)

if DI>=1.:
    if DII>=1.:
        if DI>DII: # case1
            if first_dI_info[0]==False:
                first_dI_info = [True, time_current, zz.id(), 'D_I_mix']
                # common_set(tmpzz, coh, fric, ten, group):
                common_set(zz, c_res, f_res_I, 0,0, 'D_I_mix')
                zz.set_extra(2,1.0)
                zz.set_extra(3,1.0)
                n_mode1 += 1
                if n_mode1 == 1:
                    ft_info = [ts1, ts2, ts3, sI, sIc]
        else:      # case2
            if first_dII_info[0]==False:
                first_dII_info = [True, time_current, zz.id(), 'D_II_mix']
                common_set(zz, c_res, f_res_II, 0,0, 'D_II_mix')
                zz.set_extra(2,1.0)
                zz.set_extra(3,1.0)
                n_mode2 += 1
    else:          # case3
        if first_dI_info[0]==False:
            first_dI_info = [True, time_current, zz.id(), 'D_I']
            common_set(zz, c_res, f_res_I, 0,0, 'D_I')
            zz.set_extra(2,1.0)
            n_mode1 += 1
            if n_mode1 == 1:
                ft_info = [ts1, ts2, ts3, sI, sIc]
else:
    if DII>=1.:    # case4
        if first_dII_info[0]==False:
                first_dII_info = [True, time_current, zz.id(), 'D_II']
                common_set(zz, c_res, f_res_II, 0,0, 'D_II')
                zz.set_extra(3,1.0)
                n_mode2 += 1
    # else:          # pass
        
zyield = [z for z in it.zone.list() if ((z.extra(1) == 0) and (z.state(True)!=0))]    
for id, tmp in enumerate(zyield):    
    tmp.set_group('YIELD')    
tt_count = len(zyield)               
d_I_sum = za.extra(2).sum()    
d_II_sum = za.extra(3).sum()           