#%%
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#%% Measurementes
ms_path = 'examples/bs2019/measurements.csv'
ms = pd.read_csv(ms_path, index_col=0)
ms.index = pd.to_datetime(ms.index)

t0 = pd.to_datetime('2018-04-05 00:00:00')
t1 = pd.to_datetime('2018-04-08 00:00:00')

ms = ms.loc[t0:t1]

solrad = ms['solrad'].values
Tout = ms['Tout'].values
occ = ms['occ'].values
t = (ms.index - ms.index[0]).total_seconds() / 3600.

fig, ax = plt.subplots(3, 1, figsize=(5, 3), sharex=True)
fig.set_dpi(120)

ax[0].plot(t, Tout, 'b-')
ax[0].set_ylabel("$T_{out}$ [$^\circ$C]")

ax[1].plot(t, solrad, 'b-')
ax[1].set_ylabel("$q_{sol}$ [W/m$^2$]")

ax[2].plot(t, occ, 'b-')
ax[2].set_ylabel("$n_{occ}$ [-]")

ax[2].set_xticks(np.arange(0, 73, 24))
ax[2].set_xlabel("$t$ [h]")

plt.subplots_adjust(0.13, 0.15, 0.98, 0.98)
fig.savefig('examples/bs2019/figs/inputs_mpc.pdf')


### Case 1 ####################################################################
#%% Compare estimates between FMUs

est_dir = 'examples/bs2019/case1/results/est/'

tols = [
    '1e-4',
    '1e-6',
    '1e-7',
    '1e-9',
    '1e-11'
    ]

cols = pd.read_csv(est_dir + 'r1c1_dymola_' + tols[0] +
                   '/parameters_rel.csv').columns

parameters = pd.DataFrame(index=pd.Index(tols, name='tol'),
                          columns=cols)

for t in tols:
    for p in cols:
        parameters.loc[t, p] = pd.read_csv(est_dir + 'r1c1_dymola_'
                      + t + '/parameters_rel.csv')[p].iloc[0]
    
parameters.T.plot(kind='bar')

#%% Parameter estimation: validation

est_dir = 'examples/bs2019/case1/results/est/'

tols = [
    '1e-4',
    '1e-6',
    '1e-7',
    '1e-9',
    '1e-11'
    ]

idl = pd.read_csv(est_dir + 'ideal.csv', index_col=0)
idl

vld = pd.DataFrame()

for t in tols:
    res = pd.read_csv(est_dir + 'r1c1_dymola_' + t +
                      '/vld_res.csv', index_col=0)
    vld[t] = res['T']

idl = idl.loc[:vld.index[-1]]

idl.index /= 3600.
vld.index /= 3600.

# Plot
fig, ax = plt.subplots(1, 1, figsize=(5, 3))
fig.set_dpi(130)

ax.plot(idl['T'], ls='-', label='Measurement')
ax.plot(vld['1e-11'], ls='-.', label='R1C1')
ax.legend(loc='lower right')
ax.set_xticks(np.arange(0, vld.index[-1] + 1, 24))
ax.set_ylabel('$T$ [$^\circ$C]')
ax.set_xlabel('$t$ [h]')

ax.vlines(5*24., ymin=19.5, ymax=26.9, linestyles='--', lw=0.75, color='k')
ax.set_ylim(19.5, 26.9)
ax.text(80, 26.4, "Training")
ax.text(128, 26.4, "Validation")

plt.subplots_adjust(0.1, 0.18, 0.99, 0.99)
fig.savefig('examples/bs2019/figs/validation_T.pdf')

#%% Result overview
fmu = 'r1c1_dymola_1e-11'
hrz = 4
outdir = 'examples/bs2019/case1/results/mpc/{}/h{}/'.format(fmu, hrz)

# Constraints
constr = pd.read_csv(outdir + 'constr.csv', index_col=0)
constr.index /= 3600.
constr -= 273.15

# Emulation states
xemu = pd.read_csv(outdir + '/xemu.csv')\
    .set_index('time')
xemu.index /= 3600.
xemu -= 273.15

# Control states
xctr = pd.read_csv(outdir + 'xctr.csv')\
    .set_index('time')
xctr.index /= 3600.
xctr -= 273.15

# Control inputs
u = pd.read_csv(outdir + 'u.csv')\
    .set_index('time')
u.index /= 3600.

# Optimized inputs
fig, ax = plt.subplots(2, 1, sharex=True, sharey=False,
    figsize=(5, 4))
fig.set_dpi(130)

# ax[0]
ax[0].plot(u['vpos'], 'k-', lw=2)
ax[0].set_ylim(-100, 100)
ax[0].set_ylabel('$q$ [%]')

# ax[1]
# ax[1].plot(xctr['x0'], label='Control')
ax[1].plot(xemu['cair.T'], 'r-', label=fmu)
ax[1].legend(loc='upper right')
ax[1].plot(constr['Tmin'], 'k--', lw=0.5)
ax[1].plot(constr['Tmax'], 'k--', lw=0.5)

ax[1].set_xticks(np.arange(0, u.index.values[-1] + 1, 24))
plt.minorticks_off()
ax[1].set_yticks(np.arange(19, 25, 1))

ax[1].set_xlabel('$t$ [h]')
ax[1].set_ylabel('$T_i$ [$^\circ$C]')

# ax[0] - subinterval solutions
files = os.listdir(outdir)
ufiles = list()
for f in files:
    fname = f.split('.')[0]
    if fname[0] == 'u' and len(fname) > 1:
        ufiles.append(f)

udfs = list()

for i in range(len(ufiles)):
    df = pd.read_csv(outdir + 'u{}.csv'.format(i), index_col=0)
    df.index /= 3600.
    ax[0].plot(df['vpos'], ls='--', lw=1.)

# plt.show()

#%% Compare horizons
    
fmu = 'r1c1_dymola_1e-11'
horizons = [2, 4, 6, 8, 10]

fig, ax = plt.subplots(2, 1, sharex=True, sharey=False,
    figsize=(5, 4))
fig.set_dpi(120)

Qrc = dict()

i = 0
for hrz in horizons:
    outdir = 'examples/bs2019/case1/results/mpc/{}/h{}/'.format(fmu, hrz)
    
    # Constraints
    constr = pd.read_csv(outdir + 'constr.csv', index_col=0)
    constr.index /= 3600.
    constr -= 273.15
    
    # Emulation states
    xemu = pd.read_csv(outdir + '/xemu.csv')\
        .set_index('time')
    xemu.index /= 3600.
    xemu -= 273.15
    
    # Control states
    xctr = pd.read_csv(outdir + 'xctr.csv')\
        .set_index('time')
    xctr.index /= 3600.
    xctr -= 273.15
    
    # Control inputs
    u = pd.read_csv(outdir + 'u.csv')\
        .set_index('time')
    u.index /= 3600.
    u['vpos'] *= 20  # [%] -> [W]

    Qrc[hrz] = u['vpos'].abs().sum() / 1000.

    # Actual horizon string, e.g. "6h"
    ahrz = "{}h".format(hrz)
    
    # Color map
    lspace = np.linspace(0, 1, len(horizons))
    colors = [plt.cm.winter(x) for x in lspace]
    
    # ax[0]
    ax[0].plot(u['vpos'], c=colors[i], label=ahrz)
        
    # ax[1]
    ax[1].plot(xemu['cair.T'], c=colors[i], label=ahrz)

    i += 1

ax[1].legend(loc='center', bbox_to_anchor=(0.5,-0.5), ncol=5)

ax[1].plot(constr['Tmin'], 'k--', lw=0.5)
ax[1].plot(constr['Tmax'], 'k--', lw=0.5)

ax[0].set_ylim(-2200, 2200)
    
ax[1].set_xticks(np.arange(0, u.index.values[-1] + 1, 24))
plt.minorticks_off()
ax[1].set_yticks(np.arange(19, 25, 1))

ax[0].set_ylabel('$q$ [W]')

ax[1].set_xlabel('$t$ [h]')
ax[1].set_ylabel('$T$ [$^\circ$C]')

ax[0].set_title('(a)')
ax[1].set_title('(b)')

plt.subplots_adjust(left=0.16, right=0.99, top=0.93, bottom=0.24)
fig.tight_layout()
fig.savefig('examples/bs2019/figs/case1_horizon_tol_1e-11.pdf')

#%% Computational time

# FMU
wd1 = 'examples/bs2019/case1/results/mpc/r1c1_dymola_1e-11/'

# SVM
wd2 = 'examples/bs2019/case2/results/mpc-lin/'

hdirs1 = [x[0].split('/')[-1] for x in os.walk(wd1)][1:]
hdirs2 = [x[0].split('/')[-1] for x in os.walk(wd2)][1:]

hix = [int(x[1:]) for x in hdirs1]
hix = sorted(hix)

ct1 = list()
ct2 = list()

# Number of optimization variables
nv = [x * 2 for x in hix]

# Optimization horizon [h]
oh = [x for x in hix]

for h in hix:
    with open(wd1 + "h" + str(h) + '/cputime.txt') as f:
        s = f.read().split(' ')
        x = int(s[-2])
        ct1.append(x / 60.)

    with open(wd2 + "h" + str(h) + '/cputime.txt') as f:
        s = f.read().split(' ')
        x = int(s[-2])
        ct2.append(x / 60.)

fig, ax = plt.subplots(1, 1, figsize=(5,3))
fig.set_dpi(120)

plt.plot(oh, ct1, marker='s', c='k', ls=':', lw=1., label='R1C1 FMU (tol=1e-11)')
plt.plot(oh, ct2, marker='v', c='r', ls=':', lw=1., label='SVR')

ax.set_xlabel('Optimization horizon [h]')
ax.set_ylabel('Total CPU time [min]')

ax2 = ax.twiny()
ax2.set_xticks(ax.get_xticks())
ax2.set_xlim(ax.get_xlim())
ax2.set_xticklabels([int(x * 2) for x in ax.get_xticks()])
ax2.set_xlabel('Number of optimization variables')

ax.legend()
ax.grid()

plt.subplots_adjust(0.1, 0.18, 0.99, 0.85)
fig.savefig('examples/bs2019/figs/cputime.pdf')

plt.show()

#%% Solution quality - omit CVode FMUs, they seem not working correctly

# Read all inputs and states
wd = 'examples/bs2019/case1/results/mpc/'
fmus = os.listdir(wd)
hz = '/h10/'

new_names = [y[5:].replace('_', ' ') for y in fmus]
for i in range(len(new_names)):
        new_names[i] = new_names[i].replace('dymola ', 'tol=')

cdirs = [wd + x + hz for x in fmus]
cmap = {x:y  for x, y in zip(cdirs, new_names)}

uall = pd.DataFrame()
xall = pd.DataFrame()

for c, f in zip(cdirs, fmus):
    u = pd.read_csv(c + 'u.csv', index_col=0)
    x = pd.read_csv(c + 'xemu.csv', index_col=0)
    
    uall[c] = u['vpos']
    xall[c] = x['cair.T']
    
uall = uall.rename(columns=cmap)  # Inputs
xall = xall.rename(columns=cmap)  # States

# Energy consumption
q = uall * 20.
Q = q.abs().sum() / 1000. # [kWh]

# Constraint violation
cstr = pd.read_csv(wd + 'r1c1_dymola_1e-9/h2/constr.csv')
cstr['time'] = cstr['time'].astype(int)
cstr = cstr.set_index('time')

vup = xall.copy()
vlo = xall.copy()

for c in xall:
    vup[c] = xall[c] - cstr['Tmax']
    vup[c].loc[vup[c] < 0] = 0
    
    vlo[c] = cstr['Tmin'] - xall[c]
    vlo[c].loc[vlo[c] < 0] = 0

vtot = vup + vlo
vtot = vtot.sum()

# Case order for plots
cord = ['tol=1e-4', 'tol=1e-6', 'tol=1e-7', 'tol=1e-9', 'tol=1e-11']

# Ordered results
Qord = [Q.loc[x] for x in cord]
vord = [vtot.loc[x] for x in cord]

# Show both on scatter plot
n_horizons = 5
lspace = np.linspace(0, 1, n_horizons)
colors = [plt.cm.jet(x) for x in lspace]
markers = ['o', 's', 'D', 'v', '^']

fig, ax = plt.subplots(figsize=(5, 3))
fig.set_dpi(120)

for q, v, l, c, m in zip(Qord, vord, cord, colors, markers):
    plt.scatter(q, v, label=l, c=c, s=100, marker=m)

ax.set_xlabel('Total energy consumption $Q$ [kWh]')
ax.set_ylabel('Temperature violation $v_T$ [Kh]')
ax.legend(loc='center', ncol=3, bbox_to_anchor=(0.45,-0.4))
ax.grid()

plt.subplots_adjust(0.18, 0.35, 0.97, 0.95)
fig.savefig('examples/bs2019/figs/solution_quality.pdf')

# Case 2 ######################################################################
#%% Model validation

svr_x = pd.read_csv('examples/bs2019/case2/results/mpc-lin/vld_xctr.csv', index_col=0)
svr_x = svr_x.rename(columns={'cair.T':'T'})
svr_x.index /= 3600.
svr_x['T'] -= 273.15

rc_x = pd.read_csv('examples/bs2019/case2/results/mpc-lin/vld_xemu.csv', index_col=0)
rc_x = rc_x.rename(columns={'cair.T':'T'})
rc_x.index /= 3600.
rc_x['T'] -= 273.15

fig, ax = plt.subplots(1, 1, figsize=(5, 3))
fig.set_dpi(130)

ax.plot(rc_x.index, rc_x['T'].values, label='R1C1')
ax.plot(svr_x.index, svr_x['T'].values, label='SVR', ls='--')

ax.legend()
ax.set_xlabel('$t$ [h]')
ax.set_ylabel('$T$ [$^\circ$C]')

ax.set_xticks(np.arange(0, 97, 24))
plt.subplots_adjust(0.13, 0.15, 0.98, 0.98)
fig.savefig('examples/bs2019/figs/svr_validation.pdf')

#%% Overview

hrz = 6
outdir = 'examples/bs2019/case2/results/mpc-lin/h{}/'.format(hrz)

# Constraints
constr = pd.read_csv(outdir + 'constr.csv', index_col=0)
constr.index /= 3600.
constr -= 273.15

# Emulation states
xemu = pd.read_csv(outdir + '/xemu.csv')\
    .set_index('time')
xemu.index /= 3600.
xemu -= 273.15

# Control states
xctr = pd.read_csv(outdir + 'xctr.csv')\
    .set_index('time')
xctr.index /= 3600.
xctr -= 273.15

# Control inputs
u = pd.read_csv(outdir + 'u.csv')\
    .set_index('time')
u.index /= 3600.

# Optimized inputs
fig, ax = plt.subplots(2, 1, sharex=True, sharey=False,
    figsize=(5, 4))
fig.set_dpi(130)

# ax[0]
ax[0].plot(u['vpos'], 'k-', lw=2)
ax[0].set_ylim(-100, 100)
ax[0].set_ylabel('$q$ [%]')

# ax[1]
ax[1].plot(xctr['x0'], label='Control')
ax[1].plot(xemu['cair.T'], 'r--', label='Emulation')
ax[1].legend(loc='upper left')
ax[1].plot(constr['Tmin'], 'k--', lw=0.5)
ax[1].plot(constr['Tmax'], 'k--', lw=0.5)

ax[1].set_xticks(np.arange(0, u.index.values[-1] + 1, 24))
plt.minorticks_off()
ax[1].set_yticks(np.arange(19, 25, 1))

ax[1].set_xlabel('$t$ [h]')
ax[1].set_ylabel('$T_i$ [$^\circ$C]')

# ax[0] - subinterval solutions
files = os.listdir(outdir)
ufiles = list()
for f in files:
    fname = f.split('.')[0]
    if fname[0] == 'u' and len(fname) > 1:
        ufiles.append(f)

udfs = list()

for i in range(len(ufiles)):
    df = pd.read_csv(outdir + 'u{}.csv'.format(i), index_col=0)
    df.index /= 3600.
    ax[0].plot(df['vpos'], ls='--', lw=1.)

#%%
horizons = [2, 4, 6, 8, 10]

fig, ax = plt.subplots(2, 1, sharex=True, sharey=False,
    figsize=(5, 4))
fig.set_dpi(120)

Qsvr = dict()

i = 0
for hrz in horizons:
    outdir = 'examples/bs2019/case2/results/mpc-lin/h{}/'.format(hrz)
    
    # Constraints
    constr = pd.read_csv(outdir + 'constr.csv', index_col=0)
    constr.index /= 3600.
    constr -= 273.15
    
    # Emulation states
    xemu = pd.read_csv(outdir + '/xemu.csv')\
        .set_index('time')
    xemu.index /= 3600.
    xemu -= 273.15
    
    # Control states
    xctr = pd.read_csv(outdir + 'xctr.csv')\
        .set_index('time')
    xctr.index /= 3600.
    xctr -= 273.15
    
    # Control inputs
    u = pd.read_csv(outdir + 'u.csv')\
        .set_index('time')
    u.index /= 3600.
    u['vpos'] *= 20.  # [%] -> [W]
    
    Qsvr[hrz] = u['vpos'].abs().sum() / 1000.
        
    # Actual horizon string, e.g. "6h"
    ahrz = "{}h".format(hrz)
    
    # Color map
    lspace = np.linspace(0, 1, len(horizons))
    colors = [plt.cm.winter(x) for x in lspace]
    
    # ax[0]
    ax[0].plot(u['vpos'], c=colors[i], label=ahrz)
        
    # ax[1]
    ax[1].plot(xemu['cair.T'], c=colors[i], label=ahrz)

    i += 1

ax[1].legend(loc='center', bbox_to_anchor=(0.5,-0.5), ncol=5)

ax[1].plot(constr['Tmin'], 'k--', lw=0.5)
ax[1].plot(constr['Tmax'], 'k--', lw=0.5)

ax[0].set_ylim(-2200, 2200)
    
ax[1].set_xticks(np.arange(0, u.index.values[-1] + 1, 24))
plt.minorticks_off()
ax[1].set_yticks(np.arange(19, 25, 1))

ax[0].set_ylabel('$q$ [W]')

ax[1].set_xlabel('$t$ [h]')
ax[1].set_ylabel('$T$ [$^\circ$C]')

ax[0].set_title('(a)')
ax[1].set_title('(b)')

plt.subplots_adjust(left=0.16, right=0.99, top=0.93, bottom=0.24)
fig.tight_layout()
fig.savefig('examples/bs2019/figs/case2_horizon.pdf')

### Case 3 ####################################################################
#%% Result vs. horizon
horizons = [9]#, 4, 6, 8, 10]

fig, ax = plt.subplots(2, 2, sharex=True, sharey=False,
    figsize=(6, 4))
fig.set_dpi(120)

i = 0
for hrz in horizons:
    outdir = 'examples/bs2019/case3/results/mpc/h{}/'.format(hrz)
    
    # Constraints
    constr = pd.read_csv(outdir + 'constr.csv', index_col=0)
    constr.index /= 3600.
    constr['Tmin'] -= 273.15
    constr['Tmax'] -= 273.15
    
    # Emulation states
    xemu = pd.read_csv(outdir + '/xemu.csv')\
        .set_index('time')
    xemu.index /= 3600.
    xemu['cair.T'] -= 273.15
    
    # Control states
    xctr = pd.read_csv(outdir + 'xctr.csv')\
        .set_index('time')
    xctr.index /= 3600.
    xctr['x0'] -= 273.15
    
    # Control inputs
    u = pd.read_csv(outdir + 'u.csv')\
        .set_index('time')
    u.index /= 3600.
        
    # Actual horizon string, e.g. "6h"
    ahrz = "{}h".format(hrz)
    
    # Color map
    lspace = np.linspace(0, 1, len(horizons))
    colors = [plt.cm.winter(x) for x in lspace]
    
    # ax[0]
    ax[0][0].plot(u['vpos'], c=colors[i], label=ahrz)
    ax[0][1].plot(u['dpos'], c=colors[i], label=ahrz)
        
    # ax[1]
    ax[1][0].plot(xemu['cair.T'], c=colors[i], label=ahrz)
    ax[1][1].plot(xemu['co2.balance.CO2ppmv_i'], c=colors[i], label=ahrz)

    i += 1

#ax[1][0].legend(loc='center', bbox_to_anchor=(0.5,-0.5), ncol=5)

ax[1][0].plot(constr['Tmin'], 'k--', lw=0.5)
ax[1][0].plot(constr['Tmax'], 'k--', lw=0.5)

ax[1][1].plot(constr['CO2min'], 'k--', lw=0.5)
ax[1][1].plot(constr['CO2max'], 'k--', lw=0.5)

ax[0][0].set_ylim(0, 105)
    
ax[1][0].set_xticks(np.arange(0, u.index.values[-1] + 1, 24))
plt.minorticks_off()
ax[1][0].set_yticks(np.arange(19, 25, 1))

ax[0][0].set_ylabel('$v_{p}$ [%]')
ax[1][0].set_xlabel('$t$ [h]')
ax[1][0].set_ylabel('$T_i$ [$^\circ$C]')

ax[0][1].set_ylabel('$d_{p}$ [%]')
ax[1][1].set_xlabel('$t$ [h]')
ax[1][1].set_ylabel('$C_i$ [ppm]')

fig.tight_layout()

#%% MPC vs PID - using emulation model results

hrz = 3

# Outputs
ydf_pid = pd.read_csv('examples/bs2019/case3/results/pid/ydf.csv',
                      index_col=0)
ydf_pid.index /= 3600.
y_pid = ydf_pid.drop(['vpos', 'dpos'], axis=1)

y_mpc = pd.read_csv('examples/bs2019/case3/results/mpc/h{}/yemu.csv'\
                      .format(hrz), index_col=0)
y_mpc.index /= 3600.

# Control inputs
u_pid = ydf_pid[['vpos', 'dpos']]
y_pid.index /= 3600.

u_mpc = pd.read_csv('examples/bs2019/case3/results/mpc/h{}/u.csv'.format(hrz),
                    index_col=0)
u_mpc.index /= 3600.

# States
x_pid = pd.read_csv('examples/bs2019/case3/results/pid/xdf.csv',
                    index_col=0)
x_pid.index /= 3600.

x_mpc = pd.read_csv('examples/bs2019/case3/results/mpc/h{}/xemu.csv'\
                      .format(hrz), index_col=0)
x_mpc.index /= 3600.

# Constraints
constr = pd.read_csv('examples/bs2019/case3/results/mpc/h{}/constr.csv'\
                     .format(hrz), index_col=0)
constr.index /= 3600.

# Plot
fig, axes = plt.subplots(2, 2, figsize=(5,3), sharex=True)
fig.set_dpi(120)

ax = axes[0][0]
ax.plot(x_pid['cair.T'] - 273.15, label='PID')
ax.plot(x_mpc['cair.T'] - 273.15, label='MPC')
ax.plot(constr['Tmin'] - 273.15, 'k--', lw=0.5)
ax.plot(constr['Tmax'] - 273.15, 'k--', lw=0.5)
ax.set_ylabel('$T_i$ [$^\circ$C]')
ax.set_yticks(np.arange(19, 26, 2))

ax = axes[0][1]
ax.plot(x_pid['co2.balance.CO2ppmv_i'], label='PID')
ax.plot(x_mpc['co2.balance.CO2ppmv_i'], label='MPC')
ax.plot(constr['CO2min'], 'k--', lw=0.5)
ax.plot(constr['CO2max'], 'k--', lw=0.5)
ax.set_ylabel('$C_i$ [ppm]')
ax.set_yticks(np.arange(400, 1001, 200))

ax = axes[1][0]
ax.plot(u_pid['vpos'], label='PID')
ax.plot(u_mpc['vpos'], label='MPC')
ax.set_ylim(0, 100)
ax.set_xlabel('$t$ [h]')
ax.set_ylabel('$v_{p}$ [%]')

ax = axes[1][1]
ax.plot(u_pid['dpos'], label='PID')
ax.plot(u_mpc['dpos'], label='MPC')
ax.set_ylim(0, 100)
ax.set_xlabel('$t$ [h]')
ax.set_ylabel('$d_{p}$ [%]')
ax.set_xticks(np.arange(0, u_pid.index.values[-1] + 1, 24))

fig.tight_layout()

axes[1][0].legend(loc='center', bbox_to_anchor=(1.15, -0.6), ncol=2)

plt.subplots_adjust(0.15, 0.25, 0.99, 0.98)
fig.savefig('examples/bs2019/figs/case3_mpc_pid.pdf')

#%% MPC vs PID - using control model results

hrz = 10

# Outputs
ydf_pid = pd.read_csv('examples/bs2019/case3/results/pid/ydf.csv',
                      index_col=0)
ydf_pid.index /= 3600.
y_pid = ydf_pid.drop(['vpos', 'dpos'], axis=1)

y_mpc = pd.read_csv('examples/bs2019/case3/results/mpc/h{}/yemu.csv'\
                      .format(hrz), index_col=0)
y_mpc.index /= 3600.

# Control inputs
u_pid = ydf_pid[['vpos', 'dpos']]
y_pid.index /= 3600.

u_mpc = pd.read_csv('examples/bs2019/case3/results/mpc/h{}/u.csv'.format(hrz),
                    index_col=0)
u_mpc.index /= 3600.

# States
x_pid = pd.read_csv('examples/bs2019/case3/results/pid/xdf.csv',
                    index_col=0)
x_pid.index /= 3600.

x_mpc = pd.read_csv('examples/bs2019/case3/results/mpc/h{}/xctr.csv'\
                      .format(hrz), index_col=0)
x_mpc.index /= 3600.

# Constraints
constr = pd.read_csv('examples/bs2019/case3/results/mpc/h{}/constr.csv'\
                     .format(hrz), index_col=0)
constr.index /= 3600.

# Plot
fig, axes = plt.subplots(2, 2, figsize=(5,3), sharex=True)
fig.set_dpi(120)

ax = axes[0][0]
ax.plot(x_pid['cair.T'] - 273.15, label='PID')
ax.plot(x_mpc['x0'] - 273.15, label='MPC')
ax.plot(constr['Tmin'] - 273.15, 'k--', lw=0.5)
ax.plot(constr['Tmax'] - 273.15, 'k--', lw=0.5)
ax.set_ylabel('$T_i$ [$^\circ$C]')
ax.set_yticks(np.arange(19, 26, 2))

ax = axes[0][1]
ax.plot(x_pid['co2.balance.CO2ppmv_i'], label='PID')
ax.plot(x_mpc['x1'], label='MPC')
ax.plot(constr['CO2min'], 'k--', lw=0.5)
ax.plot(constr['CO2max'], 'k--', lw=0.5)
ax.set_ylabel('$C_i$ [ppm]')
ax.set_yticks(np.arange(400, 1001, 200))

ax = axes[1][0]
ax.plot(u_pid['vpos'], label='PID')
ax.plot(u_mpc['vpos'], label='MPC')
ax.set_ylim(0, 100)
ax.set_xlabel('$t$ [h]')
ax.set_ylabel('$v_{p}$ [%]')

ax = axes[1][1]
ax.plot(u_pid['dpos'], label='PID')
ax.plot(u_mpc['dpos'], label='MPC')
ax.set_ylim(0, 100)
ax.set_xlabel('$t$ [h]')
ax.set_ylabel('$d_{p}$ [%]')
ax.set_xticks(np.arange(0, u_pid.index.values[-1] + 1, 24))

fig.tight_layout()

axes[1][0].legend(loc='center', bbox_to_anchor=(1.15, -0.6), ncol=2)

plt.subplots_adjust(0.15, 0.25, 0.99, 0.98)
fig.savefig('examples/bs2019/figs/case3_mpc_pid.pdf')