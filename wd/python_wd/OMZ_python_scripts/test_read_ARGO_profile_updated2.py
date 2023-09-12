# -*- coding: utf-8 -*-
"""
Created on Mon Aug 21 09:13:06 2023

@author: Sergey
"""

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy.ma as ma
import cmocean
from scipy import interpolate
from datetime import datetime, timedelta
import gsw
import csv

# this is a function to fill nan values in 1D array
def fill_nan(A):
    '''
    interpolate to fill nan values
    '''
    inds = np.arange(A.shape[0])
    good = np.where(np.isfinite(A))
    f = interpolate.interp1d(inds[good], A[good],bounds_error=False, fill_value="extrapolate")
    B = np.where(np.isfinite(A),A,f(inds))
    return B

# opening our sprof file
path='D:\\ARGO_workshop\\6901183_Sprof.nc'
fin = Dataset(path)

#getting our data from netcdf
T=np.asarray(fin.variables['TEMP'][:])
P=np.asarray(fin.variables['PRES'][:])
J=np.asarray(fin.variables['JULD'][:])
#CHL=np.asarray(fin.variables['CHLA_ADJUSTED'][:])
CHL=np.asarray(fin.variables['CHLA'][:])
BBP532=np.asarray(fin.variables['BBP532'][:])
BBP700=np.asarray(fin.variables['BBP700'][:])
S=np.asarray(fin.variables['PSAL'][:])

#Here we create a depth array from 1m to 1000m depth
Depth=np.linspace(1,1000,1000)

#we replace all missing values with nans (missing value is 999999 or something like that)
P[P>50000]=np.nan
J[J>50000]=np.nan
T[T>50000]=np.nan
CHL[CHL>50000]=np.nan
BBP532[BBP532>50000]=np.nan
BBP700[BBP700>50000]=np.nan

#this is to convert Julian dayys to datetime
start_datetime = datetime(1950,1,1,0,0,0) 

#here we create empty array to fill them with interpolated values
data_T=[]
data_J=[]
data_S=[]
data_CHL=[]
data_BBP532=[]
data_BBP700=[]

# we interpolate all profiles to Depth array which we created (from 1m to 1000m)
for i in range(len(T)):
    # f = interpolate.interp1d(T[i,:], P[i,:],bounds_error=False)
    f = interpolate.interp1d( P[i,:],T[i,:],bounds_error=False)
    B = f(Depth)
    data_T.append(B)
    f = interpolate.interp1d( P[i,:],S[i,:],bounds_error=False)
    B = f(Depth)
    data_S.append(B)
    data_J.append(start_datetime + timedelta(days = J[i]))
    
    CHL_ar=fill_nan(CHL[i,:])
    
    f = interpolate.interp1d( P[i,:],CHL_ar,bounds_error=False)
    B = f(Depth)
    data_CHL.append(B)
    
    BBP532_ar=fill_nan(BBP532[i,:])
    f = interpolate.interp1d( P[i,:],BBP532_ar,bounds_error=False)
    B = f(Depth)
    data_BBP532.append(B)
    
    BBP700_ar=fill_nan(BBP700[i,:])
    f = interpolate.interp1d( P[i,:],BBP700_ar,bounds_error=False)
    B = f(Depth)
    data_BBP700.append(B)


#trasnpose for plotting
data_CHL=np.transpose(data_CHL)
data_BBP532=np.transpose(data_BBP532)
data_BBP700=np.transpose(data_BBP700)
data_T=np.transpose(np.asarray(data_T))
data_S=np.transpose(np.asarray(data_S))


# from 1d Depth arrray Icreate 2d depth array to match our ARGO data
data_P=np.tile(Depth, (len(T), 1))
data_P=np.transpose(data_P)

#calulating absolute saliniuty conservative temperature and density
SA=gsw.conversions.SA_from_SP(data_S, data_P, -80, 20)
CT=gsw.conversions.CT_from_pt(SA, data_T)
Dens=gsw.density.sigma0(SA, CT)

plt.contourf(data_J,Depth,data_T,10,cmap=cmocean.cm.thermal)
plt.colorbar()
plt.xticks(rotation=70)
plt.gca().invert_yaxis()
plt.plot([data_J[147],data_J[147]],[0,1000], lw=1, c="r", linestyle='--')
plt.plot([data_J[160],data_J[160]],[0,1000], lw=1, c="r", linestyle='--')

fig=plt.figure()

plt.contourf(data_J,Depth,Dens,10,cmap=cmocean.cm.dense)
plt.colorbar()
plt.xticks(rotation=70)
plt.gca().invert_yaxis()
plt.plot([data_J[147],data_J[147]],[0,1000], lw=1, c="r", linestyle='--')
plt.plot([data_J[160],data_J[160]],[0,1000], lw=1, c="r", linestyle='--')

#function to calculate MLD
def mld(SA, CT, p, criterion="pdvar"):
    r"""
    Compute the mixed layer depth.

    Parameters
    ----------
    SA : array_like
         Absolute Salinity  [g/kg]
    CT : array_like
         Conservative Temperature [:math:`^\circ` C (ITS-90)]
    p : array_like
        sea pressure [dbar]
    criterion : str, optional
               MLD Criteria

    Mixed layer depth criteria are:

    'temperature' : Computed based on constant temperature difference
    criterion, CT(0) - T[mld] = 0.5 degree C.

    'density' : computed based on the constant potential density difference
    criterion, pd[0] - pd[mld] = 0.125 in sigma units.

    'pdvar' : computed based on variable potential density criterion
    pd[0] - pd[mld] = var(T[0], S[0]), where var is a variable potential
    density difference which corresponds to constant temperature difference of
    0.5 degree C.

    Returns
    -------
    MLD : array_like
          Mixed layer depth
    idx_mld : bool array
              Boolean array in the shape of p with MLD index.

    References
    ----------
    Monterey, G., and S. Levitus, 1997: Seasonal variability of mixed
    layer depth for the World Ocean. NOAA Atlas, NESDIS 14, 100 pp.
    Washington, D.C.

    """

    SA, CT, p = list(map(np.asanyarray, (SA, CT, p)))
    SA, CT, p = np.broadcast_arrays(SA, CT, p)
    SA, CT, p = list(map(ma.masked_invalid, (SA, CT, p)))

    p_min, idx = p.min(), p.argmin()

    sigma = gsw.rho(SA, CT, p_min) - 1000.0

    # Temperature and Salinity at the surface,
    T0, S0, Sig0 = CT[idx], SA[idx], sigma[idx]

    # NOTE: The temperature difference criterion for MLD
    Tdiff = T0 - 0.5  # 0.8 on the matlab original

    if criterion == "temperature":
        idx_mld = CT > Tdiff
    elif criterion == "pdvar":
        pdvar_diff = gsw.rho(S0, Tdiff, p_min) - 1000.0
        idx_mld = sigma <= pdvar_diff
    elif criterion == "density":
        # sig_diff = Sig0 + 0.125
        sig_diff = Sig0 + 0.03
        idx_mld = sigma <= sig_diff
    else:
        raise NameError(f"Unknown criteria {criterion}")

    MLD = ma.masked_all_like(p)
    MLD[idx_mld] = p[idx_mld]

    return MLD.max(axis=0), idx_mld

data_MLD=[]
for i in range(len(T)):
    M,MLD=mld(SA[:,i], CT[:,i], Depth, criterion="density")
    data_MLD.append(M)

data_MLD=np.asarray(data_MLD)


# here we make a map with hurricane track and argo track
Lat=np.asarray(fin.variables['LATITUDE'][:])
Lon=np.asarray(fin.variables['LONGITUDE'][:])


Lat_H=[]
Lon_H=[]
Time_H=[]
i=0
with open('maria_2017_ibtracs_data.csv', newline='') as f:
    reader = csv.reader(f)
    for row in reader:
        i+=1
        if i==1:continue
        else:
            Lat_H.append(float(row[9]))
            Lon_H.append(float(row[10]))
            Time_H.append(datetime.strptime(row[7], '%Y-%m-%d %H:%M:%S'))

Lat_H=np.asarray(Lat_H)
Lon_H=np.asarray(Lon_H)

fig = plt.figure()


ax = fig.add_axes([0.1,0.1,0.85,0.85])

m = Basemap(llcrnrlon=-71.,llcrnrlat=21.,urcrnrlon=-69.5,urcrnrlat=22.,
            projection='merc', resolution ='i')

m.drawcoastlines(linewidth=1.25)
m.fillcontinents(color='0.8')
m.drawmapboundary(fill_color='0.9')




x, y = m(Lon[150:154], Lat[150:154])


x1, y1 = m(Lon_H, Lat_H)


m.scatter(x, y, 10, marker='o',color='m')
n=data_J[150:154]
for i, txt in enumerate(n):
    ax.annotate(txt, (x[i], y[i]))

m.plot(x1, y1, marker='.',color='r')
m.scatter(x1, y1, marker='o',color='r')

n=Time_H
for i, txt in enumerate(n):
    ax.annotate(txt, (x1[i], y1[i]), color='b')



m.drawparallels(np.arange(-90.,99.,.5))
m.drawmeridians(np.arange(0.,360.,.5))

plt.show()


#plotting the data

fig=plt.figure()

levels=np.linspace(0,0.6, 20)
plt.contourf(data_J,Depth,data_CHL,levels)
plt.colorbar()
plt.xticks(rotation=70)
plt.gca().invert_yaxis()
plt.plot([data_J[147],data_J[147]],[0,1000], lw=1, c="r", linestyle='--')
plt.plot([data_J[160],data_J[160]],[0,1000], lw=1, c="r", linestyle='--')


fig=plt.figure()

levels=np.linspace(0,0.001, 20)

plt.contourf(data_J,Depth,data_BBP700,levels)
plt.colorbar()
plt.xticks(rotation=70)
plt.gca().invert_yaxis()

plt.plot([data_J[147],data_J[147]],[0,1000], lw=1, c="r", linestyle='--')
plt.plot([data_J[160],data_J[160]],[0,1000], lw=1, c="r", linestyle='--')


fig=plt.figure()
levels=np.linspace(0,0.003, 20)
plt.contourf(data_J,Depth,data_BBP532,levels)
plt.colorbar()
plt.xticks(rotation=70)
plt.gca().invert_yaxis()
plt.plot([data_J[147],data_J[147]],[0,1000], lw=1, c="r", linestyle='--')
plt.plot([data_J[160],data_J[160]],[0,1000], lw=1, c="r", linestyle='--')

int_CHL=np.nansum(data_CHL[:250,:],axis=0)

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()


ax1.plot(data_J,data_MLD,c='r')
ax2.plot(data_J,int_CHL,c='b')

ax1.set_ylabel('MLD', color='r')
ax2.set_ylabel('integrated CHL (top250m)', color='b')
ax1.tick_params(axis='y', colors='r')
ax2.tick_params(axis='y', colors='b')

ax1.plot([data_J[147],data_J[147]],[0,140], lw=1, c="k", linestyle='--')
ax1.plot([data_J[160],data_J[160]],[0,140], lw=1, c="k", linestyle='--')
plt.xticks(rotation=70)

# lower= 147
# upper= 160


fig=plt.figure()
levels=np.linspace(0,0.6, 20)
plt.contourf(data_J[147:160],Depth[:300],data_CHL[0:300,147:160],levels)
plt.colorbar()
plt.xticks(rotation=70)
plt.gca().invert_yaxis()

plt.xticks(rotation=70)

fig=plt.figure()
levels=np.linspace(0,0.003, 20)
plt.contourf(data_J[147:160],Depth[:300],data_BBP532[0:300,147:160],levels)
plt.colorbar()
plt.xticks(rotation=70)
plt.gca().invert_yaxis()
plt.xticks(rotation=70)


fig=plt.figure()
levels=np.linspace(0,0.001, 20)
plt.contourf(data_J[147:160],Depth[:300],data_BBP700[0:300,147:160],levels)
plt.colorbar()
plt.xticks(rotation=70)
plt.gca().invert_yaxis()

plt.xticks(rotation=70)


fig=plt.figure()

plt.contourf(data_J[147:160],Depth[:300],data_T[0:300,147:160],10, cmap=cmocean.cm.thermal)
plt.colorbar()
plt.xticks(rotation=70)
plt.gca().invert_yaxis()

plt.xticks(rotation=70)


fig=plt.figure()

plt.contourf(data_J[147:160],Depth[:300],Dens[0:300,147:160],20, cmap=cmocean.cm.dense)
plt.colorbar()
plt.xticks(rotation=70)
plt.gca().invert_yaxis()

plt.xticks(rotation=70)



fig, ax1 = plt.subplots()
ax2 = ax1.twinx()


ax1.plot(data_J[147:160],data_MLD[147:160],c='r')
ax2.plot(data_J[147:160],int_CHL[147:160],c='b')
plt.xticks(rotation=70)

ax1.set_ylabel('MLD', color='r')
ax2.set_ylabel('integrated CHL (top250m)', color='b')
ax1.tick_params(axis='y', colors='r')
ax2.tick_params(axis='y', colors='b')


def runningMeanFast(x, N):
    return np.convolve(x, np.ones((N,))/N)[(N-1):]



f, axs = plt.subplots(1, 3, figsize=(18, 10))

#plt.subplot(1,3,1)

axs[0].plot(Dens[0:300,150],Depth[:300], 'b--')
axs[0].plot(Dens[0:300,151],Depth[:300], 'b')
axs[0].plot(Dens[0:300,152],Depth[:300], 'r')
axs[0].plot(Dens[0:300,153],Depth[:300], 'r--')
axs[0].set_ylim([0,300])
axs[0].invert_yaxis()
axs[0].set_xlabel('Potential Density (kg m^-^3)')
axs[0].set_ylabel('Depth (m)')

BBP700_pr1=runningMeanFast(data_BBP700[:,150],10)

axs[1].plot(BBP700_pr1[0:300],Depth[:300], 'b--')

BBP700_pr2=runningMeanFast(data_BBP700[:,151],10)

axs[1].plot(BBP700_pr2[0:300],Depth[:300], 'b')

BBP700_pr3=runningMeanFast(data_BBP700[:,152],10)

axs[1].plot(BBP700_pr3[0:300],Depth[:300], 'r')

BBP700_pr4=runningMeanFast(data_BBP700[:,153],10)

axs[1].plot(BBP700_pr4[0:300],Depth[:300], 'r--')
axs[1].set_ylim([0,300])
axs[1].set_xlabel('BBP700 (m^-^1))')
axs[1].invert_yaxis()


CHL_pr1=runningMeanFast(data_CHL[:,150],10)

axs[2].plot(CHL_pr1[0:300],Depth[:300], 'b--')

CHL_pr2=runningMeanFast(data_CHL[:,151],10)

axs[2].plot(CHL_pr2[0:300],Depth[:300], 'b')

CHL_pr3=runningMeanFast(data_CHL[:,152],10)

axs[2].plot(CHL_pr3[0:300],Depth[:300], 'r')

CHL_pr4=runningMeanFast(data_CHL[:,153],10)

axs[2].plot(CHL_pr4[0:300],Depth[:300], 'r--')

axs[2].set_ylim([0,300])

axs[2].set_xlabel('Chlorophyll (mg m^-^3))')

plt.gca().legend((data_J[150],data_J[151],data_J[152],data_J[153]))
axs[2].invert_yaxis()
