#!/usr/bin/env python3.8
"""
Created on 29/10/2020
@author: Jiacheng Wu, jcwu@pku.edu.cn
"""
import numpy as np

# Common parameters of Earth and Math
pi = np.pi                  # pi =3.14...
re = 6.371e6                # Radius of Earth, unit: m
year = 3600 * 24 * 365      # one year on Earth

# atmospheric parameters
cp_a = 1004.0               # The special heat capacity of air, unit: J/kg/K
c = 5300.0 * cp_a
"""
Radiative flux: R = S_short - IR
Shortwave radiation: S = Sol * (1 - albedo)
Longwave radiation: IR = A + B*Ta
"""
Sol = np.array([320.0, 390.0, 270.0])       # Solar insolation, unit: W/m**2
albedo = np.array([0.4, 0.25, 0.42])        # Albedo
S_short = Sol * (1.0 - albedo)              # Shortwave radiation, unit: W/m**2

A, B = 213.35, 2.22                         # IR = A + B*T (Budyko 1969)

SurfFrac = np.array([0.16267, 0.22222, 0.21069])    # Surface fraction

Ks = 2.5*1.0e13             # parameter for sensible heat flux
Kl = 1.5*5.1e17             # parameter for latent heat flux
RH = 0.8                    # Global-mean relative humidity
lat_lsa = -90               # lower latitude of southern Atmosphere, unit: degree
lat_usa = -30               # upper latitude of southern Atmosphere, unit: degree
lat_lna = 45                # lower latitude of northern Atmosphere, unit: degree
lat_una = 90                # upper latitude of northern Atmosphere, unit: degree

clat_s = np.rad2deg(np.arcsin((np.sin(np.deg2rad(lat_lsa)) + np.sin(np.deg2rad(lat_usa))) / 2))
clat_m = np.rad2deg(np.arcsin((np.sin(np.deg2rad(lat_usa)) + np.sin(np.deg2rad(lat_lna))) / 2))
clat_n = np.rad2deg(np.arcsin((np.sin(np.deg2rad(lat_lna)) + np.sin(np.deg2rad(lat_una))) / 2))
# centroid latitude of southern, tropic, and northern atmosphere

dy_s = re * (np.deg2rad(clat_m) - np.deg2rad(clat_s))
dy_n = re * (np.deg2rad(clat_n) - np.deg2rad(clat_m))

# oceanic parameters
phi = 1.5264e10            # parameter for Phi
m, n = 8.0e-4, 1.5e-4      # Phi=phi*(m*dS-n*dT), dS and dT is the different of salinity and temperature, respectively
rho_sw = 1025              # The density of sea water, unit: kg/m**3
cp_o = 4200.0              # The special heat capacity of ocean, unit: J/kg/K

lat_lso = -60               # lower latitude of southern Atlantic, unit: degree
lat_uso = -30               # upper latitude of southern Atlantic, unit: degree
lat_lno = 45                # lower latitude of northern Atlantic, unit: degree
lat_uno = 80                # upper latitude of northern Atlantic, unit: degree

lon_atl = 80                # longitude of Atlantic, unit: degree

"""
Math: delta_Area = dx*dy= (r*cos(lat)*delta_lon)*(r*delta_lat) 
"""
area_s = (re**2)*(np.abs(np.sin(np.deg2rad(lat_uso))-np.sin(np.deg2rad(lat_lso))))*(2*pi*lon_atl/360)
area_m = (re**2)*(np.abs(np.sin(np.deg2rad(lat_lno))-np.sin(np.deg2rad(lat_uso))))*(2*pi*lon_atl/360)
area_n = (re**2)*(np.abs(np.sin(np.deg2rad(lat_uno))-np.sin(np.deg2rad(lat_lno))))*(2*pi*lon_atl/360)
area_d = area_m
# Area of southern, mixing layer in tropical, northern, and deep layer in tropical Atlantic

dz1 = 600.0                # The depth of mixing layer, unit: m
dz2 = 4000.0                # The depth of ocean, unit: m
dz = dz2 - dz1             # The depth of deep layer, unit: m

vs = dz2 * area_s
vm = dz1 * area_m
vn = dz2 * area_n
vd = dz * area_d
# Volume of southern, mixing layer in tropical, northern, and deep layer in tropical Atlantic

"""
Heat flux: H = Q1-Q2*(To-Ta)
P-E: PE = (2pi*lon_atl/360)*re*Lr*cos(lat_u)*Fl = pe * Fl
"""
Q1 = np.array([10.0, 70.0, 20.0])
Q2 = np.array([50.0, 50.0, 50.0])

S_ref = 34.9                    # The global-mean salinity
rho_fw = 1000.0                 # The density of freshwater, unit: kg/m**3
L_r = 1 / (rho_fw * 2.5e6)      # The reverse of latent heat of freshwater

pes = (2*pi*lon_atl/360) * re * L_r * np.cos(np.deg2rad(lat_uso))       # parameter for PES
pen = (2*pi*2.5*lon_atl/360) * re * L_r * np.cos(np.deg2rad(lat_lno))   # parameter for PEN
