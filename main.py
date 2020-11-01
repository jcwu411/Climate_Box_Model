#!/usr/bin/env python3.8
"""
Created on 29/10/2020
@author: Jiacheng Wu, jcwu@pku.edu.cn
"""

import numpy as np
import ode_solver  # solve ode
import constant_par as con  # input constant
import plot_data  # input plot function

"""
    tso: The temperature of Southern Atlantic
    tmo: The temperature of Mixing layer in the tropic Atlantic 
    tno: The temperature of Northern Atlantic
    tdo: The temperature of Deep layer in the tropical Atlantic
    sso: The salinity of Southern Atlantic
    smo: The salinity of Mixing layer in the tropical Atlantic 
    sno: The salinity of Northern Atlantic
    sdo: The salinity of Deep layer in the tropical Atlantic
    tsa: The temperature of Southern Atmosphere
    tma: The temperature of tropical Atmosphere 
    tna: The temperature of Northern Atmosphere
"""


def boxmodel(x_0):
    tso, tmo, tno, tdo, sso, smo, sno, sdo, tsa, tma, tna = x_0
    # atmospheric term
    """
    Radiative flux: R = S_short - IR
    Longwave radiation: IR = A + B*Ta
    """
    irs = con.A + con.B * tsa
    irm = con.A + con.B * tma
    irn = con.A + con.B * tna

    rs = con.S_short[0] - irs
    rm = con.S_short[1] - irm
    rn = con.S_short[2] - irn

    """
    Heat flux in Atlantic basin: H = Q1-Q2*(To-Ta)
    Heat flux of the atmosphere and ocean interface: h = surfFrac X H
    """
    hs = con.Q1[0] - con.Q2[0] * (tso - tsa)
    hm = con.Q1[1] - con.Q2[1] * (tmo - tma)
    hn = con.Q1[2] - con.Q2[2] * (tno - tna)
    his = con.SurfFrac[0] * hs
    him = con.SurfFrac[1] * hm
    hin = con.SurfFrac[2] * hn

    """
    Meridional Heat Flux: F = Fs + Fl
    Meridional Sensible Heat Flux: Fs = Ks*dTa/dy
    Meridional Latent Heat Flux: Fl = Kl*RH*dQSdT/dy
    """
    fss = con.Ks * (tma - tsa) / con.dy_s
    fns = -1 * con.Ks * (tna - tma) / con.dy_n

    tc_sm = (tsa * (con.clat_m - con.lat_usa) + tma * (con.lat_usa - con.clat_s)) / (con.clat_m - con.clat_s)
    tc_mn = (tma * (con.clat_n - con.lat_lna)
             + tna * (con.lat_lna - con.clat_m)) / (con.clat_n - con.clat_m)

    sat_s = 6.112 * np.exp(17.67 * tc_sm / (tc_sm + 243.5))
    sat_n = 6.112 * np.exp(17.67 * tc_mn / (tc_mn + 243.5))

    dqsdt_s = 243.5 * 17.67 * (0.622 * sat_s / 1000 / ((tc_sm + 243.5) ** 2))
    dqsdt_n = 243.5 * 17.67 * (0.622 * sat_n / 1000 / ((tc_mn + 243.5) ** 2))

    fsl = con.Kl * (con.RH * dqsdt_s) / con.dy_s
    fnl = con.Kl * (con.RH * dqsdt_n) / con.dy_n

    fs = fss + fsl
    fn = fns + fnl

    """
    dF/dy = (-F1*cos(lat1)+(-F2)*cos(lat2))/(re*(sin(lat2)-sin(lat1)))
    """
    dfsdy = fs * np.cos(np.deg2rad(con.lat_usa)) / (con.re * (np.sin(np.deg2rad(con.lat_usa))
                                                              - np.sin(np.deg2rad(con.lat_lsa))))
    dfmdy = (-1 * fn * np.cos(np.deg2rad(con.lat_lna))
             + -1 * fs * np.cos(np.deg2rad(con.lat_usa))) / (con.re * (np.sin(np.deg2rad(con.lat_lna))
                                                                       - np.sin(np.deg2rad(con.lat_usa))))
    dfndy = fn * np.cos(np.deg2rad(con.lat_lna)) / (con.re * (np.sin(np.deg2rad(con.lat_una))
                                                              - np.sin(np.deg2rad(con.lat_lna))))

    # oceanic term
    phi = con.phi * (con.m * (sno - sso) - con.n * (tno - tso))  # phi must >=0
    if phi < 0:
        phi = 0

    """
    P-E: (2pi*lon_atl/360) * re * Lr * cos(lat_u) * Fl
    """
    pes = con.pes * fsl
    pen = con.pen * fnl
    pem = pes + pen

    # tendency
    fx = [-1 * (tso - tdo) * phi / con.vs + hs / (con.rho_sw * con.cp_o * con.dz2),
          -1 * (tmo - tso) * phi / con.vm + hm / (con.rho_sw * con.cp_o * con.dz1),
          -1 * (tno - tmo) * phi / con.vn + hn / (con.rho_sw * con.cp_o * con.dz2),
          -1 * (tdo - tno) * phi / con.vd,
          -1 * (sso - sdo) * phi / con.vs - con.S_ref * pes / con.vs,
          -1 * (smo - sso) * phi / con.vm + con.S_ref * pem / con.vm,
          -1 * (sno - smo) * phi / con.vn - con.S_ref * pen / con.vn,
          -1 * (sdo - sno) * phi / con.vd,
          1 / con.c * (dfsdy + rs - his),
          1 / con.c * (dfmdy + rm - him),
          1 / con.c * (dfndy + rn - hin)]
    return np.array(fx)


if __name__ == "__main__":
    debug = False  # debug
    freshwater = True  # Freshwater experiment
    plot = True  # plot or not
    save = True  # save file or not

    # Initialization
    plot_dt = 0.01
    dt = plot_dt * con.year
    cal_year = 3000
    fresh_year = 5000
    steps_fresh = 100 * fresh_year

    if debug:
        steps = 0
    else:
        steps = 100 * cal_year
    scheme = "rk4"

    x0 = [4.777404031, 24.42876625, 2.66810894, 2.67598915,
          34.40753555, 35.62585068, 34.92513657, 34.91130066,
          4.67439556, 23.30437851, 0.94061828]

    if freshwater:
        result_f = ode_solver.Ode(iv=x0, function=boxmodel, dt=dt, steps=steps_fresh, debug=1)
        result_f.integrate(scheme)
        # record the final state, and inject freshwater
        xf0 = result_f.trajectory[len(result_f.trajectory) - 1, :]
        xf0[6] = xf0[6] - 0.7  # change the sno by reducing 0.7
        # rerun for balance
        result_b = ode_solver.Ode(iv=xf0, function=boxmodel, dt=dt, steps=steps, debug=1)
        result_b.integrate(scheme)
        # Data merging
        data = np.concatenate((result_f.trajectory, result_b.trajectory), axis=0)
    else:
        result = ode_solver.Ode(iv=x0, function=boxmodel, dt=dt, steps=steps, debug=1)
        result.integrate(scheme)
        data = result.trajectory

    # ---------------------------------------
    # Plotting
    # ---------------------------------------

    if plot:
        if freshwater:
            plot_data.plot_2d_atmos_temperature_sub(data, dt=plot_dt, steps=steps_fresh + 1 + steps + 1,
                                                    lw=1.0, ti="Atmospheric Temperature", fn="AT_fresh.pdf", sa=save)
            plot_data.plot_2d_oce_temperature_sub(data, dt=plot_dt, steps=steps_fresh + 1 + steps + 1,
                                                  lw=1.0, ti="Oceanic Temperature", fn="OT_fresh.pdf", sa=save)
            plot_data.plot_2d_oce_salinity_sub(data, dt=plot_dt, steps=steps_fresh + 1 + steps + 1,
                                               lw=1.0, ti="Oceanic Salinity", fn="OS_fresh.pdf", sa=save)
            plot_data.plot_2d_phi(data, dt=plot_dt, steps=steps_fresh + 1 + steps + 1,
                                  lw=1.0, ti="Phi", fn="Phi_fresh.pdf", sa=save)
        else:
            plot_data.plot_2d_atmos_temperature_sub(data, dt=plot_dt, steps=steps + 1,
                                                    lw=1.0, ti="Atmospheric Temperature", fn="AT.pdf", sa=save)
            plot_data.plot_2d_oce_temperature_sub(data, dt=plot_dt, steps=steps + 1,
                                                  lw=1.0, ti="Oceanic Temperature", fn="OT.pdf", sa=save)
            plot_data.plot_2d_oce_salinity_sub(data, dt=plot_dt, steps=steps + 1,
                                               lw=1.0, ti="Oceanic Salinity", fn="OS.pdf", sa=save)
            plot_data.plot_2d_phi(data, dt=plot_dt, steps=steps + 1,
                                  lw=1.0, ti="Phi", fn="Phi.pdf", sa=save)
