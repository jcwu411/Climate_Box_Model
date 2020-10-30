#!/usr/bin/env python3.8
"""
Created on 29/10/2020
@author: Jiacheng Wu, jcwu@pku.edu.cn
"""

import numpy as np
import constant_par as con
import matplotlib.pyplot as plt


def plot_2d_ocean_temperature(data, dt, steps, lw=1.0, ti="Plot", fn="test2d.pdf", sa=True):
    tso, tmo, tno, tdo = data[:, 0], data[:, 1], data[:, 2], data[:, 3]

    t = np.linspace(0, steps * dt, steps)
    plt.figure()
    plt.plot(t, tso, color='red', linewidth=lw, label='S')
    plt.plot(t, tmo, color='blue', linewidth=lw, label='M')
    plt.plot(t, tno, color='green', linewidth=lw, label='N')
    plt.plot(t, tdo, color='black', linewidth=lw, label='D')
    plt.ylabel('Oceanic temperature')
    plt.xlabel('t (a)')
    plt.legend(loc='lower right')
    ax = plt.gca()
    # ax.yaxis.set_major_locator(MultipleLocator(10))
    ax.tick_params(axis='y', labelsize=10)
    ax.set_xlim(0, steps * dt)

    if sa:
        plt.savefig(fn)
    plt.show()
    plt.close()


def plot_2d_ocean_salinity(data, dt, steps, lw=1.0, ti="Plot", fn="test2d.pdf", sa=True):
    sso, smo, sno, sdo = data[:, 4], data[:, 5], data[:, 6], data[:, 7]

    t = np.linspace(0, steps * dt, steps)
    plt.figure()
    plt.plot(t, sso, color='red', linewidth=lw, label='S')
    plt.plot(t, smo, color='blue', linewidth=lw, label='M')
    plt.plot(t, sno, color='green', linewidth=lw, label='N')
    plt.plot(t, sdo, color='black', linewidth=lw, label='D')
    plt.ylabel('Salinity')
    plt.xlabel('t (a)')
    plt.legend(loc='lower right')
    ax = plt.gca()
    ax.tick_params(axis='y', labelsize=10)
    ax.set_xlim(0, steps * dt)

    if sa:
        plt.savefig(fn)
    plt.show()
    plt.close()


def plot_2d_atmosphere_temperature(data, dt, steps, lw=1.0, ti="Plot", fn="test2d.pdf", sa=True):
    tsa, tma, tna = data[:, 8], data[:, 9], data[:, 10]
    tmean = (0.5 * tsa + 1.207 * tma + 0.293 * tna) / 2
    t = np.linspace(0, steps * dt, steps)
    plt.figure()
    plt.plot(t, tsa, color='red', linewidth=lw, label='S')
    plt.plot(t, tma, color='blue', linewidth=lw, label='M')
    plt.plot(t, tna, color='green', linewidth=lw, label='N')
    plt.plot(t, tmean, color='black', linewidth=lw, label='G')
    plt.ylabel('Atmosphere temperature')
    plt.xlabel('t (a)')
    plt.legend(loc='lower right')
    ax = plt.gca()
    ax.tick_params(axis='y', labelsize=10)
    ax.set_xlim(0, steps * dt)

    if sa:
        plt.savefig(fn)
    plt.show()
    plt.close()


def plot_2d_phi(data, dt, steps, lw=1.0, ti="Plot", fn="test2d.pdf", sa=True):
    tso, tno, sso, sno = data[:, 0], data[:, 2], data[:, 4],  data[:, 6]
    phi = con.phi * (con.m * (sno - sso) - con.n * (tno - tso)) * 1e-6
    t = np.linspace(0, steps * dt, steps)
    plt.figure()
    plt.plot(t, phi, color='blue', linewidth=lw)
    plt.ylabel('Phi')
    plt.xlabel('t (a)')
    ax = plt.gca()
    ax.tick_params(axis='y', labelsize=10)
    ax.set_xlim(0, steps * dt)

    if sa:
        plt.savefig(fn)
    plt.show()
    plt.close()
