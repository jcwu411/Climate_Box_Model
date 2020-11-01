#!/usr/bin/env python3.8
"""
Created on 29/10/2020
@author: Jiacheng Wu, jcwu@pku.edu.cn
"""

import numpy as np
import constant_par as con
import matplotlib.pyplot as plt

ylabelsize = 8


def plot_2d_phi(data, dt, steps, lw=1.0, ti="Plot", fn="test2d.pdf", sa=True):
    tso, tno, sso, sno = data[:, 0], data[:, 2], data[:, 4],  data[:, 6]
    phi = con.phi * (con.m * (sno - sso) - con.n * (tno - tso)) * 1e-6
    t = np.linspace(0, steps * dt, steps)
    plt.figure()
    color = 'blue'  # The color of the line
    plt.plot(t, phi, color=color, linewidth=lw)
    plt.ylabel(r'$\Phi\ (Sv)$')
    plt.xlabel('t (a)')
    ax = plt.gca()
    ax.tick_params(axis='y', labelsize=ylabelsize)
    ax.set_xlim(0, steps * dt)
    # Mark the final value
    axr = plt.twinx()
    axr.plot(t, phi, color=color, linewidth=lw * 0.1)
    axr.set_yticks([phi[len(t) - 1]])
    axr.set_yticklabels([str('%.2f' % phi[len(t) - 1])], color=color)

    if sa:
        plt.savefig(fn)
    plt.show()
    plt.close()

# ---------------------------------------
# Subplot
# ---------------------------------------


def plot_2d_oce_temperature_sub(data, dt, steps, lw=1.0, ti="Plot", fn="test2d.pdf", sa=True):
    # Read data
    tso, tmo, tno, tdo = data[:, 0], data[:, 1], data[:, 2], data[:, 3]
    # Creat time series
    t = np.linspace(0, steps * dt, steps)
    # Plot
    plt.figure()
    plt.subplot(411)
    color = 'red'
    plt.plot(t, tso, color=color, linewidth=lw)
    plt.ylabel(r'$T_S^o$')
    ax = plt.gca()
    ax.tick_params(axis='y', labelsize=ylabelsize)
    ax.set_xlim(0, steps * dt)
    ax.set_xticklabels([])
    # Mark the final value
    axr = plt.twinx()
    axr.plot(t, tso, color=color, linewidth=lw * 0.1)
    axr.set_yticks([tso[len(t) - 1]])
    axr.set_yticklabels([str('%.2f' % tso[len(t) - 1])], color=color)
    axr.tick_params(axis='y', labelsize=ylabelsize)

    plt.subplot(412)
    color = 'blue'  # The color of this plot
    plt.plot(t, tmo, color=color, linewidth=lw)
    plt.ylabel(r'$T_M^o$')
    ax = plt.gca()
    ax.tick_params(axis='y', labelsize=ylabelsize)
    ax.set_xlim(0, steps * dt)
    ax.set_xticklabels([])
    # Mark the final value
    axr = plt.twinx()
    axr.plot(t, tmo, color=color, linewidth=lw * 0.1)
    axr.set_yticks([tmo[len(t) - 1]])
    axr.set_yticklabels([str('%.2f' % tmo[len(t) - 1])], color=color)
    axr.tick_params(axis='y', labelsize=ylabelsize)

    plt.subplot(413)
    color = 'green'
    plt.plot(t, tno, color=color, linewidth=lw)
    plt.ylabel(r'$T_N^o$')
    ax = plt.gca()
    ax.tick_params(axis='y', labelsize=ylabelsize)
    ax.set_xlim(0, steps * dt)
    ax.set_xticklabels([])
    # Mark the final value
    axr = plt.twinx()
    axr.plot(t, tno, color=color, linewidth=lw * 0.1)
    axr.set_yticks([tno[len(t) - 1]])
    axr.set_yticklabels([str('%.2f' % tno[len(t) - 1])], color=color)
    axr.tick_params(axis='y', labelsize=ylabelsize)

    plt.subplot(414)
    color = 'black'
    plt.plot(t, tdo, color=color, linewidth=lw)
    plt.ylabel(r'$T_D^o$')
    plt.xlabel('t (a)')
    ax = plt.gca()
    ax.tick_params(axis='y', labelsize=ylabelsize)
    ax.set_xlim(0, steps * dt)
    # Mark the final value
    axr = plt.twinx()
    axr.plot(t, tdo, color=color, linewidth=lw * 0.1)
    axr.set_yticks([tdo[len(t) - 1]])
    axr.set_yticklabels([str('%.2f' % tdo[len(t) - 1])], color=color)
    axr.tick_params(axis='y', labelsize=ylabelsize)

    plt.suptitle(ti + r'$({}^{\circ}C)$')

    if sa:
        plt.savefig(fn)
    plt.show()
    plt.close()


def plot_2d_oce_salinity_sub(data, dt, steps, lw=1.0, ti="Plot", fn="test2d.pdf", sa=True):
    # Read data
    sso, smo, sno, sdo = data[:, 4], data[:, 5], data[:, 6], data[:, 7]
    # Creat time series
    t = np.linspace(0, steps * dt, steps)
    # Plot
    plt.figure()
    plt.subplot(411)
    color = 'red'
    plt.plot(t, sso, color=color, linewidth=lw)
    plt.ylabel(r'$S_S^o$')
    ax = plt.gca()
    ax.tick_params(axis='y', labelsize=ylabelsize)
    ax.set_xlim(0, steps * dt)
    ax.set_xticklabels([])
    # Mark the final value
    axr = plt.twinx()
    axr.plot(t, sso, color=color, linewidth=lw * 0.1)
    axr.set_yticks([sso[len(t) - 1]])
    axr.set_yticklabels([str('%.2f' % sso[len(t) - 1])], color=color)
    axr.tick_params(axis='y', labelsize=ylabelsize)

    plt.subplot(412)
    color = 'blue'
    plt.plot(t, smo, color=color, linewidth=lw)
    plt.ylabel(r'$S_M^o$')
    ax = plt.gca()
    ax.tick_params(axis='y', labelsize=ylabelsize)
    ax.set_xlim(0, steps * dt)
    ax.set_xticklabels([])
    # Mark the final value
    axr = plt.twinx()
    axr.plot(t, smo, color=color, linewidth=lw * 0.1)
    axr.set_yticks([smo[len(t) - 1]])
    axr.set_yticklabels([str('%.2f' % smo[len(t) - 1])], color=color)
    axr.tick_params(axis='y', labelsize=ylabelsize)

    plt.subplot(413)
    color = 'green'
    plt.plot(t, sno, color=color, linewidth=lw)
    plt.ylabel(r'$S_N^o$')
    ax = plt.gca()
    ax.tick_params(axis='y', labelsize=ylabelsize)
    ax.set_xlim(0, steps * dt)
    ax.set_xticklabels([])
    # Mark the final value
    axr = plt.twinx()
    axr.plot(t, sno, color=color, linewidth=lw * 0.1)
    axr.set_yticks([sno[len(t) - 1]])
    axr.set_yticklabels([str('%.2f' % sno[len(t) - 1])], color=color)
    axr.tick_params(axis='y', labelsize=ylabelsize)

    plt.subplot(414)
    color = 'black'
    plt.plot(t, sdo, color=color, linewidth=lw)
    plt.ylabel(r'$S_D^o$')
    plt.xlabel('t (a)')
    ax = plt.gca()
    ax.tick_params(axis='y', labelsize=ylabelsize)
    ax.set_xlim(0, steps * dt)
    # Mark the final value
    axr = plt.twinx()
    axr.plot(t, sdo, color=color, linewidth=lw * 0.1)
    axr.set_yticks([sdo[len(t) - 1]])
    axr.set_yticklabels([str('%.2f' % sdo[len(t) - 1])], color=color)
    axr.tick_params(axis='y', labelsize=ylabelsize)

    plt.suptitle(ti + "(PSU)")

    if sa:
        plt.savefig(fn)
    plt.show()
    plt.close()


def plot_2d_atmos_temperature_sub(data, dt, steps, lw=1.0, ti="Plot", fn="test2d.pdf", sa=True):
    # Read data
    tsa, tma, tna = data[:, 8], data[:, 9], data[:, 10]
    # Calculate Global Atmospheric Temperature
    tmean = (0.5 * tsa + 1.207 * tma + 0.293 * tna) / 2
    # Creat time series
    t = np.linspace(0, steps * dt, steps)
    # Plot
    plt.figure()
    plt.subplot(411)
    color = 'red'
    plt.plot(t, tsa, color=color, linewidth=lw)
    plt.ylabel(r'$T_S^a$')
    ax = plt.gca()
    ax.tick_params(axis='y', labelsize=ylabelsize)
    ax.set_xlim(0, steps * dt)
    ax.set_xticklabels([])
    # Mark the final value
    axr = plt.twinx()
    axr.plot(t, tsa, color=color, linewidth=lw * 0.1)
    axr.set_yticks([tsa[len(t) - 1]])
    axr.set_yticklabels([str('%.2f' % tsa[len(t)-1])], color=color)
    axr.tick_params(axis='y', labelsize=ylabelsize)

    plt.subplot(412)
    color = 'blue'
    plt.plot(t, tma, color=color, linewidth=lw)
    plt.ylabel(r'$T_M^a$')
    ax = plt.gca()
    ax.tick_params(axis='y', labelsize=ylabelsize)
    ax.set_xlim(0, steps * dt)
    ax.set_xticklabels([])
    # Mark the final value
    axr = plt.twinx()
    axr.plot(t, tma, color=color, linewidth=lw * 0.1)
    axr.set_yticks([tma[len(t) - 1]])
    axr.set_yticklabels([str('%.2f' % tma[len(t) - 1])], color=color)
    axr.tick_params(axis='y', labelsize=ylabelsize)

    plt.subplot(413)
    color = 'green'
    plt.plot(t, tna, color=color, linewidth=lw)
    plt.ylabel(r'$T_N^a$')
    ax = plt.gca()
    ax.tick_params(axis='y', labelsize=ylabelsize)
    ax.set_xlim(0, steps * dt)
    ax.set_xticklabels([])
    # Mark the final value
    axr = plt.twinx()
    axr.plot(t, tna, color=color, linewidth=lw * 0.1)
    axr.set_yticks([tna[len(t) - 1]])
    axr.set_yticklabels([str('%.2f' % tna[len(t) - 1])], color=color)
    axr.tick_params(axis='y', labelsize=ylabelsize)

    plt.subplot(414)
    color = 'black'
    plt.plot(t, tmean, color=color, linewidth=lw)
    plt.ylabel(r'$T_{mean}$')
    plt.xlabel('t (a)')
    ax = plt.gca()
    ax.tick_params(axis='y', labelsize=ylabelsize)
    ax.set_xlim(0, steps * dt)
    # Mark the final value
    axr = plt.twinx()
    axr.plot(t, tmean, color=color, linewidth=lw * 0.1)
    axr.set_yticks([tmean[len(t) - 1]])
    axr.set_yticklabels([str('%.2f' % tmean[len(t) - 1])], color=color)
    axr.tick_params(axis='y', labelsize=ylabelsize)

    plt.suptitle(ti + r'$({}^{\circ}C)$')

    if sa:
        plt.savefig(fn)
    plt.show()
    plt.close()

# ---------------------------------------
# Testing plot function
# ---------------------------------------


def plot_2d_oce_temperature(data, dt, steps, lw=1.0, ti="Plot", fn="test2d.pdf", sa=True):
    # Read data
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
    ax.tick_params(axis='y', labelsize=10)
    ax.set_xlim(0, steps * dt)

    if sa:
        plt.savefig(fn)
    plt.show()
    plt.close()


def plot_2d_oce_salinity(data, dt, steps, lw=1.0, ti="Plot", fn="test2d.pdf", sa=True):
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


def plot_2d_atmos_temperature(data, dt, steps, lw=1.0, ti="Plot", fn="test2d.pdf", sa=True):
    tsa, tma, tna = data[:, 8], data[:, 9], data[:, 10]
    tmean = (0.5 * tsa + 1.207 * tma + 0.293 * tna) / 2
    t = np.linspace(0, steps * dt, steps)
    plt.figure()
    plt.plot(t, tsa, color='red', linewidth=lw, label='S')
    plt.plot(t, tma, color='blue', linewidth=lw, label='M')
    plt.plot(t, tna, color='green', linewidth=lw, label='N')
    plt.plot(t, tmean, color='black', linewidth=lw, label='G')
    plt.ylabel('Atmospheric temperature')
    plt.xlabel('t (a)')
    plt.legend(loc='lower right')
    ax = plt.gca()
    ax.tick_params(axis='y', labelsize=10)
    ax.set_xlim(0, steps * dt)

    if sa:
        plt.savefig(fn)
    plt.show()
    plt.close()

