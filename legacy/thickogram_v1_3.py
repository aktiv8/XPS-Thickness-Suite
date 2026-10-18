# -*- coding: utf-8 -*-
"""
Created on Sun Jun 22 12:38:40 2025

To do...
1. Modify so redline plots correctly
2. Entry for RSFs rather than the ratio hence calculate the ratio = DONE
3. Add the reference Peter J. Cumpson, Surf. Interface Anal. 29, 403-406 (2000).

@author: david
"""

import numpy as np
import matplotlib.pyplot as plt
from tkinter import *
from tkinter import ttk

# === Math Functions ===
def sinh(x):
    return (np.exp(x) - np.exp(-x)) / 2

def func(x):
    return x * ((4.2 - 0.6*x)**0.75 - 0.5)

def fAB(A, B, x):
    return np.log(A/2) - (B**0.75 - 0.5)*x - np.log(sinh(x/2))

def solve_fAB_0(A, B):
    x = 0.001
    deltax = 0.1
    eps = 1.0e-12
    while deltax > eps:
        fx = fAB(A, B, x)
        fxd = fAB(A, B, x + deltax)
        if fx * fxd >= 0.0:
            x += deltax
        else:
            if deltax <= eps:
                x += deltax
            else:
                deltax /= 10
    return x

# === Plotting Function ===
def plot_results(xthicko, Ioisro, Eoes):
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.set_facecolor('white')

    x_vals = np.linspace(0.01, 6, 600)
    y_vals = np.log(sinh(x_vals/2))
    ax.plot(2*x_vals, y_vals, color='green', linewidth=1)

    axp = np.arange(0.1, 6.1, 0.1)
    axd = [0.1, 0.2, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    ayp = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1,
           0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5, 6,
           7, 8, 9, 10, 20, 30]
    bx = np.arange(0.4, 3.2, 0.2)

    for x in axp:
        y = np.log(sinh(x/2))
        ax.plot(2*x, y, 'go', markersize=4)
    for x in axd:
        y = np.log(sinh(x/2))
        ax.plot(2*x, y, 'go', markersize=7)

    for k in range(4, len(ayp)-1):
        b_vals = np.arange(0.4, 3.1, 0.1)
        x_line = (4.2 - b_vals) / 0.6
        y_line = np.log(ayp[k]/2) - func(x_line)
        ax.plot(2*x_line, y_line, color='blue', linewidth=1)

    for b in bx:
        x = (4.2 - b) / 0.6
        y0 = np.log(ayp[4]/2) - func(x)
        y1 = np.log(ayp[-2]/2) - func(x)
        ax.plot([2*x, 2*x], [y0, y1], color='blue', linewidth=1)

    for y in ayp:
        ax.plot([0, 0.18], [np.log(y/2)]*2, color='black', linewidth=1)

    ayd = [0.01, 0.1, 1.0, 10.0]
    for y in ayd:
        ax.text(-0.4, np.log(y/2), f"{y}", verticalalignment='center', color='black')

    for x in axd:
        y = np.log(sinh(x/2))
        ax.text(2*x + 0.15, y - 0.25, f"{x}", color='green')

    for i in range(1, len(ayd)):
        xl = (4.2 - 3)/0.6
        xr = (4.2 - 0.4)/0.6
        yl = np.log(ayd[i]/2) - func(xl)
        yr = np.log(ayd[i]/2) - func(xr)
        ax.text(2*xl - 0.4, yl, f"{ayd[i]:.1f}", color='blue')
        ax.text(2*xr + 0.2, yr, f"{ayd[i]:.1f}", color='blue')

    for i in range(len(bx)):
        if i % 2 == 0:
            x = (4.2 - bx[i]) / 0.6
            y = np.log(ayp[4]/2) - func(x)
            ax.text(2*x + 0.05, y - 0.3, f"{bx[i]:.1f}", color='blue')

    x = xthicko
    y0 = np.log(sinh(x/2)) - 0.3
    y1 = np.log(sinh(x/2)) + 0.3
    ax.plot([2*x, 2*x], [y0, y1], color='red', linewidth=2)
    ax.plot([0, 2*x], [np.log(Ioisro/2), np.log(Ioisro/2) - func((4.2 - Eoes)/0.6)], color='red', linewidth=2)

    ax.text(2*x + 0.6, 2.4, f"{x:.6f}", color='red')

    ax.set_title("Thickogram XPS Thickness Calculation")
    ax.set_xlabel("t/L*cos(theta) (2x)")
    ax.set_ylabel("log(Intensity Ratio)")
    plt.grid(True, linestyle='--', alpha=0.5)
    # Remove outer axes
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

    plt.tight_layout()
    plt.show()
 
# === Calculation Trigger ===
def calculate():
    try:
        Eo = float(entry_vars['Eo'].get())
        Es = float(entry_vars['Es'].get())
        Io = float(entry_vars['Io'].get())
        Is = float(entry_vars['Is'].get())
        So = float(entry_vars['So'].get())
        Ss = float(entry_vars['Ss'].get())
        La = float(entry_vars['Lambda'].get())
        theta = float(entry_vars['Theta'].get())
        Esource = float(entry_vars['Esource'].get())

        Eo = Esource - Eo
        Es = Esource - Es
        Eoes = Eo / Es
        Ro = So / Ss
        Ioisro = Io / (Is * Ro)
        rad = np.radians(theta)
        xthicko = solve_fAB_0(Ioisro, Eoes)
        t = xthicko * La * np.cos(rad)

        result_label.config(text=f"THICKNESS: {t:.2f} nm")

        plot_results(xthicko, Ioisro, Eoes)

    except Exception as e:
        result_label.config(text=f"Error: {e}")

# === Tkinter GUI Setup ===
root = Tk()
root.title("Thickogram Calculator by DaveXPS v0.1")

mainframe = ttk.Frame(root, padding="10")
mainframe.grid(row=0, column=0, sticky=(N, W, E, S))

# Input fields along with instructional text
entry_vars = {}
fields = [
    ("Eo", "Binding Energy of Overlayer (eV):", "     e.g. 707 for Fe 2p"),
    ("Es", "Binding Energy of Substrate (eV):", "     e.g. 100 for Si 2p"),
    ("Io", "Intensity of Overlayer (Io):", "     Peak area as measured in software for overlayer)"),
    ("Is", "Intensity of Substrate (Is):", "     Peak area as measured in software for substrate"),
    ("So", "RSF of Overlayer (So):", "     Sensitivity factor of overlayer peak"),
    ("Ss", "RSF of Substrate (Ss):", "     Sensitivity factor of substrate peak"),
    ("Lambda", "Lambda (Lo) in nm:", "     Attenuation length of photoelectrons (from the overlayer) in the overlayer	"),
    ("Theta", "Theta (degrees):", "     Angle of emission. Use 0 for normal emission"),
    ("Esource", "X-ray Source Energy (eV):", "     e.g., 1486.6")
]

for i, (var, label, tip) in enumerate(fields):
    ttk.Label(mainframe, text=label).grid(row=i, column=0, sticky=W)
    entry_vars[var] = StringVar()
    ttk.Entry(mainframe, textvariable=entry_vars[var], width=15).grid(row=i, column=1)
    ttk.Label(mainframe, text=tip, foreground='green').grid(row=i, column=2, sticky=W)

# Calculate button
ttk.Button(mainframe, text="Calculate & Plot", command=calculate).grid(row=len(fields), column=0, columnspan=3, pady=10)

# Result label
result_label = ttk.Label(mainframe, text="")
result_label.grid(row=len(fields)+1, column=0, columnspan=3)

root.mainloop()
