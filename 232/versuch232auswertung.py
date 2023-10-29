#Versuch 232 Auswertung

# Aufgabe a: Geradenfit            U_I Diagramm

import easyparse
import numpy as np
from scipy.stats import linregress
from matplotlib import pyplot as plt

def main():
    data = easyparse.parse("232_a_messwerte.csv")
    u = data['U']
    i = data['I\n']
    result = linregress(i, u)
    m = result.slope
    n = result.intercept
    dm = result.stderr
    dn = result.intercept_stderr
    print(m, n, dm, dn)

    x = np.linspace(0, 75, 75)
    plt.errorbar(i, u, xerr = 3, yerr = 0.05, linestyle ='', label = 'Messwerte', color = 'blue')
    plt.plot(x, m*x+n, label = 'Geraden fit', color = 'blue')

    rx = 21.55*10**-3
    plt.plot(x, rx*x, label = 'Rx-Gerade', color = 'orange')

    plt.legend()
    plt.grid()

    plt.show()

if __name__ == '__main__': 
    main()