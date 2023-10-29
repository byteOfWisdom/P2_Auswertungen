# Aufgabe 232.e 

import easyparse
import numpy as np
from scipy.stats import linregress
from matplotlib import pyplot as plt

def main():

    # Messwerte in arrays umwandeln
    data = easyparse.parse("232_messwerte_e.csv")
    u = data['U']
    i = data['I\n']
    result = linregress(i, u)
    m = result.slope
    n = result.intercept
    dm = result.stderr
    dn = result.intercept_stderr
    print(m, n, dm, dn)

    plt.xlabel('I [mA]')
    plt.ylabel('U [V]')
    x = np.linspace(0, 150, 150)
    plt.errorbar(i, u, xerr = 3, yerr = 0.05, linestyle ='', label = 'Messwerte', color = 'blue')
    plt.plot(x, m*x+n, label = 'Geraden fit', color = 'blue')

    # rx = - 100/7 *10**-3
    # n = 3* 5/7
    # plt.plot(x, rx*x + n , label = 'Rx-Gerade', color = 'orange')

    plt.legend()
    plt.grid()

    plt.savefig('232e_plot.png')
#    plt.show()


if __name__ == '__main__': 
    main()