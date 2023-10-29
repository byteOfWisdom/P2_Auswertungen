import easyparse
import numpy as np
from scipy.stats import linregress
from matplotlib import pyplot as plt

def main():

    data = easyparse.parse("temp1.csv", ',')
    m = data['m']
    n = data ['n']
    dm = data['dm']
    dn = data['dn']

    data = easyparse.parse("232f_20.csv")
    x = data['x']
    p = m * x**2 / 20
    dp = ((5.32 * 10**-19 * x**4)+ (2.64 * 10**-16 * x**2))**1/2

    plt.xlabel('x in Skt')
    plt.ylabel('P [W]')
    x0 = np.linspace(0, 900, 900)
    plt.errorbar(x, p, xerr = 2, yerr = dp, linestyle ='', fmt = 'x', label = r'Messwerte R = 20 $\Omega$', color = 'orange')
    plt.plot(x0, m * x0**2 /20, label = r'Leistung R = 20 $\Omega$', color = 'orange')

    data = easyparse.parse("temp2.csv", ',')
    m = data['m']
    n = data ['n']
    dm = data['dm']
    dn = data['dn']

    data = easyparse.parse("232f_50.csv")
    x = data['x']
    p = m * x**2 / 50
    dp = ((9.71 * 10**-19 * x**4)+ (8.79 * 10**-16 * x**2))**1/2

    plt.xlabel('x in Skt')
    plt.ylabel('P [W]')
    x0 = np.linspace(0, 900, 900)
    plt.errorbar(x, p, xerr = 2, yerr = dp, linestyle ='', fmt = 'x', label = r'Messwerte R = 50 $\Omega$', color = 'purple')
    plt.plot(x0, m * x0**2 /50, label = r'Leistung R = 50 $\Omega$', color = 'purple')

    plt.legend()
    plt.grid()

    # plt.show()
    plt.savefig('232g_plot.png')

if __name__ == '__main__': 
    main()