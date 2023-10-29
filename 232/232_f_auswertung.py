
import easyparse
import numpy as np
from scipy.stats import linregress
from matplotlib import pyplot as plt

def main():

    # ohne Last
    data = easyparse.parse("232f_0.csv")
    u = data['U']
    x = data['x']
    result = linregress(x, u)
    m = result.slope
    n = result.intercept
    dm = result.stderr
    dn = result.intercept_stderr
    print(m, n, dm, dn)

    plt.xlabel('x in Skt')
    plt.ylabel('U [V]')
    x0 = np.linspace(0, 900, 900)
    plt.errorbar(x, u, xerr = 2, yerr = 0.05, linestyle ='', fmt = 'x', label = r'Messwerte ohne Last', color = 'blue')
    plt.plot(x0, m*x0+n, label = r'Geraden fit ohne Last', color = 'blue')

    # R = 20 omega
    
    data = easyparse.parse("232f_20.csv")
    u = data['U']
    x = data['x']
    result = linregress(x, u)
    m = result.slope
    n = result.intercept
    dm = result.stderr
    dn = result.intercept_stderr
    print(m, n, dm, dn)

    output = {'m': [m], 'n':[n], 'dm':[dm], 'dn':[dn]}
    easyparse.write_printable(output, 'temp1.csv')

    plt.xlabel('x in Skt')
    plt.ylabel('U [V]')
    x0 = np.linspace(0, 900, 900)
    plt.errorbar(x, u, xerr = 2, yerr = 0.05, linestyle ='', fmt = 'x', label = r'Messwerte R = 20 $\Omega$', color = 'orange')
    plt.plot(x0, m*x0+n, label = r'Geraden fit R = 20 $\Omega$', color = 'orange')

    # R = 50 omega

    data = easyparse.parse("232f_50.csv")
    u = data['U']
    x = data['x']
    result = linregress(x, u)
    m = result.slope
    n = result.intercept
    dm = result.stderr
    dn = result.intercept_stderr
    print(m, n, dm, dn)

    output = {'m': [m], 'n':[n], 'dm':[dm], 'dn':[dn]}
    easyparse.write_printable(output, 'temp2.csv')

    plt.xlabel('x in Skt')
    plt.ylabel('U [V]')
    x0 = np.linspace(0, 900, 900)
    plt.errorbar(x, u, xerr = 2, yerr = 0.05, linestyle ='', fmt = 'x', label = 'Messwerte R = 50 $\Omega$', color = 'green')
    plt.plot(x0, m*x0+n, label = r'Geraden fit R = 50 $\Omega$', color = 'green')



    plt.legend()
    plt.grid()
    plt.show()

    # plt.savefig('232f_plot.png')



if __name__ == '__main__': 
    main()