#!python3

from sys import argv

import easyparse
import numpy as np
from scipy.stats import linregress
from matplotlib import pyplot as plt

misc_values = ''


def note(s):
    global misc_values
    misc_values += s + "\n"


def note_var(var, value, err=None, unit=''):
    if err == None:
        note('{} = {} {}'.forma(var, value, unit))
    else:
        note('{} = ({} +- {}) {}'.format(var, value, err, unit))

def write_notes(file):
    global misc_values
    with open(file, 'w') as handle:
        handle.write(misc_values)


def save(file):
    plt.savefig(file)
    plt.clf()


def sq(x): return x * x


def amp_meter_resistance(max_amplitude):
    ri_base = 50 # aus der Anleitung
    base_max_I = 2e-3

    shunt = ri_base * base_max_I / (max_amplitude - base_max_I)
    return (ri_base * shunt) / (ri_base + shunt)



def voltmeter_resistance(max_amplitude):
    ri_base = 50
    I = 2e-3

    Rv = (max_amplitude / I) - ri_base
    return Rv + ri_base



def a(preview):
    data = easyparse.parse("daten/232a.csv")
    u = data['U'] * (data['U_max'] / data['U_skt'])
    i = data['I'] * (data['I_max'] / data['I_skt'])
    rx = data['Rx']
    result = linregress(i, u)

    m = result.slope
    n = result.intercept
    dm = result.stderr
    dn = result.intercept_stderr

    note_var('Rb', round(m, 3), round(dm, 3), 'Ohm')
    note_var('Rx', rx, round(rx * 1e-2, 1), 'Ohm')

    #print(m, n, dm, dn)
    #print(amp_meter_resistance(data['I_max']))
    #print(voltmeter_resistance(data['U_max']))

    x = np.linspace(min(i), max(i), 1000)
    plt.errorbar(i, u, xerr =  (data['I_max'] / data['I_skt']), yerr = (data['U_max'] / data['U_skt']), fmt='x', linestyle ='', label = 'Messwerte', color = 'blue')
    plt.plot(x, m*x+n, label = 'Geraden fit', color = 'blue')

    plt.plot(x, (rx)*x, label = 'Rx-Gerade', color = 'orange')

    plt.xlabel("I [A]")
    plt.ylabel("U [V]")
    plt.legend()
    plt.grid()

    if preview:
        plt.show()
    else:
        save('results/232a.png')


def e(preview):
    data = easyparse.parse("daten/232d.csv")
    u = data['U'] * (data['U_max'] / data['U_skt'])
    i = data['I'] * (data['I_max'] / data['I_skt'])
    result = linregress(i, u)
    m = result.slope
    n = result.intercept
    dm = result.stderr
    dn = result.intercept_stderr
    note('Gradenfit mit:')
    note_var('a', m, dm, 'Ohm')
    note_var('b', n, dn, 'Ohm')


    plt.xlabel('I [mA]')
    plt.ylabel('U [V]')
    x = np.linspace(min(i), max(i), 1000)
    plt.errorbar(i, u, xerr = (data['I_max'] / data['I_skt']), yerr = (data['U_max'] / data['U_skt']), linestyle ='', label = 'Messwerte')
    plt.plot(x, m*x+n, label = 'Geraden fit')

    # rx = - 100/7 *10**-3
    # n = 3* 5/7
    # plt.plot(x, rx*x + n , label = 'Rx-Gerade', color = 'orange')

    plt.legend()
    plt.grid()

    if preview:
        plt.show()
    else:
        save('results/232e.png')


def f(preview):
    # ohne Last
    data = easyparse.parse("daten/232f_inf.csv")
    u = data['U'] * (data['U_max'] / data['U_skt'])
    x = data['R']
    result = linregress(x, u)
    m = result.slope
    n = result.intercept
    dm = result.stderr
    dn = result.intercept_stderr
    print(m, n, dm, dn)

    plt.xlabel('x in Skt')
    plt.ylabel('U [V]')
    x0 = np.linspace(min(x), max(x), 1000)
    plt.errorbar(x, u, xerr = 2, yerr = 0.05, linestyle ='', fmt = 'x', label = r'Messwerte ohne Last', color = 'blue')
    plt.plot(x0, m*x0+n, label = r'Geraden fit ohne Last', color = 'blue')

    # R = 20 ohm
    
    data = easyparse.parse("daten/232f_20.csv")
    u = data['U'] * (data['U_max'] / data['U_skt'])
    x = data['R']
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
    x0 = np.linspace(min(x), max(x), 1000)
    plt.errorbar(x, u, xerr = 2, yerr = 0.05, linestyle ='', fmt = 'x', label = r'Messwerte R = 20 $\Omega$', color = 'orange')
    plt.plot(x0, m*x0+n, label = r'Geraden fit R = 20 $\Omega$', color = 'orange')

    # R = 50 ohm

    data = easyparse.parse("daten/232f_50.csv")
    u = data['U'] * (data['U_max'] / data['U_skt'])
    x = data['R']
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
    x0 = np.linspace(min(x), max(x), 1000)
    plt.errorbar(x, u, xerr = 2, yerr = 0.05, linestyle ='', fmt = 'x', label = 'Messwerte R = 50 $\Omega$', color = 'green')
    plt.plot(x0, m*x0+n, label = r'Geraden fit R = 50 $\Omega$', color = 'green')



    plt.legend()
    plt.grid()

    if preview:
        plt.show()
    else:
        save('results/232f.png')



def g(preview):
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



def fit_and_plot(x, y, f, preview):
    #x = x + 273.15
    fit = linregress(x, y)
    a = fit.slope
    b = fit.intercept
    da = fit.stderr
    db = fit.intercept_stderr
    plt.plot(x, y, 'x', label='Messdaten')
    x_ = np.linspace(min(x), max(x), 1000)
    plt.plot(x_, a * x_ + b, label='funktionsfit')
    plt.legend()
    plt.grid()
    if preview:
        plt.show()
    return [a, b], [da, db]


def n(preview):
    data = easyparse.parse('daten/232n.csv')
    t_ntp = 'T1'
    r_ntp = 'R1'
    t_konstantan = 'T2'
    r_konstantan = 'R2'
    t_ptc = 'T3'
    r_ptc = 'R3'
    t_platin = 'T4'
    r_platin = 'R4'
    t_carbon = 'T5'
    r_carbon = 'R4'

    metal = lambda t, R, a : R * (1 + a * t)
    semi_conductor = lambda t, R, Ek: R * np.exp(Ek / (2 * t))

    plt.title('Platin')
    plt.xlabel(r'T $[^\circ C]$')
    plt.ylabel(r'R $[\Omega]$')
    fit, err = fit_and_plot(data[t_platin], data[r_platin], metal, preview)
    print('Platin:')
    print('b = ' + str(fit[1]) + '+-' + str(err[1]))
    print('a = ' + str(fit[0]) + '+-' + str(err[0]))
    print('b/a = ' + str(fit[0] / fit[1]))
    print('d(a/b) = ' + str(np.sqrt(sq(err[0] / fit[1]) + sq(fit[0] * err[1] / sq(fit[1])))))
    if not preview: save('results/232platin.png')


    plt.title('Konstantan')
    plt.xlabel(r'T $[^\circ C]$')
    plt.ylabel(r'R $[\Omega]$')
    plt.ylim([min(data[r_konstantan]) * (1/1.5), max(data[r_konstantan] * 1.5)])
    fit, err = fit_and_plot(data[t_konstantan], data[r_konstantan], metal, preview)
    print('Konstantan:')
    print('b = ' + str(fit[1]) + '+-' + str(err[1]))
    print('a = ' + str(fit[0]) + '+-' + str(err[0]))
    print('a/b = ' + str(fit[0] / fit[1]))
    print('d(a/b) = ' + str(np.sqrt(sq(err[0] / fit[1]) + sq(fit[0] * err[1] / sq(fit[1])))))
    if not preview: save('results/232konstantan.png')


    plt.title('Kohleschicht')
    plt.xlabel(r'T $[^\circ C]$')
    plt.ylabel(r'R $[\Omega]$')
    plt.ylim([min(data[r_carbon]) * (1/1.5), max(data[r_carbon] * 1.5)])
    fit, err = fit_and_plot(data[t_carbon], data[r_carbon], metal, preview)
    print('Kohleschichtwiderstand:')
    print('b = ' + str(fit[1]) + '+-' + str(err[1]))
    print('a = ' + str(fit[0]) + '+-' + str(err[0]))
    print('a/b = ' + str(fit[0] / fit[1]))
    print('d(a/b) = ' + str(np.sqrt(sq(err[0] / fit[1]) + sq(fit[0] * err[1] / sq(fit[1])))))
    if not preview: save('results/232kohleschicht.png')


    plt.title('NTC')
    plt.xlabel(r'$\frac{1}{T} [K^{-1}]$')
    plt.ylabel(r'$ln{R}$')
    fit, err = fit_and_plot(1/(data[t_ntp] + 273.15), np.log(data[r_ntp]), semi_conductor, preview)
    print('NTC:')
    print('b = ' + str(fit[1]) + '+-' + str(err[1]))
    print('a = ' + str(fit[0]) + '+-' + str(err[0]))
    print('a/b = ' + str(fit[0] / fit[1]))
    print('d(a/b) = ' + str(np.sqrt(sq(err[0] / fit[1]) + sq(fit[0] * err[1] / sq(fit[1])))))
    if not preview: save('results/232ntc.png')
    else: plt.show()

    plt.title('PTC')
    #plt.yscale('log')
    plt.plot(data[t_ptc] + 273.15, np.log(data[r_ptc]), 'x', label='Messdaten')
    plt.vlines([330], 4, 13, linestyle='--', color='yellow', label='Curie-Temperatur')
    plt.xlabel(r'$T [K]$')
    plt.ylabel(r'$ln{R}$')
    plt.legend()
    plt.grid()
    if not preview: save('results/232ptc.png')
    else: plt.show()



def main():
    preview = False

    if '-p' in argv or 'pv' in argv:
        preview = True
    
    if '-a' in argv[1:]:
        a(preview)
        note('')
    elif '-e' in argv[1:]:
        e(preview)
        note('')
    elif '-f' in argv[1:]:
        f(preview)
        note('')
    elif '-g' in argv[1:]:
        g(preview)
        note('')
    elif '-n' in argv[1:]:
        n(preview)
        note('')
    else:
        a(preview)
        note('')
        e(preview)
        note('')
        f(preview)
        note('')
        g(preview)
        note('')
        n(preview)
        note('')

    if not preview:
        write_notes('results/notes.txt')
    else:
        global misc_values
        print(misc_values)




if __name__ == '__main__':
    main()