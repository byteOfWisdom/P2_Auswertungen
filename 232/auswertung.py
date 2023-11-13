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
        note('{} = {} {}'.format(var, value, unit))
    else:
        note('{} = ({} +- {}) {}'.format(var, value, err, unit))

def write_notes(file):
    global misc_values
    with open(file, 'w') as handle:
        handle.write(misc_values)

plot_n = 0
def plot_num():
    global plot_n
    plot_n += 1
    return plot_n


def save(file):
    plt.savefig(file, dpi=250)
    plt.clf()


def sq(x): return x * x


def amp_meter_resistance(max_amplitude):
    v = 50 * 2e-3
    return v / max_amplitude


def voltmeter_resistance(max_amplitude):
    ri_base = 50
    I = 2e-3

    Rv = (max_amplitude / I) - ri_base
    return Rv + ri_base



def a(preview):
    data = easyparse.parse("daten/232a.csv")
    u = data['U'] * (data['U_max'] / data['U_skt'])
    i = data['I'] * (data['I_max'] / data['I_skt']) #* 0.6
    result = linregress(i, u)

    m = result.slope
    n = result.intercept
    dm = result.stderr
    dn = result.intercept_stderr

    note('U: {}'.format(str(list(map(lambda x: round(x, 3), u)))))
    note('I: {}'.format(str(list(map(lambda x: round(x, 3), i)))))

    rx = m - amp_meter_resistance(0.05)

    note_var('a', m, dm, 'Ohm')
    note_var('b', n, dn, 'Volt')

    note_var('Ri', amp_meter_resistance(0.05), unit='Ohm')
    note_var('Rb', m, round(dm, 3), 'Ohm')
    note_var('Rx', round(m - amp_meter_resistance(0.5), 3), round(dm, 3), 'Ohm')

    x = np.linspace(min(i), max(i), 1000)
    plt.errorbar(
        i, u, 
        xerr =  (data['I_max'] / data['I_skt']), 
        yerr = (data['U_max'] / data['U_skt']), 
        fmt='x', linestyle ='', 
        label = 'Messwerte', 
        color = 'blue')

    plt.plot(x, m*x+n, label = 'Geraden fit', color = 'blue')

    plt.plot(x, (rx)*x, label = 'Rx-Gerade', color = 'orange')

    plt.title(r"Fig {}: U-I Diagramm mit fit und $R_x$".format(plot_num()))
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
    i = data['I'] * (data['I_max'] / data['I_skt']) #* 0.6
    result = linregress(i, u)
    m = result.slope
    n = result.intercept
    dm = result.stderr
    dn = result.intercept_stderr
    note('Gradenfit mit:')
    note_var('a', m, dm, 'Ohm')
    note_var('b', n, dn, 'V')


    plt.xlabel('I [mA]')
    plt.ylabel('U [V]')
    x = np.linspace(min(i), max(i), 1000)
    plt.errorbar(i, u, xerr = (data['I_max'] / data['I_skt']), yerr = (data['U_max'] / data['U_skt']), linestyle ='', label = 'Messwerte')
    plt.plot(x, m*x+n, label = 'Geraden fit')


    plt.title('Fig {}: Leerlaufspannung und Innenwiderstand'.format(plot_num()))
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
    x = 100 - data['R']
    result = linregress(x, u)
    m = result.slope
    n = result.intercept
    dm = result.stderr
    dn = result.intercept_stderr
    note('Aufgabe f:')
    note('R = inf')
    note_var('a', m, dm)
    note_var('b', n, dn)

    plt.xlabel('x in Skt')
    plt.ylabel('U [V]')
    x0 = np.linspace(min(x), max(x), 1000)
    plt.errorbar(x, u, xerr = 2, yerr = 0.05, linestyle ='', fmt = 'x', label = r'Messwerte ohne Last', color = 'blue')
    plt.plot(x0, m*x0+n, label = r'Geraden fit ohne Last', color = 'blue')

    # R = 20 ohm
    
    data = easyparse.parse("daten/232f_20.csv")
    u = data['U'] * (data['U_max'] / data['U_skt'])
    x = 100 - data['R']
    result = linregress(x, u)
    m = result.slope
    n = result.intercept
    dm = result.stderr
    dn = result.intercept_stderr
    note('Aufgabe f:')
    note('R = 20 Ohm')
    note_var('a', m, dm)
    note_var('b', n, dn)

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
    x = 100 - data['R']
    result = linregress(x, u)
    m = result.slope
    n = result.intercept
    dm = result.stderr
    dn = result.intercept_stderr
    note('Aufgabe f:')
    note('R = 50 Ohm')
    note_var('a', m, dm)
    note_var('b', n, dn)

    output = {'m': [m], 'n':[n], 'dm':[dm], 'dn':[dn]}
    easyparse.write_printable(output, 'temp2.csv')



    plt.xlabel('x in Skt')
    plt.ylabel('U [V]')
    x0 = np.linspace(min(x), max(x), 1000)
    plt.errorbar(x, u, xerr = 2, yerr = 0.05, linestyle ='', fmt = 'x', label = 'Messwerte R = 50 $\Omega$', color = 'green')
    plt.plot(x0, m*x0+n, label = r'Geraden fit R = 50 $\Omega$', color = 'green')

    plt.title('Fig {}: Belasteter Spannungsteiler'.format(plot_num()))

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

    data = easyparse.parse("daten/232f_20.csv")
    x = 100 - data['R']
    p = m * x**2 / 20
    dp = ((5.32 * 10**-19 * x**4)+ (2.64 * 10**-16 * x**2))**1/2

    plt.xlabel('x in Skt')
    plt.ylabel('P [W]')
    x0 = np.linspace(0, max(x), 1000)
    raw_data = easyparse.parse("daten/232f_20.csv")

    dp = np.sqrt(sq(raw_data['U'] * raw_data['U_max'] * raw_data['I_max'] / (raw_data['I_skt'] * raw_data['U_skt'])) + sq(raw_data['I'] * raw_data['U_max'] * raw_data['I_max'] / (raw_data['I_skt'] * raw_data['U_skt'])))

    #plt.errorbar(x, p, xerr = 2, yerr = dp, linestyle ='', fmt = 'x', label = r'Messwerte R = 20 $\Omega$', color = 'orange')
    plt.errorbar(100 - raw_data['R'], raw_data['U'] * raw_data['U_max'] * raw_data['I'] * raw_data['I_max'] / (raw_data['U_skt'] * raw_data['I_skt']), xerr = 2, yerr = dp, linestyle ='', fmt = 'x', label = r'Messwerte R = 20 $\Omega$', color = 'orange')
    plt.plot(x0, (m * x0)**2 /20, label = r'Leistung R = 20 $\Omega$', color = 'orange')

    data = easyparse.parse("temp2.csv", ',')
    m = data['m']
    n = data ['n']
    dm = data['dm']
    dn = data['dn']

    data = easyparse.parse("daten/232f_50.csv")
    x = 100 - data['R']
    p = m * x**2 / 50
    dp = ((9.71 * 10**-19 * x**4)+ (8.79 * 10**-16 * x**2))**1/2

    raw_data = easyparse.parse("daten/232f_50.csv")

    dp = np.sqrt(sq(raw_data['U'] * raw_data['U_max'] * raw_data['I_max'] / (raw_data['I_skt'] * raw_data['U_skt'])) + sq(raw_data['I'] * raw_data['U_max'] * raw_data['I_max'] / (raw_data['I_skt'] * raw_data['U_skt'])))

    plt.xlabel('x in Skt')
    plt.ylabel('P [W]')
    x0 = np.linspace(0, max(x), 1000)
    #plt.errorbar(x, p, xerr = 2, yerr = dp, linestyle ='', fmt = 'x', label = r'Messwerte R = 50 $\Omega$', color = 'purple')
    plt.errorbar(100 - raw_data['R'], raw_data['U'] * raw_data['U_max'] * raw_data['I'] * raw_data['I_max'] / (raw_data['U_skt'] * raw_data['I_skt']), xerr = 2, yerr = dp, linestyle ='', fmt = 'x', label = r'Messwerte R = 50 $\Omega$', color = 'purple')

    plt.plot(x0, (m * x0)**2 / 50, label = r'Leistung R = 50 $\Omega$', color = 'purple')

    plt.legend()
    plt.grid()
    plt.title('Fig {}: Leistungskurve pro Lastwiderstand'.format(plot_num()))

    # plt.show()
    if preview: plt.show()
    else: save('results/232g_plot.png')



def fit_and_plot(x, y, f, err, preview):
    #x = x + 273.15
    fit = linregress(x, y)
    a = fit.slope
    b = fit.intercept
    da = fit.stderr
    db = fit.intercept_stderr
    plt.errorbar(x, y, xerr=err, yerr=(y * 1e-2), fmt='x', label='Messdaten')
    #plt.plot(x, y, 'x', label='Messdaten')
    x_ = np.linspace(min(x), max(x), 1000)
    plt.plot(x_, a * x_ + b, label='funktionsfit')
    plt.legend()
    plt.grid()
    if preview:
        plt.show()
    return [a, b], [da, db]


def n(preview):
    note('Aufgabe n:')
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
    r_carbon = 'R5'

    metal = lambda t, R, a : R * (1 + a * t)
    semi_conductor = lambda t, R, Ek: R * np.exp(Ek / (2 * t))

    plt.title('Fig {}: Platinwiderstand'.format(plot_num()))
    plt.xlabel(r'T $[^\circ C]$')
    plt.ylabel(r'R $[\Omega]$')
    fit, err = fit_and_plot(data[t_platin], data[r_platin], metal, 0.5, preview)
    note('Platin:')
    note('b = ' + str(fit[1]) + '+-' + str(err[1]))
    note('a = ' + str(fit[0]) + '+-' + str(err[0]))
    note('b/a = ' + str(fit[0] / fit[1]))
    note('d(a/b) = ' + str(np.sqrt(sq(err[0] / fit[1]) + sq(fit[0] * err[1] / sq(fit[1])))))
    if not preview: save('results/232platin.png')


    plt.title('Fig {}: Konstantanwiderstand'.format(plot_num()))
    plt.xlabel(r'T $[^\circ C]$')
    plt.ylabel(r'R $[\Omega]$')
    plt.ylim([min(data[r_konstantan]) * (1/1.5), max(data[r_konstantan] * 1.5)])
    fit, err = fit_and_plot(data[t_konstantan], data[r_konstantan], metal, 0.5, preview)
    note('Konstantan:')
    note('b = ' + str(fit[1]) + '+-' + str(err[1]))
    note('a = ' + str(fit[0]) + '+-' + str(err[0]))
    note('a/b = ' + str(fit[0] / fit[1]))
    note('d(a/b) = ' + str(np.sqrt(sq(err[0] / fit[1]) + sq(fit[0] * err[1] / sq(fit[1])))))
    if not preview: save('results/232konstantan.png')


    plt.title('Fig {}:Kohleschichtwiderstand'.format(plot_num()))
    plt.xlabel(r'T $[^\circ C]$')
    plt.ylabel(r'R $[\Omega]$')
    plt.ylim([min(data[r_carbon]) * (1/1.5), max(data[r_carbon] * 1.5)])
    fit, err = fit_and_plot(data[t_carbon], data[r_carbon], metal, 0.5, preview)
    note('Kohleschichtwiderstand:')
    note('b = ' + str(fit[1]) + '+-' + str(err[1]))
    note('a = ' + str(fit[0]) + '+-' + str(err[0]))
    note('a/b = ' + str(fit[0] / fit[1]))
    note('d(a/b) = ' + str(np.sqrt(sq(err[0] / fit[1]) + sq(fit[0] * err[1] / sq(fit[1])))))
    if not preview: save('results/232kohleschicht.png')


    plt.title('Fig {}: NTC-Widerstand'.format(plot_num()))
    plt.xlabel(r'$\frac{1}{T} [K^{-1}]$')
    plt.ylabel(r'$ln{R}$')
    fit, err = fit_and_plot(1/(data[t_ntp] + 273.15), np.log(data[r_ntp]), semi_conductor, 0.5 / sq(data[t_ntp] + 273.15), preview)
    note('NTC:')
    note('b = ' + str(fit[1]) + '+-' + str(err[1]))
    note('a = ' + str(fit[0]) + '+-' + str(err[0]))
    note('a/b = ' + str(fit[0] / fit[1]))
    note('d(a/b) = ' + str(np.sqrt(sq(err[0] / fit[1]) + sq(fit[0] * err[1] / sq(fit[1])))))
    if not preview: save('results/232ntc.png')
    else: plt.show()

    #plt.yscale('log')
    #plt.errorbar(data[t_ptc] + 273.15, np.log(data[r_ptc]), xerr=0.5, yerr=1e-2 ,fmt='x', label='Messdaten')
    plt.errorbar(data[t_ptc] + 273.15, np.log(data[r_ptc]), xerr=0.5, yerr=1 / data[r_ptc] ,fmt='x', label='Messdaten')
    plt.vlines([330], 4, 13, linestyle='--', color='yellow', label='Curie-Temperatur')
    plt.xlabel(r'$T [K]$')
    plt.ylabel(r'$ln{R}$')
    plt.legend()
    plt.grid()
    plt.title('Fig {}: PTC-Widerstand'.format(plot_num()))
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