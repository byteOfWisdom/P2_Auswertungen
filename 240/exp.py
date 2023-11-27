from labtools.defaults import *
from scipy.optimize import curve_fit
from numba import njit


@njit
def approx_mu_max(x, y):
    tolerance = 0.001
    possible_tangents = np.linspace(10, 0, 1_000_000)
    for tangent in possible_tangents:
        points = x * tangent
        if len(list(filter(lambda x: x, np.abs(points - y) <= tolerance))) > 0:
            return tangent


def c(preview, data=None):
    table = data['b']['cassy']['data']

    I = table[1] * data['b']['mirror']
    dI = I * 1e-2 + 3 * 5e-2 # 1 % plus 5% des messbereichs
    I = ev(I, dI) # I = (I +- dI) A

    B = ev(table[2], table [2] * 3e-2) * 1e-3
    N = data['b']['N']
    d = ev(2e-3, 0.05e-3)
    l = ev(477e-3, 2e-3)

    mu0 = 4 * np.pi * 1e-7

    H = (2 * N * I / l) - B * (d / (mu0 * l))

    value_table = {
        'I [A]': I,
        'B [T]': B,
        'H [A/m]': H,
    }

    if not preview:
        write_printable(value_table, 'results/240.csv')


    start = np.where(np.abs(I) > 0.1)[0][0]
    end = np.where(np.abs(I) > 2)[0][0]


    x_mu_A = np.where(np.abs(B) > (0.002 + abs(B[0])))[0][0]
    mu_A = (B[x_mu_A] / H[x_mu_A]).value - 1e-6

    mu_max = approx_mu_max(value(H[start:end]), value(B[start:end]))

    note_var('mu_max', mu_max, unit='N/A^2')
    note_var('mu_A', mu_A, unit='N/A^2')
    note_var('mu_max', mu_max / mu0, unit='mu_0')
    note_var('mu_A', mu_A / mu0, unit='mu_0')
    note_var('d_mu', 1e-5 / mu0, unit='mu_0')


    if not some(mu_max):
        print('failed to find mu_max')
        return


    plot = Plot(r'$H [\frac{A}{m}]$', 'B [T]')
    plot.title = 'Hysterese Kurve'
    plot.add_element(ev(value(H), 0), ev(value(B), 0), 'Messdaten')
    plot.add_element(lambda x: x * mu_A, r'$\mu_A$', color='red')
    plot.add_element(lambda x: x * mu_max, r'$\mu_{max}$', color='green')

    plot.finish(preview, 'results/240c.png')


    plot = Plot(r'$H [\frac{A}{m}]$', 'B [T]')
    plot.title = 'Hysterse Kurve mit Fehlerbalken'
    plot.add_element(H, B, 'Messdaten')
    plot.add_element(lambda x: x * mu_A, r'$\mu_A$', color='red')
    plot.add_element(lambda x: x * mu_max, r'$\mu_{max}$', color='green')

    plot.finish(preview, 'results/240c_fehlerbalken.png')