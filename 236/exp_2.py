import labtools
from labtools.defaults import *

from scipy.stats import linregress


def h(preview, data_file=None):
    labtools.notes.note('h:')
    data = data_file['h']
    #data = labtools.ep.parse('data/236h.csv')
    time = data['t']
    time_err = data['dt']
    angle = data['phi'] / 2 # skala wert ist 2 * die auslenkung
    angle_err = data['dphi'] / 2
    capacitance = data['c']

    write_printable({'phi_m': ev(angle, angle_err), 'ln(phi)': ev(np.log(angle), angle_err / angle)}, 'results/236h.csv')

    plot = labtools.Plot('t [s]', r'$\ln{\varphi}$')
    res = linregress(time, labtools.np.log(angle))

    plot.add_element(
        time, 
        labtools.np.log(angle), 
        time_err,
        angle_err / angle,
        'Messdaten')


    plot.add_element(lambda x: res.slope * x + res.intercept, 'Fit Funktion')


    a = labtools.p.make_errval(res.slope, res.stderr)
    b = labtools.p.make_errval(res.intercept, res.intercept_stderr)
    labtools.notes.note_var('a', a)
    labtools.notes.note_var('b', b)

    R = - 1 / (a * capacitance)
    labtools.notes.note_var('R', R)


    plot.title = 'Fig 2: Widerstandsmessung'
    plot.finish(preview, 'results/236h.png')