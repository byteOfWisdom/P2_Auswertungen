import labtools as tools
from labtools.defaults import *
from labtools import math


temp = 0.0


def viscocity():
    return np.interp(temp, [0, 20, 40], [17.2e-6, 18.9e-6, 19.12e-6])

#   return 18.9e-6


density_air = 1.225
density_oil = 886
grav_const = 9.81
plate_seperation = None
delta_E_field = 0

class Drop():
    def __init__(self, data, s, dt):
        v_data = t_to_v(data, s, dt)
        self.id = v_data[0]
        self.v0 = v_data[1]
        self.v_up = v_data[2]
        self.v_down = v_data[3]
        self.E_field = ev(data[4] / plate_seperation, delta_E_field)
        #print('drop: {}, {}, {}'.format(self.v0.value, self.v_down.value, self.v_up.value))
        #print('diff: {}'.format(2 * self.v0.value - (self.v_down.value - self.v_up.value)))


def is_valid(droplet):
    # checks if the droplet fullfills 2v_0 = v_down - v_up
    # where v_up has a negative sign
    return tools.perror.within_deviation(2 * droplet.v0, droplet.v_down - droplet.v_up)


def t_to_v(droplet_data, s, dt):
    t_0 = ev(droplet_data[1], dt)
    t_up = ev(droplet_data[2], dt)
    t_down = ev(droplet_data[3], dt)
    v_drop = [droplet_data[0], ev(s, 0.01e-3) / t_0, ev(s, 0.01e-3) / t_up, ev(s, 0.01e-3) / t_down]
    return v_drop


def radius(drop):
    return tools.np.sqrt((9 * viscocity() * (drop.v_down - drop.v_up)) / (4 * grav_const * (density_oil - density_air)))


def charge(drop):
    return 3 * tools.np.pi * viscocity() * radius(drop) * (drop.v_down + drop.v_up) / drop.E_field


def average(list_of_drops):
    averaged_drops = []

    current_id = list_of_drops[0].id
    acc = [list_of_drops[0]]
    i = 1

    while i < len(list_of_drops):
        if list_of_drops[i].id == current_id:
            acc.append(list_of_drops[i])
        else:
            avd = Drop([1, 1, 1, 1, -1], 1, 1)

            n = len(acc)
            avd.id = acc[0].id
            avd.v0 = sum(d.v0 for d in acc) / n
            avd.v_down = sum(d.v_down for d in acc) / n
            avd.v_up = sum(d.v_up for d in acc) / n
            avd.E_field = sum(d.E_field for d in acc) / n

            averaged_drops.append(avd)
            current_id = list_of_drops[i].id
            acc = [list_of_drops[i]]

        i += 1

    return averaged_drops


def find_elemental_charge(charges, preview):
    e = math.agcd(value(charges), 0.000001)
    ns = np.round(value(charges) / e)
    return e, charges / ns


def g(preview, data=None):
    global temp, plate_seperation, delta_E_field
    data = data['c']

    plate_seperation = data['plate_seperation']
    delta_E_field = data['dU'] / plate_seperation
    #E_field = data['U'] / data['plate_seperation']
    temp = data['temp']

    note_var('viscocity', viscocity())

    drops = [Drop(drop, data['distance'], data['dt']) for drop in data['droplets']]

    if not preview: write_printable({
        'Tropfen': [d.id for d in drops],
        r'v0 [$\frac{m}{s}$]': [d.v0 for d in drops],
        r'v_up [$\frac{m}{s}$]': [d.v_up for d in drops],
        r'v_down [$\frac{m}{s}$]': [d.v_down for d in drops],},
        'results/242g_rohdaten.csv')

    drops = average(drops)
    drops = list(filter(is_valid, drops))
    charges = tools.np.array([charge(drop) for drop in drops])
    r = np.array([radius(d) for d in drops])
    if not preview: write_printable({
        r'Tropfen': [d.id for d in drops],
        r'r [$m$]': r,
        r'v0 [$\frac{m}{s}$]': [d.v0 for d in drops],
        r'v_up [$\frac{m}{s}$]': [d.v_up for d in drops],
        r'v_down [$\frac{m}{s}$]': [d.v_down for d in drops],
        r'q [$C$]': charges,
    }, 'results/242g_ladungen.csv')


    avg_err_charges = sum([c.error for c in charges]) / len(charges)


    elemental_charge, esi = find_elemental_charge(charges, preview)
    #elemental_charge = tools.math.agcd([c.value for c in charges], avg_err_charges / 100000) # well fuck this ain't working
    note_var('e_si', elemental_charge, unit='C')

    plot = tools.Plot('Tropfen', 'Ladung [C]')
    plot.add_element(tools.np.array([n for n in range(len(charges))]), charges)
    plot.finish(preview, 'results/242g.png')

    # aufgabe i

    if not preview: write_printable({
        'Tropfen': [d.id for d in drops],
        r'r [$m$]': r,
        r'v0 [$\frac{m}{s}$]': [d.v0 for d in drops],
        r'v_up [$\frac{m}{s}$]': [d.v_up for d in drops],
        r'v_down [$\frac{m}{s}$]': [d.v_down for d in drops],
        r'q [$C$]': charges,
        r'e_si [$C$]': esi,
    }, 'results/242g_ladungen.csv')


    y = (esi ** (2 / 3))
    x = (1 / r)

    plot = Plot(r'$\frac{1}{r} [m^{-1}]$', r'$e_{s,i}^{2/3} [C^{3/2}]$')
    points = plot.add_element(x, y, r'$e_{s,i}^{3/2}$')
    func, params = fit_func(lambda x, a, b: a * x + b, x, y)
    plot.add_element(func, 'Gradenfit')
    plot.title = 'Cunningham - Korrektur'
    plot.finish(preview, 'results/242i.png')

    y = value(esi ** (2 / 3))
    x = value(1 / r)

    plot = Plot(r'$\frac{1}{r} [m^{-1}]$', r'$e_{s,i}^{2/3} [C^{3/2}]$')
    points = plot.add_element(x, y, r'$e_{s,i}^{3/2}$')
    func, params = fit_func(lambda x, a, b: a * x + b, x, y)
    plot.add_element(func, 'Gradenfit')
    plot.title = 'Cunningham - Korrektur'
    plot.finish(preview, 'results/242i_keine_fehlerbalken.png')


    note_var('a', params[0])
    note_var('b', params[1])
    note_var('e0', params[1] ** (3/2))
