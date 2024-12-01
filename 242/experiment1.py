from labtools.defaults import *

u0 = 1.256e-6
n = 130
R = 150e-3
c = 0.716 * u0 * n / R

def Bs(I):
	return c * I


def b(preview, data=None):
	data = data['a']
	r = 0.5 * ev((data['d_base'] - data['d']) * 1e-2, 2e-3)
	U = ev(data['U'], 1)
	Ia = ev(data['I_a'], 0.01)
	Ib = ev(data['I_b'], 0.01)

	I = (Ia + Ib) / 2

	B = Bs((Ia + Ib) / 2)

	plot = Plot(r'$U [V]$', r'$(r * I)^2 [(m * A)^2]$')
	plot.title = 'Spezifische Ladung'
	p1 = plot.add_element(U, sq(r * I), 'Messdaten')
	func, params = fit_func(lambda x, a, b: a * x + b, U, sq(r * I))
	plot.add_element(func, 'Gradenfit')
	plot.finish(preview, 'results/242b.png')

	note_var('a', params[0], unit='(m * A)^2 / V')

	em = 2 / (params[0] * sq(c))
	note_var('e/m', em, unit='C / kg')


	Be = 0.5 * Bs(Ib - Ia)

	note_var('B_E', sum(Be) / len(Be), unit='T')

	table = {
		'r': r,
		'B': B,
		'Be': Be,
	}

	if not preview: write_printable(table, 'results/242b.csv')
