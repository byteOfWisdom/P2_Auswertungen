from labtools.defaults import *

def c(preview, data=None):
	data = data['c']
	U_a = ev([round(v, 1) for v in data['data']['Spannung_A']], 0.1)
	U_b = ev([round(v, 1) for v in data['data']['Spannung_B']], 0.1)
	I_a = ev([round(v, 2) for v in data['data']['Strom_A']], 0.01)
	I_b = ev([round(v, 2) for v in data['data']['Strom_B']], 0.01)
	P_a = ev([round(v, 1) for v in data['data']['P_ein']], 0.5)
	P_b = ev([round(v, 1) for v in data['data']['P_aus']], 0.5)

	Rv = data['Rv']


	Ps1 = U_a * I_a
	Ps2 = U_b * I_b
	Pv = P_a - P_b
	Pcu = data["R1"] * sq(I_a) +  data["R2"] * sq(I_b)
	Pfe = Pv - Pcu
	eta = P_b / P_a

	#short_circut = I_b == max(I_b)
	#open_circut = I_b == min(I_b)
	short_circut = -1
	open_circut = 1

	table = {
		'U_1': U_a,
		'U_2': U_b,
		'I_1': I_a,
		'I_2': I_b,
		'P_w1': P_a,
		'P_w2': P_b,
	}

	table2 = {
		'P_s1': Ps1,
		'P_s2': Ps2,
		'P_v': Pv,
		'P_cu': Pcu,
		'P_fe': Pfe,
		'eta': eta,
	}

	if not preview:
		write_printable(table, 'results/rohdaten_c.csv', sig_digits=3)
		write_printable(table2, 'results/berechnete_c.csv', sig_digits=3)





	plot = Plot(r'$I_2 [A]$', r'P [W]')
	plot.add_element(I_b, P_a, r'$P_{W, 1}$')
	plot.add_element(I_b, P_b, r'$P_{W, 2}$')
	plot.finish(preview, 'results/238c1.png')


	plot = Plot(r'$I_2 [A]$', r'P [W]')
	plot.add_element(I_b, Pcu, r'$P_{Cu}$')
	plot.add_element(I_b, Pfe, r'$P_{Fe}$')
	plot.add_element(I_b, Pv, r'$P_{V}$')
	plot.finish(preview, 'results/238c2.png')

	from math import inf, nan
	plot = Plot(r'$I_2 [A]$', r'$\eta$')
	plot.add_element(I_b[eta != 0.0], eta[eta != 0.0])
	plot.finish(preview, 'results/238c3.png')


	# aufgabe e

	wL = (U_a / I_a)[open_circut]
	note_var('wL', wL, 'Ohm')
	L = wL / (np.pi * 50)
	note_var('L', L, 'H')



	# aufgabe f

	sigma_one = 1 - sq(I_b[short_circut] / I_a[short_circut])#[0]
	note_var('first sigma', sigma_one)


	sigma_two = 1 - sq(U_b[open_circut] / U_a[open_circut])#[0]
	note_var('second sigma', sigma_two)


	sigma_three = (U_a[short_circut] * I_a[open_circut]) / (U_a[open_circut] * I_a[short_circut])
	note_var('third sigma', sigma_three)


	sigma_four = (U_a / (I_b * wL))[short_circut]
	note_var('fourth sigma', sigma_four)

	# aufgabe g


	R = U_b  / I_b
	temp = 2 * Rv + R
	temp_2 = np.array([wL * sigma_three for _ in temp])
	ml = np.array([U_b[open_circut] / U_a[open_circut] for _ in temp])
	U_frac = (R / (R + Rv)) * ml / np.sqrt(1 + (temp_2 / temp))

	table = {
		'U2 / U1 gemessen': U_b / U_a,
		'U2 / U1 berechnet': U_frac,
	}

	if not preview:
		write_printable(table, 'results/werte_g.csv')

	plot = Plot(r'$I_2 [A]$', r'U [V]')
	plot.add_element(I_b, U_b / U_a, r'Messwerte $\frac{U_2}{U_1}$')
	plot.add_element(I_b, U_frac, r'Berechnete $\frac{U_2}{U_1}$')
	plot.finish(preview, 'results/238g.png')