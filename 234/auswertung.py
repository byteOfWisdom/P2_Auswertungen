#!python3

from sys import argv

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Circle

import easyparse
from perror import make_errval as ev
import pointer_diagrams

misc_values = ''

def save(file):
	plt.savefig(file, bbox_inches='tight')
	plt.clf()


def note(s):
	global misc_values
	misc_values += s + "\n"


def write_notes(file):
	global misc_values
	with open(file, 'w') as handle:
		handle.write(misc_values)


def sq(x): return x * x


def percent_error(n):
	return ev(n, n * 1e-2)


def c(preview):
	note("c:")
	data = easyparse.parse('data/234c.csv')
	r = ev(data["r"], 0.1)
	u = ev(data["u"], 1e-4)
	i = ev(data["i"], 1e-5)
	w = ev(data["f"] * 2 * np.pi, 0)

	z = u/i

	wl = (z**2 - r**2) ** 0.5

	L = wl / w

	note("L = " + str(L) + ' H')

	phi = np.arctan(float(wl) / float(r))
	d_phi = (sq(float(r) * (wl.error) / sq(float(z))) + sq(float(wl) * (r.error) / sq(float(z)))) ** 0.5
	note('phi = {} +- {}'.format(round(phi, 3), round(d_phi, 3)))

	labelR = "R=" + str(r) + r"$\Omega$"
	labelZ = "|Z|=" + str(z) + r"$\Omega$"
	labelWL = r"$\omega L$=" + str(wl) + r"$\Omega$"

	plt.axis('scaled')
	plt.grid()
	plt.ylim([-0.1, float(wl) + 0.5])
	plt.xlim([-0.1, float(r) + 0.5])
	plt.ylabel("Im(Z)")
	plt.xlabel("Re(Z)")
	plt.arrow(0, 0, float(r), float(wl), width=0.01, head_width=0.1, head_length=0.1, label=labelZ, color="blue", length_includes_head=True)
	plt.arrow(float(r), 0, 0, float(wl), width=0.01, head_width=0.1, head_length=0.1, label=labelWL, color="red", length_includes_head=True)
	plt.arrow(0, 0, float(r), 0, width=0.01, head_width=0.1, head_length=0.1, label=labelR, color="green", length_includes_head=True)
	plt.legend()
	plt.title('Fig 1: Zeigerdiagram des komplexen Widerstands einer Spule')

	if preview:
		plt.show()
	else:
		plt.savefig('results/234_c.png')
	plt.clf()



def d(preview):
	data = easyparse.parse("data/234d.csv")


	pointer_diagrams.pd_presets()

	plt.xlim([-data["ue"] - 0.01, data["ue"] + 0.01])
	plt.ylim([-0.01 , data["ue"] + 0.1])

	for pair in zip(data["ur"], data["uc"]):
		pointer_diagrams.pointer(0, 0, *pair)

	sc = lambda x, r: (r**2 - x**2) ** 0.5
	xr = np.linspace(-data["ue"], data["ue"], 10000)
	plt.plot(xr, sc(xr, data["ue"]), label=r"Halbkreis mit $r = U_{r1} = U_{r2}$")
	plt.xlabel("Re(U) [V]")
	plt.ylabel("Im(U) [V]")
	plt.legend()
	plt.title('Fig 2: Zeigerdiagram des Phasenschiebers')

	if preview:
		plt.show()
	else:
		plt.savefig("results/234_d.png")
	plt.clf()


def crosses(a, b, c):
	up = a <= c and b >= c
	down = a >= c and b <= c
	return up or down


def interpolate_intersect(x, y, y_cutoff):
	i = 0
	if y[0] > y_cutoff:
		x = np.flip(x)
		y = np.flip(y)

	while y[i] < y_cutoff: i += 1
	if y[i] == y_cutoff: return x[i]

	ydiff = y[i + 1] - y[i]
	xdiff = x[i + 1] - x[i]
	ycd = y_cutoff - y[i]
	return x[i] + xdiff * (ycd / ydiff)


def paired(arr):
	i = 0
	while i < len(arr) - 1:
		i += 1
		yield arr[i - 1], arr[i]


def all_intersects(x, y, y_cutoff):
	intersects =[]
	for pair in paired(y):
		if crosses(*pair, y_cutoff):
			id0 = list(y).index(pair[0])
			
			ydiff = y[id0 + 1] - y[id0]
			xdiff = x[id0 + 1] - x[id0]
			ycd = y_cutoff - y[id0]
			intersects.append(x[id0] + xdiff * (ycd / ydiff))

	return intersects


def db(v, ref):
	val = 20 * np.log10(v / ref)
	return val


def db_err(v, ref):
	err = np.sqrt(sq(20 / (v * np.log(10))) * sq(v * 1e-2) + sq(20 / (ref * np.log(10))) * sq(ref * 1e-2))
	return err


def filter_graph(data, suffix, preview):
	ue = data['ue'] # from khz to hz
	ua = data['ua']
	f = data['f']


	ua_err = lambda x: ev(x, 0.01)
	ua_ = np.array(list(map(ua_err, ua)))
	ue_ = ev(ue, 0.01)
	a_err = ua_ / ue_

	R = 100
	C = 1.5e-6

	forward = lambda x: 20 * np.log10(x/ue)
	inverse = lambda y: (10 ** (y / 20)) * ue

	#plt.loglog(base=10)
	db_values = forward(ua)
	db_err = np.sqrt(sq(20 / (ua * np.log(10))) * sq(ua * 1e-2) + sq(20 / (ue * np.log(10))) * sq(ue * 1e-2))

	plt.xscale('log')
	#plt.yscale('function', functions=(inverse, forward))
	plt.ylim([min(db_values) - 1, max([max(db_values), 0.0]) + 0.5])

	plt.ylabel(r'$dB(A) = 20 \times \log_{10}{(\frac{U_a}{U_e})}$')
	plt.xlabel(r'$\Omega = \frac{\nu}{\nu_{gr}}$')

	if suffix == 'b':
		note('Sperrfilter:')
		f0_index = np.where(ua == min(ua))[0][0]
		f0 = f[f0_index]
		f0_ = ev(f0, f0 * 5e-2)

		q_theo = 1 / (2 * np.pi * f0_ * 1e3 * R * C)
		note('Q_theo = ' + str(q_theo))


		intersects = all_intersects(f, ua, ue / np.sqrt(2))
		df1 = intersects[0]
		df2 = intersects[1]
		df1_ = ev(df1, 25)
		df2_ = ev(df2, 25)

		note('fgr1 = ' + str(df1_) + ' kHz')
		note('fgr2 = ' + str(df2_) + ' kHz')
		note('f0 = ' + str(f0_) + ' kHz')
		q_exp = f0_ / (df2_ - df1_)
		note('Q_exp = ' + str(q_exp))

		plt.vlines([df1 / f0, df2 / f0], min(db_values) - 2, 1, color="green", linestyle='--', label=r'$\nu_{gr}$')
		plt.vlines([1], min(db_values) - 2, 1, color="blue", linestyle='--', label=r'$\nu_{0}$')
		plt.hlines(forward(ue/np.sqrt(2)), -1, max(f / f0) + 1, color='yellow', linestyle='--', label=r'$\frac{U_e}{\sqrt{2}}$')

		plt.errorbar(f / f0, db_values, db_err, fmt='x', label='Messwerte')


	else:	
		if suffix == 'lp': note('Tiefpass Filter:')
		else: note('Hochpass Filter:')
		fgr_theo = 1e-3 / (2 * np.pi * R * C) #1e-3 hz -> khz
		fgr = all_intersects(f, ua, ue / np.sqrt(2))[0]
		fgr_label = r'$ \nu_{gr} =' + str(round(fgr, 3)) + ' kHz$'
		fgr_theo_label = r'$ \nu_{gr, theo} =' + str(round(fgr_theo, 3)) + ' kHz$'

		note("f_gr_theo_" + suffix + " = " + str(round(fgr_theo, 3)) + ' kHz')
		note("f_gr_exp_" + suffix + " = " + str(round(fgr, 3)) + ' kHz')

		plt.errorbar(f / fgr, db_values, db_err,  fmt='x', label='Messwerte')

		plt.hlines(forward(ue/np.sqrt(2)), -1, max(f / fgr) + 1, color='yellow', linestyle='--', label=r'$\frac{U_e}{\sqrt{2}}$')
		plt.vlines(1, min(db_values) - 2, 1, color="green", linestyle='--', label=fgr_label)
		plt.vlines(fgr_theo / fgr, min(db_values) - 2, 1, color="blue", linestyle='--', label=fgr_theo_label)

		plt.xlim([min(f / fgr) - 0.05, max(f / fgr) + 0.05])
		#plt.ylim([0, 1.05])


	plt.legend()

	if suffix == 'lp':
		plt.title(r'Fig 3a: Tiefpass Filter')
	if suffix == 'hp':
		plt.title(r'Fig 3b: Hochpass Filter')
	if suffix == 'b':
		plt.title(r'Fig 3c: Sperrfilter')

	if preview:
		plt.show()
	else:
		plt.savefig('results/234_e_' + suffix + '.png', bbox_inches='tight')
	plt.clf()

	return db_values, db_err


def e(preview):
	note('e:')
	lowpass = easyparse.parse('data/234eLP.csv')
	highpass = easyparse.parse('data/234eHP.csv')
	blocking = easyparse.parse('data/234eB.csv')

	try:
		lp, lpe = filter_graph(lowpass, 'lp', preview)
	except:
		lpe = db_err(lowpass['ua'], lowpass['ue'])
		lp = db(lowpass['ua'], lowpass['ue'])
		plt.plot(lowpass['f'], lp, 'x')
		plt.title('Tiefpass Filter, Daten fehlerhaft')

		if preview:	plt.show()
		else: save('results/234_e_lp.png')

	try:
		hp, hpe = filter_graph(highpass, 'hp', preview)
	except:
		hpe = db_err(highpass['ua'], highpass['ue']) 
		hp = db(highpass['ua'], highpass['ue'])
		plt.plot(highpass['f'], hp, 'x')
		plt.title('Hochpass Filter, Daten fehlerhaft')

		if preview:	plt.show()
		else: save('results/234_e_hp.png')

	try:
		b, be = filter_graph(blocking, 'b', preview)
	except:
		be = db_err(blocking['ua'], blocking['ue'])
		b = db(blocking['ua'], blocking['ue'])
		plt.plot(blocking['f'], b, 'x')
		plt.title('Sperrfilter, Daten fehlerhaft')

		if preview:	plt.show()
		else: save('results/234_e_b.png')

	easyparse.write_printable({
		'dB_Tiefpass': lp, 
		'delta_dB_Hochpass': lpe,
	}, 'results/234e_lp.csv')

	easyparse.write_printable({
		'dB_Hochpass': hp, 
		'delta_dB_highpass': hpe,
	}, 'results/234e_hp.csv')

	easyparse.write_printable({
		'dB_Sperrfilter': b, 
		'delta_dB_Sperrfilter': be,
	}, 'results/234e_b.csv')


def i(preview):
	data = easyparse.parse('data/234i.csv')
	f = data['f']
	u = data['ua']
	#ue = data['ue']
	C = data['c']
	R = data['rl']
	L = data['l']

	f_max_index = np.where(u == max(u))[0][0]
	f_max = f[f_max_index]

	note('i:')
	note('f_max_exp = ' + str(f_max))
	note('U_max = ' + str(ev(max(u), 0.01)) + ' V')
	note('U(f=0) = ' + str(ev(u[1], 0.01)) + ' V')

	intersects = all_intersects(f, u, max(u) / np.sqrt(2))
	f1, f2 = intersects[0], intersects[1]

	f1_ = ev(f1, 25)
	f2_ = ev(f2, 25)
	#ue_ = percent_error(ue)
	u_max = ev(max(u), 0.01)

	delta_f = f2_ - f1_

	#Q aus resonanzueberhoeung:
	Q_u = u_max / ev(u[1], 0.01)
	note('Q aus der Resonanzueberhoehung:')
	note('Q_u = ' + str(Q_u))

	#Q aus der Resonanzbreite
	Q_delta_w = ((f_max / delta_f)**2 - 0.5) ** 0.5
	note('Q aus der Resonanzbreite:')
	note('Q_delta_f = ' + str(Q_delta_w))


	w_max = f_max / (2 * np.pi)
	#Q aus Q = w0 * L / Rl
	w0 = (1 / np.sqrt(L * C))
	#w0 = w_max / ((1 - (1 / (2 * sq(Q_delta_w)))) ** 0.5)
	note('w_0=' + str(round(w0, 3)))
	note('f_0=' + str(round(w0 / (2 * np.pi), 3)))
	Q_3 = w0 * L / R
	note('Q aus der der letzten Formel:')
	note('Q_form = ' + str(round(Q_3, 3)))

	#ev(w_max, 25 * 2 * np.pi)
	w0_q = w_max / (np.sqrt(1 - 1 / (2 * sq(Q_3))))
	note('w0_theo = {} s^-1'.format(str(w0_q)))

	L = 1 / (w0_q * C)
	note('L = {}'.format(str(L)))


	plt.grid()
	plt.ylim([-0.1, max(u) + 0.2])
	plt.xlim([-10, max(f) + 10])
	plt.plot(f, u, 'x', label='Messwerte')
	plt.xlabel(r'$\nu$ [Hz]')
	plt.ylabel(r'$U$ [V]')

	plt.vlines([f_max], 0, max(u) + 1, color='green', linestyle='--', label=r'$\nu_{max} =' + str(f_max) +'Hz$')
	plt.vlines(f1, 0, max(u) + 1, color='orange', linestyle='--', label=r'$\nu_{gr, 1} =' + str(round(f1, 2)) +'Hz$')
	plt.vlines(f2, 0, max(u) + 1, color='red', linestyle='--', label=r'$\nu_{gr, 2} =' + str(round(f2, 2)) +'Hz$')
	plt.hlines([max(u) / (np.sqrt(2))], 0, max(f) + 10, color='yellow', linestyle='--', label=r'$\frac{U_{max}}{\sqrt{2}}$')

	plt.title('Fig 4: Resonanzkurve des Schwingkreises')
	plt.legend()
	
	if preview:
		plt.show()
	else:
		plt.savefig('results/234i.png')
	plt.clf()


def main():
	preview = False

	if '-p' in argv or 'pv' in argv:
		preview = True

	if '-c' in argv[1:]:
		c(preview)
		note('')
	elif '-d' in argv[1:]:
		d(preview)
		note('')
	elif '-e' in argv[1:]:
		e(preview)
		note('')
	elif '-i' in argv[1:]:
		i(preview)
	else:
		c(preview)
		note('')
		d(preview)
		note('')
		e(preview)
		note('')
		i(preview)
		note('')

	if not preview:
		write_notes('results/notes.txt')
	else:
		global misc_values
		print(misc_values)


if __name__ == '__main__':
	main()