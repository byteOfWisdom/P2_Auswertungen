import labtools as tools
from scipy.stats import linregress

def c(preview, data_file=None):
	tools.notes.note('c-f:')
	#data = tools.ep.parse('data/236c.csv')
	data = data_file['c']
	angle = data['skt'] / 2
	angle_err = data['dskt'] / 2
	resistance = data['R']
	res_err = data['dR']

	U0 = tools.p.make_errval(data['U0'], data['dU0'])
	R1 = tools.p.make_errval(data['R1'], data['dR1'])
	R2 = tools.p.make_errval(data['R2'], data['dR2'])


	plot = tools.Plot(r'R [$\Omega$]', r'$\frac{1}{\varphi}$')

	res = linregress(resistance, 1 / angle)
	a = tools.p.make_errval(res.slope, res.stderr)
	b = tools.p.make_errval(res.intercept, res.intercept_stderr)

	tools.notes.note_var('a', a)
	tools.notes.note_var('b', b)
	c1 = (R1 + R2) / (a * U0 * R2)
	Rg = b * c1 * U0 * R2 / (R1 + R2)

	tools.notes.note_var('cI', c1)
	tools.notes.note_var('Rg', Rg)


	plot.add_element(
		resistance, 
		1 / angle, 
		xerr = res_err,
		yerr = angle_err / tools.misc.sq(angle), 
		label='Messdaten')
	
	temp1 = a.value # i don't fucking know why this is needed
	temp2 = b.value
	plot.add_element(lambda x: temp1 * x + temp2, 'Geradenfit', color='red')

	a_tad_modified = False # do not even ask
	if a_tad_modified:
		res = linregress(resistance[2:-2], 1 / angle[2:-2])
		tools.notes.note_var('a', res.slope, res.stderr)
		tools.notes.note_var('b', res.intercept, res.intercept_stderr)
		a = tools.p.make_errval(res.slope, res.stderr)
		b = tools.p.make_errval(res.intercept, res.intercept_stderr)

		c1 = (R1 + R2) / (a * U0 * R2)
		tools.notes.note_var('cI', c1)

		Rg = b * c1 * U0 * R2 / (R1 + R2)
		tools.notes.note_var('Rg', Rg)
		plot.add_element(lambda x: res.slope * x + res.intercept, 'Geradenfit ohne die ersten und letzten Werte')




	plot.title = 'Fig 1: Stromempfindlichkeit'
	plot.finish(preview, 'results/236c.png')