from labtools.defaults import *

def c(preview, data=None):
    if data['accommodate_weirdo']:
        set_setting('localization', 'de_De')

    data = data['c']

    U_a_index = data['data']['headers']['"Spannung" U_B1 / V']
    U_b_index = data['data']['headers']['"Spannung" U_B2 / V']
    I_a_index = data['data']['headers']['"Strom" I_A1 / A']
    I_b_index = data['data']['headers']['"Strom" I_A2 / A']
    P_a_index = data['data']['headers']['"P_ein" P_1 / W']
    P_b_index = data['data']['headers']['"P_aus" P_2 / W']

    # the error assumptions are taken from the manual
    dU_a = data['data']['data'][U_a_index] * 1e-2 + 70 * 0.5e-2
    U_a = ev(data['data']['data'][U_a_index], dU_a)

    dU_b = data['data']['data'][U_b_index] * 1e-2 + 70 * 0.5e-2
    U_b = ev(data['data']['data'][U_b_index], dU_b)

    dI_a = data['data']['data'][I_a_index] * 2e-2 + 21.5 * 0.5e-3
    I_a = ev(data['data']['data'][I_a_index], dI_a)

    dI_b = data['data']['data'][I_b_index] * 2e-2 + 21.5 * 0.5e-3
    I_b = ev(data['data']['data'][I_b_index], dI_b)

    dP_a = error(I_a * U_a)
    dP_b = error(I_b * U_b)

    P_a = ev(data['data']['data'][P_a_index]  * data['pa_correction'], dP_a)
    P_b = ev(data['data']['data'][P_b_index] * data['pb_correction'], dP_b)


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
    open_circut = 0

    table = {
        'U_1 [V]': U_a,
        'U_2 [V]': U_b,
        'I_1 [A]': I_a,
        'I_2 [A]': I_b,
        r'P_{w, 1} [W]': P_a,
        r'P_{w, 2} [W]': P_b,
    }

    table2 = {
        r'P_{s, 1} [W]': Ps1,
        r'P_{s, 2} [W]': Ps2,
        'P_v [W]': Pv,
        r'P_{cu} [W]': Pcu,
        r'P_{fe} [W]': Pfe,
        r'\eta': eta,
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


    # task e

    note_var('U_open', U_a[open_circut], unit='V')
    note_var('I_open', I_a[open_circut], unit='A')

    wL = (U_a / I_a)[open_circut]
    note_var('wL', wL, 'Ohm')
    L = wL / (np.pi * 50)
    note_var('L', L, 'H')



    # task f

    sigma_one = 1 - sq(I_b[short_circut] / I_a[short_circut])#[0]
    note_var('first sigma', sigma_one)


    sigma_two = 1 - sq(U_b[open_circut] / U_a[open_circut])#[0]
    note_var('second sigma', sigma_two)


    sigma_three = (U_a[short_circut] * I_a[open_circut]) / (U_a[open_circut] * I_a[short_circut])
    note_var('third sigma', sigma_three)


    sigma_four = (U_a / (I_b * wL))[short_circut]
    note_var('fourth sigma', sigma_four)

    # task g


    R = U_b  / I_b
    temp = 2 * Rv + R
    temp_2 = np.array([wL * sigma_three for _ in temp])
    ml = np.array([U_b[open_circut] / U_a[open_circut] for _ in temp])
    U_frac = (R / (R + Rv)) * ml / np.sqrt(1 + (temp_2 / temp))

    table = {
        r'\frac{U_2}{U_1} gemessen': U_b / U_a,
        r'\frac{U_2}{U_1} berechnet': U_frac,
    }

    if not preview:
        write_printable(table, 'results/werte_g.csv')

    plot = Plot(r'$I_2 [A]$', r'U [V]')
    not_nonsene = error(U_frac) < 1.5
    plot.add_element(I_b, U_b / U_a, r'Messwerte $\frac{U_2}{U_1}$')
    plot.add_element(I_b[not_nonsene], U_frac[not_nonsene], r'Berechnete $\frac{U_2}{U_1}$')
    plot.finish(preview, 'results/238g.png')