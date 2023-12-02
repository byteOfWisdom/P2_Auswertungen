import labtools as tools
from labtools.defaults import *

def b(preview, data=None):
    if data['accommodate_weirdo']:
        set_setting('localization', 'de_De')

    correction = data['a']['pa_correction']
    data = data['a']['data']

    U_a_index = data['headers']['"Spannung" U_B1 / V']
    U_r_index = data['headers']['"Spannung" U_B2 / V']
    I_a_index = data['headers']['"Strom" I_A1 / A']
    I_b_index = data['headers']['"Strom" I_A2 / A']
    P_a_index = data['headers']['"P_ein" P_1 / W']
    P_b_index = data['headers']['"P_aus" P_2 / W']


    dU_a = data['data'][U_a_index] * 1e-2 + 70 * 0.5e-2
    U_a = ev(data['data'][U_a_index], dU_a)
    dU_r = data['data'][U_r_index] * 1e-2 + 70 * 0.5e-2
    U_r = ev(data['data'][U_r_index], dU_r)
    I = data['data'][I_a_index]
    dI = I * 2e-2 + 21.5 * 0.5e-3
    I = ev(I, dI)
    P_w = ev(data['data'][P_a_index], error(I * U_a)) * correction


    R = U_r / I
    P_s = U_a * I
    cos_phi = U_r / U_a

    U_eff = sum(U_a) / (len(U_a))
    pmax = sq(U_eff) * np.pi * 50 * 80e-6

    rmax = 1 / (2 * np.pi * 50 * 80e-6)


    note_var('P_max', pmax, unit="W")
    note_var('R_max', rmax, unit="Ohm")



    plot = tools.Plot(r'R $[\Omega]$', r'P $[W]$')
    plot.mark_x(rmax, r'$R_{W, max}$', 'red')
    plot.mark_y(pmax, r'$P_{W, max}$', 'violet')
    plot.add_element(R, P_w, r'$P_W$')
    plot.add_element(R, P_s, r'$P_S$')
    plot.add_element(R, P_s * cos_phi, r'$P_S \times \cos{\phi}$')


    plot.finish(preview, 'results/238b.png')


    if not preview:
        tools.easyparse.write_printable(
            {'U_1 [V]': U_a,
            'I [A]': I,
            'U_R [V]': U_r,
            'P_W [W]': P_w,
            }, 
            'results/rohdaten_a.csv', 3)

        tools.easyparse.write_printable(
            {
            'R [Ohm]': R,
            'P_S [W]': P_s,
            'cos(phi)': cos_phi,
            'P_S * cos_phi [W]': P_s * cos_phi,
            },
            'results/berechnete_a.csv', 3)