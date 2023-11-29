import labtools as tools

def b(preview, data=None):
    data = data['a']['data']

    U_a_index = data['headers']['"Spannung" U_B1 / V']
    U_r_index = data['headers']['"Spannung" U_B2 / V']
    I_a_index = data['headers']['"Strom" I_A1 / A']
    I_b_index = data['headers']['"Strom" I_A2 / A']
    P_a_index = data['headers']['"P_ein" P_1 / W']
    P_b_index = data['headers']['"P_aus" P_2 / W']



    U_a = tools.p.ev([round(v, 1) for v in data['data'][U_a_index]], 0.1)
    U_r = tools.p.ev([round(v, 1) for v in data['data'][U_r_index]], 0.1)
    I = tools.p.ev([round(v, 2) for v in data['data'][I_a_index]], 0.01)
    P_w = tools.p.ev([round(v, 1) for v in data['data'][P_a_index]], 0.5)

    P_w = P_w * -1

    R = U_r / I
    P_s = U_a * I
    cos_phi = U_r / U_a

    plot = tools.Plot(r'R $[\Omega]$', r'P $[W]$')
    plot.mark_x(39.79, r'$R_{W, max}$', 'red')
    plot.mark_y(27.76, r'$P_{W, max}$', 'violet')
    plot.add_element(R, P_w, r'$P_W$')
    plot.add_element(R, P_s, r'$P_S$')
    plot.add_element(R, P_s * cos_phi, r'$P_S \times \cos{\phi}$')


    plot.finish(preview, 'results/238b.png')


    if not preview:
        tools.easyparse.write_printable(
            {'U_1': U_a,
            'I': I,
            'U_R': U_r,
            'P_W': P_w,
            }, 
            'results/rohdaten_a.csv', 3)

        tools.easyparse.write_printable(
            {
            'R': R,
            'P_S': P_s,
            'cos(phi)': cos_phi,
            'P_S * cos_phi': P_s * cos_phi,
            },
            'results/berechnete_a.csv', 3)