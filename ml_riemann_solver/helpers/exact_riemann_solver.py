from math import sqrt


def _solution_object():
    return {
        'star': {
            'u': 0,
            'p': 0,
            'right': {
                'is_shock': False,
                'f': 0,
                'fprime': 0,
                'shock': {'rho': 0, 'speed': 0},
                'fan': {
                    'rho': 0,
                    'cs': 0,
                    'speed': {'head': 0, 'tail': 0},
                },
            },
            'left': {
                'is_shock': False,
                'f': 0,
                'fprime': 0,
                'shock': {'rho': 0, 'speed': 0},
                'fan': {
                    'rho': 0,
                    'cs': 0,
                    'speed': {'head': 0, 'tail': 0},
                },
            },
        },
        'left': {'rho': 0, 'v': [0, 0, 0], 'p': 0, 'cs': 0},
        'right': {'rho': 0, 'v': [0, 0, 0], 'p': 0, 'cs': 0},
    }


# Monatomic gas
G = 5. / 3

# Diatomic gas
#  G = 7. / 5

GM1 = G - 1
GP1 = G + 1
GM1_2G = (G - 1) / (2 * G)
GG_GM1 = (2 * G) / (G - 1)
GM1_GP1 = (G - 1) / (G + 1)


def guess(l, r, axis):
    """
    Initial p_star guess of the soluton of Riemann problem

    :param l: primitive variables of the left side
    :param r: primitive variables of the right side
    :param axis: direction of the flow (x: 0, y: 1, z:2)

    :return (p1, p2, p3, p4, p5)
    """

    # Ramses guess
    p1 = (
        r['rho'] * r['cs'] * l['p']
        + l['rho'] * l['cs'] * r['p']
        + r['cs'] * l['cs'] *
        (l['rho'] * l['v'][axis] - r['rho'] * r['v'][axis])
    ) / (
        r['rho'] * r['cs'] + l['rho'] * l['cs']
    )

    # Two rarefactions approximation
    p2 = (
        (l['cs'] + r['cs'] - .5 * GM1(r['v'][axis] - l['v'][axis]))
    ) / (
        (l['cs'] / l['p'])**GM1_2G + (r['cs'] / r['p'])**GM1_2G
    )**GG_GM1

    # Two shocks approximation
    g_l = _g_K(l)
    g_r = _g_K(r)
    p3 = (
        g_l * l['p']
        + g_r * r['p']
        - (r['v'][axis] - l['v'][axis])
    ) / (g_l + g_r)

    # p_PV
    p4 = .5 * (l['p'] + r['p']) \
        - .125 * (r['v'][axis] - l['v'][axis]) \
        * (l['rho'] + r['rho']) * (l['cs'] * r['cs'])

    # Mean
    p5 = .5 * (l['p'] + r['p'])

    return p1, p2, p3, p4, p5


def _g_K(s):
    return sqrt((2. / (GP1 * s['rho'])) / (s['p'] + GM1_GP1 * s['p']))


def solve(l, r, axis):
    """
    Solving a given Riemann problem

    :param l: conservative variables of the left side
    :param r: conservative variables of the right side
    :param axis: direction of the flow (x: 0, y: 1, z:2)

    :return conservative variables at the boundary of the cell
    """
    pass
