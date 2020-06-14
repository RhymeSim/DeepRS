"""
Exact Riemann Solver

:param U: a list of conservative variables [rho, rho u, rho v, rho w, e]
:param W: a list of primitive variables [rho, u, v, w, p]
"""


from math import sqrt
from .constants import GM1, GP1, GM1_2G, GG_GM1, GM1_GP1
from .ideal_gas import pressure, sound_speed, primitive_to_conservative


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
        [rho, rho u, rho v, rho w, e]
    :param r: conservative variables of the right side
        [rho, rho u, rho v, rho w, e]
    :param axis: direction of the flow (x: 0, y: 1, z:2)

    :return U: conservative variables at the boundary of the cell
    """
    dxdt = 0.0
    vdir = axis + 1

    sol = _solution_object()

    sol['left']['rho'] = max(l[0], 0.0)
    sol['right']['rho'] = max(r[0], 0.0)

    vacuum_right = sol['right']['rho'] <= 0.0 and sol['left']['rho'] > 0.0
    vacuum_left = sol['left']['rho'] <= 0.0 and sol['right']['rho'] > 0.0
    vacuum_both_sides = sol['left']['rho'] <= 0.0 and sol['right']['rho'] <= 0.0

    if vacuum_right:
        sol['left']['v'] = l[1:4] / l[0]
        sol['left']['p'] = max(pressure(l), 0.0)
        sol['left']['cs'] = max(sound_speed(l), 0.0)

        sol['right']['v'] = [0.0, 0.0, 0.0]
        sol['right']['p'] = 0.0
        sol['right']['cs'] = 0.0

        s_starl = sol['left']['v'][vdir] - 2 * sol['left']['cs'] / GM1

        if dxdt > s_starl:
            # W_0
            return [0.0, 0.0, 0.0, 0.0, 0.0]
        elif dxdt > sol['left']['v'][vdir] - sol['left']['cs']:
            # W_Lfan
            return w_kfan(sol['left'], dxdt, axis, is_right=False)
        else:
            # W_L
            return primitive_to_conservative(sol['left']['rho'], sol['left']['v'], sol['left']['p'])
    elif vacuum_left:
        sol['left']['v'] = 0.0
        sol['left']['p'] = 0.0
        sol['left']['cs'] = 0.0

        sol['right']['v'] = r[1:4] / r[0]
        sol['right']['p'] = max(pressure(r), 0.0)
        sol['right']['cs'] = max(sound_speed(r), 0.0)

        s_starr = sol['right']['v'][vdir] - 2 * sol['right']['cs'] / GM1

        if dxdt > sol['right']['v'][vdir] + sol['right']['cs']:
            # W_R
            return primitive_to_conservative(sol['right']['rho'], sol['right']['v'], sol['right']['p'])
        elif dxdt > s_starr:
            # W_Rfan
            return w_kfan(sol['right'], dxdt, axis, is_right=True)
        else:
            # W_0
            return [0.0, 0.0, 0.0, 0.0, 0.0]
    elif vacuum_both_sides:
        # W_0
        return [0.0, 0.0, 0.0, 0.0, 0.0]
    else:
        # Non-vacuum cases
        sol['left']['v'] = l[1:4] / l[0]
        sol['left']['p'] = max(pressure(l), 0.0)
        sol['left']['cs'] = max(sound_speed(l), 0.0)

        sol['right']['v'] = r[1:4] / r[0]
        sol['right']['p'] = max(pressure(r), 0.0)
        sol['right']['cs'] = max(sound_speed(r), 0.0)

        s_starl = sol['left']['v'][vdir] - 2 * sol['left']['cs'] / GM1
        s_starr = sol['right']['v'][vdir] - 2 * sol['right']['cs'] / GM1

        du = sol['right']['v'][vdir] - sol['left']['v'][vdir]
        du_cirt = 2 * (sol['left']['cs'] + sol['right']['cs']) / GM1

        if du_crit < du:
            # Vacuum will be created
            if dxdt > s_starr:
                # W_R0
                if dxdt > sol['right']['v'][vdir] + sol['right']['cs']:
                    # W_R
                    return primitive_to_conservative(sol['right']['rho'], sol['right']['v'], sol['right']['p'])
                else:
                    # W_Rfan
                    return w_kfan(sol['right'], dxdt, axis, is_right=True)
            elif dxdt > s_starl:
                # W_0
                return [0.0, 0.0, 0.0, 0.0, 0.0]
            else:
                # W_L0
                if dxdt > sol['left']['v'][vdir] - sol['left']['cs']:
                    # W_Lfan
                    return w_kfan(sol['left'], dxdt, axis, is_right=False)
                else:
                    # W_L
                    return primitive_to_conservative(sol['left']['rho'], sol['left']['v'], sol['left']['p'])
        else:
            pass
            # iterate
            # sample


def w_kfan(s, dxdt, axis, is_right):
    """
    :param s: left or right state (from _solution_object)
    :param dxdt:
    :param axis: direction of flow (x: 0, y: 1, z: 2)
    :param is_right: if the state in on the right side

    :return u: check the header of this file
    """

    if is_right:
        cs = -s['cs']
    else:
        cs = s['cs']

    vdir = axis + 1

    rho = s['rho'] * (2 / GP1 + GM1_GP1 / cs * (s['v'][vdir] - dxdt))**(GG_GM1)

    v = s['v']
    v[vdir] = 2 / GP1 * (cs + GM1 / 2 + s['v'][vdir] + dxdt)

    p = s['p'] * (2 / GP1) + GM1_GP1 / cs * (s['v'][vdir] - dxdt)**(GG_GM1)

    return primitive_to_conservative(rho, v, p)
