"""
Exact Riemann Solver

:param U: a list of conservative variables [rho, rho u, rho v, rho w, e]
:param W: a list of primitive variables [rho, u, v, w, p]
"""


from .ideal_gas import pressure, sound_speed, primitive_to_conservative
from math import sqrt
from .constants import G, GM1, GP1, GM1_2G, GP1_2G, GG_GM1, GM1_GP1, G_INV, N_ITERATIONS, DESIRED_ACCURACY
import sys


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


def nonlinear_wave_function(s, p):
    Ak = 2.0 / (GP1 * s['rho'])
    Bk = GM1 / GP1 * s['p']

    if p > s['p']:
        factor = sqrt(Ak / (Bk + p))
        f = factor * (p - s['p'])
        fprime = factor * (1.0 - (p - s['p']) / (2.0 * (Bk + p)))
    else:
        f = 2.0 * s['cs'] / GM1 * ((p / s['p'])**(GM1_2G) - 1.0)
        fprime = 1.0 / (s['rho'] * s['cs']) * (p / s['p'])**(-GP1_2G)

    return f, fprime


def iterate(sol, axis):
    """
    Finding the root of the pressure function

    :param sol: filled _solution_object
    :param axis: direction of the flow (x: 0, y: 1, z:2)

    :return sol: updated _solution_object
    """
    ps = guess(sol['left'], sol['right'], axis)

    found = False

    for p in ps:
        p_prev = 0.0
        sol['star']['p'] = p

        for i in range(N_ITERATIONS):
            sol['star']['left']['f'], sol['star']['left']['fprime'] = \
                nonlinear_wave_function(sol['left'], sol['star']['p'])
            sol['star']['right']['f'], sol['star']['right']['fprime'] = \
                nonlinear_wave_function(sol['right'], sol['star']['p'])

            sol['star']['p'] -= (
                sol['star']['left']['f'] + sol['star']['right']['f']
                + (sol['right']['v'][axis] - sol['left']['v'][axis])
            ) / (sol['star']['left']['fprime'] + sol['star']['right']['fprime'])

            if sol['star']['p'] < 0:
                break

            if abs(sol['star']['p'] - p_prev) / abs(sol['star']['p'] + p_prev) < DESIRED_ACCURACY:
                found = True
                break

            p_prev = sol['star']['p']

        if found:
            break

    if sol['star']['p'] < 0:
        print('Not converged')
        sys.exit(0)

    sol['star']['u'] = 0.5 * (
        (sol['right']['v'][axis] + sol['left']['v'][axis])
        + (sol['star']['right']['f'] - sol['star']['left']['f'])
    )

    ps_pl = sol['star']['p'] / sol['left']['p']
    ps_pr = sol['star']['p'] / sol['right']['p']

    if sol['star']['p'] > sol['left']['p']:
        sol['star']['left']['is_shock'] = True
        sol['star']['left']['shock']['rho'] = sol['left']['rho'] * \
            (GM1_GP1 + ps_pl) / (GM1_GP1 * ps_pl + 1)
        sol['star']['left']['shock']['speed'] = sol['left']['v'][axis] - \
            sol['left']['cs'] * sqrt(GP1_2G * ps_pl + GM1_2G)
    else:
        sol['star']['left']['is_shock'] = False
        sol['star']['left']['fan']['rho'] = sol['left']['rho'] * ps_pl**G_INV
        sol['star']['left']['fan']['cs'] = sol['left']['cs'] * ps_pl**GM1_2G
        sol['star']['left']['fan']['speed']['head'] = sol['left']['v'][axis] * \
            sol['left']['cs']
        sol['star']['left']['fan']['speed']['tail'] = sol['star']['u'] - \
            sol['star']['left']['fan']['cs']

    if sol['star']['p'] > sol['right']['p']:
        sol['star']['right']['is_shock'] = True
        sol['star']['right']['shock']['rho'] = sol['right']['rho'] * \
            (GM1_GP1 + ps_pr) / (GM1_GP1 * ps_pr + 1)
        sol['star']['right']['shock']['speed'] = sol['right']['v'][axis] + \
            sol['right']['cs'] * sqrt(GP1_2G * ps_pr + GM1_2G)
    else:
        sol['star']['right']['is_shock'] = False
        sol['star']['right']['fan']['rho'] = sol['right']['rho'] * ps_pr**G_INV
        sol['star']['right']['fan']['cs'] = sol['right']['cs'] * ps_pr**GM1_2G
        sol['star']['right']['fan']['speed']['head'] = sol['right']['v'][axis] + \
            sol['right']['cs']
        sol['star']['right']['fan']['speed']['tail'] = sol['star']['u'] + \
            sol['star']['right']['fan']['cs']


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
