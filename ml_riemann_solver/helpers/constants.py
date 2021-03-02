from astropy.constants import m_p, k_B, u
from astropy.units import cm, Myr, Mpc, K


# Iteration constants
N_ITERATIONS = 10000
DESIRED_ACCURACY = 1e-10


# Monatomic gas
G = 5. / 3

# Diatomic gas
#  G = 7. / 5

G_INV = 1.0 / G
GM1 = G - 1
GP1 = G + 1
GM1_2G = (G - 1) / (2 * G)
GP1_2G = (G + 1) / (2 * G)
GG_GM1 = (2 * G) / (G - 1)
GM1_GP1 = (G - 1) / (G + 1)


kB_amu = k_B.to(m_p / cm**3 * Mpc**5 / Myr**2 / K).value / \
    u.to(m_p / cm**3 * Mpc**3).value
