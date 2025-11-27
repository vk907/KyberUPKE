import math
from math import exp, sqrt, erf, ceil
from math import factorial as fac

def gaussian_center_weight(sigma, t):
    """ Weight of the gaussian of std deviation s, on the interval [-t, t]
    :param x: (float)
    :param y: (float)
    :returns: erf( t / (sigma*\sqrt 2) )
    """
    return erf(t / (sigma * sqrt(2.)))


def binomial(x, y):
    """ Binomial coefficient
    :param x: (integer)
    :param y: (integer)
    :returns: y choose x
    """
    try:
        binom = fac(x) // fac(y) // fac(x - y)
    except ValueError:
        binom = 0
    return binom


def centered_binomial_pdf(k, x):
    """ Probability density function of the centered binomial law of param k at x
    :param k: (integer)
    :param x: (integer)
    :returns: p_k(x)
    """
    return binomial(2*k, x+k) / 2.**(2*k)


def build_centered_binomial_law(k):
    """ Construct the binomial law as a dictionnary
    :param k: (integer)
    :param x: (integer)
    :returns: A dictionnary {x:p_k(x) for x in {-k..k}}
    """
    D = {}
    for i in range(-k, k+1):
        D[i] = centered_binomial_pdf(k, i)
    return D

def build_centered_normal(sigma):
    """ Construct the normal law as a dictionary
    :param k: (integer)
    :param x: (integer)
    :returns: A dictionnary {x:p_k(x) for x in {-k..k}}
    """

    kmax = 1
    while True:
        partial_weight = 2*exp(-kmax**2/(2*sigma**2))
        if partial_weight < 2**(-400):
            break
        kmax += 1
    
    weight = 1
    for k in range(kmax, 0, -1):
        partial_weight = 2*exp(-k**2/(2*sigma**2))
        weight += partial_weight
        k += 1

    D = {}
    for i in range(-kmax, kmax+1):
        D[i] = exp(-i**2/(2*sigma**2))/weight

    return D


def mod_switch(x, q, rq):
    """ Modulus switching (rounding to a different discretization of the Torus)
    :param x: value to round (integer)
    :param q: input modulus (integer)
    :param rq: output modulus (integer)
    """
    return int(round(1.* rq * x / q) % rq)


def mod_centered(x, q):
    """ reduction mod q, centered (ie represented in -q/2 .. q/2)
    :param x: value to round (integer)
    :param q: input modulus (integer)
    """
    a = x % q
    if a < q/2:
        return a
    return a - q


def build_mod_switching_error_law(q, rq):
    """ Construct Error law: law of the difference introduced by switching from and back a uniform value mod q
    :param q: original modulus (integer)
    :param rq: intermediate modulus (integer)
    """
    D = {}
    V = {}
    for x in range(q):
        y = mod_switch(x, q, rq)
        z = mod_switch(y, rq, q)
        d = mod_centered(x - z, q)
        D[d] = D.get(d, 0) + 1./q
        V[y] = V.get(y, 0) + 1

    return D


def law_convolution(A, B):
    """ Construct the convolution of two laws (sum of independent variables from two input laws)
    :param A: first input law (dictionnary)
    :param B: second input law (dictionnary)
    """

    C = {}
    for a in A:
        for b in B:
            c = a+b
            C[c] = C.get(c, 0) + A[a] * B[b]
    return C


def law_product(A, B):
    """ Construct the law of the product of independent variables from two input laws
    :param A: first input law (dictionnary)
    :param B: second input law (dictionnary)
    """
    C = {}
    for a in A:
        for b in B:
            c = a*b
            C[c] = C.get(c, 0) + A[a] * B[b]
    return C

def clean_dist(A):
    """ Clean a distribution to accelerate further computation (drop element of the support with proba less than 2^-300)
    :param A: input law (dictionnary)
    """
    B = {}
    for (x, y) in A.items():
        if y>2**(-400):
            B[x] = y
    return B


def iter_law_convolution(A, i):
    """ compute the -ith forld convolution of a distribution (using double-and-add)
    :param A: first input law (dictionnary)
    :param i: (integer)
    """
    D = {0: 1.0}
    i_bin = bin(i)[2:]  # binary representation of n
    for ch in i_bin:
        print("iter", ch)
        D = law_convolution(D, D)
        D = clean_dist(D)
        if ch == '1':
            D = law_convolution(D, A)
            D = clean_dist(D)

    return D


def tail_probability(D, t):
    '''
    Probability that an drawn from D is strictly greater than t in absolute value
    :param D: Law (Dictionnary)
    :param t: tail parameter (integer)
    '''
    s = 0
    ma = max(D.keys())
    if t >= ma:
        return 0
    for i in reversed(range(int(ceil(t)), ma)):  # Summing in reverse for better numerical precision (assuming tails are decreasing)
        s += D.get(i, 0) + D.get(-i, 0)
    return s


def find_tail_for_probability(D, p):
    s = 0
    ma = max(D.keys())
    for i in reversed(range(1, ma)):  # Summing in reverse for better numerical precision (assuming tails are decreasing)
        s += D.get(i, 0) + D.get(-i, 0)
        if s >= p:
            return i-1
    return s

import math

# quick check eta=2 in CBD
bitsec = 128
n, k = 256, 4            # k changed from 3 to 4
NK = n * k               # = 1024
q = 44620801
# For k=3, q_max=455681, k1_max=46202
# For k=4, q_max=44620801, k1_max=332272506
target = 2**(-128) / n   # Common standard for DFR per-poly in Kyber/MLWE

def DFR(k1, q):
    # Single coefficient variance
    var1 = 2 * (k1 + 1) + 1
    # Total variance - changed to 1024
    sigma = math.sqrt(NK * var1)
    z = (q / 4) / sigma
    # Normal tail probability (two-tailed)
    tail = math.erfc(z / math.sqrt(2))
    return n * tail   # n=256 (per-poly failure)

# Binary search for maximum k1
lo, hi = 0, 20000000000
while lo < hi:
    mid = (lo + hi + 1) // 2
    if DFR(mid, q) < 2**(-128):
        lo = mid
    else:
        hi = mid - 1

print("max k1 =", lo)
print("log2(max k1) =", math.log(lo, 2))
# Test example
print("DFR(log2):", math.log(DFR(lo, q), 2))

# Size calculations (uncompressed public key and ciphertext, fixed A/rho)
pk_size = n * k * q.bit_length() / 8 + 32
ct_size = (n * k * q.bit_length() + n * q.bit_length()) / 8
print("public key size (bytes):", pk_size)
print("ciphertext size (bytes):", ct_size)


import math
from math import log
# simulation check tail-bound

def build_centered_binomial(param):
    """Build centered binomial distribution ψ_param"""
    D = {}
    for i in range(-param, param+1):
        if abs(i) <= param:
            binom_coeff = math.comb(2*param, param + i)
            prob = binom_coeff / (4 ** param)
            if prob > 1e-12:
                D[i] = prob
    return D

def multiply_distribution(D, factor):
    """Multiply distribution D by factor"""
    new_D = {}
    for x, prob in D.items():
        new_x = factor * x
        new_D[new_x] = new_D.get(new_x, 0) + prob
    return new_D

# Build distributions (keep your original naming)
CBD3 = build_centered_binomial(2)  # for s and r
CBD2 = build_centered_binomial(2)  # for e and e'

k1=lo 

# ---- Construct S = s0 + s1 + s2... ----
S = build_centered_binomial(2)
for _ in range(k1):
    S = law_convolution(S, CBD3)
    S = clean_dist(S)

# ---- Construct E_pub = e0 + e1 + e2...----
E_pub = build_centered_binomial(2)
for _ in range(k1):
    E_pub = law_convolution(E_pub, CBD3)
    E_pub = clean_dist(E_pub)

# ---- Construct R----
R = CBD3

# ---- Construct Eprime----
Eprime = CBD2

# ---- Additional separate e (the + e in the expression) ----
E_extra = CBD2

# T1 = e * r (corresponding to sum(ei) * r
T1 = law_product(E_pub, R)
T1 = clean_dist(T1)

# T2 = e_extra (direct additional noise)
T2 = E_extra
T2 = clean_dist(T2)

# T3 = e' * s (corresponding to e_r * sum(si))
T3 = law_product(Eprime, S)
T3 = clean_dist(T3)

# noise_coef = T1 + T2 - T3
neg_T3 = { -x: p for x, p in T3.items() }

noise_coef = law_convolution(T1, T2)        # T1 + T2
noise_coef = clean_dist(noise_coef)
noise_coef = law_convolution(noise_coef, neg_T3)  # - T3
noise_coef = clean_dist(noise_coef)

# Now accumulate single coefficient noise over n*k independent coefficients (each coefficient independent)
tot_terms = n * k
noise_total = iter_law_convolution(noise_coef, tot_terms)
noise_total = clean_dist(noise_total)

# Target probability (consistent with your original script)
target_prob = 2**(-128) / n

# Calculate final tail bound
f_final = find_tail_for_probability(noise_total, target_prob)
print(f"Kyber-CPA tail-bound: {f_final}")

# Verify if modulus q satisfies q > 4 * f_final (check from your original script)

required_condition = q > 4 * f_final
print(f"is q={q} enough: {required_condition}")
if not required_condition:
    print(f"The requested q: {4 * f_final}")
    print(f"Recommendation: 12289, 18433, 40961")
else:
    print(f"Rest: {q/4 - f_final}")

print("log2 of n * tail probability at (q/4):", math.log(n * tail_probability(noise_total, q/4), 2))
