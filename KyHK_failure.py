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

def DFR(k1, Q, P, N, K, eta=2):
    """
    Adjusted DFR calculation for 1-dimensional v (single polynomial)
    
    Args:
        k1: number of additional public keys summed
        Q: modulus
        P: compression parameter  
        N: ring dimension
        eta: noise parameter for CBD
    """
    
    # 1. sum(pk_i) * r ≈ (k1+1) * (η²/4) * N
    term1_per_coeff_variance = (k1 + 1) * (eta**2 / 4) * N*K
    
    # 2. sum(s_i) * u  
    term2_per_coeff_variance = (k1 + 1) * (eta**2 / 4) * N*K
    
    # 3. e' 
    term3_per_coeff_variance = eta / 2
    
    # 4. delta_v 
    compression_per_coeff_variance = (Q / (2.0 * P))**2 / 3.0
    
    # compress u
    compress_u=0 # No compression
#     compress_u=(Q / (2.0 * P*2**8))**2 / 3 * N*K*(k1 + 1)

    
    # all
    variance_per_coefficient = (term1_per_coeff_variance + 
                              term2_per_coeff_variance + 
                              term3_per_coeff_variance + 
                              compression_per_coeff_variance+compress_u)
    
    total_variance = variance_per_coefficient
    
    sigma = math.sqrt(total_variance)
    z = (Q / 4.0) / sigma
    tail_prob = math.erfc(z / math.sqrt(2)) / 2
    
    return 2*N * tail_prob 

def quick_check_1(N, K, Q,P, eta=2, target_bitsec=128):
    """
    Comprehensive parameter check with adjustable eta
    New condition: DFR × k1 ≤ 2^(-target_bitsec)
    """
    NK = N * K
    target = 2**(-target_bitsec)
    
    print(f"Parameter Analysis:")
    print(f"n={N}, k={K}, q={Q}, eta={eta}, log2Q={math.log(Q,2):.2f}")
    print(f"Target: DFR × k1 ≤ 2^{-target_bitsec}")
    
    # Binary search for maximum k1 satisfying DFR × k1 ≤ 2^(-128)
    lo, hi = 1, 10**30  # Start from 1 since k1 >= 1
    while lo < hi:
        mid = (lo + hi + 1) // 2
        current_dfr = DFR(mid, Q,P, N, K, eta)
        current_product = current_dfr * mid
        
        if current_product < target:
            lo = mid
        else:
            hi = mid - 1
    
    # Final verification
    final_k1 = lo
    final_dfr = DFR(final_k1, Q,P, N, K, eta)
    final_product = final_dfr * final_k1
    
    # Size calculations
    Q_bits = math.ceil(math.log2(Q))
    P_bits = math.ceil(math.log2(P))

    pk_size = N * K * Q_bits / 8 + 32
    ct_size = (N * K * Q_bits + N * P_bits) / 8
    print(f"Results:")
    print(f"Max k1 = {final_k1}")
    print(f"log2(Max k1) = {math.log(final_k1, 2):.6f}")
    print(f"Final DFR = {final_dfr:.2e}")
    print(f"Final DFR (log2) = {math.log(final_dfr, 2):.6f}")
    print(f"DFR × k1 = {final_product:.2e}")
    print(f"DFR × k1 (log2) = {math.log(final_product, 2):.6f}")
    print(f"Public key size = {pk_size:.1f} bytes")
    print(f"Ciphertext size = {ct_size:.1f} bytes")
    
    return final_k1


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
n, k = 256, 3            # k changed from 3 to 4
q = 130561
p=2**5
lo=quick_check_1(n,k,q,p)
# lo=64
# quick_check_1(256,5,4294966273)
# quick_check_1(256,6,274877905921) 
# quick_check_1(256,7,35184372088321) 

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

k1 = lo 

R = CBD3
Eprime = CBD2
E_extra = CBD2

E_v=build_mod_switching_error_law(q,2**5)

# T1 = sum(ei) * r (n×k terms)
T1 = law_product(E_pub, R)
T1 = clean_dist(T1)

# T2 = e_extra (n terms - one per coefficient)
T2 = E_extra
T2 = clean_dist(T2)

# T3 = e' * sum(si) (n×k terms)  
T3 = law_product(Eprime, S)
T3 = clean_dist(T3)

# noise_subtotal = T1 - T3 (n×k terms)
neg_T3 = {-x: p for x, p in T3.items()}
noise_subtotal = law_convolution(T1, neg_T3)
noise_subtotal = clean_dist(noise_subtotal)

# Apply n×k convolution to (T1 - T3)
tot_terms_t1_t3 = n * k
noise_subtotal_expanded = iter_law_convolution(noise_subtotal, tot_terms_t1_t3)
noise_subtotal_expanded = clean_dist(noise_subtotal_expanded)

# Finally combine: (T1 - T3) expanded over n×k + T2 expanded over n
noise_total = law_convolution(noise_subtotal_expanded, T2)
noise_total = clean_dist(noise_total)
noise_total = law_convolution(noise_total, E_v)
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

print("log2 of k1*n * tail probability at (q/4):", math.log(lo*n * tail_probability(noise_total, q/4), 2))


"""
Parameter Analysis:
n=256, k=3, q=130561, eta=2, log2Q=16.99
Target: DFR × k1 ≤ 2^-128
Results:
Max k1 = 2590
log2(Max k1) = 11.338736
Final DFR = 1.13e-42
Final DFR (log2) = -139.341873
DFR × k1 = 2.93e-39
DFR × k1 (log2) = -128.003137
Public key size = 1664.0 bytes
Ciphertext size = 1792.0 bytes
iter 1
iter 1
iter 0
iter 0
iter 0
iter 0
iter 0
iter 0
iter 0
iter 0
Kyber-CPA tail-bound: 29066
is q=130561 enough: True
Rest: 3574.25
log2 of k1*n * tail probability at (q/4): -151.2139471988895
"""
