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
    
    # 1. Variance of sum(pk_i) * r (polynomial multiplication)
    # Each coefficient is a sum of N products,
    # variance ≈ (k1+1) * (η²/4) * N
    term1_per_coeff_variance = (k1 + 1) * (eta**2 / 4) * N*K
    
    # 2. Variance of sum(s_i) * u (polynomial multiplication)
    term2_per_coeff_variance = (k1 + 1) * (eta**2 / 4) * N*K
    
    # 3. Variance of e' (a single noise polynomial)
    term3_per_coeff_variance = eta / 2
    
    # 4. Compression error delta_v – per coefficient variance
    compression_per_coeff_variance = (Q / (2.0 * P))**2 / 3.0
    
    # compress u
    compress_u = 0  # No compression
#     compress_u = (Q / (2.0 * P*2**8))**2 / 3 * N*K*(k1 + 1)

    # Total variance per coefficient (all noise sources)
    variance_per_coefficient = (term1_per_coeff_variance + 
                              term2_per_coeff_variance + 
                              term3_per_coeff_variance + 
                              compression_per_coeff_variance +
                              compress_u)
    
    # Since v is 1-dimensional, total variance equals per-coefficient variance
    # No need to multiply by K (because K=1)
    total_variance = variance_per_coefficient
    
    sigma = math.sqrt(total_variance)
    z = (Q / 4.0) / sigma
    tail_prob = math.erfc(z / math.sqrt(2)) / 2
    
    return 2*N * tail_prob  # polynomial failure rate

def quick_check_1(N, K, Q, P, eta=2, target_bitsec=128):
    """
    Comprehensive parameter check with adjustable eta.
    New condition: DFR × k1 ≤ 2^(-target_bitsec)
    """
    NK = N * K
    target = 2**(-target_bitsec)
    
    print(f"Parameter Analysis:")
    print(f"n={N}, k={K}, q={Q}, eta={eta}, log2Q={math.log(Q,2):.2f}")
    print(f"Target: DFR × k1 ≤ 2^{-target_bitsec}")
    
    # Binary search for max k1 satisfying DFR × k1 ≤ 2^-128
    lo, hi = 1, 10**30
    while lo < hi:
        mid = (lo + hi + 1) // 2
        current_dfr = DFR(mid, Q, P, N, K, eta)
        current_product = current_dfr * mid
        
        if current_product < target:
            lo = mid
        else:
            hi = mid - 1
    
    # Final verification
    final_k1 = lo
    final_dfr = DFR(final_k1, Q, P, N, K, eta)
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

def test_condition():
    """Test various parameter configurations"""
    test_cases = [
        (256, 3, 130561, 2**5),
        (256, 4, 2**21+2**19+2**18, 2**6),
        (256, 5, 2**27, 2**5),
        (256, 6, 2**30, 2**5),
        (256, 7, 2**34, 2**5),
        (256, 8, 2**38, 2**5),
        (256, 9, 2**41, 2**5),
        (256, 10, 2**44+2**43, 2**5),
    ]
    
    for N, K, Q, P in test_cases:
        print("=" * 60)
        k1_max = quick_check_1(N, K, Q, P)
        print(f"For (n={N}, k={K}, q={Q}), p={P}, the maximum k1 satisfying DFR×k1 ≤ 2^-128 is {k1_max}")

if __name__ == "__main__":
    test_condition()

"""
============================================================
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
For (n=256, k=3, q=130561), p=32, the maximum k1 satisfying DFR×k1 ≤ 2^-128 is 2590
============================================================
Parameter Analysis:
n=256, k=4, q=2883584, eta=2, log2Q=21.46
Target: DFR × k1 ≤ 2^-128
Results:
Max k1 = 1122386
log2(Max k1) = 20.098137
Final DFR = 2.62e-45
Final DFR (log2) = -148.098206
DFR × k1 = 2.94e-39
DFR × k1 (log2) = -128.000069
Public key size = 2848.0 bytes
Ciphertext size = 3008.0 bytes
For (n=256, k=4, q=2883584), p=64, the maximum k1 satisfying DFR×k1 ≤ 2^-128 is 1122386
============================================================
Parameter Analysis:
n=256, k=5, q=134217728, eta=2, log2Q=27.00
Target: DFR × k1 ≤ 2^-128
Results:
Max k1 = 1384135668
log2(Max k1) = 30.366338
Final DFR = 2.12e-48
Final DFR (log2) = -158.366338
DFR × k1 = 2.94e-39
DFR × k1 (log2) = -128.000000
Public key size = 4352.0 bytes
Ciphertext size = 4480.0 bytes
For (n=256, k=5, q=134217728), p=32, the maximum k1 satisfying DFR×k1 ≤ 2^-128 is 1384135668
============================================================
Parameter Analysis:
n=256, k=6, q=1073741824, eta=2, log2Q=30.00
Target: DFR × k1 ≤ 2^-128
Results:
Max k1 = 70310789919
log2(Max k1) = 36.033027
Final DFR = 4.18e-50
Final DFR (log2) = -164.033027
DFR × k1 = 2.94e-39
DFR × k1 (log2) = -128.000000
Public key size = 5792.0 bytes
Ciphertext size = 5920.0 bytes
For (n=256, k=6, q=1073741824), p=32, the maximum k1 satisfying DFR×k1 ≤ 2^-128 is 70310789919
============================================================
Parameter Analysis:
n=256, k=7, q=17179869184, eta=2, log2Q=34.00
Target: DFR × k1 ≤ 2^-128
Results:
Max k1 = 14462819057903
log2(Max k1) = 43.717414
Final DFR = 2.03e-52
Final DFR (log2) = -171.717414
DFR × k1 = 2.94e-39
DFR × k1 (log2) = -128.000000
Public key size = 7648.0 bytes
Ciphertext size = 7776.0 bytes
For (n=256, k=7, q=17179869184), p=32, the maximum k1 satisfying DFR×k1 ≤ 2^-128 is 14462819057903
============================================================
Parameter Analysis:
n=256, k=8, q=274877906944, eta=2, log2Q=38.00
Target: DFR × k1 ≤ 2^-128
Results:
Max k1 = 3040692512117796
log2(Max k1) = 51.433321
Final DFR = 9.66e-55
Final DFR (log2) = -179.433321
DFR × k1 = 2.94e-39
DFR × k1 (log2) = -128.000000
Public key size = 9760.0 bytes
Ciphertext size = 9888.0 bytes
For (n=256, k=8, q=274877906944), p=32, the maximum k1 satisfying DFR×k1 ≤ 2^-128 is 3040692512117796
============================================================
Parameter Analysis:
n=256, k=9, q=2199023255552, eta=2, log2Q=41.00
Target: DFR × k1 ≤ 2^-128
Results:
Max k1 = 165126874596713871
log2(Max k1) = 57.196353
Final DFR = 1.78e-56
Final DFR (log2) = -185.196353
DFR × k1 = 2.94e-39
DFR × k1 (log2) = -128.000000
Public key size = 11840.0 bytes
Ciphertext size = 11968.0 bytes
For (n=256, k=9, q=2199023255552), p=32, the maximum k1 satisfying DFR×k1 ≤ 2^-128 is 165126874596713871
============================================================
Parameter Analysis:
n=256, k=10, q=26388279066624, eta=2, log2Q=44.58
Target: DFR × k1 ≤ 2^-128
Results:
Max k1 = 20253915276037617662
log2(Max k1) = 64.134835
Final DFR = 1.45e-58
Final DFR (log2) = -192.134835
DFR × k1 = 2.94e-39
DFR × k1 (log2) = -128.000000
Public key size = 14432.0 bytes
Ciphertext size = 14560.0 bytes
For (n=256, k=10, q=26388279066624), p=32, the maximum k1 satisfying DFR×k1 ≤ 2^-128 is 20253915276037617662
"""
