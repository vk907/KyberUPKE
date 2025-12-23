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
    lo, hi = 1, 10**100
    while lo < hi:
        mid = (lo + hi + 1) // 2
        current_dfr = DFR(mid, Q, P, N, K, eta)
        current_product = current_dfr * mid
        
        if current_product < target:
            lo = mid
        else:
            hi = mid - 1
    
    # Final verification
    if 2*(lo+1)<Q/4:
        final_k1 = lo
    else:
        final_k1 = Q/8-1
#     final_k1 = lo
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
        (256, 4, 2**22-2**19, 2**5),
#         (256, 5, 2**27, 2**5),
#         (256, 6, 2**30, 2**5),
#         (256, 7, 2**34, 2**5),
#         (256, 8, 2**38, 2**5),
#         (256, 9, 2**41, 2**5),
#         (256, 10, 2**45, 2**5),
#         (256, 11, 2**49, 2**5),
#         (256, 12, 2**53, 2**5),
#         (256, 13, 2**57, 2**5),
#         (256, 14, 2**61, 2**5),
#         (256, 15, 2**65, 2**5),
        (256, 16, 2**67, 2**5),
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
n=256, k=4, q=3670016, eta=2, log2Q=21.81
Target: DFR × k1 ≤ 2^-128
Results:
Max k1 = 458751.0
log2(Max k1) = 18.807352
Final DFR = 1.59e-89
Final DFR (log2) = -294.980098
DFR × k1 = 7.31e-84
DFR × k1 (log2) = -276.172747
Public key size = 2848.0 bytes
Ciphertext size = 2976.0 bytes
For (n=256, k=4, q=3670016), p=32, the maximum k1 satisfying DFR×k1 ≤ 2^-128 is 458751.0
============================================================
Parameter Analysis:
n=256, k=16, q=147573952589676412928, eta=2, log2Q=67.00
Target: DFR × k1 ≤ 2^-128
Results:
Max k1 = 1.8446744073709552e+19
log2(Max k1) = 64.000000
Final DFR = 1.25e-166
Final DFR (log2) = -551.114997
DFR × k1 = 2.31e-147
DFR × k1 (log2) = -487.114997
Public key size = 34336.0 bytes
Ciphertext size = 34464.0 bytes
For (n=256, k=16, q=147573952589676412928), p=32, the maximum k1 satisfying DFR×k1 ≤ 2^-128 is 1.8446744073709552e+19
"""
