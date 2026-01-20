# Subgroup Membership Testing Methods for BLS12-381 G1

This document explains two methods for testing whether a point belongs to the correct subgroup G₁ on the BLS12-381 elliptic curve.

## Background

BLS12-381 is a pairing-friendly elliptic curve defined by `y² = x³ + 4` over a prime field 𝔽p. The curve has:
- Field characteristic: p = 381 bits
- Subgroup order: r = 255 bits (prime)
- Cofactor: h₁ = (p + 1 - t)/r where t is the trace

The group E(𝔽p) of rational points has structure:
```
E(𝔽p) ≅ ℤₑ₁ ⊕ ℤₑ₂·r
```
where e₁ and e₂ are uniquely determined with e₁ | e₂.

For BLS12-381:
- e₁ = |z-1|/3 where z is the curve seed parameter
- e₂ = |z-1|
- m = e₂/e₁ = 3

The prime-order subgroup G₁ = E(𝔽p)[r] is the r-torsion subgroup, which is the cryptographically relevant group for pairing-based protocols.

**Why subgroup membership testing matters:** In pairing-based cryptography, malicious inputs that lie outside the correct subgroup can lead to small subgroup attacks, potentially leaking information about secret keys. Therefore, verifying that points belong to G₁ before using them in cryptographic operations is crucial for security.

## Method 1: Scott's GLV-based Method (gnark-crypto IsInSubGroup)

### Overview
This is the current state-of-the-art method, proposed by Michael Scott in 2021 (ePrint 2021/1130). It exploits an efficiently computable endomorphism on curves with j-invariant 0.

### Mathematical Foundation

For BLS12-381, there exists a GLV (Gallant-Lambert-Vanstone) endomorphism φ:
```
φ: (x, y) ↦ (ω·x, y)
```
where ω is a primitive cube root of unity in 𝔽ₚ*: `ω³ = 1 (mod p)`.

The dual endomorphism is:
```
φ̂: (x, y) ↦ (ω²·x, y)
```

These endomorphisms satisfy the characteristic polynomial `X² + X + 1 = 0`.

### The Test

A point P ∈ E(𝔽p) belongs to G₁ if and only if:
```
u₀·P + u₁·φ(P) = O
```
where (u₀, u₁) is a short vector in the 2D lattice:
```
Lφ = {(a, b) ∈ ℤ² | a + b·λ₁ ≡ 0 (mod r)}
```
and λ₁ is the eigenvalue of φ acting on G₁.

The key condition is:
```
gcd(#E(𝔽p), u₀² - u₀u₁ + u₁²) = r
```

### Implementation Details

For BLS12-381, the test simplifies to checking:
```
φ̂²(P) = [-z²]P
```
This is equivalent to:
```
φ²(P) + [z²]P = O
```

**Cost:** ~log₂(r)/2 ≈ 127 group operations (scalar multiplication by z²)

**Advantages:**
- No precomputation required
- No additional storage needed
- Simple to implement
- ω is already available for scalar multiplication optimizations

**Disadvantages:**
- Requires about half of a full scalar multiplication
- Specific to curves with efficiently computable endomorphisms

## Method 2: Tate Pairing-based Method (Dai et al. 2024)

### Overview
This method, presented in "Revisiting subgroup membership testing on pairing-friendly curves via the Tate pairing" (ePrint 2024/1790), uses two small Tate pairings with a shared Miller loop.

### Mathematical Foundation

The **reduced Tate pairing** of order n on E(𝔽ₚ) is:
```
Tn: E(𝔽p)[n] × E(𝔽p)/nE(𝔽p) → μn
    (P, R) ↦ f_{n,P}(R)^((p-1)/n)
```
where:
- E(𝔽p)[n] is the n-torsion subgroup
- μn is the group of n-th roots of unity in 𝔽p*
- fn,P is the Miller function (computed via Miller's algorithm)

**Key property:** The Tate pairing is non-degenerate, meaning if Tn(P, R) = 1 for all P ∈ E(𝔽p)[n], then R ∈ nE(𝔽p).

### The Test (Theorem 3 from the paper)

Let Q be a point of order e₂, and let m̃ be an integer with gcd(e₁, m̃) = 1.

A point P ∈ E(𝔽p) belongs to G₁ if and only if:
```
Tₑ₂(Q, φ̂(P))^m̃ = 1  AND  Tₑ₂(Q, P) = 1
```

For BLS12-381, we can choose m̃ = m = 3 (since gcd(e₁, 3) = 1), which simplifies the final exponentiation.

### Algorithm Structure

**Phase 1: Shared Miller Loop** (Algorithm 4 in the paper)
Both pairings share the same first argument Q, so we can:
1. Compute the Miller function fₑ₂,Q at both P and φ̂(P) simultaneously
2. This reduces the Miller loop iterations from 2×log₂(e₂) to ~log₂(e₂)

**Phase 2: Final Exponentiations**
Compute:
- f₁ = fₑ₂,Q(P)^((p-1)/e₂)
- f₂ = fₑ₂,Q(φ̂(P))^((p-1)/e₁)

Check that f₁ = f₂ = 1.

### Precomputation Optimization

Since Q is a fixed system parameter (a generator of E(𝔽p)[e₂]), we can precompute all line function coefficients needed during the Miller loop:

**Precomputed values** (Algorithm 3):
- Slopes λ_T for each doubling step at T = iQ
- x-coordinates and y-coordinates of intermediate points
- Combined values for doubling-addition steps

This reduces the online computation to:
- Line evaluations at P and φ̂(P)
- Field multiplications and squarings
- No point additions or doublings during the test

**Storage:** ~250 field elements (~12 KB) for the lookup table

### Implementation Details for BLS12-381

**Miller loop length:** log₂(e₂) = log₂(|z-1|) ≈ 64 bits

**NAF representation:** The value e₂-1 is represented in Non-Adjacent Form (NAF) for efficient computation with fewer non-zero digits.

**Final exponentiations:**
- exp₁ = (p-1)/e₂: ~381 bits, computed using sliding window method
- exp₂ = (p-1)/e₁ = |z⁵ - z⁴ - z³ + z² + z + 2|: optimized using addition chains

**Cost (without precomputation):**
- Shared Miller loop: ~64 iterations with dual evaluations
- Two final exponentiations: ~311 squarings + 70 multiplications (exp₁) + ~5 scalar exponentiations by z (exp₂)

**Cost (with precomputation):**
- Shared Miller loop: ~64 iterations, only line evaluations
- Two final exponentiations: same as above
- Total: significantly faster, approaching parity with Method 1

### Advantages
- With precomputation: competitive performance (2.7% faster than Scott's method per paper)
- Applicable to a broader class of curves
- Theoretical interest: novel use of pairings for membership testing

### Disadvantages
- Without precomputation: significantly slower (~50% slower)
- Requires storage for precomputed table (~12 KB)
- More complex implementation
- Requires fixed generator Q to be a system parameter

## Comparison

| Aspect | Scott's Method | Tate Pairing Method |
|--------|---------------|---------------------|
| **Computation** | ~log₂(r)/2 ≈ 127 group ops | ~log₂(e₂) ≈ 64 Miller iterations + 2 exponentiations |
| **Precomputation** | None | Optional, ~250 field elements |
| **Storage** | 0 bytes | 0 bytes (without) / ~12 KB (with) |
| **Complexity** | Simple | Moderate |
| **Performance (BLS12-381)** | Baseline | -50% (no precomp) / +2.7% (with precomp) |
| **Applicability** | Curves with cheap endomorphisms | Pairing-friendly curves with e₂ \| (p-1) |

## Performance Notes

Based on benchmarks in this repository:
- **gnark-crypto (Go)**: Scott's method runs in ~41-49 µs
- **RELIC (C)**: Scott's method runs in ~52 µs (at 3.7 GHz)
- **RELIC (C)**: Tate method with precomputation runs in ~49 µs (at 3.7 GHz)

The performance gap between Go and C implementations suggests:
1. RELIC uses more generic scalar multiplication for Scott's method
2. gnark-crypto uses optimized short addition chains for the seed z

## Conclusion

For BLS12-381:
- **Scott's method** remains the practical choice for general use due to simplicity and zero storage requirements
- **Tate pairing method with precomputation** is theoretically interesting and becomes competitive when:
  - The fixed generator Q can be standardized
  - Storage for precomputation is acceptable
  - The additional implementation complexity is justified

The Tate method becomes more attractive for curves with smaller ρ-values (ρ = log₂(p)/log₂(r)), such as BW13-310 (ρ ≈ 1.17) where it shows 62-110% speedup over Scott's method.

## References

1. Scott, M. "A note on group membership tests for G₁, G₂ and G_T on BLS pairing-friendly curves" (2021), https://eprint.iacr.org/2021/1130
2. Dai, Y., He, D., Koshelev, D., Peng, C., Yang, Z. "Revisiting subgroup membership testing on pairing-friendly curves via the Tate pairing" (2024), https://eprint.iacr.org/2024/1790
3. Dai, Y., Lin, K., Zhao, C.A., Zhou, Z. "Fast subgroup membership testings for G₁, G₂ and G_T on pairing-friendly curves" (2023)
4. Koshelev, D. "Subgroup membership testing on elliptic curves via the Tate pairing" (2023)
