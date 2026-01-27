# Optimizations Applied to Cubical Tate Pairing

## 1. Karatsuba-like Optimization for cADD (Main Optimization)

**Original cADD formula** required 4 multiplications for products:
```
U = X_P * X_Q       (1M)
V = Z_P * Z_Q       (1M)
t1 = X_P * Z_Q      (1M)
t2 = Z_P * X_Q      (1M)
Vdiff = t1 - t2
X_out = (U - V)² * inv
Z_out = Vdiff²
```
**Cost: 4M + 1M(inv) + 2S = 5M + 2S**

**Optimized formula** uses a Karatsuba-like identity:
```
(X_T + Z_T)(X_R - Z_R) = X_T*X_R - X_T*Z_R + Z_T*X_R - Z_T*Z_R
                       = (X_T*X_R - Z_T*Z_R) - (X_T*Z_R - Z_T*X_R)
                       = (A - B) - Vdiff
```

So we compute:
```
A = X_T * X_R           (1M)
B = Z_T * Z_R           (1M)
sumT = X_T + Z_T        (1A)
diffR = X_R - Z_R       (1A)
P3 = sumT * diffR       (1M)
diff1 = A - B           (1A)
Vdiff = diff1 - P3      (1A)  // This gives X_T*Z_R - Z_T*X_R
X_out = diff1² * inv    (1S + 1M)
Z_out = Vdiff²          (1S)
```
**Cost: 3M + 1M(inv) + 2S = 4M + 2S**

**Savings: 1M per cADD = 20% reduction in multiplications**

## 2. Shared `diffR` Between T1 and T2

Since both T1 and T2 additions use the same R point at each iteration:
- When `bit == 0`: both use R0, so `diffR0 = xR0 - zR0` is computed once
- When `bit == 1`: both use R1, so `diffR1 = xR1 - zR1` is computed once

This saves 1 subtraction per iteration (minor, but free).

## 3. Precomputed `diffR` Values in Ladder State

The `ladderState` struct now stores:
```go
type ladderState struct {
    xR0, zR0   fp.Element  // R0 coordinates
    xR1, zR1   fp.Element  // R1 coordinates
    diffR0     fp.Element  // xR0 - zR0 (precomputed)
    diffR1     fp.Element  // xR1 - zR1 (precomputed)
    bit        uint8
}
```

These are computed once during table initialization, eliminating runtime subtractions.

## 4. Applied to All Code Paths

The Karatsuba optimization was applied to:
- `membershipTestCubicalPrecomputed` (main benchmark path)
- `IsInSubGroupCubicalMontgomery` (Montgomery-native path)
- `cubicalLadderCombined` (used by other functions)

## Performance Summary

| Method | Before | After | Improvement |
|--------|--------|-------|-------------|
| Tate-Cubical | ~72µs | ~60µs | ~17% |
| Tate-Cubical-MontgomeryNative | ~58µs | ~54µs | ~7% |

The cubical approach is still slower than Miller (~45µs) because it fundamentally requires computing two separate accumulator chains (T1 and T2) through 63 iterations, whereas Miller can share more computation between the two evaluation points.

## Benchmark Results (Apple M1)

```
BenchmarkG1IsInSubGroupComparison/GLV-8                             ~42µs
BenchmarkG1IsInSubGroupComparison/Tate-Miller-8                     ~45µs
BenchmarkG1IsInSubGroupComparison/Tate-Miller-Precomputed-8         ~46µs
BenchmarkG1IsInSubGroupComparison/Tate-Cubical-8                    ~60µs
BenchmarkG1IsInSubGroupComparison/Tate-Cubical-MontgomeryNative-8   ~54µs
```

## 5. Differential Addition Chain / NAF Analysis (NOT Applicable)

We investigated using the Non-Adjacent Form (NAF) representation to reduce T additions, as suggested by Costello-Smith (Section 6 of "Montgomery curves and their arithmetic").

**NAF of e₂ - 1:**
```
e₂ - 1 = 0x8508bfffffffffff
NAF:     0+0000+0+0000+00+0-000...000-
Non-zero digits: 6 (vs 51 ones in binary)
```

**Decomposition:**
```
e₂ - 1 = 545932 × 2^44 - 1
```

**Why NAF/Addition Chains don't help for cubical pairing:**

The T updates in the cubical pairing are NOT simply accumulating a scalar - they're accumulating the **pairing Z-coordinate**. Each cADD operation incorporates a factor into Z that corresponds to a line evaluation in the equivalent Miller loop.

Testing confirmed that skipping T updates at NAF=0 positions **breaks the pairing computation**:
- Standard approach: f1IsOne = true ✓
- NAF-skipped approach: f1IsOne = false ✗

**Conclusion:** The Montgomery ladder's T updates are essential for the pairing and cannot be reduced via NAF or other addition chain optimizations. The structure of the scalar (many consecutive 1s giving only 6 NAF digits) is mathematically interesting but doesn't translate to fewer operations for the cubical Tate pairing.

The fundamental issue is that the cubical pairing formula:
```
e_c(Q, P) = Z_{[ℓ]Q+P} / Z_{[ℓ]Q}
```
requires the Z-coordinate to accumulate through ALL ladder steps, not just non-zero NAF positions.

## 6. Why Cubical Tate Cannot Beat Miller Tate for Subgroup Testing

The paper (ePrint 2025/672) shows faster **general pairings**, but the **subgroup membership test** has a specific structure that favors Miller:

### The Subgroup Test Structure
- Requires **two** pairing evaluations: `e(Q, P)` and `e(Q, φ(P))`
- **Same Q**, different evaluation points (P and φ(P))
- GLV endomorphism: φ(P) = (ωx, y)

### Why Miller Wins Here

**Miller (shared loop):**
```
T ← 2T              ← SHARED (expensive curve op, done ONCE)
ℓ ← tangent at T    ← SHARED (computed ONCE)
f₁ ← f₁² · ℓ(P)     ← cheap field evaluation
f₂ ← f₂² · ℓ(φP)    ← cheap field evaluation
```
- Curve operations (doubling T) are **expensive** → done once
- Field operations (evaluating ℓ at P, φP) are **cheap** → done twice

**Cubical:**
```
R ← ladder step     ← Q trajectory
T1 ← T1 + R         ← INDEPENDENT (expensive cADD)
T2 ← T2 + R         ← INDEPENDENT (expensive cADD)
```
- T1 = [k]Q + P and T2 = [k]Q + φ(P) are **different points**
- The Z-coordinate accumulation (which IS the pairing value) depends on P
- **Cannot merge** T1 and T2 operations because they accumulate different Z values

### The Fundamental Mismatch

| Scenario | Cubical Excels | Miller Excels |
|----------|---------------|---------------|
| Single pairing e(Q,P) | ✓ | - |
| Same P, different Q | ✓ | - |
| **Same Q, different P** | ✗ | ✓ ← subgroup test |

The cubical pairing's advantage is in x-only arithmetic for a single computation. But for "same Q, different P", the Miller loop's ability to **share the curve point trajectory** while cheaply evaluating at multiple points is superior.

### Bottom Line

**63 shared doublings + 126 cheap evals** (Miller) beats **63 ladder steps + 126 independent cADDs** (Cubical)
