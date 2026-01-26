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
