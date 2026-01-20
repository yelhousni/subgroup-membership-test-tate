timization Ideas for Tate-based Subgroup Membership Testing

This document outlines potential optimizations to make the Tate-based method significantly faster than the GLV-based method for BLS12-381 G1 subgroup membership testing.

## Current Performance Gap

- **gnark-crypto (Go) GLV method**: ~41-49 µs
- **Current Tate implementation (Go)**: ~49 µs (with precomputation)
- **Goal**: Achieve 20-40% speedup over GLV method

## Priority 1: Critical Path Optimizations

### 1. Optimize exp₁ Final Exponentiation ⭐⭐⭐ (HIGHEST IMPACT)

**Current Status:** The exp₁ exponentiation computes `x^((p-1)/e₂)` using a generic sliding window method:
- ~311 squarings
- ~70 multiplications
- This is likely the **single biggest bottleneck** (~40-50% of total time)

**Problem:** The exponent `(p-1)/e₂` is treated as an arbitrary 381-bit integer.

**Proposed Solution:** Express `(p-1)/e₂` as a polynomial in the seed parameter `z` with small coefficients.

For BLS12-381:
```
p - 1 = (z - 1)² · (z⁴ - z² + 1) + z
e₂ = |z - 1|
Therefore: (p-1)/e₂ = |z - 1| · (z⁴ - z² + 1) + z/(z-1)
```

**Mathematical Analysis Needed:**
1. Compute the exact closed form of `(p-1)/e₂` in terms of z
2. Find the optimal addition chain
3. Compare against the current 311 squarings

**Expected Impact:**
- If expressible with coefficients ≤ 10: reduce to ~80-120 operations
- **Estimated speedup: 40-50% on final exponentiation**
- **Overall speedup: 25-30%**

**Implementation:**
```go
// Instead of generic exponentiation:
func f1IsOne(x *fp.Element) bool {
    // Current: 311 squarings + 70 mults

    // Optimized: express as polynomial in z
    // (p-1)/e₂ = c₀ + c₁·z + c₂·z² + ...
    var z, z2, z3, z4 fp.Element
    expBySeed(&z, x)      // x^z
    expBySeed(&z2, &z)    // x^(z²)
    expBySeed(&z3, &z2)   // x^(z³)
    expBySeed(&z4, &z3)   // x^(z⁴)

    // Combine according to polynomial structure
    // ... (specific formula once derived)
}
```

### 2. Lazy Reduction in Miller Loop ⭐⭐

**Current Status:** Modular reduction is performed after nearly every field operation.

**Proposed Solution:** Delay reductions using "lazy reduction" technique.

BLS12-381 uses 6 limbs of 64 bits. After operations:
- Addition: result ≤ 2p (1 bit slack)
- Squaring: result ≤ p² (381 bits slack)
- Multiplication: result ≤ p² (381 bits slack)

**Strategy:**
```go
// Track "slack" - how many times over p values can be
type lazyElement struct {
    limbs [6]uint64
    slack uint8  // 0 = fully reduced, 1 = ≤ 2p, 2 = ≤ 4p, etc.
}

// In Miller loop:
n1.lazySquare()   // slack += 1
n1.lazySquare()   // slack += 1
n1.lazyMul(&f1)   // slack = slack * f1.slack
if n1.slack > threshold {
    n1.reduce()   // Only reduce when necessary
}
```

**Expected Impact:**
- Save 20-30% of reductions in Miller loop
- Miller loop is ~30% of total time
- **Overall speedup: 6-9%**

**Implementation Complexity:** Medium (requires careful tracking of slack)

### 3. Exploit Shared Structure Between P and φ̂(P) ⭐⭐

**Observation:** Since `φ̂(P) = (ω²·x, y)`, the y-coordinates are identical.

**Current Code:**
```go
// Two separate subtraction and multiplication sequences
f1.Sub(&p.X, &tab[k+1])       // x - a
f2.Sub(&p2.X, &tab[k+1])      // ω²·x - a
g1.Mul(&f1, &tab[k])          // (x - a) · b
g2.Mul(&f2, &tab[k])          // (ω²·x - a) · b
```

**Optimization 1 - Factor out common terms:**
```go
// Precompute during table generation:
// tab_scaled[k] = tab[k] / ω²
// Then: f2 = ω² · f1 - adjustment

f1.Sub(&p.X, &tab[k+1])           // x - a
f1_scaled.Mul(&f1, &omega_squared) // ω²·(x - a)
// Adjust for the ω²·a vs a difference
f2.Sub(&f1_scaled, &tab_adjustment[k])
```

**Optimization 2 - Shared y-coordinate operations:**
```go
// Any operation involving only y can be computed once
y_term.Sub(&p.Y, &tab[k+2])  // Computed once, used for both paths
```

**Expected Impact:**
- Save 10-15% of field operations in Miller loop
- **Overall speedup: 3-5%**

## Priority 2: Architectural Optimizations

### 4. SIMD Vectorization of Miller Loop ⭐⭐⭐

**Observation:** The Miller loop computes field operations for both P and φ̂(P) independently.

**Proposed Solution:** Use SIMD instructions to process both paths simultaneously.

**For x86-64 with AVX2:**
```go
// Conceptually process [n1, n2] and [d1, d2] as vectors
type Vec2Fp struct {
    v0, v1 fp.Element
}

func (z *Vec2Fp) Square(x *Vec2Fp) {
    // Use AVX2 to compute both squares in parallel
    // __m256i for 256-bit operations
}

// In Miller loop:
var n, d Vec2Fp
n.Square(&n)  // Squares n1 and n2 simultaneously
d.Square(&d)  // Squares d1 and d2 simultaneously
```

**Expected Impact:**
- 1.3-1.8x speedup on Miller loop portion
- Miller loop is ~30% of total time
- **Overall speedup: 10-24%**

**Implementation Complexity:** High (requires assembly/intrinsics, platform-specific)

**Alternatives:**
- Start with ARM NEON for better Go assembly support
- Use compiler auto-vectorization with careful code structuring

### 5. Specialized Multiplication by ω² ⭐

**Observation:** Multiplication by ω² appears frequently in the computation of φ̂(P).

For BLS12-381:
```
ω² = 0x5f19672fdf76ce51ba69c6076a0f77eaddb3a93be6f89688de17d813620a00022e01fffffffefffe
```

**Proposed Solution:** Create a specialized multiplication function exploiting the structure of ω².

```go
// Instead of generic multiplication
func (z *fp.Element) MulByOmegaSquare(x *fp.Element) *fp.Element {
    // Exploit specific bit pattern of ω²
    // ω² has form 2^381 - small_value (approximately)

    // Method 1: Use Montgomery form tricks
    // Method 2: Decompose ω² into sum/difference of powers of 2
    // Example: if ω² ≈ 2^k - c where c is small
    // x·ω² = x·2^k - x·c (shift and subtract)
}
```

**Expected Impact:**
- 2-3x faster than generic multiplication by ω²
- Used ~64 times per test
- **Overall speedup: 3-5%**

**Implementation Complexity:** Medium

### 6. Optimize Field Arithmetic Assembly ⭐⭐

**Current Status:** gnark-crypto has assembly for basic field operations, but not for the specific patterns in the Miller loop.

**Proposed Solution:** Write hand-optimized assembly for hot paths:

1. **Miller loop iteration** - the entire inner loop body
2. **expBySeed** - currently Go code, could be assembly
3. **Combined operations** - fused multiply-add, square-multiply chains

```asm
; Example: Fused square-square-multiply
; Computes x² · x² · y without intermediate reductions
; Input: x in registers, y in memory
; Output: result in registers
miller_loop_iteration:
    ; Square x
    ; Square result
    ; Multiply by y
    ; Single reduction at end
    ret
```

**Expected Impact:**
- 10-20% faster field operations in critical paths
- **Overall speedup: 5-10%**

**Implementation Complexity:** High (assembly expertise required)

## Priority 3: Algorithmic Refinements

### 7. Batch Final Exponentiation Checks ⭐

**Current Approach:** Check `f₁^exp₁ == 1` and `f₂^exp₂ == 1` separately.

**Proposed Solution:** Use combined checking strategies.

**Option 1 - Merged Computation:**
```go
// Instead of:
// f1^exp1 == 1 AND f2^exp2 == 1

// Compute:
// (f1^exp1) · (f2^exp2) == 1
// Then verify one component

// This allows sharing some squarings if gcd(exp1, exp2) > 1
```

**Option 2 - Early Abort:**
```go
// Check the "cheaper" exponentiation first
if !f2IsOne(&n2) {  // exp2 is faster (5 expBySeed calls)
    return false     // Early abort saves exp1 computation
}
return f1IsOne(&n1)
```

**Expected Impact:**
- Option 1: 5-10% if significant sharing possible
- Option 2: 25-30% on rejections (assuming 50% rejection rate = 12-15% average)
- **Overall speedup: 5-15%** (depends on workload)

### 8. Cache-Optimized Precomputation Table ⭐

**Current Status:** Lookup table is ~12 KB (~250 field elements × 48 bytes).

**Analysis:**
- L1 cache: typically 32 KB per core
- L2 cache: typically 256 KB per core
- 12 KB fits comfortably in L1

**Proposed Optimizations:**

**Option 1 - Table Compression:**
```go
// Store only necessary bits
// BLS12-381 field elements are 381 bits, but stored in 384 bits (48 bytes)
// Could compress to 48 bytes → 44 bytes (saving 8%)

type compressedFp struct {
    limbs [5]uint64  // 320 bits
    top   uint64     // 61 bits used, top 3 bits unused
}
```

**Option 2 - On-the-fly Reconstruction:**
```go
// For infrequently-used table entries:
// - Store a seed value
// - Recompute when needed
// Trade: 1-2 operations vs cache miss

if isRareEntry(index) {
    reconstructEntry(seed, index)
} else {
    return precomputed[index]
}
```

**Expected Impact:**
- Marginal for current 12 KB table (already cache-friendly)
- More important if scaling to other curves with larger tables
- **Overall speedup: 1-3%**

## Priority 4: Implementation Details

### 9. Improve NAF Representation ⭐

**Current Status:** NAF (Non-Adjacent Form) has 65 digits for e₂-1.

**Proposed Solution:** Use wNAF (width-w NAF) for denser representation.

```go
// Current NAF: digits in {-1, 0, 1}
// Average density: ~33% non-zero

// wNAF with w=4: digits in {0, ±1, ±3, ±5, ..., ±15}
// Average density: ~20% non-zero
// Requires precomputing multiples of Q: [1Q, 3Q, 5Q, ..., 15Q]
```

**Trade-off:**
- Fewer iterations (20% reduction)
- More table entries (8× larger)
- More complex iteration logic

**Expected Impact:**
- 10-15% fewer Miller loop iterations
- But: table increases from 12 KB to ~96 KB
- May hurt cache performance
- **Net impact: unclear, needs profiling**

### 10. Multi-Exponentiation for Final Exponentiation ⭐

**If exp₁ cannot be optimized to polynomial form:**

Use multi-exponentiation techniques:

```go
// Instead of computing f1^exp1 and f2^exp2 separately,
// use simultaneous multi-exponentiation if there's structure

// Straus algorithm or Pippenger's algorithm for
// f1^exp1 · f2^exp2 (if we reformulate the check)
```

**Expected Impact:**
- 15-25% faster than two separate exponentiations
- **Overall speedup: 10-15%** (if exp₁ optimization fails)

## Summary of Expected Speedups

| Optimization | Difficulty | Estimated Speedup | Priority |
|--------------|-----------|-------------------|----------|
| 1. Optimize exp₁ polynomial | Medium-High | 25-30% | ⭐⭐⭐ Critical |
| 2. Lazy reduction | Medium | 6-9% | ⭐⭐ High |
| 3. Shared P/φ̂(P) structure | Low-Medium | 3-5% | ⭐⭐ High |
| 4. SIMD vectorization | High | 10-24% | ⭐⭐⭐ High |
| 5. Specialized ω² mult | Medium | 3-5% | ⭐ Medium |
| 6. Assembly hot paths | High | 5-10% | ⭐⭐ Medium |
| 7. Batch exponentiations | Low | 5-15% | ⭐ Medium |
| 8. Cache optimization | Low | 1-3% | ⭐ Low |
| 9. wNAF representation | Medium | Unclear | ⭐ Low |
| 10. Multi-exponentiation | Medium-High | 10-15% | ⭐ Fallback |

**Realistic Combined Speedup:**
- Conservative: 1.5-1.7× faster (implementing #1, #2, #3)
- Optimistic: 2.0-2.5× faster (implementing #1-#6)
- **Target: 20-40% faster than GLV method**

## Implementation Roadmap

### Phase 1: Mathematical Analysis (Week 1)
1. Derive exact polynomial form of `(p-1)/e₂` in terms of z
2. Compute optimal addition chain
3. Verify correctness with test vectors

### Phase 2: Quick Wins (Week 2)
1. Implement polynomial exp₁ optimization
2. Add shared P/φ̂(P) structure exploitation
3. Implement early abort in exponentiation checks
4. **Expected: 30-40% speedup**

### Phase 3: Advanced Optimizations (Weeks 3-4)
1. Implement lazy reduction framework
2. Add specialized ω² multiplication
3. Profile and identify remaining bottlenecks
4. **Expected: Additional 15-25% speedup**

### Phase 4: Platform-Specific (Weeks 5-6)
1. SIMD vectorization for Miller loop
2. Hand-optimized assembly for hot paths
3. Platform-specific tuning (x86 vs ARM)
4. **Expected: Additional 10-20% speedup**

## Validation Strategy

For each optimization:
1. **Correctness:** Run full test suite with random inputs
2. **Performance:** Benchmark against baseline and GLV method
3. **Regression:** Ensure no performance loss on other curves
4. **Cross-platform:** Test on x86-64, ARM64, and WASM

## Open Questions

1. **Can `(p-1)/e₂` be expressed with small coefficients in z?**
   - This is the make-or-break optimization
   - Needs algebraic analysis

2. **What is the optimal w for wNAF on BLS12-381?**
   - Depends on cache characteristics
   - Needs profiling on target hardware

3. **Is SIMD worth the portability cost?**
   - Could maintain two implementations: portable and SIMD
   - Similar to how gnark-crypto handles field arithmetic

4. **Should we optimize for the positive case (P ∈ G₁) or negative case?**
   - Early abort helps for rejections
   - Might want different strategies for different workloads

## References

1. Dai et al., "Revisiting subgroup membership testing on pairing-friendly curves via the Tate pairing" (2024)
2. Hamburg, "Fast and compact elliptic-curve cryptography" (2012) - lazy reduction techniques
3. Lim & Lee, "More flexible exponentiation with precomputation" (1994) - multi-exponentiation
4. Bernstein & Lange, "Analysis and optimization of elliptic-curve single-scalar multiplication" (2007)
