# Implementation Differences: BLS12-381 vs BLS12-377

This document describes the key differences between the Tate-based G1 subgroup membership test implementations for BLS12-381 and BLS12-377 curves.

## Overview

Both implementations follow the algorithm from ["Revisiting subgroup membership testing on pairing-friendly curves via the Tate pairing"](https://eprint.iacr.org/2024/1790.pdf) by Dai et al. However, the sign of the curve parameter `z` leads to significant structural differences.

## Curve Parameters

| Parameter | BLS12-381 | BLS12-377 |
|-----------|-----------|-----------|
| z | `-0xd201000000010000` (negative) | `0x8508c00000000001` (positive) |
| e2 = \|z-1\| | `0xd201000000010001` (= \|z\|+1) | `0x8508c00000000000` (= z-1) |
| e2-1 (Miller loop scalar) | `0xd201000000010000` (= \|z\|) | `0x8508bfffffffffff` (= z-2) |

## NAF Representation

The Non-Adjacent Form (NAF) of the Miller loop scalar `e2-1` differs significantly:

### BLS12-381: NAF of |z| = `0xd201000000010000`
```go
var naf = [65]int8{
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // bits 0-15
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // bits 16-31
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // bits 32-47
    1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, -1, 0, 1,  // bits 48-64
}
```
- **MSB**: `naf[64] = 1`
- **LSB**: `naf[0] = 0`
- Loop starts at `i = 63`

### BLS12-377: NAF of z-2 = `0x8508bfffffffffff`
```go
var naf = [65]int8{
    -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 0,
    0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0,
}
```
- **MSB**: `naf[63] = 1`
- **LSB**: `naf[0] = -1` ← Key difference!
- Loop starts at `i = 62`

## Miller Loop Start Index

| Curve | MSB Position | Loop Start |
|-------|--------------|------------|
| BLS12-381 | `naf[64] = 1` | `i = 63` |
| BLS12-377 | `naf[63] = 1` | `i = 62` |

## Special Case: `naf[0] < 0` Handling

This is the most significant implementation difference.

### BLS12-381
Since `naf[0] = 0`, the `naf[0] < 0` code path is **never executed**. The main loop processes all bits down to `i = 0`.

```go
// BLS12-381: This code is never reached
if naf[0] < 0 {
    g1.Sub(&p.X, &tab[k])
    g2.Sub(&p2.X, &tab[k])
    n1.Square(n1)
    d1.Square(d1)
    d1.Mul(d1, &g1)
    n2.Square(n2)
    d2.Square(d2)
    d2.Mul(d2, &g2)
}
```

### BLS12-377
Since `naf[0] = -1`, the main loop stops at `i = 1` (setting `j = 1`), and special handling is required.

**Critical Issue**: At the end of the main loop, the running point `T = [e2/2]Q = (-1, 0)` is a **2-torsion point** (Y = 0). This happens because:
- NAF sum without bit 0 = `z-1 = e2`
- But we've only done half the doublings, so `T = [e2/2]Q`
- For BLS12-377, `[e2/2]Q` is a point with Y = 0

The tangent line at a 2-torsion point is **vertical**, requiring special handling:

```go
// BLS12-377: Handle 2-torsion point case
if naf[0] < 0 {
    // T = [e2/2]Q = (-1, 0) is a 2-torsion point
    // Tangent is vertical: P.X - T.X
    // Doubled point 2T = O (infinity)
    // Subtraction O - Q = -Q contributes nothing (cancels)

    g1.Sub(&p.X, &tab[k-2])  // P.X - T.X where T.X = -1
    g2.Sub(&p2.X, &tab[k-2])

    n1.Square(n1)
    d1.Square(d1)
    d1.Mul(d1, &g1)

    n2.Square(n2)
    d2.Square(d2)
    d2.Mul(d2, &g2)
}
```

## Precomputation Table

### BLS12-381
```go
if naf[0] < 0 {
    tab = append(tab, t0.X)  // Never executed
}
```

### BLS12-377
No additional table entries are needed. The coordinates of the 2-torsion point `(-1, 0)` are already stored by the last iteration of the main loop at `tab[k-2]` and `tab[k-1]`.

```go
if naf[0] < 0 {
    // No extra entries needed - T.X is already at tab[k-2]
}
```

## Final Exponentiation (f2IsOne)

The second final exponentiation formula differs based on the sign of z.

### BLS12-381 (z < 0)
```
exp2 = |z|^5 + |z|^4 - |z|^3 - |z|^2 + |z| - 2
```
Check: `x^{|z|^5 + |z|^4} == x^{|z|^3 + |z|^2 + 2} * x^|z|`

### BLS12-377 (z > 0)
```
exp2 = z^5 - z^4 - z^3 + z^2 + z + 2
```
Check: `x^{z^5 + z^2 + z + 2} == x^{z^4 + z^3}`

```go
// BLS12-377 f2IsOne implementation
func f2IsOne(x *fp.Element) bool {
    var u0, u1, u2, u3, u4, u5 fp.Element
    u0.Square(x)        // x^2
    expBySeed(&u1, x)   // x^z
    expBySeed(&u2, &u1) // x^{z^2}
    expBySeed(&u3, &u2) // x^{z^3}
    expBySeed(&u4, &u3) // x^{z^4}
    expBySeed(&u5, &u4) // x^{z^5}

    // Left side: x^{2 + z + z^2 + z^5}
    u0.Mul(&u0, &u1)
    u0.Mul(&u0, &u2)
    u0.Mul(&u0, &u5)

    // Right side: x^{z^3 + z^4}
    u3.Mul(&u3, &u4)

    return u0.Equal(&u3)
}
```

## Cube Root of Unity

Both curves use a primitive cube root of unity ω for the endomorphism φ(P) = (ωx, y), but the values differ:

| Curve | ω |
|-------|---|
| BLS12-381 | `4002409555221667392624310435006688643935503118305586438271171395842971157480381377015405980053539358417135540939436` |
| BLS12-377 | `80949648264912719408558363140637477264845294720710499478137287262712535938301461879813459410945` |

## Summary Table

| Aspect | BLS12-381 | BLS12-377 |
|--------|-----------|-----------|
| z sign | Negative | Positive |
| NAF LSB | 0 | -1 |
| Loop start | i = 63 | i = 62 |
| naf[0] < 0 executed | No | Yes |
| 2-torsion handling | Not needed | Required |
| Table entries for naf[0]<0 | 1 (unused) | 0 |
| f2IsOne formula | Different signs | Different signs |
