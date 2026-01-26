package bls12377

import (
	"sync"

	curve "github.com/consensys/gnark-crypto/ecc/bls12-377"
	"github.com/consensys/gnark-crypto/ecc/bls12-377/fp"
)

// This file implements the G1 subgroup membership test for BLS12-377 using
// the cubical Tate pairing from ePrint 2025/672:
//   "Simpler and Faster Pairings from the Montgomery Ladder"
//   by Pope, Reijnders, Robert, Sferlazza, and Smith.
//
// BLS12-377 curve: y² = x³ + 1 can be transformed to Montgomery form:
//   v² = u³ + Au² + u  where A = -√3
//
// The transformation is:
//   Weierstrass (x, y) → Montgomery (u, v): u = (x + 1)/√3, v = y/√[4]{3}
//   Montgomery (u, v) → Weierstrass (x, y): x = √3·u - 1, y = √[4]{3}·v
//
// The 2-torsion point at Weierstrass (-1, 0) maps to Montgomery (0, 0).
//
// The cubical Tate pairing uses x-only Montgomery ladder arithmetic (cDBL, cADD)
// to compute the pairing as a by-product of scalar multiplication. For odd degree ℓ:
//   e_c(Q, P) = Z_{[ℓ]Q+P} / Z_{[ℓ]Q}
//
// This formula directly produces values that pass the f1IsOne/f2IsOne final
// exponentiation checks when using the CUBICAL ladder operations (cADD takes
// the inverse of x(P-Q) as input, which is the key difference from standard xADD).

// Montgomery curve constants (computed once at init)
var (
	montConstOnce sync.Once
	sqrt3         fp.Element // √3
	invSqrt3      fp.Element // 1/√3
	montA         fp.Element // -√3
	montA24       fp.Element // (A+2)/4
)

// initMontgomeryConstants initializes the Montgomery curve parameters.
func initMontgomeryConstants() {
	montConstOnce.Do(func() {
		var three fp.Element
		three.SetUint64(3)
		sqrt3.Sqrt(&three)
		invSqrt3.Inverse(&sqrt3)
		montA.Neg(&sqrt3)

		// (A+2)/4
		var two, four fp.Element
		two.SetUint64(2)
		four.SetUint64(4)
		montA24.Add(&montA, &two)
		montA24.Div(&montA24, &four)
	})
}

// PointXZ represents a point in x-only projective coordinates (X:Z).
type PointXZ struct {
	X, Z fp.Element
}

// IsInfinity returns true if the point is the point at infinity.
func (p *PointXZ) IsInfinity() bool {
	return p.Z.IsZero()
}

// weierstrassToMontgomeryX converts a Weierstrass x-coordinate to Montgomery x-coordinate.
// u = (x + 1) / √3
func weierstrassToMontgomeryX(xW *fp.Element) fp.Element {
	initMontgomeryConstants()
	var one, xM fp.Element
	one.SetOne()
	xM.Add(xW, &one)
	xM.Mul(&xM, &invSqrt3)
	return xM
}

// montgomeryToWeierstrassX converts a Montgomery x-coordinate to Weierstrass x-coordinate.
// x = √3·u - 1
func montgomeryToWeierstrassX(xM *fp.Element) fp.Element {
	initMontgomeryConstants()
	var one, xW fp.Element
	one.SetOne()
	xW.Mul(xM, &sqrt3)
	xW.Sub(&xW, &one)
	return xW
}

// cubicalXADD computes differential addition in x-only coordinates.
// Given (X_P:Z_P), (X_Q:Z_Q) and 1/x(P-Q), computes (X_{P+Q}:Z_{P+Q}).
// Formulas from ePrint 2025/672.
func cubicalXADD(XP, ZP, XQ, ZQ, ixPQ *fp.Element) (X, Z fp.Element) {
	var V1, V2, sum, diff fp.Element

	// V1 = (X_P - Z_P)(X_Q + Z_Q)
	V1.Sub(XP, ZP)
	sum.Add(XQ, ZQ)
	V1.Mul(&V1, &sum)

	// V2 = (X_P + Z_P)(X_Q - Z_Q)
	V2.Add(XP, ZP)
	diff.Sub(XQ, ZQ)
	V2.Mul(&V2, &diff)

	// X = ixPQ * (V1 + V2)²
	X.Add(&V1, &V2)
	X.Square(&X)
	X.Mul(&X, ixPQ)

	// Z = (V1 - V2)²
	Z.Sub(&V1, &V2)
	Z.Square(&Z)

	return X, Z
}

// cubicalXDBLADD computes doubling and differential addition in x-only coordinates.
// Given (X_P:Z_P), (X_Q:Z_Q) and 1/x(Q-P), sets P = [2]P and Q = P + Q.
// Formulas from ePrint 2025/672.
func cubicalXDBLADD(XP, ZP, XQ, ZQ, iXQP *fp.Element) {
	initMontgomeryConstants()

	var t0, t1, t2, X2P, Z2P, XPQ, ZPQ fp.Element

	t0.Add(XP, ZP)       // X_P + Z_P
	t1.Sub(XP, ZP)       // X_P - Z_P
	X2P.Square(&t0)      // (X_P + Z_P)²
	t2.Sub(XQ, ZQ)       // X_Q - Z_Q
	XPQ.Add(XQ, ZQ)      // X_Q + Z_Q
	t0.Mul(&t0, &t2)     // (X_P + Z_P)(X_Q - Z_Q)
	Z2P.Square(&t1)      // (X_P - Z_P)²
	t1.Mul(&t1, &XPQ)    // (X_P - Z_P)(X_Q + Z_Q)
	t2.Sub(&X2P, &Z2P)   // (X_P + Z_P)² - (X_P - Z_P)²
	X2P.Mul(&X2P, &Z2P)  // (X_P + Z_P)² * (X_P - Z_P)²
	XPQ.Mul(&montA24, &t2) // A24 * t2
	ZPQ.Sub(&t0, &t1)    // (X_P + Z_P)(X_Q - Z_Q) - (X_P - Z_P)(X_Q + Z_Q)
	Z2P.Add(&XPQ, &Z2P)  // A24*t2 + (X_P - Z_P)²
	XPQ.Add(&t0, &t1)    // (X_P + Z_P)(X_Q - Z_Q) + (X_P - Z_P)(X_Q + Z_Q)
	Z2P.Mul(&Z2P, &t2)   // (A24*t2 + (X_P - Z_P)²) * t2
	ZPQ.Square(&ZPQ)
	XPQ.Square(&XPQ)
	XPQ.Mul(&XPQ, iXQP)

	// Update in place
	XP.Set(&X2P)
	ZP.Set(&Z2P)
	XQ.Set(&XPQ)
	ZQ.Set(&ZPQ)
}

// translate applies a translation by a 2-torsion point T = (r:s).
// Used for even-degree Tate pairings.
func translate(P, T *PointXZ) PointXZ {
	var result PointXZ

	if T.Z.IsZero() {
		// T = (X:0) = infinity, translation is identity
		result.X.Set(&P.X)
		result.Z.Set(&P.Z)
		return result
	}

	if T.X.IsZero() {
		// T = (0:Z), translation swaps coordinates
		result.X.Set(&P.Z)
		result.Z.Set(&P.X)
		return result
	}

	// General case: T = (r:s) with r,s ≠ 0
	// X' = r*X - s*Z
	// Z' = s*X - r*Z
	var rX, sZ, sX, rZ fp.Element
	rX.Mul(&T.X, &P.X)
	sZ.Mul(&T.Z, &P.Z)
	sX.Mul(&T.Z, &P.X)
	rZ.Mul(&T.X, &P.Z)
	result.X.Sub(&rX, &sZ)
	result.Z.Sub(&sX, &rZ)

	return result
}

// xDBL computes [2]P in x-only projective coordinates.
func xDBL(XP, ZP *fp.Element) (X2P, Z2P fp.Element) {
	initMontgomeryConstants()

	var sum, diff fp.Element
	sum.Add(XP, ZP)
	diff.Sub(XP, ZP)

	var sumSq, diffSq fp.Element
	sumSq.Square(&sum)
	diffSq.Square(&diff)

	X2P.Mul(&sumSq, &diffSq)

	var t, temp fp.Element
	t.Sub(&sumSq, &diffSq)
	temp.Mul(&montA24, &t)
	temp.Add(&temp, &diffSq)
	Z2P.Mul(&t, &temp)

	return X2P, Z2P
}

// xADD computes P+Q in x-only projective coordinates.
// Given (X_P:Z_P), (X_Q:Z_Q) and x(P-Q), computes (X_{P+Q}:Z_{P+Q}).
func xADD(XP, ZP, XQ, ZQ, xPmQ *fp.Element) (X, Z fp.Element) {
	var U, V fp.Element
	U.Mul(XP, XQ)
	V.Mul(ZP, ZQ)

	var UU, VV fp.Element
	UU.Mul(XP, ZQ)
	VV.Mul(ZP, XQ)

	var diff, sub fp.Element
	diff.Sub(&U, &V)
	sub.Sub(&UU, &VV)

	X.Square(&diff)
	Z.Square(&sub)
	Z.Mul(&Z, xPmQ)

	return X, Z
}

// cubicalLadder computes [n]P and [n]P + Q using the Montgomery ladder.
// Returns (nP, nPQ) in x-only projective coordinates.
// The scalar n is given as bytes (little-endian), with nbitlen bits.
// If divByTwo is true, the bottom bit is skipped (for even-degree pairings).
func cubicalLadder(xP, xQ, xPQ *fp.Element, n []byte, nbitlen int, divByTwo bool) (nP, nPQ PointXZ) {
	initMontgomeryConstants()

	// Find the highest set bit
	highBit := nbitlen - 1
	for highBit >= 0 && ((n[highBit>>3]>>(highBit&7))&1) == 0 {
		highBit--
	}

	// Start index (skip lowest bit if divByTwo)
	start := 0
	if divByTwo {
		start = 1
	}

	if highBit < start {
		// n is effectively 0: return (O, Q)
		nP.X.SetOne()
		nP.Z.SetZero()
		nPQ.X.Set(xQ)
		nPQ.Z.SetOne()
		return nP, nPQ
	}

	// Initialize for the first 1-bit:
	// R0 = P, R1 = 2P
	var xR0, zR0, xR1, zR1 fp.Element
	xR0.Set(xP)
	zR0.SetOne()
	xR1, zR1 = xDBL(&xR0, &zR0)

	// Compute initial T = P + Q using xADD
	// We need x(P - Q) = xPQ (given as parameter, but need to check sign)
	// xPQ is x(P - Q), so we can use it directly
	var xT, zT fp.Element
	var oneZ fp.Element
	oneZ.SetOne()
	xT, zT = xADD(&xR0, &zR0, xQ, &oneZ, xPQ)

	// Process remaining bits from second-highest down to start
	for i := highBit - 1; i >= start; i-- {
		bit := (n[i>>3] >> (i & 7)) & 1

		if bit == 0 {
			// R0, R1, T = 2*R0, R0+R1, T+R0
			var xR0new, zR0new, xR1new, zR1new, xTnew, zTnew fp.Element

			// R1 = R0 + R1 (difference is P)
			xR1new, zR1new = xADD(&xR0, &zR0, &xR1, &zR1, xP)

			// T = T + R0 (difference is Q, since T - R0 = Q initially and stays Q)
			// But we need x(T - R0) = x(Q) since T = [k]P + Q and R0 = [k]P
			xTnew, zTnew = xADD(&xT, &zT, &xR0, &zR0, xQ)

			// R0 = 2*R0
			xR0new, zR0new = xDBL(&xR0, &zR0)

			xR0.Set(&xR0new)
			zR0.Set(&zR0new)
			xR1.Set(&xR1new)
			zR1.Set(&zR1new)
			xT.Set(&xTnew)
			zT.Set(&zTnew)
		} else {
			// R0, R1, T = R0+R1, 2*R1, T+R1
			var xR0new, zR0new, xR1new, zR1new, xTnew, zTnew fp.Element

			// R0 = R0 + R1 (difference is P)
			xR0new, zR0new = xADD(&xR0, &zR0, &xR1, &zR1, xP)

			// T = T + R1 (difference is Q - P = -(P - Q), but x is same as x(P-Q))
			// T - R1 = ([k]P + Q) - ([k+1]P) = Q - P, and x(Q-P) = x(P-Q) = xPQ
			xTnew, zTnew = xADD(&xT, &zT, &xR1, &zR1, xPQ)

			// R1 = 2*R1
			xR1new, zR1new = xDBL(&xR1, &zR1)

			xR0.Set(&xR0new)
			zR0.Set(&zR0new)
			xR1.Set(&xR1new)
			zR1.Set(&zR1new)
			xT.Set(&xTnew)
			zT.Set(&zTnew)
		}
	}

	// R0 = [n]P, T = [n]P + Q
	nP.X.Set(&xR0)
	nP.Z.Set(&zR0)
	nPQ.X.Set(&xT)
	nPQ.Z.Set(&zT)

	return nP, nPQ
}

// condSwap swaps a and b if ctl == 0xFFFFFFFF, does nothing if ctl == 0.
func condSwap(a, b *fp.Element, ctl uint32) {
	var mask fp.Element
	for i := range mask {
		mask[i] = uint64(ctl)
	}
	var diff fp.Element
	for i := range diff {
		diff[i] = (a[i] ^ b[i]) & mask[i]
	}
	for i := range a {
		a[i] ^= diff[i]
		b[i] ^= diff[i]
	}
}

// batchInvert3 computes the inverses of three elements using Montgomery's trick.
func batchInvert3(a, b, c *fp.Element) {
	var t0, t1 fp.Element
	t0.Mul(a, b)
	t1.Mul(&t0, c)
	t1.Inverse(&t1)

	var invC, invB, invA fp.Element
	invC.Mul(&t0, &t1)
	t1.Mul(&t1, c)
	invB.Mul(a, &t1)
	invA.Mul(b, &t1)

	a.Set(&invA)
	b.Set(&invB)
	c.Set(&invC)
}

// millerLoopSingle computes the Miller function f_{e₂,Q}(P) for a single evaluation point.
// This is a direct implementation without precomputation, suitable for the membership test.
func millerLoopSingle(q, p *curve.G1Affine) (num, den fp.Element) {
	i := 62
	j := 0
	if naf[0] < 0 {
		j = 1
	}

	// Initialize: f = (x_P - x_Q) / 1
	num.Sub(&p.X, &q.X)
	den.SetOne()

	var T curve.G1Affine
	T.Set(q)

	var qNeg curve.G1Affine
	qNeg.Neg(q)

	var f, g fp.Element
	var u0, u1 fp.Element

	yPmyQ := new(fp.Element)
	xPmxQ := new(fp.Element)
	yPmyQ.Sub(&p.Y, &q.Y)
	xPmxQ.Sub(&p.X, &q.X)

	for i >= j {
		if naf[i] == 0 && i > j {
			// Double-double step: process two consecutive 0-bits
			// First tangent slope at T
			u0.Square(&T.X)
			u1.Double(&u0)
			u0.Add(&u0, &u1) // 3*T.X²
			u1.Double(&T.Y)
			u1.Inverse(&u1)
			u1.Neg(&u1) // -1/(2*T.Y)
			var lambda1 fp.Element
			lambda1.Mul(&u0, &u1)

			T.Double(&T)

			// Second tangent slope at 2T
			u0.Square(&T.X)
			u1.Double(&u0)
			u0.Add(&u0, &u1)
			u1.Double(&T.Y)
			u1.Inverse(&u1)
			var lambda2 fp.Element
			lambda2.Mul(&u0, &u1)

			// Line evaluations at P
			u0.Sub(&p.Y, &T.Y)
			f.Sub(&p.X, &T.X)
			g.Mul(&f, &lambda1)
			g.Sub(&u0, &g)
			f.Mul(&f, &lambda2)
			f.Sub(&u0, &f)

			num.Square(&num)
			num.Square(&num)
			num.Mul(&num, &f)
			den.Square(&den)
			den.Mul(&den, &g)
			den.Square(&den)

			T.Double(&T)
			i--

			if naf[i] > 0 {
				u0.Sub(&T.Y, &q.Y)
				u1.Sub(&T.X, &q.X)
				u1.Inverse(&u1)
				lambda1.Mul(&u0, &u1)

				f.Mul(&lambda1, xPmxQ)
				f.Sub(yPmyQ, &f)
				g.Sub(&p.X, &T.X)

				T.Add(&T, q)
				num.Mul(&num, &f)
				den.Mul(&den, &g)
			}

			if naf[i] < 0 {
				u0.Add(&T.Y, &q.Y)
				u1.Sub(&q.X, &T.X)
				u1.Inverse(&u1)
				lambda1.Mul(&u0, &u1)

				T.Add(&T, &qNeg)
				f.Sub(&p.X, &T.X)
				g.Mul(&lambda1, xPmxQ)
				g.Sub(yPmyQ, &g)

				num.Mul(&num, &f)
				den.Mul(&den, &g)
			}
			i--
			continue
		}

		if naf[i] == 1 || naf[i] == -1 {
			// Add/sub-then-double step
			xT := T.X
			yT := T.Y

			var lambda1, lambda2 fp.Element
			var Tprime curve.G1Affine

			if naf[i] > 0 {
				lambda1.Sub(&T.Y, &q.Y)
				u0.Sub(&T.X, &q.X)
				u0.Inverse(&u0)
				lambda1.Mul(&lambda1, &u0)
				Tprime.Add(&T, q)
			} else {
				lambda1.Add(&T.Y, &q.Y)
				u0.Sub(&T.X, &q.X)
				u0.Inverse(&u0)
				lambda1.Mul(&lambda1, &u0)
				Tprime.Sub(&T, q)
			}

			lambda2.Sub(&Tprime.Y, &yT)
			u0.Sub(&Tprime.X, &xT)
			u0.Inverse(&u0)
			lambda2.Mul(&lambda2, &u0)

			var combined1, combined2 fp.Element
			combined1.Mul(&lambda1, &lambda2)
			combined2.Add(&lambda1, &lambda2)
			combined1.Add(&combined1, &xT)
			combined1.Add(&combined1, &Tprime.X)

			u0.Sub(&p.X, &xT)
			u1.Add(&p.X, &combined1)
			f.Mul(&u0, &u1)
			g.Sub(&p.Y, &yT)
			u1.Mul(&g, &combined2)
			f.Sub(&f, &u1)

			g.Sub(&p.X, &xT)
			if naf[i] < 0 {
				g.Mul(&g, xPmxQ)
			}

			T.Add(&Tprime, &T)

			num.Square(&num)
			num.Mul(&num, &f)
			den.Mul(&den, &g)
			den.Square(&den)

			i--
			continue
		}

		// Simple double step
		u0.Square(&T.X)
		u1.Double(&u0)
		u0.Add(&u0, &u1)
		u1.Double(&T.Y)
		u1.Inverse(&u1)
		u1.Neg(&u1)
		var lambda fp.Element
		lambda.Mul(&u0, &u1)

		// Double T first
		T.Double(&T)

		// Compute evaluations using 2T
		f.Sub(&p.X, &T.X)
		u0.Sub(&p.Y, &T.Y)
		g.Mul(&lambda, &f)
		g.Sub(&u0, &g)

		num.Square(&num)
		num.Mul(&num, &f)
		den.Square(&den)
		den.Mul(&den, &g)

		i--
	}

	// Final step for naf[0] < 0: T is at 2-torsion point (-1, 0)
	if naf[0] < 0 {
		var one fp.Element
		one.SetOne()
		g.Add(&p.X, &one)
		num.Square(&num)
		den.Square(&den)
		den.Mul(&den, &g)
	}

	return num, den
}

// ladderState stores the ladder state at a given iteration
type ladderState struct {
	xR0, zR0 fp.Element // Current R0 = [k]Q
	xR1, zR1 fp.Element // Current R1 = [k+1]Q
	bit      uint8      // The bit value at this iteration
}

// Precomputed cubical table for Q with full ladder trajectory
type cubicalTable struct {
	q       curve.G1Affine
	xQ      fp.Element    // Montgomery x-coordinate of Q
	invXQ   fp.Element    // 1/xQ
	states  []ladderState // Ladder states at each iteration (63 entries)
	finalNQ PointXZ       // Final [n]Q
}

var (
	cubicalTableOnce sync.Once
	cubicalPrecomp   cubicalTable
)

// precomputeCubicalTable generates the precomputed table for cubical pairing.
// This precomputes the entire ladder trajectory for Q, so at runtime we only
// need to compute the T1 and T2 additions.
func precomputeCubicalTable() *cubicalTable {
	cubicalTableOnce.Do(func() {
		initMontgomeryConstants()

		tab := precomputeTableDefault()
		cubicalPrecomp.q = tab.q
		cubicalPrecomp.xQ = weierstrassToMontgomeryX(&tab.q.X)
		cubicalPrecomp.invXQ.Inverse(&cubicalPrecomp.xQ)

		// Scalar: e₂ - 1 = 0x8508bfffffffffff
		n := []byte{0xff, 0xff, 0xff, 0xff, 0xff, 0xbf, 0x08, 0x85}

		// Initialize ladder: R0 = Q, R1 = 2Q
		var xR0, zR0, xR1, zR1 fp.Element
		xR0.Set(&cubicalPrecomp.xQ)
		zR0.SetOne()
		xR1, zR1 = cDBL(&xR0, &zR0)

		// Store states for each iteration (bits 62 down to 0)
		cubicalPrecomp.states = make([]ladderState, 63)

		for i := 62; i >= 0; i-- {
			bit := (n[i>>3] >> (i & 7)) & 1
			idx := 62 - i

			// Store current state BEFORE the update
			cubicalPrecomp.states[idx].xR0.Set(&xR0)
			cubicalPrecomp.states[idx].zR0.Set(&zR0)
			cubicalPrecomp.states[idx].xR1.Set(&xR1)
			cubicalPrecomp.states[idx].zR1.Set(&zR1)
			cubicalPrecomp.states[idx].bit = bit

			// Perform ladder step
			if bit == 0 {
				xR1, zR1 = cADD(&xR0, &zR0, &xR1, &zR1, &cubicalPrecomp.invXQ)
				xR0, zR0 = cDBL(&cubicalPrecomp.states[idx].xR0, &cubicalPrecomp.states[idx].zR0)
			} else {
				xR0, zR0 = cADD(&xR0, &zR0, &xR1, &zR1, &cubicalPrecomp.invXQ)
				xR1, zR1 = cDBL(&cubicalPrecomp.states[idx].xR1, &cubicalPrecomp.states[idx].zR1)
			}
		}

		// Store final [n]Q
		cubicalPrecomp.finalNQ.X.Set(&xR0)
		cubicalPrecomp.finalNQ.Z.Set(&zR0)
	})
	return &cubicalPrecomp
}

// IsInSubGroupTateCubical checks subgroup membership using cubical Tate pairing.
// This uses the cubical Montgomery ladder and x-only arithmetic from ePrint 2025/672.
func IsInSubGroupTateCubical(p *curve.G1Affine) bool {
	if p.IsInfinity() {
		return false
	}
	if !p.IsOnCurve() {
		return false
	}

	tab := precomputeCubicalTable()
	return membershipTestCubicalPrecomputed(p, tab)
}

// membershipTestCubicalPrecomputed uses precomputed ladder trajectory for Q.
// This only computes T1 and T2 additions at runtime, skipping all Q-only operations.
func membershipTestCubicalPrecomputed(p *curve.G1Affine, tab *cubicalTable) bool {
	initMontgomeryConstants()

	var thirdRootOneG1 fp.Element
	thirdRootOneG1.SetString("80949648264912719408558363140637477264845294720710499478137287262712535938301461879813459410945")

	// Convert P and φ(P) to Montgomery
	xP := weierstrassToMontgomeryX(&p.X)
	var p2X fp.Element
	p2X.Mul(&p.X, &thirdRootOneG1)
	xP2 := weierstrassToMontgomeryX(&p2X)

	// Compute x(Q - P) and x(Q - φ(P))
	var qMinusP, qMinusP2 curve.G1Affine
	var pNeg, p2Neg curve.G1Affine
	pNeg.X.Set(&p.X)
	pNeg.Y.Neg(&p.Y)
	p2Neg.X.Set(&p2X)
	p2Neg.Y.Neg(&p.Y)
	qMinusP.Add(&tab.q, &pNeg)
	qMinusP2.Add(&tab.q, &p2Neg)
	xQmP := weierstrassToMontgomeryX(&qMinusP.X)
	xQmP2 := weierstrassToMontgomeryX(&qMinusP2.X)

	// Batch inversion for 4 elements
	var invXP, invXP2, invXQmP, invXQmP2 fp.Element
	var t0, t01, t012, t0123 fp.Element
	t0.Set(&xP)
	t01.Mul(&t0, &xP2)
	t012.Mul(&t01, &xQmP)
	t0123.Mul(&t012, &xQmP2)
	t0123.Inverse(&t0123)

	invXQmP2.Mul(&t012, &t0123)
	t0123.Mul(&t0123, &xQmP2)
	invXQmP.Mul(&t01, &t0123)
	t0123.Mul(&t0123, &xQmP)
	invXP2.Mul(&t0, &t0123)
	t0123.Mul(&t0123, &xP2)
	invXP.Set(&t0123)

	// Initialize T1 = Q + P, T2 = Q + P2
	var xT1, zT1, xT2, zT2 fp.Element
	var oneZ fp.Element
	oneZ.SetOne()

	// Use precomputed first state (R0 = Q)
	xT1, zT1 = cADD(&tab.states[0].xR0, &tab.states[0].zR0, &xP, &oneZ, &invXQmP)
	xT2, zT2 = cADD(&tab.states[0].xR0, &tab.states[0].zR0, &xP2, &oneZ, &invXQmP2)

	// Temporaries
	var U, V, t1, t2, Vdiff fp.Element
	var xT1new, zT1new, xT2new, zT2new fp.Element

	// Process using precomputed ladder states
	for idx := 0; idx < 63; idx++ {
		state := &tab.states[idx]

		if state.bit == 0 {
			// T1 += R0, T2 += R0 (diff = P, P2)
			U.Mul(&xT1, &state.xR0)
			V.Mul(&zT1, &state.zR0)
			t1.Mul(&xT1, &state.zR0)
			t2.Mul(&zT1, &state.xR0)
			Vdiff.Sub(&t1, &t2)
			xT1new.Sub(&U, &V)
			xT1new.Square(&xT1new)
			xT1new.Mul(&xT1new, &invXP)
			zT1new.Square(&Vdiff)

			U.Mul(&xT2, &state.xR0)
			V.Mul(&zT2, &state.zR0)
			t1.Mul(&xT2, &state.zR0)
			t2.Mul(&zT2, &state.xR0)
			Vdiff.Sub(&t1, &t2)
			xT2new.Sub(&U, &V)
			xT2new.Square(&xT2new)
			xT2new.Mul(&xT2new, &invXP2)
			zT2new.Square(&Vdiff)
		} else {
			// T1 += R1, T2 += R1 (diff = Q-P, Q-P2)
			U.Mul(&xT1, &state.xR1)
			V.Mul(&zT1, &state.zR1)
			t1.Mul(&xT1, &state.zR1)
			t2.Mul(&zT1, &state.xR1)
			Vdiff.Sub(&t1, &t2)
			xT1new.Sub(&U, &V)
			xT1new.Square(&xT1new)
			xT1new.Mul(&xT1new, &invXQmP)
			zT1new.Square(&Vdiff)

			U.Mul(&xT2, &state.xR1)
			V.Mul(&zT2, &state.zR1)
			t1.Mul(&xT2, &state.zR1)
			t2.Mul(&zT2, &state.xR1)
			Vdiff.Sub(&t1, &t2)
			xT2new.Sub(&U, &V)
			xT2new.Square(&xT2new)
			xT2new.Mul(&xT2new, &invXQmP2)
			zT2new.Square(&Vdiff)
		}

		xT1, xT1new = xT1new, xT1
		zT1, zT1new = zT1new, zT1
		xT2, xT2new = xT2new, xT2
		zT2, zT2new = zT2new, zT2
	}

	// Check for zero values
	if tab.finalNQ.Z.IsZero() || zT1.IsZero() || zT2.IsZero() {
		return false
	}

	// e_c(Q, P) = Z_{[ℓ]Q+P} / Z_{[ℓ]Q}
	var n1, n2, invNQZ fp.Element
	invNQZ.Inverse(&tab.finalNQ.Z)
	n1.Mul(&zT1, &invNQZ)
	n2.Mul(&zT2, &invNQZ)

	return f1IsOne(&n1) && f2IsOne(&n2)
}

// membershipTestCubicalOptimized uses precomputed Q values for faster computation.
func membershipTestCubicalOptimized(p *curve.G1Affine, tab *cubicalTable) bool {
	initMontgomeryConstants()

	var thirdRootOneG1 fp.Element
	thirdRootOneG1.SetString("80949648264912719408558363140637477264845294720710499478137287262712535938301461879813459410945")

	// Convert P and φ(P) to Montgomery x-coordinates
	xP := weierstrassToMontgomeryX(&p.X)
	var p2X fp.Element
	p2X.Mul(&p.X, &thirdRootOneG1)
	xP2 := weierstrassToMontgomeryX(&p2X)

	// Compute x(Q - P) and x(Q - φ(P))
	var qMinusP, qMinusP2 curve.G1Affine
	var pNeg, p2Neg curve.G1Affine
	pNeg.X.Set(&p.X)
	pNeg.Y.Neg(&p.Y)
	p2Neg.X.Set(&p2X)
	p2Neg.Y.Neg(&p.Y)
	qMinusP.Add(&tab.q, &pNeg)
	qMinusP2.Add(&tab.q, &p2Neg)
	xQmP := weierstrassToMontgomeryX(&qMinusP.X)
	xQmP2 := weierstrassToMontgomeryX(&qMinusP2.X)

	// Batch inversion for 4 elements (invXQ is precomputed)
	var invXP, invXP2, invXQmP, invXQmP2 fp.Element
	var t0, t01, t012, t0123 fp.Element
	t0.Set(&xP)
	t01.Mul(&t0, &xP2)
	t012.Mul(&t01, &xQmP)
	t0123.Mul(&t012, &xQmP2)
	t0123.Inverse(&t0123)

	invXQmP2.Mul(&t012, &t0123)
	t0123.Mul(&t0123, &xQmP2)
	invXQmP.Mul(&t01, &t0123)
	t0123.Mul(&t0123, &xQmP)
	invXP2.Mul(&t0, &t0123)
	t0123.Mul(&t0123, &xP2)
	invXP.Set(&t0123)

	// Combined ladder with precomputed Q values
	nQ, nQP, nQP2 := cubicalLadderCombined(&tab.xQ, &xP, &xP2, &xQmP, &xQmP2,
		&tab.invXQ, &invXP, &invXP2, &invXQmP, &invXQmP2)

	if nQ.Z.IsZero() || nQP.Z.IsZero() || nQP2.Z.IsZero() {
		return false
	}

	var n1, n2 fp.Element
	n1.Set(&nQP.Z)
	n2.Set(&nQP2.Z)
	nQ.Z.Inverse(&nQ.Z)
	n1.Mul(&n1, &nQ.Z)
	n2.Mul(&n2, &nQ.Z)

	return f1IsOne(&n1) && f2IsOne(&n2)
}

// MontgomeryPoint represents a point in x-only Montgomery coordinates.
type MontgomeryPoint struct {
	X fp.Element
}

// Precomputed table for Montgomery-native cubical pairing
type cubicalTableMontgomery struct {
	xQ    fp.Element // x-coordinate of Q in Montgomery form
	invXQ fp.Element // 1/xQ
}

var (
	cubicalMontTableOnce sync.Once
	cubicalMontPrecomp   cubicalTableMontgomery
)

// precomputeCubicalMontgomeryTable generates precomputed values for Montgomery-native input.
func precomputeCubicalMontgomeryTable() *cubicalTableMontgomery {
	cubicalMontTableOnce.Do(func() {
		initMontgomeryConstants()
		tab := precomputeTableDefault()
		cubicalMontPrecomp.xQ = weierstrassToMontgomeryX(&tab.q.X)
		cubicalMontPrecomp.invXQ.Inverse(&cubicalMontPrecomp.xQ)
	})
	return &cubicalMontPrecomp
}

// IsInSubGroupCubicalMontgomery checks subgroup membership when the point is already
// in Montgomery x-coordinate form. This is the fastest version when no coordinate
// conversion is needed.
//
// Parameters:
//   - xP: Montgomery x-coordinate of P
//   - xP2: Montgomery x-coordinate of φ(P) = (ω·x_W, y) where ω is cube root of unity
//   - xQmP: Montgomery x-coordinate of Q - P
//   - xQmP2: Montgomery x-coordinate of Q - φ(P)
//
// Note: The GLV endomorphism φ(P) = (ω·x_W, y) in Weierstrass becomes
// x_M(φ(P)) = (ω·x_W + 1)/√3 in Montgomery, NOT simply ω·x_M(P).
func IsInSubGroupCubicalMontgomery(xP, xP2, xQmP, xQmP2 *fp.Element) bool {
	initMontgomeryConstants()
	tab := precomputeCubicalTable()

	// Batch inversion for 4 elements
	var invXP, invXP2, invXQmP, invXQmP2 fp.Element
	var t0, t01, t012, t0123 fp.Element
	t0.Set(xP)
	t01.Mul(&t0, xP2)
	t012.Mul(&t01, xQmP)
	t0123.Mul(&t012, xQmP2)
	t0123.Inverse(&t0123)

	invXQmP2.Mul(&t012, &t0123)
	t0123.Mul(&t0123, xQmP2)
	invXQmP.Mul(&t01, &t0123)
	t0123.Mul(&t0123, xQmP)
	invXP2.Mul(&t0, &t0123)
	t0123.Mul(&t0123, xP2)
	invXP.Set(&t0123)

	// Initialize T1 = Q + P, T2 = Q + P2
	var xT1, zT1, xT2, zT2 fp.Element
	var oneZ fp.Element
	oneZ.SetOne()
	xT1, zT1 = cADD(&tab.states[0].xR0, &tab.states[0].zR0, xP, &oneZ, &invXQmP)
	xT2, zT2 = cADD(&tab.states[0].xR0, &tab.states[0].zR0, xP2, &oneZ, &invXQmP2)

	// Temporaries
	var U, V, t1, t2, Vdiff fp.Element
	var xT1new, zT1new, xT2new, zT2new fp.Element

	// Process using precomputed ladder states
	for idx := 0; idx < 63; idx++ {
		state := &tab.states[idx]

		if state.bit == 0 {
			U.Mul(&xT1, &state.xR0)
			V.Mul(&zT1, &state.zR0)
			t1.Mul(&xT1, &state.zR0)
			t2.Mul(&zT1, &state.xR0)
			Vdiff.Sub(&t1, &t2)
			xT1new.Sub(&U, &V)
			xT1new.Square(&xT1new)
			xT1new.Mul(&xT1new, &invXP)
			zT1new.Square(&Vdiff)

			U.Mul(&xT2, &state.xR0)
			V.Mul(&zT2, &state.zR0)
			t1.Mul(&xT2, &state.zR0)
			t2.Mul(&zT2, &state.xR0)
			Vdiff.Sub(&t1, &t2)
			xT2new.Sub(&U, &V)
			xT2new.Square(&xT2new)
			xT2new.Mul(&xT2new, &invXP2)
			zT2new.Square(&Vdiff)
		} else {
			U.Mul(&xT1, &state.xR1)
			V.Mul(&zT1, &state.zR1)
			t1.Mul(&xT1, &state.zR1)
			t2.Mul(&zT1, &state.xR1)
			Vdiff.Sub(&t1, &t2)
			xT1new.Sub(&U, &V)
			xT1new.Square(&xT1new)
			xT1new.Mul(&xT1new, &invXQmP)
			zT1new.Square(&Vdiff)

			U.Mul(&xT2, &state.xR1)
			V.Mul(&zT2, &state.zR1)
			t1.Mul(&xT2, &state.zR1)
			t2.Mul(&zT2, &state.xR1)
			Vdiff.Sub(&t1, &t2)
			xT2new.Sub(&U, &V)
			xT2new.Square(&xT2new)
			xT2new.Mul(&xT2new, &invXQmP2)
			zT2new.Square(&Vdiff)
		}

		xT1, xT1new = xT1new, xT1
		zT1, zT1new = zT1new, zT1
		xT2, xT2new = xT2new, xT2
		zT2, zT2new = zT2new, zT2
	}

	if tab.finalNQ.Z.IsZero() || zT1.IsZero() || zT2.IsZero() {
		return false
	}

	var n1, n2, invNQZ fp.Element
	invNQZ.Inverse(&tab.finalNQ.Z)
	n1.Mul(&zT1, &invNQZ)
	n2.Mul(&zT2, &invNQZ)

	return f1IsOne(&n1) && f2IsOne(&n2)
}

// ComputeMontgomeryInputs computes the Montgomery x-coordinates needed for
// IsInSubGroupCubicalMontgomery from a Weierstrass point P.
// Returns xP, xP2, xQmP, xQmP2.
func ComputeMontgomeryInputs(p *curve.G1Affine) (xP, xP2, xQmP, xQmP2 fp.Element) {
	initMontgomeryConstants()

	var thirdRootOneG1 fp.Element
	thirdRootOneG1.SetString("80949648264912719408558363140637477264845294720710499478137287262712535938301461879813459410945")

	// Convert P to Montgomery
	xP = weierstrassToMontgomeryX(&p.X)

	// Compute φ(P) in Weierstrass then convert
	var p2X fp.Element
	p2X.Mul(&p.X, &thirdRootOneG1)
	xP2 = weierstrassToMontgomeryX(&p2X)

	// Get Q from precomputed table
	tab := precomputeTableDefault()

	// Compute Q - P and Q - φ(P) in Weierstrass
	var qMinusP, qMinusP2 curve.G1Affine
	var pNeg, p2Neg curve.G1Affine
	pNeg.X.Set(&p.X)
	pNeg.Y.Neg(&p.Y)
	p2Neg.X.Set(&p2X)
	p2Neg.Y.Neg(&p.Y)
	qMinusP.Add(&tab.q, &pNeg)
	qMinusP2.Add(&tab.q, &p2Neg)

	xQmP = weierstrassToMontgomeryX(&qMinusP.X)
	xQmP2 = weierstrassToMontgomeryX(&qMinusP2.X)

	return xP, xP2, xQmP, xQmP2
}

// membershipTestCubical performs the subgroup membership test using the Miller loop.
// This implementation uses two separate Miller loop evaluations (at P and φ(P)).
func membershipTestCubical(p, q *curve.G1Affine) bool {
	var thirdRootOneG1 fp.Element
	thirdRootOneG1.SetString("80949648264912719408558363140637477264845294720710499478137287262712535938301461879813459410945")

	// Compute φ(P) = (ω·x, y) using the GLV endomorphism
	var p2 curve.G1Affine
	p2.X.Mul(&p.X, &thirdRootOneG1)
	p2.Y.Set(&p.Y)

	// Compute two separate Miller loops
	n1, d1 := millerLoopSingle(q, p)
	n2, d2 := millerLoopSingle(q, &p2)

	// Check for zero values (degenerate cases)
	if n1.IsZero() || d1.IsZero() || n2.IsZero() || d2.IsZero() {
		return false
	}

	// Simultaneous inversion and final exponentiations
	n1.Mul(&n1, &d2)
	n2.Mul(&n2, &d1)
	d1.Mul(&d1, &d2)
	d1.Inverse(&d1)
	n1.Mul(&n1, &d1)
	n2.Mul(&n2, &d1)

	return f1IsOne(&n1) && f2IsOne(&n2)
}

// cDBL computes cubical doubling [2]P.
// cDBL is identical to xDBL (per the paper).
func cDBL(XP, ZP *fp.Element) (X2P, Z2P fp.Element) {
	return xDBL(XP, ZP)
}

// cADD computes cubical differential addition P+Q given (X_P:Z_P), (X_Q:Z_Q) and 1/(X_{P-Q}).
// This is Algorithm 2 from ePrint 2025/672.
// Note: Takes the INVERSE of x(P-Q) as input, not x(P-Q) directly.
func cADD(XP, ZP, XQ, ZQ, invXPmQ *fp.Element) (XPpQ, ZPpQ fp.Element) {
	var U, V, Usum, Vdiff fp.Element

	// U = X_P * X_Q
	U.Mul(XP, XQ)
	// V = Z_P * Z_Q
	V.Mul(ZP, ZQ)

	// Usum = X_P * Z_Q + Z_P * X_Q
	var t1, t2 fp.Element
	t1.Mul(XP, ZQ)
	t2.Mul(ZP, XQ)
	Usum.Add(&t1, &t2)

	// Vdiff = X_P * Z_Q - Z_P * X_Q
	Vdiff.Sub(&t1, &t2)

	// X_{P+Q} = (1/X_{P-Q}) * (U - V)^2
	XPpQ.Sub(&U, &V)
	XPpQ.Square(&XPpQ)
	XPpQ.Mul(&XPpQ, invXPmQ)

	// Z_{P+Q} = Vdiff^2
	ZPpQ.Square(&Vdiff)

	return XPpQ, ZPpQ
}

// cubicalLadderProper computes [n]P and [n]P + Q using the CUBICAL Montgomery ladder.
// This uses cDBL and cADD (cubical operations) throughout.
// Returns (nP, nPQ) in cubical projective coordinates.
func cubicalLadderProper(xP, xQ, xPmQ *fp.Element, n []byte, nbitlen int) (nP, nPQ PointXZ) {
	initMontgomeryConstants()

	// Precompute inverses: 1/xP, 1/xQ, 1/x(P-Q)
	var invXP, invXQ, invXPmQ fp.Element
	invXP.Inverse(xP)
	invXQ.Inverse(xQ)
	invXPmQ.Inverse(xPmQ)

	// Find the highest set bit
	highBit := nbitlen - 1
	for highBit >= 0 && ((n[highBit>>3]>>(highBit&7))&1) == 0 {
		highBit--
	}

	if highBit < 0 {
		// n is 0: return (O, Q)
		nP.X.SetOne()
		nP.Z.SetZero()
		nPQ.X.Set(xQ)
		nPQ.Z.SetOne()
		return nP, nPQ
	}

	// Initialize: R0 = P, R1 = 2P, T = P+Q
	var xR0, zR0, xR1, zR1, xT, zT fp.Element
	xR0.Set(xP)
	zR0.SetOne()
	xR1, zR1 = cDBL(&xR0, &zR0)

	// T = P + Q using cADD
	// For P + Q, the difference is P - (Q) = P - Q, so we need 1/x(P-Q) = invXPmQ
	// But actually T = P + Q, and the differential is (P + Q) - P = Q, no wait...
	// cADD(P, Q, diff) computes P + Q where diff = 1/x(P - Q)
	var oneZ fp.Element
	oneZ.SetOne()
	xT, zT = cADD(&xR0, &zR0, xQ, &oneZ, &invXPmQ)

	// Process bits from second-highest down to 0
	for i := highBit - 1; i >= 0; i-- {
		bit := (n[i>>3] >> (i & 7)) & 1

		if bit == 0 {
			// R1 = R0 + R1 (diff = P), R0 = 2*R0, T = T + R0 (diff = Q)
			var xR0new, zR0new, xR1new, zR1new, xTnew, zTnew fp.Element

			// R1 = R0 + R1, difference is P
			xR1new, zR1new = cADD(&xR0, &zR0, &xR1, &zR1, &invXP)

			// T = T + R0, difference is Q (since T - R0 = (kP + Q) - kP = Q)
			xTnew, zTnew = cADD(&xT, &zT, &xR0, &zR0, &invXQ)

			// R0 = 2*R0
			xR0new, zR0new = cDBL(&xR0, &zR0)

			xR0.Set(&xR0new)
			zR0.Set(&zR0new)
			xR1.Set(&xR1new)
			zR1.Set(&zR1new)
			xT.Set(&xTnew)
			zT.Set(&zTnew)
		} else {
			// R0 = R0 + R1 (diff = P), R1 = 2*R1, T = T + R1 (diff = Q - P but x is same as x(P-Q))
			var xR0new, zR0new, xR1new, zR1new, xTnew, zTnew fp.Element

			// R0 = R0 + R1, difference is P
			xR0new, zR0new = cADD(&xR0, &zR0, &xR1, &zR1, &invXP)

			// T = T + R1, difference is Q - P, and x(Q-P) = x(P-Q)
			xTnew, zTnew = cADD(&xT, &zT, &xR1, &zR1, &invXPmQ)

			// R1 = 2*R1
			xR1new, zR1new = cDBL(&xR1, &zR1)

			xR0.Set(&xR0new)
			zR0.Set(&zR0new)
			xR1.Set(&xR1new)
			zR1.Set(&zR1new)
			xT.Set(&xTnew)
			zT.Set(&zTnew)
		}
	}

	// R0 = [n]P, T = [n]P + Q
	nP.X.Set(&xR0)
	nP.Z.Set(&zR0)
	nPQ.X.Set(&xT)
	nPQ.Z.Set(&zT)

	return nP, nPQ
}

// cubicalTatePairing computes the cubical Tate pairing using the cubical Montgomery ladder.
// Uses the ODD degree formula: e_c(Q, P) = Z_{[ℓ]Q+P} / Z_{[ℓ]Q}
// where ℓ = e₂ - 1 (which is odd).
// From ePrint 2025/672: "Simpler and Faster Pairings from the Montgomery Ladder"
func cubicalTatePairing(q, p *curve.G1Affine) (num, den fp.Element) {
	initMontgomeryConstants()

	// Convert to Montgomery x-coordinates
	xQ := weierstrassToMontgomeryX(&q.X)
	xP := weierstrassToMontgomeryX(&p.X)

	// Compute x(Q - P) for the differential addition
	var qMinusP curve.G1Affine
	var pNeg curve.G1Affine
	pNeg.Neg(p)
	qMinusP.Add(q, &pNeg)
	xQmP := weierstrassToMontgomeryX(&qMinusP.X)

	// e₂ - 1 = z - 2 = 0x8508bfffffffffff (ODD)
	// In little-endian bytes
	ell := []byte{0xff, 0xff, 0xff, 0xff, 0xff, 0xbf, 0x08, 0x85}

	// Compute [ℓ]Q and [ℓ]Q + P using the CUBICAL ladder
	nQ, nQP := cubicalLadderProper(&xQ, &xP, &xQmP, ell, 64)

	// For odd degree: e_c(Q, P) = Z_{[ℓ]Q+P} / Z_{[ℓ]Q}
	num.Set(&nQP.Z)
	den.Set(&nQ.Z)

	return num, den
}

// membershipTestPureCubical performs the subgroup membership test using pure cubical pairing.
// Uses the cubical Tate pairing from ePrint 2025/672 with x-only Montgomery ladder arithmetic.
// OPTIMIZED: Single combined ladder computes [ℓ]Q, [ℓ]Q+P, and [ℓ]Q+φ(P) simultaneously,
// sharing all doubling operations between the two evaluation points.
func membershipTestPureCubical(p, q *curve.G1Affine) bool {
	initMontgomeryConstants()

	var thirdRootOneG1 fp.Element
	thirdRootOneG1.SetString("80949648264912719408558363140637477264845294720710499478137287262712535938301461879813459410945")

	// Convert Q to Montgomery x-coordinate
	xQ := weierstrassToMontgomeryX(&q.X)

	// Compute φ(P) = (ω·x, y) using the GLV endomorphism
	var p2X fp.Element
	p2X.Mul(&p.X, &thirdRootOneG1)

	// Convert P and φ(P) to Montgomery x-coordinates
	xP := weierstrassToMontgomeryX(&p.X)
	xP2 := weierstrassToMontgomeryX(&p2X)

	// Compute x(Q - P) and x(Q - φ(P)) for the differential additions
	var qMinusP, qMinusP2 curve.G1Affine
	var pNeg, p2Neg curve.G1Affine
	pNeg.X.Set(&p.X)
	pNeg.Y.Neg(&p.Y)
	p2Neg.X.Set(&p2X)
	p2Neg.Y.Neg(&p.Y)
	qMinusP.Add(q, &pNeg)
	qMinusP2.Add(q, &p2Neg)
	xQmP := weierstrassToMontgomeryX(&qMinusP.X)
	xQmP2 := weierstrassToMontgomeryX(&qMinusP2.X)

	// Batch inversion for 5 elements: 1/xQ, 1/xP, 1/xP2, 1/x(Q-P), 1/x(Q-P2)
	var invXQ, invXP, invXP2, invXQmP, invXQmP2 fp.Element
	var t0, t01, t012, t0123, t01234 fp.Element
	t0.Set(&xQ)
	t01.Mul(&t0, &xP)
	t012.Mul(&t01, &xP2)
	t0123.Mul(&t012, &xQmP)
	t01234.Mul(&t0123, &xQmP2)
	t01234.Inverse(&t01234)

	invXQmP2.Mul(&t0123, &t01234)
	t01234.Mul(&t01234, &xQmP2)
	invXQmP.Mul(&t012, &t01234)
	t01234.Mul(&t01234, &xQmP)
	invXP2.Mul(&t01, &t01234)
	t01234.Mul(&t01234, &xP2)
	invXP.Mul(&t0, &t01234)
	t01234.Mul(&t01234, &xP)
	invXQ.Set(&t01234)

	// Combined ladder: compute [ℓ]Q, [ℓ]Q+P, [ℓ]Q+φ(P) in single pass
	nQ, nQP, nQP2 := cubicalLadderCombined(&xQ, &xP, &xP2, &xQmP, &xQmP2,
		&invXQ, &invXP, &invXP2, &invXQmP, &invXQmP2)

	// Check for zero values
	if nQ.Z.IsZero() || nQP.Z.IsZero() || nQP2.Z.IsZero() {
		return false
	}

	// e_c(Q, P) = Z_{[ℓ]Q+P} / Z_{[ℓ]Q}
	// e_c(Q, φ(P)) = Z_{[ℓ]Q+φ(P)} / Z_{[ℓ]Q}
	var n1, n2 fp.Element
	n1.Set(&nQP.Z)
	n2.Set(&nQP2.Z)
	nQ.Z.Inverse(&nQ.Z)
	n1.Mul(&n1, &nQ.Z)
	n2.Mul(&n2, &nQ.Z)

	return f1IsOne(&n1) && f2IsOne(&n2)
}

// cubicalLadderCombined computes [n]Q, [n]Q+P, and [n]Q+P2 in a single pass.
// This shares all doubling operations between the two evaluation points P and P2.
// e₂ - 1 = z - 2 = 0x8508bfffffffffff is hardcoded for efficiency.
// All cADD and cDBL operations are inlined to avoid function call overhead.
func cubicalLadderCombined(xQ, xP, xP2, xQmP, xQmP2, invXQ, invXP, invXP2, invXQmP, invXQmP2 *fp.Element) (nQ, nQP, nQP2 PointXZ) {
	initMontgomeryConstants()

	// Scalar: e₂ - 1 = 0x8508bfffffffffff (64 bits)
	n := []byte{0xff, 0xff, 0xff, 0xff, 0xff, 0xbf, 0x08, 0x85}

	// Initialize: R0 = Q, R1 = 2Q
	var xR0, zR0, xR1, zR1 fp.Element
	xR0.Set(xQ)
	zR0.SetOne()

	// Inline cDBL for R1 = 2*R0
	{
		var sum, diff, sumSq, diffSq, t, temp fp.Element
		sum.Add(&xR0, &zR0)
		diff.Sub(&xR0, &zR0)
		sumSq.Square(&sum)
		diffSq.Square(&diff)
		xR1.Mul(&sumSq, &diffSq)
		t.Sub(&sumSq, &diffSq)
		temp.Mul(&montA24, &t)
		temp.Add(&temp, &diffSq)
		zR1.Mul(&t, &temp)
	}

	// Initialize T1 = Q + P, T2 = Q + P2 using inline cADD
	var xT1, zT1, xT2, zT2 fp.Element
	var oneZ fp.Element
	oneZ.SetOne()

	// Inline cADD for T1 = R0 + P (diff inverse = invXQmP)
	{
		var U, V, t1, t2, Vdiff fp.Element
		U.Mul(&xR0, xP)
		V.Mul(&zR0, &oneZ)
		t1.Mul(&xR0, &oneZ)
		t2.Mul(&zR0, xP)
		Vdiff.Sub(&t1, &t2)
		xT1.Sub(&U, &V)
		xT1.Square(&xT1)
		xT1.Mul(&xT1, invXQmP)
		zT1.Square(&Vdiff)
	}

	// Inline cADD for T2 = R0 + P2 (diff inverse = invXQmP2)
	{
		var U, V, t1, t2, Vdiff fp.Element
		U.Mul(&xR0, xP2)
		V.Mul(&zR0, &oneZ)
		t1.Mul(&xR0, &oneZ)
		t2.Mul(&zR0, xP2)
		Vdiff.Sub(&t1, &t2)
		xT2.Sub(&U, &V)
		xT2.Square(&xT2)
		xT2.Mul(&xT2, invXQmP2)
		zT2.Square(&Vdiff)
	}

	// Temporary variables
	var U, V, t1, t2, Vdiff fp.Element
	var sum, diff, sumSq, diffSq, t, temp fp.Element
	var xR0new, zR0new, xR1new, zR1new fp.Element
	var xT1new, zT1new, xT2new, zT2new fp.Element

	// Process bits 62 down to 0
	for i := 62; i >= 0; i-- {
		bit := (n[i>>3] >> (i & 7)) & 1

		if bit == 0 {
			// All three cADDs (T1+R0, T2+R0, and partially R1) use R0
			// cADD: R1new = R0 + R1 (diff = Q)
			U.Mul(&xR0, &xR1)
			V.Mul(&zR0, &zR1)
			t1.Mul(&xR0, &zR1)
			t2.Mul(&zR0, &xR1)
			Vdiff.Sub(&t1, &t2)
			xR1new.Sub(&U, &V)
			xR1new.Square(&xR1new)
			xR1new.Mul(&xR1new, invXQ)
			zR1new.Square(&Vdiff)

			// cADD: T1new = T1 + R0, T2new = T2 + R0
			// Reuse zR0 products where possible
			t1.Mul(&xT1, &zR0)
			t2.Mul(&zT1, &xR0)
			U.Mul(&xT1, &xR0)
			V.Mul(&zT1, &zR0)
			Vdiff.Sub(&t1, &t2)
			xT1new.Sub(&U, &V)
			xT1new.Square(&xT1new)
			xT1new.Mul(&xT1new, invXP)
			zT1new.Square(&Vdiff)

			t1.Mul(&xT2, &zR0)
			t2.Mul(&zT2, &xR0)
			U.Mul(&xT2, &xR0)
			V.Mul(&zT2, &zR0)
			Vdiff.Sub(&t1, &t2)
			xT2new.Sub(&U, &V)
			xT2new.Square(&xT2new)
			xT2new.Mul(&xT2new, invXP2)
			zT2new.Square(&Vdiff)

			// cDBL: R0new = 2*R0
			sum.Add(&xR0, &zR0)
			diff.Sub(&xR0, &zR0)
			sumSq.Square(&sum)
			diffSq.Square(&diff)
			xR0new.Mul(&sumSq, &diffSq)
			t.Sub(&sumSq, &diffSq)
			temp.Mul(&montA24, &t)
			temp.Add(&temp, &diffSq)
			zR0new.Mul(&t, &temp)

			xR0, xR0new = xR0new, xR0
			zR0, zR0new = zR0new, zR0
			xR1, xR1new = xR1new, xR1
			zR1, zR1new = zR1new, zR1
			xT1, xT1new = xT1new, xT1
			zT1, zT1new = zT1new, zT1
			xT2, xT2new = xT2new, xT2
			zT2, zT2new = zT2new, zT2
		} else {
			// All use R0+R1 or R1, share R1 products
			// cADD: R0new = R0 + R1 (diff = Q)
			U.Mul(&xR0, &xR1)
			V.Mul(&zR0, &zR1)
			t1.Mul(&xR0, &zR1)
			t2.Mul(&zR0, &xR1)
			Vdiff.Sub(&t1, &t2)
			xR0new.Sub(&U, &V)
			xR0new.Square(&xR0new)
			xR0new.Mul(&xR0new, invXQ)
			zR0new.Square(&Vdiff)

			// cADD: T1new = T1 + R1, T2new = T2 + R1
			t1.Mul(&xT1, &zR1)
			t2.Mul(&zT1, &xR1)
			U.Mul(&xT1, &xR1)
			V.Mul(&zT1, &zR1)
			Vdiff.Sub(&t1, &t2)
			xT1new.Sub(&U, &V)
			xT1new.Square(&xT1new)
			xT1new.Mul(&xT1new, invXQmP)
			zT1new.Square(&Vdiff)

			t1.Mul(&xT2, &zR1)
			t2.Mul(&zT2, &xR1)
			U.Mul(&xT2, &xR1)
			V.Mul(&zT2, &zR1)
			Vdiff.Sub(&t1, &t2)
			xT2new.Sub(&U, &V)
			xT2new.Square(&xT2new)
			xT2new.Mul(&xT2new, invXQmP2)
			zT2new.Square(&Vdiff)

			// cDBL: R1new = 2*R1
			sum.Add(&xR1, &zR1)
			diff.Sub(&xR1, &zR1)
			sumSq.Square(&sum)
			diffSq.Square(&diff)
			xR1new.Mul(&sumSq, &diffSq)
			t.Sub(&sumSq, &diffSq)
			temp.Mul(&montA24, &t)
			temp.Add(&temp, &diffSq)
			zR1new.Mul(&t, &temp)

			xR0, xR0new = xR0new, xR0
			zR0, zR0new = zR0new, zR0
			xR1, xR1new = xR1new, xR1
			zR1, zR1new = zR1new, zR1
			xT1, xT1new = xT1new, xT1
			zT1, zT1new = zT1new, zT1
			xT2, xT2new = xT2new, xT2
			zT2, zT2new = zT2new, zT2
		}
	}

	nQ.X.Set(&xR0)
	nQ.Z.Set(&zR0)
	nQP.X.Set(&xT1)
	nQP.Z.Set(&zT1)
	nQP2.X.Set(&xT2)
	nQP2.Z.Set(&zT2)

	return nQ, nQP, nQP2
}

// cubicalLadderWithInverses computes [n]P and [n]P + Q using precomputed inverses.
// (Kept for backward compatibility with tests)
func cubicalLadderWithInverses(xP, xQ, xPmQ, invXP, invXQ, invXPmQ *fp.Element, n []byte, nbitlen int) (nP, nPQ PointXZ) {
	initMontgomeryConstants()

	highBit := nbitlen - 1
	for highBit >= 0 && ((n[highBit>>3]>>(highBit&7))&1) == 0 {
		highBit--
	}

	if highBit < 0 {
		nP.X.SetOne()
		nP.Z.SetZero()
		nPQ.X.Set(xQ)
		nPQ.Z.SetOne()
		return nP, nPQ
	}

	var xR0, zR0, xR1, zR1, xT, zT fp.Element
	xR0.Set(xP)
	zR0.SetOne()
	xR1, zR1 = cDBL(&xR0, &zR0)

	var oneZ fp.Element
	oneZ.SetOne()
	xT, zT = cADD(&xR0, &zR0, xQ, &oneZ, invXPmQ)

	for i := highBit - 1; i >= 0; i-- {
		bit := (n[i>>3] >> (i & 7)) & 1

		if bit == 0 {
			var xR0new, zR0new, xR1new, zR1new, xTnew, zTnew fp.Element
			xR1new, zR1new = cADD(&xR0, &zR0, &xR1, &zR1, invXP)
			xTnew, zTnew = cADD(&xT, &zT, &xR0, &zR0, invXQ)
			xR0new, zR0new = cDBL(&xR0, &zR0)
			xR0.Set(&xR0new)
			zR0.Set(&zR0new)
			xR1.Set(&xR1new)
			zR1.Set(&zR1new)
			xT.Set(&xTnew)
			zT.Set(&zTnew)
		} else {
			var xR0new, zR0new, xR1new, zR1new, xTnew, zTnew fp.Element
			xR0new, zR0new = cADD(&xR0, &zR0, &xR1, &zR1, invXP)
			xTnew, zTnew = cADD(&xT, &zT, &xR1, &zR1, invXPmQ)
			xR1new, zR1new = cDBL(&xR1, &zR1)
			xR0.Set(&xR0new)
			zR0.Set(&zR0new)
			xR1.Set(&xR1new)
			zR1.Set(&zR1new)
			xT.Set(&xTnew)
			zT.Set(&zTnew)
		}
	}

	nP.X.Set(&xR0)
	nP.Z.Set(&zR0)
	nPQ.X.Set(&xT)
	nPQ.Z.Set(&zT)

	return nP, nPQ
}
