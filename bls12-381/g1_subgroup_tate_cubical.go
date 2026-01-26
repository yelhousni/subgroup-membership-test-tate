package bls12381

import (
	curve "github.com/consensys/gnark-crypto/ecc/bls12-381"
	"github.com/consensys/gnark-crypto/ecc/bls12-381/fp"
)

// This file implements a Tate pairing for the G1 subgroup membership test using
// two separate Miller loops instead of the shared Miller loop approach.
//
// The key difference from the original g1_subgroup_tate.go is that this version
// evaluates the Miller function at P and φ(P) using two independent loop executions,
// rather than interleaving the computations in a single shared loop.
//
// Performance comparison (Apple M1):
//   - GLV (standard):           ~40 µs
//   - Tate shared loop:         ~46 µs
//   - Tate two loops:           ~48 µs  (this implementation)
//
// The two-loops approach has ~2% overhead compared to the shared loop, but provides
// a cleaner separation of concerns and serves as a foundation for further optimizations.

// IsInSubGroupTateTwoLoops checks subgroup membership using two separate
// precomputed Miller loop evaluations (no shared state between P and φ(P)).
func IsInSubGroupTateTwoLoops(tab loopkupTable, p *curve.G1Affine) bool {
	if p.IsInfinity() {
		return false
	}
	if !p.IsOnCurve() {
		return false
	}
	return membershipTestTwoLoops(p, &tab.q, tab.tab)
}

// membershipTestTwoLoops performs subgroup membership test using
// two separate evaluations of the precomputed Miller loop.
func membershipTestTwoLoops(p, q *curve.G1Affine, tab []fp.Element) bool {
	var thirdRootOneG1 fp.Element
	thirdRootOneG1.SetString("4002409555221667392624310435006688643935503118305586438271171395842971157480381377015405980053539358417135540939436")

	// Compute φ(P) = (ω·x, y) using the GLV endomorphism
	var p2 curve.G1Affine
	p2.X.Mul(&p.X, &thirdRootOneG1)
	p2.Y.Set(&p.Y)

	// Run two separate Miller loops using precomputation
	n1, d1 := millerLoopPrecomputed(tab, q, p)
	n2, d2 := millerLoopPrecomputed(tab, q, &p2)

	// Check for zero values (degenerate cases)
	if n1.IsZero() || d1.IsZero() || n2.IsZero() || d2.IsZero() {
		return false
	}

	// Simultaneous inversion: compute 1/(d1*d2) with one inversion
	// f1 = n1/d1, f2 = n2/d2 → n1*d2/(d1*d2), n2*d1/(d1*d2)
	n1.Mul(&n1, &d2)
	n2.Mul(&n2, &d1)
	d1.Mul(&d1, &d2)
	d1.Inverse(&d1)
	n1.Mul(&n1, &d1)
	n2.Mul(&n2, &d1)

	// Final exponentiations to check membership
	return f1IsOne(&n1) && f2IsOne(&n2)
}

// millerLoopPrecomputed evaluates the Miller function f_{e₂-1,Q}(P) at a single
// point using precomputed values. This uses the same precomputation table as
// the shared loop but for a single evaluation point.
func millerLoopPrecomputed(tab []fp.Element, q, p *curve.G1Affine) (num, den fp.Element) {
	i := 63
	j := 0
	k := 0
	if naf[0] < 0 {
		j = 1
	}

	// Initialize accumulators
	num.Sub(&p.X, &q.X)
	den.SetOne()

	var f, g fp.Element
	var u0, u1, u2 fp.Element

	// Pre-compute constants used in line evaluations
	u1.Sub(&p.Y, &q.Y) // y_P - y_Q
	u2.Sub(&p.X, &q.X) // x_P - x_Q

	for i >= j {
		if naf[i] == 0 && i > j {
			// Double-double step: process two consecutive 0-bits
			// tab[k]:   λ₁ = -3x_T²/(2y_T) (negated tangent at T)
			// tab[k+1]: x_{2T}
			// tab[k+2]: y_{2T}
			// tab[k+3]: λ₂ = 3x_{2T}²/(2y_{2T}) (tangent at 2T)

			u0.Sub(&p.Y, &tab[k+2]) // y_P - y_{2T}
			f.Sub(&p.X, &tab[k+1])  // x_P - x_{2T}

			// g = (y_P - y_{2T}) - λ₁ * (x_P - x_{2T})
			g.Mul(&f, &tab[k])
			g.Sub(&u0, &g)

			// f = (y_P - y_{2T}) - λ₂ * (x_P - x_{2T})
			f.Mul(&f, &tab[k+3])
			f.Sub(&u0, &f)

			// Accumulate: num⁴ * f, den² * g then den²
			num.Square(&num)
			num.Square(&num)
			num.Mul(&num, &f)
			den.Square(&den)
			den.Mul(&den, &g)
			den.Square(&den)

			k += 4
			i--

			// Addition step if naf[i] > 0
			if naf[i] > 0 {
				// tab[k]:   λ = (y_T - y_Q)/(x_T - x_Q)
				// tab[k+1]: x_T (BEFORE addition)
				f.Mul(&tab[k], &u2)
				f.Sub(&u1, &f)
				g.Sub(&p.X, &tab[k+1])

				num.Mul(&num, &f)
				den.Mul(&den, &g)
				k += 2
			}

			// Addition step if naf[i] < 0
			if naf[i] < 0 {
				// tab[k]:   λ = (y_T + y_Q)/(x_Q - x_T)
				// tab[k+1]: x_{T-Q} (AFTER subtraction)
				f.Sub(&p.X, &tab[k+1])
				g.Mul(&tab[k], &u2)
				g.Sub(&u1, &g)

				num.Mul(&num, &f)
				den.Mul(&den, &g)
				k += 2
			}
			i--
			continue
		}

		if naf[i] == 1 || naf[i] == -1 {
			// Add-then-double step using combined line formula
			// tab[k]:   x_T
			// tab[k+1]: y_T
			// tab[k+2]: λ₁*λ₂ + x_T + x_{T'}
			// tab[k+3]: λ₁ + λ₂

			u0.Sub(&p.X, &tab[k])    // x_P - x_T
			g.Add(&p.X, &tab[k+2])   // x_P + (λ₁*λ₂ + x_T + x_{T'})
			f.Mul(&u0, &g)           // (x_P - x_T)(x_P + combined1)
			g.Sub(&p.Y, &tab[k+1])   // y_P - y_T
			u0.Mul(&g, &tab[k+3])    // (y_P - y_T) * (λ₁ + λ₂)
			f.Sub(&f, &u0)           // Combined line evaluation

			// Vertical line factor
			u0.Sub(&p.X, &tab[k])
			num.Square(&num)
			num.Mul(&num, &f)
			den.Mul(&den, &u0)
			den.Square(&den)

			// For negative NAF, multiply by additional vertical line
			if naf[i] < 0 {
				den.Mul(&den, &u2)
			}

			k += 4
			i--
			continue
		}

		// Simple double step (fallback for isolated naf[i] == 0 at boundary)
		// tab[k]:   λ = -3x_T²/(2y_T)
		// tab[k+1]: x_{2T}
		// tab[k+2]: y_{2T}
		f.Sub(&p.X, &tab[k+1])  // x_P - x_{2T}
		u0.Sub(&p.Y, &tab[k+2]) // y_P - y_{2T}
		g.Mul(&tab[k], &f)      // λ * (x_P - x_{2T})
		g.Sub(&u0, &g)          // Line evaluation

		num.Square(&num)
		num.Mul(&num, &f)
		den.Square(&den)
		den.Mul(&den, &g)

		k += 3
		i--
	}

	// Final step if naf[0] < 0: multiply by vertical line
	if naf[0] < 0 {
		g.Sub(&p.X, &tab[k])
		num.Square(&num)
		den.Square(&den)
		den.Mul(&den, &g)
	}

	return num, den
}

// ============================================================================
// Reference implementations below (for testing and comparison)
// ============================================================================

// torsionPointCubical is the fixed torsion point Q of order e₂ = |z|+1.
var torsionPointCubical curve.G1Affine

func init() {
	torsionPointCubical.X.SetString("0xD82B23C3EE86C6B55930A7755FEB499A697AAE08D97E677F61EBF6894E57EC7434DA198FE1FBF0EF1C7004640A74203")
	torsionPointCubical.Y.SetString("0x75868854578CF684F73F747280EF3F0A86CD94B3FB5954BC8B6FA4888BE7B2FB766E6DAF6F4F0AB9FE3E757B4BE8404")
}

// IsInSubGroupTateCubical is the reference implementation using on-the-fly
// Miller loop computation (without precomputation). This is slower but useful
// for testing correctness.
func IsInSubGroupTateCubical(p *curve.G1Affine) bool {
	if p.IsInfinity() {
		return false
	}
	if !p.IsOnCurve() {
		return false
	}
	return membershipTestCubical(p, &torsionPointCubical)
}

// membershipTestCubical performs the subgroup membership test using two
// separate Miller loop computations without precomputation.
func membershipTestCubical(p, q *curve.G1Affine) bool {
	var thirdRootOneG1 fp.Element
	thirdRootOneG1.SetString("4002409555221667392624310435006688643935503118305586438271171395842971157480381377015405980053539358417135540939436")

	var p2 curve.G1Affine
	p2.X.Mul(&p.X, &thirdRootOneG1)
	p2.Y.Set(&p.Y)

	n1, d1 := millerLoopSingle(q, p)
	n2, d2 := millerLoopSingle(q, &p2)

	if n1.IsZero() || d1.IsZero() || n2.IsZero() || d2.IsZero() {
		return false
	}

	n1.Mul(&n1, &d2)
	n2.Mul(&n2, &d1)
	d1.Mul(&d1, &d2)
	d1.Inverse(&d1)
	n1.Mul(&n1, &d1)
	n2.Mul(&n2, &d1)

	return f1IsOne(&n1) && f2IsOne(&n2)
}

// millerLoopSingle computes the Miller function f_{e₂-1,Q}(P) for a single point P.
// This is a reference implementation with explicit inversions (slow but correct).
func millerLoopSingle(q, p *curve.G1Affine) (num, den fp.Element) {
	i := 63
	j := 0
	if naf[0] < 0 {
		j = 1
	}

	num.Sub(&p.X, &q.X)
	den.SetOne()

	var T curve.G1Affine
	T.Set(q)

	var qNeg curve.G1Affine
	qNeg.Neg(q)

	var f, g fp.Element
	var u0, u1 fp.Element

	var yPmyQ, xPmxQ fp.Element
	yPmyQ.Sub(&p.Y, &q.Y)
	xPmxQ.Sub(&p.X, &q.X)

	for i >= j {
		if naf[i] == 0 && i > j {
			// Double-double step
			u0.Square(&T.X)
			u1.Double(&u0)
			u0.Add(&u0, &u1)
			u1.Double(&T.Y)
			u1.Inverse(&u1)
			u1.Neg(&u1)
			var lambda1 fp.Element
			lambda1.Mul(&u0, &u1)

			T.Double(&T)

			u0.Square(&T.X)
			u1.Double(&u0)
			u0.Add(&u0, &u1)
			u1.Double(&T.Y)
			u1.Inverse(&u1)
			var lambda2 fp.Element
			lambda2.Mul(&u0, &u1)

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

				f.Mul(&lambda1, &xPmxQ)
				f.Sub(&yPmyQ, &f)
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
				g.Mul(&lambda1, &xPmxQ)
				g.Sub(&yPmyQ, &g)

				num.Mul(&num, &f)
				den.Mul(&den, &g)
			}

			i--
			continue
		}

		if naf[i] == 1 || naf[i] == -1 {
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
				g.Mul(&g, &xPmxQ)
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

		f.Sub(&p.X, &T.X)
		g.Mul(&f, &lambda)
		u0.Sub(&p.Y, &T.Y)
		f.Sub(&u0, &g)

		T.Double(&T)
		g.Sub(&p.X, &T.X)

		num.Square(&num)
		num.Mul(&num, &f)
		den.Square(&den)
		den.Mul(&den, &g)

		i--
	}

	if naf[0] < 0 {
		g.Sub(&p.X, &T.X)
		num.Square(&num)
		den.Square(&den)
		den.Mul(&den, &g)
	}

	return num, den
}
