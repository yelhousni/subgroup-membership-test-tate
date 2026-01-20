package bls12381

import (
	curve "github.com/consensys/gnark-crypto/ecc/bls12-381"
	"github.com/consensys/gnark-crypto/ecc/bls12-381/fp"
)

// IsInSubGroupTate checks whether p is in the correct subgroup using the
// Tate-based test with precomputation.
//
// It follows "Revisiting subgroup membership testing on pairing-friendly
// curves via the Tate pairing" by Y.  Dai et al.
// https://eprint.iacr.org/2024/1790.pdf (Alg.4 and 5).
func IsInSubGroupTate(tab loopkupTable, p *curve.G1Affine) bool {
	if p.IsInfinity() {
		return false
	}
	if !p.IsOnCurve() {
		return false
	}
	return membershipTest(p, &tab.q, tab.tab)
}

// membershipTest performs Algorithm 5.
func membershipTest(p, q *curve.G1Affine, tab []fp.Element) bool {

	var thirdRootOneG1 fp.Element
	thirdRootOneG1.SetString("4002409555221667392624310435006688643935503118305586438271171395842971157480381377015405980053539358417135540939436")

	// Step 1: Algorithm 4.
	var p2 /* phi(p) */ curve.G1Affine
	p2.X.Mul(&p.X, &thirdRootOneG1)
	p2.Y.Set(&p.Y)

	var n1, d1, n2, d2 fp.Element
	n1.Sub(&p.X, &q.X)
	n2.Sub(&p2.X, &q.X)
	d1.SetOne()
	d2.SetOne()
	// shared Miller loop with precomputation
	sharedMillerloop(tab, &n1, &d1, &n2, &d2, q, p, &p2)

	if n1.IsZero() || d1.IsZero() || n2.IsZero() || d2.IsZero() {
		return false
	}

	// steps 5-8: simultaneous inversion and final exponentiations.
	n1.Mul(&n1, &d2)
	n2.Mul(&n2, &d1)
	d1.Mul(&d1, &d2)
	d1.Inverse(&d1)
	n1.Mul(&n1, &d1)
	n2.Mul(&n2, &d1)

	// early abort since f2IsOne is cheaper
	if !f2IsOne(&n2) {
		return false
	}

	return f1IsOne(&n1)
}

func sharedMillerloop(tab []fp.Element, n1, d1, n2, d2 *fp.Element, q, p, p2 *curve.G1Affine) {
	i := 63
	j := 0
	k := 0
	if naf[0] < 0 {
		j = 1
	}

	var f1, g1, f2, g2 fp.Element
	var u0, u1, u2, u3 fp.Element
	var v0, v1, v2 fp.Element

	u1.Sub(&p.Y, &q.Y)
	u2.Sub(&p.X, &q.X)
	u3.Sub(&p2.X, &q.X)

	for i >= j {
		if naf[i] == 0 && i > j {
			u0.Sub(&p.Y, &tab[k+2])
			f1.Sub(&p.X, &tab[k+1])
			f2.Sub(&p2.X, &tab[k+1])
			g1.Mul(&f1, &tab[k])
			g1.Sub(&u0, &g1)
			g2.Mul(&f2, &tab[k])
			g2.Sub(&u0, &g2)
			f1.Mul(&f1, &tab[k+3])
			f1.Sub(&u0, &f1)
			f2.Mul(&f2, &tab[k+3])
			f2.Sub(&u0, &f2)

			n1.Square(n1)
			n1.Square(n1)
			n1.Mul(n1, &f1)
			d1.Square(d1)
			d1.Mul(d1, &g1)
			d1.Square(d1)

			n2.Square(n2)
			n2.Square(n2)
			n2.Mul(n2, &f2)
			d2.Square(d2)
			d2.Mul(d2, &g2)
			d2.Square(d2)

			k += 4
			i--

			if naf[i] > 0 {
				f1.Mul(&tab[k], &u2)
				f1.Sub(&u1, &f1)
				f2.Mul(&tab[k], &u3)
				f2.Sub(&u1, &f2)
				g1.Sub(&p.X, &tab[k+1])
				g2.Sub(&p2.X, &tab[k+1])

				n1.Mul(n1, &f1)
				d1.Mul(d1, &g1)
				n2.Mul(n2, &f2)
				d2.Mul(d2, &g2)
				k += 2
			}

			if naf[i] < 0 {
				f1.Sub(&p.X, &tab[k+1])
				f2.Sub(&p2.X, &tab[k+1])
				g1.Mul(&tab[k], &u2)
				g1.Sub(&u1, &g1)
				g2.Mul(&tab[k], &u3)
				g2.Sub(&u1, &g2)

				n1.Mul(n1, &f1)
				d1.Mul(d1, &g1)
				n2.Mul(n2, &f2)
				d2.Mul(d2, &g2)
				k += 2
			}
			i--
			continue
		}

		if naf[i] == 1 || naf[i] == -1 {
			u0.Sub(&p.X, &tab[k])
			g1.Sub(&p2.X, &tab[k])
			g2.Add(&p.X, &tab[k+2])
			f1.Add(&p2.X, &tab[k+2])
			v0.Mul(&u0, &g2)
			v1.Mul(&g1, &f1)
			g2.Sub(&p.Y, &tab[k+1])
			v2.Mul(&g2, &tab[k+3])
			v0.Sub(&v0, &v2)
			v1.Sub(&v1, &v2)
			f1.Set(&v0)
			f2.Set(&v1)

			n1.Square(n1)
			n1.Mul(n1, &f1)
			d1.Mul(d1, &u0)
			d1.Square(d1)

			n2.Square(n2)
			n2.Mul(n2, &f2)
			d2.Mul(d2, &g1)
			d2.Square(d2)
			if naf[i] < 0 {
				d1.Mul(d1, &u2)
				d2.Mul(d2, &u3)
			}

			k += 4
			i--
			continue
		}

		f1.Sub(&p.X, &tab[k+1])
		f2.Sub(&p2.X, &tab[k+1])
		u0.Sub(&p.Y, &tab[k+2])
		g1.Mul(&tab[k], &f1)
		g1.Sub(&u0, &g1)
		g2.Mul(&tab[k], &f2)
		g2.Sub(&u0, &g2)

		n1.Square(n1)
		n1.Mul(n1, &f1)
		d1.Square(d1)
		d1.Mul(d1, &g1)
		n2.Square(n2)
		n2.Mul(n2, &f2)
		d2.Square(d2)
		d2.Mul(d2, &g2)

		k += 3
		i--
	}

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
}

// f2IsOne raises x to exp2 = |z^5-z^4-z^3+z^2+z+2| and checks if the result is 1
func f2IsOne(x *fp.Element) bool {
	var u0, u1, u2, u3 fp.Element

	// Section 4.3: exp2 final exponentiation for the Tate test.
	u0.Square(x)
	expBySeed(&u1, x)
	expBySeed(&u2, &u1)
	expBySeed(&u3, &u2)
	u0.Mul(&u0, &u2)
	u0.Mul(&u0, &u3)
	expBySeed(&u3, &u3)
	u1.Mul(&u1, &u3)
	expBySeed(&u3, &u3)
	u1.Mul(&u1, &u3)

	return u0.Equal(&u1)
}

// expBySeed computes z = x^z where z = -0x396c8c005555e1568c00aaab0000aaab
func expBySeed(z, x *fp.Element) {
	var u0 fp.Element
	u0.Square(x)
	u0.Mul(&u0, x)
	u0.Square(&u0)
	u0.Square(&u0)
	u0.Mul(&u0, x)
	u0.Square(&u0)
	u0.Square(&u0)
	u0.Square(&u0)
	u0.Mul(&u0, x)
	for i := 0; i < 9; i++ {
		u0.Square(&u0)
	}
	u0.Mul(&u0, x)
	for i := 0; i < 32; i++ {
		u0.Square(&u0)
	}
	u0.Mul(&u0, x)
	for i := 0; i < 16; i++ {
		u0.Square(&u0)
	}

	z.Set(&u0)
}

// f1IsOne raises x to exp1 = (p-1)/e2 and checks if the result is 1
//
// Optimized implementation using the closed form:
// exp1 = (|z|⁵ + |z|⁴ - |z|³ - |z|² + |z| - 2) / 3
//
// Since we only check if x^exp1 = 1, we can avoid the division by 3:
// x^exp1 = 1 ⟺ x^(3·exp1) = 1 ⟺ x^(|z|⁵ + |z|⁴ - |z|³ - |z|² + |z| - 2) = 1
func f1IsOne(x *fp.Element) bool {
	var z, z2, z3, z4, z5 fp.Element
	var result, denom, xSquared fp.Element

	// Compute powers: x^|z|, x^|z²|, x^|z³|, x^|z⁴|, x^|z⁵|
	// These are computed via repeated expBySeed which handles the seed parameter
	expBySeed(&z, x)     // x^|z|
	expBySeed(&z2, &z)   // x^|z²|
	expBySeed(&z3, &z2)  // x^|z³|
	expBySeed(&z4, &z3)  // x^|z⁴|
	expBySeed(&z5, &z4)  // x^|z⁵|

	// Compute x^(|z|⁵ + |z|⁴ - |z|³ - |z|² + |z| - 2)
	// = (x^|z⁵| · x^|z⁴| · x^|z|) / (x^|z³| · x^|z²| · x²)

	// Numerator: x^|z⁵| · x^|z⁴| · x^|z|
	result.Mul(&z5, &z4)
	result.Mul(&result, &z)

	// Denominator: x^|z³| · x^|z²| · x²
	denom.Mul(&z3, &z2)
	xSquared.Square(x)
	denom.Mul(&denom, &xSquared)

	// Divide: multiply by inverse
	denom.Inverse(&denom)
	result.Mul(&result, &denom)

	// Check if result = 1
	return result.IsOne()
}
