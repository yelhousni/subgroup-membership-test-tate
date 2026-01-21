package bls12377

import (
	"sync"

	curve "github.com/consensys/gnark-crypto/ecc/bls12-377"
	"github.com/consensys/gnark-crypto/ecc/bls12-377/fp"
)

type loopkupTable struct {
	q   curve.G1Affine
	tab []fp.Element
}

var (
	precomputeTableOnce sync.Once
	precomputedTable    loopkupTable
)

// NAF digits used by Algorithms 3/4 (Section 4.2) for BLS12-377.
// e2-1 = z-2 = 0x8508bfffffffffff (z > 0)
// NAF(z-2) with padding to match BLS12-381 structure (65 elements, naf[64]=0)
var naf = [65]int8{
	-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 0,
	0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0,
}

// precomputeTableDefault returns a cached precomputation table for torsionPoint.
func precomputeTableDefault() loopkupTable {
	precomputeTableOnce.Do(func() {
		// Point with exact order e2 = |z-1| = 0x8508c00000000000
		var torsionPoint curve.G1Affine
		torsionPoint.X.SetString("50650456740282261254444037957372148278657828508220767125535350349422402980993056080580432127086515431542182301593")
		torsionPoint.Y.SetString("211780388751464357047145090878738766994067119315669641699537871325723620439548899400758824542777430317160612183399")
		precomputedTable = generateTable(&torsionPoint)
	})
	return precomputedTable
}

// generateTable precomputes the lookup table used by the Tate-based G1 membership test.
// Algorithm 3 (Section 4.2): lookup table generation for the shared Miller loop.
// For BLS12-377: naf[63]=1 is handled by initialization, so loop starts at i=62.
func generateTable(q *curve.G1Affine) loopkupTable {
	i := 62
	j := 0
	if naf[0] < 0 {
		j = 1
	}

	tab := make([]fp.Element, 0, len(naf)*4)
	var t0, t1, qNeg curve.G1Affine
	t0.Set(q)
	qNeg.Neg(q)

	var u0, u1 fp.Element

	for i >= j {
		if naf[i] == 0 && i > j {
			u0.Square(&t0.X)
			u1.Double(&u0)
			u0.Add(&u0, &u1)
			u1.Double(&t0.Y)
			u1.Inverse(&u1)
			u1.Neg(&u1)
			u0.Mul(&u0, &u1)
			tab = append(tab, u0)

			t0.Double(&t0)
			u0.Square(&t0.X)
			u1.Double(&u0)
			u0.Add(&u0, &u1)
			u1.Double(&t0.Y)
			u1.Inverse(&u1)
			tab = append(tab, t0.X, t0.Y)
			u0.Mul(&u0, &u1)
			tab = append(tab, u0)

			t0.Double(&t0)
			i--

			if naf[i] > 0 {
				u0.Sub(&t0.Y, &q.Y)
				u1.Sub(&t0.X, &q.X)
				u1.Inverse(&u1)
				u0.Mul(&u0, &u1)
				tab = append(tab, u0, t0.X)
				t0.Add(&t0, q)
			}
			if naf[i] < 0 {
				u0.Add(&t0.Y, &q.Y)
				u1.Sub(&q.X, &t0.X)
				u1.Inverse(&u1)
				u0.Mul(&u0, &u1)
				t0.Add(&t0, &qNeg)
				tab = append(tab, u0, t0.X)
			}
			i--
			continue
		}

		if naf[i] == 1 {
			tab = append(tab, t0.X, t0.Y)

			var lambda1, lambda2 fp.Element
			lambda1.Sub(&t0.Y, &q.Y)
			lambda2.Sub(&t0.X, &q.X)
			lambda2.Inverse(&lambda2)
			lambda1.Mul(&lambda1, &lambda2)

			t1.Add(&t0, q)
			lambda2.Sub(&t1.Y, &t0.Y)
			u0.Sub(&t1.X, &t0.X)
			u0.Inverse(&u0)
			lambda2.Mul(&lambda2, &u0)

			u0.Mul(&lambda1, &lambda2)
			lambda2.Add(&lambda1, &lambda2)
			u0.Add(&u0, &t0.X)
			u0.Add(&u0, &t1.X)
			tab = append(tab, u0, lambda2)

			t0.Add(&t1, &t0)
			i--
			continue
		}

		if naf[i] == -1 {
			tab = append(tab, t0.X, t0.Y)

			var lambda1, lambda2 fp.Element
			lambda1.Add(&t0.Y, &q.Y)
			lambda2.Sub(&t0.X, &q.X)
			lambda2.Inverse(&lambda2)
			lambda1.Mul(&lambda1, &lambda2)

			t1.Sub(&t0, q)
			lambda2.Sub(&t1.Y, &t0.Y)
			u0.Sub(&t1.X, &t0.X)
			u0.Inverse(&u0)
			lambda2.Mul(&lambda2, &u0)

			u0.Mul(&lambda1, &lambda2)
			lambda2.Add(&lambda1, &lambda2)
			u0.Add(&u0, &t0.X)
			u0.Add(&u0, &t1.X)
			tab = append(tab, u0, lambda2)

			t0.Add(&t1, &t0)
			i--
			continue
		}

		u0.Square(&t0.X)
		u1.Double(&u0)
		u0.Add(&u0, &u1)
		u1.Double(&t0.Y)
		u1.Inverse(&u1)
		u1.Neg(&u1)
		u0.Mul(&u0, &u1)
		tab = append(tab, u0)

		t0.Double(&t0)
		tab = append(tab, t0.X, t0.Y)
		i--
	}

	if naf[0] < 0 {
		// Special case for BLS12-377: at this point t0 = [e2/2]Q which is a 2-torsion
		// point (has Y = 0). The last default doubling at i=1 stored the coordinates
		// of t0 = [e2/2]Q = (-1, 0) as tab[k-2], tab[k-1].
		//
		// For a 2-torsion point (Y = 0):
		// - The tangent is vertical: P.X - T.X
		// - Doubling gives O (infinity)
		//
		// For the subtraction O - Q = -Q:
		// - Line through O and -Q is vertical at Q: P.X - Q.X
		// - Vertical at -Q: P.X - Q.X (same)
		// - These cancel in n/d, contributing nothing
		//
		// We don't need extra table entries. The Miller loop will use tab[k-2] = T.X = -1
		// for the vertical tangent at the 2-torsion point.
	}

	return loopkupTable{
		q:   *q,
		tab: tab,
	}
}
