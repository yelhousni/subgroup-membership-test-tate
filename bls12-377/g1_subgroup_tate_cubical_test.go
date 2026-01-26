package bls12377

import (
	"math/big"
	"testing"

	curve "github.com/consensys/gnark-crypto/ecc/bls12-377"
	"github.com/consensys/gnark-crypto/ecc/bls12-377/fp"

	"github.com/leanovate/gopter"
	"github.com/leanovate/gopter/prop"
)

func TestMontgomeryTransformation(t *testing.T) {
	initMontgomeryConstants()

	// Test the 2-torsion point: Weierstrass (-1, 0) -> Montgomery (0, 0)
	var negOne fp.Element
	negOne.SetInt64(-1)
	xM := weierstrassToMontgomeryX(&negOne)
	if !xM.IsZero() {
		t.Errorf("Weierstrass (-1, 0) should map to Montgomery (0, 0), got x = %s", xM.String())
	}

	// Test round-trip: W -> M -> W
	_, _, g1, _ := curve.Generators()
	xM = weierstrassToMontgomeryX(&g1.X)
	xW := montgomeryToWeierstrassX(&xM)
	if !xW.Equal(&g1.X) {
		t.Errorf("Round-trip W -> M -> W failed: got %s, expected %s", xW.String(), g1.X.String())
	}
}

func TestCubicalXDBLADD(t *testing.T) {
	initMontgomeryConstants()

	// Test with generator point
	_, _, g1, _ := curve.Generators()

	// Convert to Montgomery
	xP := weierstrassToMontgomeryX(&g1.X)

	// Start with P
	var xS, zS fp.Element
	xS.Set(&xP)
	zS.SetOne()

	// Compute a simple differential point for testing
	var g1Double curve.G1Affine
	g1Double.Double(&g1)
	xQ := weierstrassToMontgomeryX(&g1Double.X)
	var xT, zT fp.Element
	xT.Set(&xQ)
	zT.SetOne()

	// Get 1/xP
	var ixP fp.Element
	ixP.Inverse(&xP)

	// Do one doubling step
	cubicalXDBLADD(&xS, &zS, &xT, &zT, &ixP)

	// xS/zS should be x([2]P) in Montgomery
	if zS.IsZero() {
		t.Error("doubling resulted in infinity")
	}

	var x2P fp.Element
	var invZS fp.Element
	invZS.Inverse(&zS)
	x2P.Mul(&xS, &invZS)

	// Convert back to Weierstrass
	x2PW := montgomeryToWeierstrassX(&x2P)

	// Compare with expected [2]G1
	if !x2PW.Equal(&g1Double.X) {
		t.Errorf("cubicalXDBLADD doubling failed: got %s, expected %s", x2PW.String(), g1Double.X.String())
	}
}

func TestG1SubGroupMembershipTateCubical(t *testing.T) {
	t.Parallel()
	parameters := gopter.DefaultTestParameters()
	if testing.Short() {
		parameters.MinSuccessfulTests = 10
	} else {
		parameters.MinSuccessfulTests = 100
	}

	properties := gopter.NewProperties(parameters)

	properties.Property("[BLS12-377] IsInSubGroupTateCubical should output true for points on G1", prop.ForAll(
		func(a fp.Element) bool {
			p := curve.MapToG1(a)
			return IsInSubGroupTateCubical(&p)
		},
		GenFp(),
	))

	properties.TestingRun(t, gopter.ConsoleReporter(false))
}

func TestPureCubicalPairing(t *testing.T) {
	initMontgomeryConstants()
	tab := precomputeTableDefault()
	q := &tab.q

	// Test with generator point
	_, _, p, _ := curve.Generators()

	// First, verify the ladder computes correct scalar multiplication
	// Compute [e₂-1]Q using standard scalar mult
	var e2m1 curve.G1Affine
	// e₂ - 1 = z - 2 = 0x8508bfffffffffff
	var scalar [8]byte
	scalar[0] = 0xff
	scalar[1] = 0xff
	scalar[2] = 0xff
	scalar[3] = 0xff
	scalar[4] = 0xff
	scalar[5] = 0xbf
	scalar[6] = 0x08
	scalar[7] = 0x85
	// gnark uses big-endian for scalars
	var scalarBE [32]byte
	for i := 0; i < 8; i++ {
		scalarBE[31-i] = scalar[i]
	}
	e2m1.ScalarMultiplication(q, new(big.Int).SetBytes(scalarBE[:]))
	t.Logf("[e₂-1]Q (standard) = (%s, %s)", e2m1.X.String(), e2m1.Y.String())

	// Now compute via ladder
	xQ := weierstrassToMontgomeryX(&q.X)
	xP := weierstrassToMontgomeryX(&p.X)
	var qMinusP curve.G1Affine
	var pNeg curve.G1Affine
	pNeg.Neg(&p)
	qMinusP.Add(q, &pNeg)
	xQmP := weierstrassToMontgomeryX(&qMinusP.X)

	ell := []byte{0xff, 0xff, 0xff, 0xff, 0xff, 0xbf, 0x08, 0x85}
	nQ, _ := cubicalLadder(&xQ, &xP, &xQmP, ell, 64, false)

	// Convert ladder result back to Weierstrass x-coordinate
	var nQxAff fp.Element
	var invZ fp.Element
	invZ.Inverse(&nQ.Z)
	nQxAff.Mul(&nQ.X, &invZ)
	nQxW := montgomeryToWeierstrassX(&nQxAff)
	t.Logf("[e₂-1]Q (ladder) x = %s", nQxW.String())
	t.Logf("Match: %v", nQxW.Equal(&e2m1.X))

	// Now test the pairing
	n, d := cubicalTatePairing(q, &p)
	t.Logf("cubicalTatePairing num = %s", n.String())
	t.Logf("cubicalTatePairing den = %s", d.String())

	// Compute ratio
	var ratio fp.Element
	d.Inverse(&d)
	ratio.Mul(&n, &d)
	t.Logf("cubicalTatePairing ratio = %s", ratio.String())

	// Compare with Miller loop
	nM, dM := millerLoopSingle(q, &p)
	var ratioM fp.Element
	dM.Inverse(&dM)
	ratioM.Mul(&nM, &dM)
	t.Logf("millerLoopSingle ratio = %s", ratioM.String())

	// Test if squaring helps
	ratio.Square(&ratio)
	t.Logf("cubical ratio^2 = %s", ratio.String())
	t.Logf("f2IsOne(ratio^2) = %v", f2IsOne(&ratio))

	ratioM.Square(&ratioM)
	t.Logf("miller ratio^2 = %s", ratioM.String())
	t.Logf("f2IsOne(miller ratio^2) = %v", f2IsOne(&ratioM))
}

func TestPureCubicalMembershipTest(t *testing.T) {
	tab := precomputeTableDefault()

	// Test with generator point
	_, _, p, _ := curve.Generators()

	result := membershipTestPureCubical(&p, &tab.q)
	t.Logf("membershipTestPureCubical(G1) = %v", result)

	// Also test with a random point
	var a fp.Element
	a.SetRandom()
	pRandom := curve.MapToG1(a)
	result2 := membershipTestPureCubical(&pRandom, &tab.q)
	t.Logf("membershipTestPureCubical(random) = %v", result2)

	// Test if squaring is still needed
	n1, d1 := cubicalTatePairing(&tab.q, &p)
	d1.Inverse(&d1)
	var ratio fp.Element
	ratio.Mul(&n1, &d1)
	t.Logf("Without squaring: f2IsOne(ratio) = %v", f2IsOne(&ratio))
	ratio.Square(&ratio)
	t.Logf("With squaring: f2IsOne(ratio^2) = %v", f2IsOne(&ratio))
}

func TestCubicalPairingRatio(t *testing.T) {
	initMontgomeryConstants()
	tab := precomputeTableDefault()
	q := &tab.q

	_, _, p, _ := curve.Generators()

	// Third root of unity for GLV endomorphism
	var thirdRootOneG1 fp.Element
	thirdRootOneG1.SetString("80949648264912719408558363140637477264845294720710499478137287262712535938301461879813459410945")

	var p2 curve.G1Affine
	p2.X.Mul(&p.X, &thirdRootOneG1)
	p2.Y.Set(&p.Y)

	// Compute cubical pairings
	n1, d1 := cubicalTatePairing(q, &p)
	n2, d2 := cubicalTatePairing(q, &p2)

	t.Logf("e_c(Q, P):    num=%s den=%s", n1.String(), d1.String())
	t.Logf("e_c(Q, φ(P)): num=%s den=%s", n2.String(), d2.String())

	// The denominator Z_{[ℓ]Q} should be the same for both!
	if d1.Equal(&d2) {
		t.Logf("Denominators are EQUAL (as expected)")
	} else {
		t.Logf("Denominators are DIFFERENT")
	}

	// Compute cubical ratio
	var cubicalRatio fp.Element
	d1.Inverse(&d1)
	cubicalRatio.Mul(&n1, &d1)
	t.Logf("Cubical e_c(Q, P) = %s", cubicalRatio.String())

	// Compute Miller ratio (unreduced)
	nM, dM := millerLoopSingle(q, &p)
	var millerRatio fp.Element
	dM.Inverse(&dM)
	millerRatio.Mul(&nM, &dM)
	t.Logf("Miller f(Q, P) = %s", millerRatio.String())

	// The paper claims: e_c^{2ℓ} = f^{2ℓ} (for the unreduced Miller function)
	// where ℓ = e₂ - 1 (the degree used in the ladder)
	// Let's verify: compute e_c^{2*(e₂-1)} and f^{2*(e₂-1)}

	// First, compute (p-1)/e₂ exponentiation (the reduced Tate pairing)
	var cubicalTate, millerTate fp.Element
	// f1IsOne checks x^{(p-1)/e₁} = 1, f2IsOne checks x^{(p-1)/e₂} = 1
	// We want the actual exponentiated values

	// For P ∈ G1: the Miller Tate pairing should equal 1
	// Let's verify by computing the full Tate test
	t.Logf("f2IsOne(cubicalRatio) = %v", f2IsOne(&cubicalRatio))
	t.Logf("f2IsOne(millerRatio) = %v", f2IsOne(&millerRatio))

	// Now let's check if (cubical/miller) is a 2*(e₂-1)-th root of unity
	var quotient fp.Element
	millerRatio.Inverse(&millerRatio)
	quotient.Mul(&cubicalRatio, &millerRatio)
	t.Logf("cubical/miller = %s", quotient.String())

	// Check if quotient^{2*(e₂-1)} = 1
	// 2*(e₂-1) = 2*(z-2) = 2*0x8508bfffffffffff = 0x10117ffffffffffe
	// This is a large exponent, let's use expBySeed equivalent
	// Actually, if cubical = miller * ζ where ζ^{2*(e₂-1)} = 1,
	// then quotient = ζ should satisfy quotient^{2*(e₂-1)} = 1

	// For now, just check small powers
	var q2, q4, q8 fp.Element
	q2.Square(&quotient)
	q4.Square(&q2)
	q8.Square(&q4)
	t.Logf("quotient^2 = %s", q2.String())
	t.Logf("quotient^4 = %s", q4.String())
	t.Logf("quotient^8 = %s", q8.String())

	// Check if quotient is a root of unity by checking quotient^{p-1} = 1
	// (should always be true by Fermat)
	var qPm1 fp.Element
	qPm1.Exp(quotient, new(big.Int).Sub(fp.Modulus(), big.NewInt(1)))
	t.Logf("quotient^{p-1} = %s (should be 1)", qPm1.String())

	// Check quotient^{e₂} and quotient^{2*e₂}
	// e₂ = z - 1 = 0x8508c00000000000
	e2 := new(big.Int)
	e2.SetString("8508c00000000000", 16)
	var qE2 fp.Element
	qE2.Exp(quotient, e2)
	t.Logf("quotient^{e₂} = %s", qE2.String())

	e2x2 := new(big.Int).Mul(e2, big.NewInt(2))
	var qE2x2 fp.Element
	qE2x2.Exp(quotient, e2x2)
	t.Logf("quotient^{2*e₂} = %s", qE2x2.String())

	// Also check cubical^{2*e₂} and miller^{2*e₂}
	n1, d1 = cubicalTatePairing(q, &p)
	d1.Inverse(&d1)
	cubicalRatio.Mul(&n1, &d1)

	nM, dM = millerLoopSingle(q, &p)
	dM.Inverse(&dM)
	millerRatio.Mul(&nM, &dM)

	cubicalTate.Exp(cubicalRatio, e2x2)
	millerTate.Exp(millerRatio, e2x2)
	t.Logf("cubical^{2*e₂} = %s", cubicalTate.String())
	t.Logf("miller^{2*e₂} = %s", millerTate.String())
	t.Logf("Are they equal? %v", cubicalTate.Equal(&millerTate))
}

// BenchmarkG1IsInSubGroupComparison compares all implementations.
func BenchmarkG1IsInSubGroupComparison(b *testing.B) {
	_, _, p, _ := curve.Generators()
	tab := precomputeTableDefault()

	b.Run("GLV", func(b *testing.B) {
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			p.IsInSubGroup()
		}
	})

	b.Run("Tate-Miller", func(b *testing.B) {
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			tab = precomputeTableDefault()
			IsInSubGroupTate(tab, &p)
		}
	})

	b.Run("Tate-Miller-Precomputed", func(b *testing.B) {
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			IsInSubGroupTate(tab, &p)
		}
	})

	b.Run("Tate-Cubical", func(b *testing.B) {
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			IsInSubGroupTateCubical(&p)
		}
	})

	// Benchmark with precomputed Montgomery inputs (simulating Montgomery-native input)
	xP, xP2, xQmP, xQmP2 := ComputeMontgomeryInputs(&p)
	b.Run("Tate-Cubical-MontgomeryNative", func(b *testing.B) {
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			IsInSubGroupCubicalMontgomery(&xP, &xP2, &xQmP, &xQmP2)
		}
	})
}
