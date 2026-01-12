// Copyright 2020-2026 Consensys Software Inc.
// Licensed under the Apache License, Version 2.0. See the LICENSE file for details.

package bls12381

import (
	"fmt"
	"testing"

	"github.com/leanovate/gopter"
	"github.com/leanovate/gopter/prop"

	curve "github.com/consensys/gnark-crypto/ecc/bls12-381"
	"github.com/consensys/gnark-crypto/ecc/bls12-381/fp"
)

func TestIsInSubGroupTatePrecomputed(t *testing.T) {
	t.Parallel()
	parameters := gopter.DefaultTestParameters()
	if testing.Short() {
		parameters.MinSuccessfulTests = 20
	} else {
		parameters.MinSuccessfulTests = 100
	}

	properties := gopter.NewProperties(parameters)

	properties.Property("[BLS12-381] IsInSubGroupTatePrecomputed and IsInSubGroup should return the same result on the prime-torsion", prop.ForAll(
		func(a fp.Element) bool {
			aff := curve.MapToG1(a)
			return IsInSubGroupTatePrecomputed(&aff) != aff.IsInSubGroup()
		},
		GenFp(),
	))
}

func BenchmarkIsInSubGroupCompare(b *testing.B) {
	var a curve.G1Affine
	a.X.SetString("3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054947676345821548104185464507")
	a.Y.SetString("1339506544944476473020471379941921221584933875938349620426543736416511423956333506472724655353366534992391756441569")
	b.Run(fmt.Sprintf("method=Tate"), func(b *testing.B) {
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			IsInSubGroupTatePrecomputed(&a)
		}
	})

	b.Run(fmt.Sprintf("method=GLV"), func(b *testing.B) {
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			a.IsInSubGroup()
		}
	})

}

// GenFp generates an Fp element
func GenFp() gopter.Gen {
	return func(genParams *gopter.GenParameters) *gopter.GenResult {
		var elmt fp.Element
		elmt.MustSetRandom()

		return gopter.NewGenResult(elmt, gopter.NoShrinker)
	}
}
