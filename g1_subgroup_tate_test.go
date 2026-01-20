package bls12381

import (
	curve "github.com/consensys/gnark-crypto/ecc/bls12-381"
	"github.com/consensys/gnark-crypto/ecc/bls12-381/fp"

	"github.com/leanovate/gopter"
	"github.com/leanovate/gopter/prop"
	"testing"
)

// GenFp generates an Fp element
func GenFp() gopter.Gen {
	return func(genParams *gopter.GenParameters) *gopter.GenResult {
		var elmt fp.Element
		elmt.MustSetRandom()

		return gopter.NewGenResult(elmt, gopter.NoShrinker)
	}
}

func TestG1SubGroupMembershipTate(t *testing.T) {
	t.Parallel()
	tab := precomputeTableDefault()
	parameters := gopter.DefaultTestParameters()
	if testing.Short() {
		parameters.MinSuccessfulTests = 10
	} else {
		parameters.MinSuccessfulTests = 100
	}

	properties := gopter.NewProperties(parameters)

	properties.Property("[BLS12-381] IsInSubGroupTate should output true for points on G1", prop.ForAll(
		func(a fp.Element) bool {
			p := curve.MapToG1(a)
			return IsInSubGroupTate(tab, &p)
		},
		GenFp(),
	))

	properties.TestingRun(t, gopter.ConsoleReporter(false))
}

// bench
func BenchmarkG1IsInSubGroup(b *testing.B) {
	_, _, p, _ := curve.Generators()
	tab := precomputeTableDefault()

	b.Run("GLV", func(b *testing.B) {
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			p.IsInSubGroup()
		}
	})

	b.Run("Tate", func(b *testing.B) {
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			IsInSubGroupTate(tab, &p)
		}
	})
}
