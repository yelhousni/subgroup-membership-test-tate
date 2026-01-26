package bls12381

import (
	curve "github.com/consensys/gnark-crypto/ecc/bls12-381"
	"github.com/consensys/gnark-crypto/ecc/bls12-381/fp"

	"github.com/leanovate/gopter"
	"github.com/leanovate/gopter/prop"
	"testing"
)

func TestG1SubGroupMembershipTateTwoLoops(t *testing.T) {
	t.Parallel()
	tab := precomputeTableDefault()
	parameters := gopter.DefaultTestParameters()
	if testing.Short() {
		parameters.MinSuccessfulTests = 10
	} else {
		parameters.MinSuccessfulTests = 100
	}

	properties := gopter.NewProperties(parameters)

	properties.Property("[BLS12-381] IsInSubGroupTateTwoLoops should output true for points on G1", prop.ForAll(
		func(a fp.Element) bool {
			p := curve.MapToG1(a)
			return IsInSubGroupTateTwoLoops(tab, &p)
		},
		GenFp(),
	))

	properties.TestingRun(t, gopter.ConsoleReporter(false))
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

	properties.Property("[BLS12-381] IsInSubGroupTateCubical should output true for points on G1", prop.ForAll(
		func(a fp.Element) bool {
			p := curve.MapToG1(a)
			return IsInSubGroupTateCubical(&p)
		},
		GenFp(),
	))

	properties.TestingRun(t, gopter.ConsoleReporter(false))
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

	b.Run("Tate-SharedLoop", func(b *testing.B) {
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			IsInSubGroupTate(tab, &p)
		}
	})

	b.Run("Tate-TwoLoops", func(b *testing.B) {
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			IsInSubGroupTateTwoLoops(tab, &p)
		}
	})

	b.Run("Tate-NoPrecompute", func(b *testing.B) {
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			IsInSubGroupTateCubical(&p)
		}
	})
}
