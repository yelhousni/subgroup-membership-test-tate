# subgroup-membership-test-tate
Implementing ePrint 2024/1790 on BLS12-381 G1 with precomputations.


## Benchmarks
On a AWS `c7a.8xlarge`:
```
goos: linux
goarch: amd64
pkg: github.com/yelhousni/subgroup-membership-test-tate
cpu: AMD EPYC 9R14
BenchmarkIsInSubGroupCompare
BenchmarkIsInSubGroupCompare/method=Tate
BenchmarkIsInSubGroupCompare/method=Tate-32                24466             49005 ns/op
BenchmarkIsInSubGroupCompare/method=GLV
BenchmarkIsInSubGroupCompare/method=GLV-32                 29142             41155 ns/op
```

The new method seems to be 19% slower compared to Scott's test, but the paper reports it to be 2.7% faster in RELIC (subsection 1.1).

The code related to the paper is available at: https://github.com/eccdaiy39/test-tate. I benchmarked on the same machine:
```
BENCH: g1_is_valid = 191852 cycles
BENCH: g1_is_valid_tate = 282430 cycles
BENCH: g1_is_valid_tate_pre = 182809 cycles
```
which translates, given the base clock ≈ 3.7 GHz on AWS c7a (EPYC Genoa), to:
```
- g1_is_valid: 191 852 cycles × 0.270 ns ≈ 51 800 ns
- g1_is_valid_tate: 282 430 cycles × 0.270 ns ≈ 76 300 ns
- g1_is_valid_tate_pre: 182 809 cycles × 0.270 ns ≈ 49 400 ns
```

The Tate-based methods with precomputations are similar (49400 ns vs. 49024 ns) but the GLV-based ones are not (51800 ns vs. 40189 ns).
Some of the reasons that may explain this are:
- RELIC uses a generic multiplication by the `seed=z` in the GLV method (`g1_mul_any`) while gnark-crypto uses a short addition chain.
- RELIC checks: `psi^2(P) == [-z^2]P` by applying `psi: (x,y)->(w*x,y)` twice instead of precomputing `w^2 mod r`.
