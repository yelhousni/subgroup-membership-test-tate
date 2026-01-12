# subgroup-membership-test-tate
Implementing ePrint 2024/1790 on BLS12-381 G1 with precomputations.


## Benchmark
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
