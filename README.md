# RLWE based Zero-Knowledge for commitments

This GitHub repository provides the implementation of the lattice-based post-quantum commitment scheme introduced by Martínez and Morillo in https://eprint.iacr.org/2019/1486.

Our implementation is very flexible and its parameters can be adjusted to obtain a trade-off between speed and memory usage, analyzing how suitable for practical use are the underlying lattice-based techniques.

**Disclaimer**: this implementation is only an academic proof of concept and is NOT ready for production.

Building the code:
----

To build the code and/or the dependencies make sure `build-essential`, `m4` and `libgmp-dev` are installed in your system.

This code depends on the following libraries (the specific tested versions are mentioned in parenthesis):

- [OpenSSL](https://www.openssl.org) (3.0.2)
- [GNU MP](https://gmplib.org) (6.2.1)
- [FLINT](https://flintlib.org) (2.9.0)

Once the dependencies are installed the code can be built running `make`.

Generating parameters:
----

In order to be able to generate new secure sets of parameters [`sagemath`](https://www.sagemath.org) has to be available and https://github.com/malb/lattice-estimator has to be added as submodule:

```
git submodule add https://github.com/malb/lattice-estimator latticeEstimator
```

Then to generate a new set of parameters with λ = 100, n = 512, q = 16381 and d = 2 just run the following line.

```
sage generate_params.sage 100 512 16381 2
```

And the full set [λ, n, q, d, k, σ, B, δOL, δM] will be computed and added to `generatedParams.csv`.

The destination file can be customized including it as an additional parameter with `-o` or `--output`.

```
sage generate_params.sage 100 512 16381 2 --output generatedParams2.csv
```


Benchmarking the protocols:
----

The general syntax for benchmarking the protocols is the following.

<pre><code>python3 run_tests.py <i>algorithm</i> <i>parameters</i></code></pre>

Where <code><i>algorithm</i></code> can be `keygen`, `commitment`, `verifier`, `opening`, `linear` or `multiplicative`.

Additionally <code><i>algorithm</i></code> can be `commit` to sequentially benchmark all commitment related protocols (`keygen`, `commitment` and `verifier`), `NIZKP` to benchmark all the Non-Interactive Zero-Knowledge Proofs (`opening`, `linear` and `multiplicative`) or `all` to run all of them.

Finally <code><i>parameters</i></code> would indicate the line number (minus two) of the set of parameters from `generatedParameters.csv` that is going to be benchmarked or `-1` if all the sets of parameters should be benchmarked sequentially.

Benchmark results indicating both running time in milliseconds and processor cycles are saved to specific files in the `./benchmarks` folder.

The number of iterations is defined by a variable `it` in the `run_tests.py` file.

If one wants to get the parameter sets from a different file then both `run_tests.py` and `./src/param.h` have to be edited modifying variables `filename` and `FILEPARAMETERS` respectively.
