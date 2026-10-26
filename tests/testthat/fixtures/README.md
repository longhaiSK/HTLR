# Legacy reference data

`legacy-3.1-1.rds` stores reference results from HTLR.old 3.1-1. Tests read
this file and do not load or install the legacy package. The tests run on CI
and CRAN. The file contains the inputs and results from the four original
sampler comparisons and the two original data generator comparisons.

The sampler cases cover binary and multiclass responses, the `t`, `neg`,
and `ghs` priors, and BCBCSF initialization. NEG and GHS use `eta = 0`, which
was the default in the original comparisons. The saved fields are coefficient
traces, log scales, log likelihoods, updated feature counts, rejection rates,
variances, and prediction probabilities. Tests allow a relative numerical
tolerance of `1e-7` for sampler results. This allows small floating-point
differences between platforms; it does not test posterior convergence.

## Source and build

The source is commit `ac8e7b34bf3a039ae2042cedcd67693ec7379c62` on the
local `legacy` branch. This commit fixes callbacks in the original legacy
source from <https://github.com/longhaiSK/HTLR>. It moves three nested callbacks
to file scope and passes their state through context structs. It keeps the
formulas and RNG call order. The fixed package passed all original live
comparisons before these fixtures were generated.

The build passed with `-Werror=trampolines`, without a linker override for
stack permissions. `readelf -W -l` showed a `GNU_STACK` entry with `RW`
permissions. Direct callback regression checks also passed.

## Rebuild the reference file

Fixture generation is a separate maintenance task. Normal tests do not need
HTLR.old. Do not regenerate expected results to make a failed test pass.
First determine the cause and review the change in behavior.

1. Check out the source commit above in a separate worktree.
2. Install the fixed HTLR.old package and its dependencies in a library
   available to R. Do not add it to `DESCRIPTION` or `renv.lock`.
3. From this repository root, run
   `Rscript -e 'source("data-raw/generate-legacy-fixtures.R")'`.
4. Review the results and run `Rscript -e 'devtools::test()'`.

The fixture metadata records the seed (1001), all three RNG settings, R and
package versions, platform, session information, and source commit.
Inputs are stored so that changes to data generation or splitting cannot
change the sampler reference cases. No test generates expected results from
the current implementation.
