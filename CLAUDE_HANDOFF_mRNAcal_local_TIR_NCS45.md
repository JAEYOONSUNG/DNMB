# Claude Handoff: DNMB mRNAcal Local TIR + NCS45 Work

## Repository Context

- Repo root: `/Users/JaeYoon/Dropbox/0.Personal folder/5. Bioinformatics/DNMB`
- Current working dir used by Codex: `/Users/JaeYoon/Dropbox/0.Personal folder/5. Bioinformatics/DNMB/data`
- Branch: `master`
- Last pushed commit before this handoff: `b163bbc Enhance mRNAcal with RNAplfold and RNAduplex`
- Current state: **uncommitted local changes** in 6 files.

## User Intent

The user wants DNMB `mRNAcal` to become a biologically useful translation-initiation scoring module, not just a generic RNA secondary-structure wrapper.

Key intent:

1. For bacterial translation, **local start-region accessibility matters more than whole-transcript/whole-ORF fold**.
2. `RNAplfold` should be the main fold/accessibility signal because it estimates local unpaired probability.
3. Runtime must improve. Do not run RNAplfold on unnecessarily long windows.
4. Include N-terminal coding sequence effects from Rongzhen Tian's NCS work, especially the first coding bases after `ATG`.
5. Lysine coding sequence signal should be explicit, not hidden.

Important papers/user references:

- Tian R. et al. 2019, *Synthetic N-terminal coding sequences for fine-tuning gene expression and metabolic engineering in Bacillus subtilis*. DOI: `10.1016/j.ymben.2019.07.001`
- Wang C. et al. 2022, *Model-driven design of synthetic N-terminal coding sequences for regulating gene expression in yeast and bacteria*. DOI: `10.1002/biot.202100655`
- RNAplfold method/manual:
  - DOI: `10.1093/bioinformatics/btk014`
  - https://viennarna.readthedocs.io/en/latest/man/RNAplfold.html

## Current Uncommitted Files

```text
M R/Output_combiner.R
M R/db_module_mrnacal.R
M R/module_runner.R
M R/run_DNMB.R
M R/run_cache.R
M README.md
```

Diff size at handoff:

```text
6 files changed, 320 insertions(+), 42 deletions(-)
```

## What Was Implemented

### 1. RNAplfold now uses a smaller local TIR window

Before:

- RNAplfold used the full extracted transcript window, normally `-60/+60` around start codon, median 120 nt.

Now:

- The full window is still extracted for context and RNAfold visualization.
- RNAplfold is run on a smaller local TIR window:

```text
start-centered local window = -35 / +30
median local length in test = 65 nt
```

This is implemented through `.dnmb_mrnacal_local_tir_window()`.

Coordinates are mapped back through `plfold_offset`, with local columns:

- `plfold_sequence`
- `plfold_offset`
- `plfold_local_upstream`
- `plfold_local_downstream`
- `rbs_plfold_start`, `rbs_plfold_end`
- `start_plfold_access_start`, `start_plfold_access_end`
- `tir_plfold_access_start`, `tir_plfold_access_end`
- `standby_plfold_start`, `standby_plfold_end`

### 2. Local accessibility now includes multiple biologically relevant regions

RNAplfold now calculates unpaired probabilities for:

- RBS region
- start-region window
- `-18..+10` TIR window relative to start codon
- upstream standby site

New/expanded result columns:

- `rbs_plfold_unpaired_probability`
- `start_plfold_unpaired_probability`
- `tir_plfold_unpaired_probability`
- `standby_plfold_unpaired_probability`
- `plfold_accessibility_score`

Composite RNAplfold score uses weighted mean:

```text
0.30 * RBS
0.25 * start
0.35 * TIR(-18..+10)
0.10 * standby
```

### 3. Upstream A/U-rich enhancer context added

Added helper `.dnmb_mrnacal_upstream_context()`.

New features:

- `upstream20_sequence`
- `upstream20_at_fraction`
- `upstream20_score`
- `upstream_enhancer_sequence`
- `upstream_enhancer_at_fraction`
- `upstream_enhancer_at_run`
- `upstream_au_score`

This captures A/U-rich upstream enhancer/standby context.

### 4. Tian-style N-terminal coding sequence / NCS45 features added

The user specifically corrected that the relevant Rongzhen Tian paper is about **N-terminal coding sequences**, not only upstream bases.

Current implementation adds:

- codon 2-8 lysine/AAA/AAG features
- `ATG` downstream `NCS45` features, meaning the 45 bp after start codon, represented as codons 2-16 when available

New columns:

- `lysine_codon_count_2_8`
- `aaa_count_2_8`
- `aag_count_2_8`
- `early_coding_at_fraction`
- `ncs45_sequence`
- `ncs45_at_fraction`
- `ncs45_lysine_codon_count`
- `ncs45_aaa_count`
- `ncs45_aag_count`
- `early_aaa_run`
- `early_poly_a_run`
- `early_poly_a_penalty`
- `skik_like`

Important implementation nuance:

- Lysine/AAA-AAG is rewarded.
- But long homopolymeric `AAA` / poly-A runs are penalized because repeated AAA can cause ribosome sliding / expression loss risk.
- This is deliberate: do not blindly reward all AAA-rich sequences.

Current early coding score logic lives in `.dnmb_mrnacal_early_context()`.

### 5. Composite score weights changed

Current `tir_score` weights:

```text
0.20 * RBS
0.18 * antiSD_duplex
0.32 * local_RNAplfold_access
0.10 * upstreamAU
0.10 * start
0.07 * earlyK/NCS coding context
0.03 * MFE
```

Rationale:

- Local accessibility is now the strongest component.
- RNAfold MFE remains only a small supporting term.
- Early coding/NCS signal is present but not allowed to dominate until more validation.

### 6. Cache invalidation updated

In `R/run_cache.R`, mRNAcal algorithm signature was changed to:

```r
mrnacal_algorithm = "rnaplfold_local_tir_ncs45_v3"
```

This prevents stale reuse of old RNAfold/RNAplfold-only results.

### 7. Plot updated

`mRNAcal_translation_efficiency.pdf` scatter plot now uses:

```text
Local Accessibility vs Score
```

instead of anti-SD binding energy as the x-axis.

Component plot now includes:

```text
UpstreamAU
```

## Verification Already Run

### R parsing / loading

These passed:

```bash
Rscript -e 'files <- c("../R/db_module_mrnacal.R", "../R/Output_combiner.R", "../R/run_cache.R"); invisible(lapply(files, parse)); cat("parse_ok\n")'
Rscript -e 'files <- list.files("../R", pattern="[.]R$", full.names=TRUE); invisible(lapply(files, parse)); cat("parse_all_ok\n")'
Rscript -e 'pkgload::load_all("..", quiet=TRUE); cat("load_all_ok\n")'
git diff --check
```

### Docker source-mount smoke test

Because local mac PATH lacks RNAfold, live mRNAcal tests were run inside the existing Docker image while mounting the modified source:

```bash
docker run --rm \
  -v "/Users/JaeYoon/Dropbox/0.Personal folder/5. Bioinformatics/DNMB:/src" \
  -v "/Users/JaeYoon/Dropbox/0.Personal folder/5. Bioinformatics/DNMB/data:/data" \
  -w /data \
  ghcr.io/jaeyoonsung/dnmbsuite:latest \
  Rscript -e 'pkgload::load_all("/src", quiet=TRUE); ...'
```

Smoke result for 10 genes:

```text
rows=10 plfold=10 earlyK=10 ncs45=10 penalty=10
```

Earlier local-TIR test result:

```text
rows=10 folds=10 plfold=10 tir=10 standby=10 upAU=10 duplex=10
summary_pdf=TRUE fold_pdf=TRUE
```

### RNAplfold speed comparison

Benchmark on 120 genes, comparing RNAplfold only:

```text
genes=120
local_elapsed=0.687
full_elapsed=1.272
median_local_len=65
median_full_len=120
```

Interpretation:

- Local TIR RNAplfold is ~46% faster than the old full-window RNAplfold run for this test.
- This supports the user's speed requirement.

## Important Caveats / Next Checks

1. **Not committed yet.**
   - Commit/push has not been done after these latest mRNAcal local-TIR/NCS45 changes.

2. **Docker image not rebuilt yet.**
   - Current `ghcr.io/jaeyoonsung/dnmbsuite:latest` includes the previous commit `b163bbc`, not these latest changes.
   - Tests used source mount (`pkgload::load_all("/src")`) to validate the current working tree.

3. **Need one final end-to-end Docker installed-package test after build.**
   - After building new Docker image, run mRNAcal without source mount to confirm the installed package contains the changes.

4. **Score calibration is heuristic.**
   - The module now reports biologically motivated features, but exact score weights are pragmatic and should be treated as a ranking heuristic.
   - NCS45/Tian-style modeling could later become a separate learned model if training data is added.

5. **Potential improvement: species-specific NCS model.**
   - Tian paper is Bacillus subtilis-focused.
   - DNMB currently applies general NCS45 features broadly.
   - Could later condition weights by taxonomy/translation domain.

## Suggested Next Steps for Claude

1. Inspect current diff:

```bash
git diff -- R/db_module_mrnacal.R R/Output_combiner.R R/run_cache.R README.md
```

2. Run final validations:

```bash
Rscript -e 'files <- list.files("../R", pattern="[.]R$", full.names=TRUE); invisible(lapply(files, parse)); cat("parse_all_ok\n")'
Rscript -e 'pkgload::load_all("..", quiet=TRUE); cat("load_all_ok\n")'
git diff --check
```

3. Run Docker source-mount smoke test again if desired.

4. If OK, commit:

```bash
git add R/Output_combiner.R R/db_module_mrnacal.R R/module_runner.R R/run_DNMB.R R/run_cache.R README.md
git commit -m "Prioritize local TIR accessibility in mRNAcal"
git push
```

5. Build and push Docker:

```bash
docker build --progress=plain -t ghcr.io/jaeyoonsung/dnmbsuite:mrnacal-local-tir .
docker tag ghcr.io/jaeyoonsung/dnmbsuite:mrnacal-local-tir ghcr.io/jaeyoonsung/dnmbsuite:latest
docker push ghcr.io/jaeyoonsung/dnmbsuite:latest
```

6. Final Docker smoke test after build:

```bash
docker run --rm \
  --user "$(id -u):$(id -g)" \
  -v "$PWD/data:/data" \
  -v "$HOME/.dnmb-cache:/opt/dnmb-cache" \
  ghcr.io/jaeyoonsung/dnmbsuite:latest \
  Rscript -e 'library(DNMB); ... run mRNAcal smoke ...'
```

Adjust volume paths depending on current working directory.

## Short Korean Summary for User

현재 작업은 mRNAcal을 whole RNAfold 중심에서 start 주변 local TIR accessibility 중심으로 바꾸는 중이다. RNAplfold 입력을 `-35/+30`으로 줄여 속도를 개선했고, RBS/start/TIR/standby unpaired probability를 따로 계산한다. 또한 Rongzhen Tian 2019 논문 흐름에 맞춰 `ATG` 뒤 `NCS45` 및 early lysine/AAA/AAG coding features를 추가했다. 현재 검증은 통과했지만 아직 commit/push 및 Docker rebuild/push는 하지 않았다.
