# REBASEfinder structural templates

This directory contains compact, chain-only experimental structures for local
REBASEfinder homology and Foldseek validation. The files cover the 11 curated
entries in `../rebasefinder_structure_refs.tsv`, add a Type I HsdR reference
that was missing from that set, and add generic helicase and DNA-repair nuclease
decoys for negative-control ranking.

## Contents

- `manifest.tsv` records role, family, source, coordinate-residue count, and
  SHA-256 checksums for every asset.
- `*.pdb.gz` contains one author-assigned protein chain. Only normalized `ATOM`,
  `TER`, and `END` records are retained.
- `*.faa` is derived from the residues that have coordinates in the paired PDB,
  not from `SEQRES`. Unresolved residues are therefore absent.

Complexes with distinct RM roles are split into role-specific templates. A
single representative is retained when an entry contains duplicate copies of
the same polymer entity. DNA, ions, waters, cofactors, ligands, unrelated
chains, hydrogen atoms, `ANISOU`, and non-selected alternate conformers are
removed. Protein selenomethionine (`MSE`) residues are normalized to standard
methionine `ATOM` records so the coordinate PDB and FASTA remain aligned.

The source PDB files were retrieved from the RCSB PDB download service on
2026-07-10:

`https://files.rcsb.org/download/<PDB_ID>.pdb.gz`

Each source entry is linked from `manifest.tsv`.

## HsdR selection

The Type I R-subunit candidates were all EcoR124I HsdR structures of the same
protein, so retaining all three would add conformational duplication rather
than family diversity.

| PDB | Method and resolution | Best author chain | Coordinate residues | Decision |
| --- | --- | --- | ---: | --- |
| 2W00 | X-ray, 2.60 A | B | 852 | C-terminal domain largely unresolved |
| 6H2J | X-ray, 2.60 A | B | 990 | Selected; best high-resolution full-length coverage |
| 7BST | cryo-EM, 4.37 A | C | 992 | Not selected; lower-resolution Ocr-bound complex |

`6H2J:B` retains the N-terminal nuclease region, tandem RecA-like ATPase motor,
and C-terminal HsdR domain. Chain A of the same entry has substantially less
coordinate coverage and was not bundled.

The `2Y7C` M/S templates are a fitted atomic model from an 18 A electron
microscopy map. They are retained to represent the existing curated reference,
but should be treated as coarse homology templates rather than high-resolution
motif-geometry evidence.

## Decoy selection

Three non-RM templates provide explicit negative controls:

| PDB | Family | Method and resolution | Coordinate residues | Purpose |
| --- | --- | --- | ---: | --- |
| 1QYR:A | KsgA rRNA adenine dimethyltransferase | X-ray, 2.10 A | 252 | RNA methylase decoy |
| 1OYW:A | E. coli RecQ SF2 helicase catalytic core | X-ray, 1.80 A | 516 | Generic helicase decoy |
| 1QTW:A | E. coli Endonuclease IV | X-ray, 1.02 A | 285 | DNA-repair nuclease decoy |

RecQ was selected instead of the 1C4O UvrB structure because 1OYW has nearly
complete coordinate coverage for its deposited catalytic-core construct (516
of 523 residues, including normalized MSE), whereas 1C4O resolves only 504 of
664 deposited residues. Endonuclease IV is a complete monomeric repair enzyme
with all 285 residues modeled at high resolution. These templates are marked
`template_class=decoy` and `enzyme_role=exclude`; they must not independently
create an RM assignment.

## License and attribution

The wwPDB states that PDB archive data files are available under the
[CC0 1.0 Universal Public Domain Dedication](https://creativecommons.org/publicdomain/zero/1.0/).
The current RCSB policy is documented at
[RCSB PDB Usage Policies](https://www.rcsb.org/pages/usage-policy).

These chain-only files are transformations of those CC0 coordinate data. RCSB
encourages attribution to the original structure authors where possible; the
PDB entry links in `manifest.tsv` provide the corresponding citation and
deposition provenance. DNMB source code licensing does not alter the CC0 status
of these coordinate-derived assets.

## Integrity

Checksums are over the distributed files exactly as stored in this directory.
The gzip streams use a zero modification time for reproducibility. Consumers
should resolve `pdb_file` and `fasta_file` relative to `manifest.tsv` and may
verify both SHA-256 columns before building a local search database.
