# Kakutani–KMS Theory of Ligand-Induced Conformational Phase Transitions

**Measure-Class Obstructions, Covariance Criticality, and Tail Allostery in Computational Pharmacology**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23042828.svg)](https://doi.org/10.5281/zenodo.23042828)

Zhengyi Chen<sup>1</sup>, Ruqing Chen<sup>2</sup>

<sup>1</sup> Guangxi Key Laboratory of Drug Discovery and Optimization, School of Pharmacy, Guilin Medical University, Guilin 541199, P. R. China — chenzhengyi@glmc.edu.cn
<sup>2</sup> GUT Geoservice Inc., Montreal, Québec, Canada — ruqing@hotmail.com

This manuscript integrates the series *Statistical Pharmacology via Kakutani Dichotomy*, Papers I–VII, into one article (`paper/kakutani_pharmacology_paper.tex`, 16 pages).

## Content

1. **Introduction.** The blindness of local and tensorized invariants. Theorem 1.1 is proved in full: identical pocket laws, identical per-residue marginals and equal binding free energies are compatible with mutually singular global ensembles, and singularity sets in exactly at decay exponent 1/2. Proposition 1.2 shows that pocket-local energies can never produce singularity.
2. **First-order index.** The marginal index DKI and the critical exponent α_c = 1/2 (Paper I).
3. **Correlated index.** CKI, the exact three-way decomposition, and 99.95 % covariance dominance, with the Szegő rate recomputed independently (Paper II).
4. **Collective modes.** Spectral Kakutani–Feldman–Hájek criticality of amplitudes and rotations (Paper III).
5. **KMS scaling and broken phases.** KMS scaling β_c = 1, tail allostery, and the chemical obstruction group (Papers IV and V).
6. **Transport cost.** The cone-geodesic Wasserstein cost W_tail and alchemical path topology (Paper VI).
7. **Pipeline and representation.** The MD estimation pipeline `kakutani_pharma` (Paper VII), and the tail-augmented representation z_local ⊕ z_symmetry ⊕ z_tail with a proposition on its information content.
8. **Discussion.** What is established, what is not, and what follows.

All numerical evidence comes from analytic models and synthetic trajectories. Validation on all-atom protein–ligand molecular dynamics has not been done.

## Verification

`code/verify_synthesis_theorems.py` checks the statements that this article proves or recomputes itself, and writes `results/synthesis_verification.txt`:

- the expansion of the local index to eighth order (the remainder over δ⁸ tends to 429/64);
- the pair divergence of Theorem 1.1(2) against a direct Bhattacharyya computation, and its partial sums for α = 1/4, 1/2 and 3/4;
- Proposition 1.2 on a Gaussian example;
- the Szegő rate s = 0.0094968218, the constant c₀ = −0.0093715854 for N = 50 to 1000, and D_mean, D_cov and CKI at N = 501.

All checks pass. The theorems of Papers I–VII are verified by the scripts archived with each paper.

## Layout

```
.
├── paper/kakutani_pharmacology_paper.{tex,pdf}
├── figures/                              # figures of Papers IV-VII
├── code/verify_synthesis_theorems.py
├── results/synthesis_verification.txt
├── CITATION.cff, LICENSE, .gitignore
└── README.md
```

## Build

```bash
pip install numpy scipy mpmath
python code/verify_synthesis_theorems.py     # prints OVERALL: PASS
cd paper
pdflatex kakutani_pharmacology_paper.tex && pdflatex kakutani_pharmacology_paper.tex && pdflatex kakutani_pharmacology_paper.tex
```

The figures in `figures/` are copied from Papers IV–VII.

## Repository

https://github.com/Ruqing1963/kakutani-statistical-pharmacology-synthesis

## Sources

Paper I 10.5281/zenodo.23005921 · II 10.5281/zenodo.23012217 · III 10.5281/zenodo.23015709 · IV 10.5281/zenodo.23018236 · V 10.5281/zenodo.23027704 · VI 10.5281/zenodo.23029360 · VII 10.5281/zenodo.23039869 · monograph 10.5281/zenodo.22883179.

## Citation

```bibtex
@article{ChenChen2026KakutaniKMS,
  author  = {Chen, Zhengyi and Chen, Ruqing},
  title   = {Kakutani--KMS Theory of Ligand-Induced Conformational Phase Transitions: Measure-Class Obstructions,
             Covariance Criticality, and Tail Allostery in Computational Pharmacology},
  year    = {2026},
  doi     = {10.5281/zenodo.23042828},
  url     = {https://doi.org/10.5281/zenodo.23042828}
}
```

## Funding

National Natural Science Foundation of China (No. 22464010), Guangxi Natural Science Foundation of China (No. 2025GXNSFAA069294), and the 2025 Bagui Youth Top Talent Project.

## License

Code under the MIT License; text and figures under CC BY 4.0 (see `LICENSE`).
