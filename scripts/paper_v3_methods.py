"""Online Methods text for manuscript v3 (Nature Methods Analysis format).

Returns (subheading, [paragraphs]) blocks. Numbers that depend on results are filled from
the report JSON files passed in by the builder.
"""

from __future__ import annotations


def methods(j: dict) -> list[tuple[str, list[str]]]:  # type: ignore[type-arg]
    s4, g3, rt, pm = j["step4_domains"], j["gate3_dbit_audit"], j["step20_runtime"], j["step25_perturbmap_lesions"]
    xen = g3["Xenium HT-29 lung (GSE319123, single cell)"]
    hw = rt["hardware"]
    return [
        ("Data sources", [
            "Discovery data were the Xenium validation section of the Perturb-DBiT study (GEO GSE319123; one formalin-fixed section of a mouse lung colonised by HT-29 cells after a genome-scale Brunello screen)[@baysoy2026,janesick2023]. "
            f"The cell-feature matrix contains {352913:,} cells and 10,029 features: 5,001 genes of the Xenium Prime 5K human panel, 27 custom guide probes and control codewords. "
            "DBiT-seq sections (GSE319277)[@liu2020] and Perturb-map Visium sections (GSE193460)[@dhainaut2022] were used for audits and the positive control. "
            "The curated lesion annotations of the Perturb-map imaging deposit (BioImage Archive S-BIAD267) were used to count lesions per knockout. "
            "A second tissue was a public Xenium mouse brain coronal subset (Zenodo record 11209134, CC BY 4.0), available only as transcript coordinates. "
            "All files were downloaded from their repositories and checked against published sizes or by decompression."
        ]),
        ("Cell filtering and normalisation", [
            f"Cells with at least 20 counts over panel genes were kept ({s4['cells_kept_min_counts']:,} cells; median {s4['median_total_counts_kept']:.0f} counts). "
            "Counts were scaled to the median total and log1p-transformed. Analysis genes were those with mean normalised expression at least 0.1 in control tumour cells "
            f"({s4['analysis_genes']} genes); the sparsity of the panel in this tissue is the reason so few genes passed."
        ]),
        ("Tumour-cell calling", [
            "A tumour score was the mean z-score of the HT-29 markers present in the panel (EPCAM, KRT20, CDX2). A two-component Gaussian mixture separated tumour cells (posterior > 0.9). "
            f"The mixture is degenerate in this sample: 69% of cells carry no marker transcript, posteriors are 0 or 1, and the call is equivalent to detecting at least one marker transcript ({s4['tumour_cells']:,} tumour cells). "
            "A stricter definition (top half of the score) was used as a sensitivity analysis."
        ]),
        ("Guide assignment", [
            "Guide probes were identified as expression features without an Ensembl identifier. A tumour cell was assigned to a knockout if exactly one guide probe had at least one UMI; cells with two or more different guides were excluded. "
            f"Negative-control probes were positive in a fraction 1.1 × 10⁻⁵ of cells, so the {len(xen['testable_upper_bound'])} knockouts with at least 60 positive cells were at least "
            f"{xen['min_enrichment_among_testable']:.0f}-fold above chance. A sensitivity analysis required at least two UMIs. Spillover was measured as the fraction of non-knockout cells with guide signal as a function of distance to the nearest knockout cell."
        ]),
        ("Niche construction", [
            "Thirty principal components were computed on control tumour cells. For every cell, niche features were the kernel-weighted mean of the components of control cells within 100 µm, "
            "with weight exp(−d²/2l²), l = 50 µm, excluding the cell itself[@hu2021]. K-means with K = 4 (the value chosen among 3-6 by silhouette on 20,000 cells in the registered Step 4) was fitted on control cells[@pedregosa2011]; this caller is used for all benchmarks and for the real-data comparison of tests, so every test sees the same niches. "
            "For pseudo-knockouts made of control cells, niche features were recomputed with the pseudo-knockout cells removed from the reference and assigned to the fixed centroids (leave-out labels), so that, like real knockout cells, they never contribute to their own niche definition. "
            "Stability was measured as the adjusted Rand index between the reference labelling and 20 refits on 80% subsamples. "
            "External niche callers were run on the same control cells: pybanksy 1.3.5 (λ = 0.8, Leiden; the first resolution with no niche above 50% of cells), "
            "SpaGCN 1.2.7 (no histology, k-means initialisation with 4 clusters because the louvain dependency could not be built on Windows; fitted on a uniform 6,000-cell subsample because its dense adjacency matrix exceeds memory for all cells, then propagated by a five-nearest-neighbour vote) "
            "and CellCharter 0.3.8 (Delaunay graph, long links removed, three-layer neighbourhood aggregation of 20 principal components in place of scVI, Gaussian mixture with 4 clusters)[@singhal2024,hu2021,varrone2024]."
        ]),
        ("Residual model", [
            "For each analysis gene g, control cells were fitted by least squares to y_ig = μ_d(i)g + β_gᵀ(c_i − c̄) + ε_ig, where d(i) is the niche of cell i and c_i holds log total counts, host-cell fraction "
            "(non-tumour cells among cells within 50 µm) and local density (cells within 50 µm). The residual of any cell is r_ig = y_ig − μ̂_d(i)g − β̂_gᵀ(c_i − c̄). Fits use the normal equations so the expression matrix stays sparse."
        ]),
        ("Tests of niche dependence", [
            "All tests compare a knockout's residuals between niches. Units: 30 µm clonal patches and 100 µm lesions are connected components of single-linkage graphs; tiles are 250 µm squares. "
            "Each unit takes its majority niche, and niches with fewer than 30 cells after this assignment are dropped (applying the minimum before assignment let a 7-cell niche drive an apparent signal).",
            "Sum-of-squares permutation: T = Σ_d n_d‖ē_d − ē‖² over genes, with niche labels permuted among cells, patches, lesions or tiles (1,000 permutations in real data, 250 in benchmarks). "
            "Studentised permutation: the sum over genes of Welch's F on unit means, permuted among units[@chung2013]. "
            "Pseudobulk ANOVA: one-way ANOVA per gene on unit means, Benjamini-Hochberg across genes[@benjamini1995]. "
            "Cluster-robust regression: per-gene OLS on niche indicators with CR1 standard errors clustered on patches, Wald F on k − 1 and G − 1 degrees of freedom. "
            "CR2: Bell-McCaffrey bias-reduced covariance with Imbens-Kolesár degrees of freedom for the contrast between the two largest niches, clustered on lesions[@bell2002,imbens2016,pustejovsky2018,cameron2015]; the implementation equals statsmodels HC2 when every cluster is a single cell (unit test). "
            "Wild cluster bootstrap: the summed per-gene CR1 Wald statistic with the restricted (no-niche) model and Webb six-point weights at the lesion level, 299-999 draws[@cameron2008,mackinnon2017,webb2023]. "
            "Within-lesion sign-flip test: for each lesion with at least three cells in each of the two largest niches, the difference of mean residuals between the niches; the statistic is the squared mean contrast summed over genes, "
            "and its null distribution flips the sign of whole lesion contrasts (exact when 2^L ≤ 2,000). Because each contrast is taken inside one lesion, the lesion's own baseline cancels. "
            "Mixed model: statsmodels MixedLM on lesion × niche pseudobulk with a lesion random intercept and a Wald χ² test per gene[@seabold2010,zimmerman2021] (Satterthwaite or Kenward-Roger degrees of freedom, and dream[@hoffman2021], need R and were not run). A test was 'applicable' to a draw when its design was estimable: lesion-unit tests need at least two lesions in at least two niches; CR2, the wild bootstrap and the sign-flip test need every compared niche to be seen in at least two lesions with three or more cells (sign-flip: the same lesion in both niches)."
        ]),
        ("Semi-synthetic benchmarks", [
            "Pseudo-knockouts were built from control tumour cells in two designs. Independent clones reproduced the patch sizes and niches of real clones, rescaled to 50-800 cells. "
            "Lesion structure placed L ∈ {1, 3, 5, 10, 20} contiguous lesions (a random seed and its nearest unused control cells) with seeds at least 400 µm apart, 400 or 800 cells in total; "
            "a realism design resampled lesion sizes from the real knockouts' heavy-tailed lesion-size distribution.",
            "Effects were added either on the log scale, as δ·sd_g·u_g to all cells plus λ·δ·sd_g·v_g in one niche (δ = 0.5; u, v random unit-RMS directions, or v along hallmark hypoxia genes[@liberzon2015]), "
            "or on raw counts, by binomial thinning (decreases) or Poisson boosting (increases) of the counts of a fraction (5%, 20% or 100%) of genes in one niche, followed by renormalisation with the original all-gene totals and median target "
            "(re-normalising unmodified counts reproduced the analysis matrix to 7 × 10⁻¹⁶). Null settings used 200-500 replicates and effect settings 100-200."
        ]),
        ("Size-adjusted power", [
            "For every test and setting, the critical value was the smaller of 0.05 and the 5th percentile of the test's p-values in the matching null run (AMENDMENT-008b: using the percentile alone would give conservative discrete tests a critical value above 0.05), and size-adjusted power was the fraction of effect replicates with p at or below it. "
            "Raw and size-adjusted power are both reported; only size-adjusted power is interpreted when a test is not calibrated."
        ]),
        ("Design model", [
            "The design-effect model n_eff = L·m/(1 + (m − 1)ρ) was evaluated with ρ, the within-lesion intraclass correlation (one-way ICC(1)) of the tested score in control pseudo-lesions, estimated both for a gene-averaged projection and per gene. "
            "The model's one free parameter (an efficiency factor) was fitted to the simulated size-adjusted power of the adopted test by least squares. Because it did not reproduce the simulated curve, SpatialPert's plan_screen_empirical() interpolates the simulated curves (log-scale and count-level scenarios) in log lesion number; plan_screen() keeps the formula for comparison."
        ]),
        ("Real-data analyses", [
            "Real knockouts were analysed with every test above on the benchmark caller's niches (Fig. 6), and, as sensitivity analyses, at 20, 30 and 50 µm patch linkage, with gene-standardised residuals, with at least two guide UMIs and with distance to the nearest host cell as a covariate. "
            f"Lesions per knockout in Perturb-map imaging were counted from {pm['lesions']:,} curated lesion polygons named by Pro-Code combination ({len(pm['mice'])} mice, {pm['lobes']} lobes)."
        ]),
        ("Lesion definition", [
            "A lesion is a connected component of the graph that joins knockout cells closer than 100 µm (single linkage). The threshold was registered before real-data testing and is about four cell diameters in this tissue, so that nearby fragments of one metastasis join while lesions seeded at least 400 µm apart (the benchmark design) stay separate. "
            "Patches (30 µm) are the same construction at a finer scale and split one lesion into several units. Because lesion counts depend on the threshold, the real-data analysis reports lesions per niche and the patch result at 20, 30 and 50 µm. "
            "Lesions are a proxy for clones: two clones that touch merge, and a clone cut by the section plane can split. Neither error can be detected without lineage barcodes, so lesion-level results are conditional on the linkage rule."
        ]),
        ("Circular and leave-out niche labels", [
            "In real screens knockout cells are excluded when niches are defined, but a knockout cell's niche still depends on its control neighbours. Pseudo-knockouts drawn from control cells were originally labelled with niches that they themselves had helped to define (circular labels). "
            "This makes each pseudo-lesion more likely to sit inside one niche than a real lesion would, and it inflates lesion-level tests (Fig. 3c). Leave-out labels remove the pseudo-knockout cells from the reference before computing neighbourhood features and assign the result to the fixed centroids; they changed the niche of 33-46% of pseudo-knockout cells. "
            "All stage-1b results use leave-out labels unless marked circular. External callers had to be run once on all control cells, so their comparison (Fig. 5b) uses circular labels for every caller, including k-means, and compares callers like for like."
        ]),
        ("Test adoption rule", [
            "Before stage 1b was run, AMENDMENT-008 fixed the rule for choosing the recommended test: among tests whose type-I error was at most 0.06 at every lesion number where they were applicable (300 null draws per level), adopt the one with the highest mean size-adjusted power over 3-20 lesions. "
            "AMENDMENT-008b corrected the size-adjusted critical value (see above) after the first analysis showed that conservative discrete tests gained power from a critical value above 0.05. Under the literal rule the within-lesion sign-flip test would have been adopted; under the corrected rule lesion-unit permutation is adopted. Both are reported. "
            "The mixed model was run on 40 null and 20 effect draws per lesion number because each fit took about 26 s (registered deviation); its type-I error was far above the threshold on this subsample, so subsampling cannot change the adoption."
        ]),
        ("Second tissue", [
            "The mouse brain subset (Zenodo 11209134) provides transcript coordinates only. Transcripts were binned into 10 µm squares, giving 81,990 bins and 50 genes; bins are pseudo-cells and are smaller and sparser than segmented cells. "
            "Niches were k-means (K = 4) on 20 principal components with the same neighbourhood kernel, and the stage-1b design (300 null and 150 effect draws per lesion number, 400 bins per pseudo-knockout, leave-out labels) was repeated with cell-level, patch, lesion, sign-flip and wild cluster bootstrap tests (299 draws). "
            "Brain niches are large and spatially smooth (one niche held 61,694 of 81,990 bins), so lesions seldom straddle niches; the sign-flip test was applicable in at most 0.7% of draws and the wild cluster bootstrap in at most 31%."
        ]),
        ("Niche-caller comparison", [
            "BANKSY, SpaGCN and CellCharter were installed in isolated virtual environments so that their dependency pins did not change the analysis environment. Each was run on the 21,773 control tumour cells, exported as AnnData with the same 5,001 genes and coordinates, and labels were joined back by cell identifier. "
            "For each caller the lesion benchmark was repeated at 3, 5, 10 and 20 lesions (200 null and 100 effect draws each) with lesion permutation, the wild cluster bootstrap (199 draws) and the sign-flip test. "
            "We also report the fraction of cells in the largest niche, because a caller that puts most cells in one niche makes most pseudo-knockouts untestable rather than miscalibrated."
        ]),
        ("Count-level realism", [
            "To test whether conclusions from log-scale injection carry over to count data, effects were applied to raw counts. For a target niche and a random fraction s of genes (s = 0.05, 0.2 or 1), each selected gene received a log fold change δ·sd_g·(±1)/√s, clipped to ±3, so that the total signal is comparable across sparsity levels (δ = 0.25, 0.5 or 1). "
            "Decreases were applied by binomial thinning and increases by adding Poisson counts, only in cells of the target niche. Lesion sizes were drawn from the 100 µm lesion sizes of the four real knockouts and rescaled to 400 cells, reproducing their heavy tail; lesions were placed with the geometry-matched sampler. "
            "Each cell was renormalised with its original all-gene total and the analysis median target, and residuals were recomputed from the control model. Null runs used 200 draws and effect settings 100 draws per lesion number (5, 10 and 20)."
        ]),
        ("Within-lesion correlation", [
            "For lesion sizes of 20, 40, 80, 133 and 400 cells, control pseudo-lesions were grown and the within-lesion intraclass correlation ICC(1) was estimated with the unbalanced one-way ANOVA estimator, both for a gene-averaged projection of standardised residuals along a random direction (the direction along which effects were injected) and for each gene separately (median over genes reported)."
        ]),
        ("Perturb-map lesion counts", [
            "The S-BIAD267 deposit contains, for each imaged lung lobe, curated QuPath annotations in which each lesion polygon is named by the Pro-Code epitope combination detected in it. We parsed the lymphoid-panel annotation files and counted lesions per code, per mouse and per lobe. "
            "We report counts per code, because the code-to-knockout mapping is not needed to estimate how many lesions a knockout receives. Only annotations, not the 15.7 GB image-segmentation tables, were downloaded."
        ]),
        ("Positive control", [
            "The effect estimator was checked on Perturb-map Visium data (GSE193460), where Tgfbr2 knockout lesions have a published, strong microenvironment phenotype[@dhainaut2022]. The Tgfbr2 versus unassigned-tumour residual profile was computed in each section and agreement between sections was compared with a null from shuffled lesion labels. "
            "This control shows that residual effects capture a real knockout effect; it does not validate niche dependence, for which no ground-truth spatial knockout is available."
        ]),
        ("Statistics", [
            "Rejection rates are reported with 95% Wilson intervals. p-values from B permutations are (1 + #null ≥ observed)/(1 + B), and values at the floor are written as p < 1/(1 + B). "
            "The minimum achievable p-value per knockout is the larger of this floor and 1/(number of distinct labellings). Multiple knockouts were combined by Benjamini-Hochberg q-values and by Fisher and Cauchy combinations."
        ]),
        ("Software and runtime", [
            f"Analyses used Python 3.12 with numpy, scipy, pandas, scikit-learn and statsmodels; the code is gated by ruff, strict mypy and pytest on every commit. On a laptop ({hw['logical_cpus']} logical CPUs, {hw['ram_gb']} GB RAM, {hw['platform']}), "
            f"loading and normalising the full {352913:,}-cell section took {rt['load_and_normalise_full_section']['seconds']:.0f} s, fitting the niche model {rt['niche_model_fit']['seconds']:.0f} s, and testing four knockouts with 1,000 permutations "
            f"{rt['location_dependence_patch_1000perm_4KO']['seconds']:.1f} s (patch units) or {rt['location_dependence_lesion_1000perm_4KO']['seconds']:.1f} s (lesion units). "
            f"A 999-draw wild cluster bootstrap of one 10-lesion, 400-cell knockout took {rt['wild_cluster_bootstrap_999_pseudoKO_10lesions_400cells']['seconds']:.0f} s and a 2,000-flip sign-flip test {rt['signflip_2000_pseudoKO_10lesions_400cells']['seconds'] * 1000:.0f} ms; peak memory was {rt['process_peak_rss_gb']} GB."
        ]),
        ("Pre-registration and amendments", [
            "The analysis plan was frozen before any knockout expression outcome was examined (git tag protocol-v1.0). Eight dated amendments record every later analysis, each registered before it was run, including the pre-set rules for adopting a test. "
            "Deviations (mixed-model subsampling for runtime; SpaGCN subsampling and initialisation) are listed in the protocol. Supplementary Note 1 gives the history and the original hypotheses."
        ]),
        ("Reproducibility", [
            "Every number in the paper is read from report files written by scripts. An independent recheck recomputes reported rates from per-replicate files and raw data through separate code paths."
        ]),
    ]
