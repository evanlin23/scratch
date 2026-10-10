## Prior art and novelty check (searched 2026-10-10)

What I searched: web searches for the Togkousidis/Stamatakis/Gascuel paper, for "IQ-TREE early stopping / stopping rule /
nstop / KH test" (2025-26), for follow-ups citing the paper, for the IQ-TREE 3 paper, and for RAxML-NG 2.0's release notes;
plus a direct read of the IQ-TREE 3.1.4 source (`utils/stoprule.cpp`, `tree/iqtree.cpp`, `utils/tools.cpp`) and
`raxml-ng --help` of the installed RAxML-NG 2.0.3.

| item | what it says / does | consequence for this idea |
|---|---|---|
| Togkousidis, Stamatakis & Gascuel, "Accelerating Maximum Likelihood Phylogenetic Inference via Early Stopping to Evade (Over-)optimization", *Syst Biol* 74(6):1020 (2025), doi:10.1093/sysbio/syaf043 (preprint "Much Ado About Nothing", bioRxiv 2024.07.04.602058) | KH test (fast normal-approx. version) between the best tree before and after each SPR round of a *simplified* RAxML-NG search (sRAxML-NG); KH-mult adds a Bonferroni correction. 300 TreeBASE MSAs. Mean speed-up of KH-mult vs RAxML-NG 1.2: 5x (DNA), 3.9x (AA) for a single parsimony start; plausibility (AU) of at least one of 10 trees 98% DNA, ~92-94% AA vs 100% for RAxML-NG 1.2. Says the criteria "can also be seamlessly integrated into other phylogenetic inference tools that use numerical convergence thresholds such as PhyML and FastTree." | Source of the idea. Note that the paper names PhyML/FastTree (hill-climbers with numeric thresholds), *not* IQ-TREE, whose stopping rule is a stochastic-search patience counter. |
| RAxML-NG 2.0 (Zenodo release 2.0.0, Mar 2026; installed 2.0.3) | ships the rule: `--stop-rule sn-rell|sn-normal|kh|kh-mult` and `--fast` = `--search --tree pars{1} --opt-topology simplified --stop-rule kh-mult` (release notes: "about 50x speedup over v1.2") | The KH stop is now a *standard baseline* in the strongest competitor. |
| IQ-TREE 3.1.4 source (git tag v3.1.4, Sep 2026) | `StopRule::meetStopCondition`: `SC_UNSUCCESS_ITERATION` = stop when `iteration > last_improved + nstop` (default 100); an improvement is *any* new-best topology, however small the lnL gain (`addTreeToCandidateSet`). Also legacy `SC_WEIBULL` (`-sr`, IQPNNI's Vinh & von Haeseler 2004 rule), `SC_REAL_TIME`, `SC_BOOTSTRAP_CORRELATION` (UFBoot). `--fast` = 2 initial trees + fixed 2 iterations. No statistical test of improvement anywhere. | No KH/statistical early stop in IQ-TREE: the specific idea is open. |
| IQ-TREE 3 paper (Wong et al., MBE 2026, msag117; EcoEvoRxiv preprint 2025) | new features: mixture models, concordance factors, dating, AliSim; nothing on stopping. | open |
| Nguyen et al. 2015 (IQ-TREE 1) | introduced the 100-unsuccessful-iterations rule, replacing IQPNNI's statistical rule; mentions users may run longer. | the old *statistical* rule (Vinh & von Haeseler 2004) predicts *more* iterations; it is about not stopping too early, not about significance of gains. |
| "A Systematic Investigation of Overfitting in ML Phylogenetic Inference" (bioRxiv 2025.10.07.680876, v2 Jan 2026; same group) | compares RAxML-NG, IQ-TREE, FastTree and RAxML-NG ES; reports topology over-fitting is rare (<1% of datasets, per the KIT thesis record). Full text not retrievable (HTTP 429). | supports "late lnL gains are mostly noise-level" but does not implement anything in IQ-TREE. |
| GitHub/web search for an IQ-TREE issue/PR on early stopping | nothing found | open (absence of evidence; the IQ-TREE team may have it unreleased). |

Verdict on novelty: "KH early stopping inside IQ-TREE's perturbation loop" has not been published or released as far as I
can find. But it is an obvious port of a 2025 method by its own authors' suggestion, and the same group's RAxML-NG 2.0
already ships it, so the novelty is incremental; a project needs a result beyond "the port works".
