# Social Decision-Making in Parkinson's
Software scripts and stan models for associated with the "*Static and dynamic decision noise abolish altruistic biases in Parkinson’s disease independent of deep brain stimulation*" manuscript (not published) from Hubschmid F., Sellitto M., Bahners H.B., Pirenne T., Schnitzler A., van Wingerden M. & Florin E.
Fully anonymized task data will be made available upon the time of publication on a dedicated OSF repository.
### Scripts
Contains necessary R code to reproduce results reported in the main text.
### Stan models 
Contains the stan softwares for the social discount models used in the main text.
### Additional disclaimers:
- Some code snippets were adapted from previous projects (Hubschmid et al., 2025; Desch et al., 2023).
- Models that map a discount function to a DDM are inspired from previous work of Peters & D'Esposito (2020).
- Dynamic noise estimation through HMM corresponds to the method of Li & al., (2024).
- Stan models are structured as the hierarchical Bayesian implementation of the hBayesDM package (Ahn et al,. 2017).
- Parameter priors for the DDM come from Matzke & Wagenmakers (2009).
### References 
- Ahn, W.-Y., Haines, N., & Zhang, L. (2017). Revealing Neurocomputational Mechanisms of Reinforcement Learning and Decision-Making With the hBayesDM Package. Computational Psychiatry (Cambridge, Mass.), 1, 24–57. https://doi.org/10.1162/CPSY_a_00002
- Desch, S., Schweinhardt, P., Seymour, B., Flor, H., & Becker, S. (2023). Evidence for dopaminergic involvement in endogenous modulation of pain relief. eLife, 12, e81436. https://doi.org/10.7554/eLife.81436
- Hubschmid, F., Desch, S., Florin, E., & Becker, S. (2026). Reward-induced endogenous pain inhibition scales with action-outcome certainty in humans. Pain, 167(9), 2026–2037. https://doi.org/10.1097/j.pain.0000000000004058
- Li, J.-J., Shi, C., Li, L., & Collins, A. G. E. (2024). Dynamic noise estimation: A generalized method for modeling noise fluctuations in decision-making. Journal of Mathematical Psychology, 119, 102842. https://doi.org/10.1016/j.jmp.2024.102842
- Matzke, D., & Wagenmakers, E.-J. (2009). Psychological interpretation of the ex-Gaussian and shifted Wald parameters: A diffusion model analysis. Psychonomic Bulletin & Review, 16(5), 798–817. https://doi.org/10.3758/PBR.16.5.798
- Peters, J., & D’Esposito, M. (2020). The drift diffusion model as the choice rule in inter-temporal and risky choice: A case study in medial orbitofrontal cortex lesion patients and controls. PLOS Computational Biology, 16(4), e1007615. https://doi.org/10.1371/journal.pcbi.1007615




