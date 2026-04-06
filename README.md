# Social Decision-Making in Parkinson's
Software scripts and stan models for associated with the "*Prosocial evidence accumulation in Parkinson's Disease*" manuscript (not published) from Hubschmid F., Sellitto M., Bahners T.B., Pirenne T., Schnitzler A., van Wingerden M. & Florin E.
Fully anonymized task data will be made available upon the time of publication on a dedicated OSF repository.
### Scripts
Contains necessary R code to reproduce results repported in the main text.
### Stan models 
Contains the stan code for the the hyperbolic discount drift diffusion model (DDM) and the empathy weighted hyperbolic discount DDM. Given there are numerous valided hierarchical implementations of the DDM both in STAN (Ahn et al,. 2017) and Python (Wiecki,. 2013) we do not here provide code for the Baseline HDDM model.
### Additional disclaimers:
- Some code snippets were adapted from a previous project (Hubschmid et al., 2025).
- Models that map a discount function to a DDM are inspired from previous work of Peters & D'Esposito (2020).
- Stan models are structured as the hierarchical Bayesian implementation of the hBayesDM package (Ahn et al,. 2017).
- Parameter priors for the DDM come from Matzke & Wagenmakers (2009).
### References 
- Ahn, W.-Y., Haines, N., & Zhang, L. (2017). Revealing Neurocomputational Mechanisms of Reinforcement Learning and Decision-Making With the hBayesDM Package. Computational Psychiatry (Cambridge, Mass.), 1, 24–57. https://doi.org/10.1162/CPSY_a_00002
- Hubschmid, F., Desch, S., Florin, E., & Becker, S. (2025). Reward-induced endogenous pain inhibition scales with action-outcome certainty in humans (p. 2025.10.31.685818). bioRxiv. https://doi.org/10.1101/2025.10.31.685818
- Matzke, D., & Wagenmakers, E.-J. (2009). Psychological interpretation of the ex-Gaussian and shifted Wald parameters: A diffusion model analysis. Psychonomic Bulletin & Review, 16(5), 798–817. https://doi.org/10.3758/PBR.16.5.798
- Peters, J., & D’Esposito, M. (2020). The drift diffusion model as the choice rule in inter-temporal and risky choice: A case study in medial orbitofrontal cortex lesion patients and controls. PLOS Computational Biology, 16(4), e1007615. https://doi.org/10.1371/journal.pcbi.1007615
