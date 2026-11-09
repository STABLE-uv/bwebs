# print methods are stable

    Code
      print(flux_parameters())
    Output
      <bwebs_flux_parameters> brown_2004
        metabolic rate per g of wet mass: exp(c + b log(M) - E / (k T)) W/g; M in g, T in K
        activation energy E = 0.63 eV (Brown et al. (2004): the value with which the rates were corrected for temperature before the intercepts were fitted (Fig. 1B))
        unicellular           c = 15.85, b = -0.29
        invertebrate          c = 17.17, b = -0.29
        fish                  c = 18.47, b = -0.29
        ectotherm_vertebrate  c = 18.18, b = -0.29
        endotherm             c = 19.50, b = -0.29
        assimilation efficiencies: animal 0.906, producer 0.77, detritus 0.158

---

    Code
      print(calc_fluxes(flux_draws(), flux_populations(), temperature = 20, engine = "internal"))
    Output
      <bwebs_fluxes> 1 sample(s) x 2 realization(s) of sampled metawebs
        biomasses: the same in every realization; parameters: brown_2004; engine: internal
        status: ok = 2
        total flux: median 0.00865 W (0.00522 to 0.0121) per unit of the sampling basis
        consumers left out without a path from a basal node: 0.5 per web on average
        use summarize_indicators() for posterior summaries; values in $values

