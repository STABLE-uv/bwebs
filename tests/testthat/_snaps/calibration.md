# the bundled calibration has the documented structure

    Code
      print(cal)
    Output
      <bwebs_calibration> 13 validation webs, 14497 links inferred, 36% observed
        Bayesian: 4 chains x 4000 draws; max R-hat 1.001, min ESS 4461
        web effect: sd = 1.58 on the logit scale; precision per web from 0.03 to 0.84
        classes ordered: a link supported by a closer relative is never less probable
           class p_average p_new logit   sd webs predicted observed
           exact     0.489 0.493 -0.05 0.41   12       549      181
           genus     0.483 0.488 -0.07 0.41   13      1725      806
       subfamily     0.318 0.367 -0.79 0.41    9      1020      134
          family     0.307 0.360 -0.84 0.41   13      4445     1410
           order     0.249 0.313 -1.15 0.41   13      5138     1564
           class     0.248 0.312 -1.15 0.41    9       772      360
         generic     0.712 0.656  0.94 0.42    3       848      716
        pairs not inferred (taxonomic decay): 15979, with 3535 recorded links; web effect sd = 0.84
          probability in a new system from 0.065 to 0.472 (in a web of 50 nodes)
          max R-hat 1.004, min ESS 4983
          consumer resolved at (logit): family -1.55, order -1.58, class -1.59, phylum -1.93, any -3.01
          resource related at (logit, added): family 0.15, order 0.11, class 0.04, phylum 0.03, any 0.00
          record of a generic group (logit, added): 1.31
          size of the web, per unit of log(nodes / 50) (logit, added): 0.03
        studies left out together: 13 + 50
        classes without validation data: beyond

