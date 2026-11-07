# the bundled calibration has the documented structure

    Code
      print(cal)
    Output
      <bwebs_calibration> 13 validation webs, 14497 links inferred, 36% observed
        Bayesian: 4 chains x 4000 draws; max R-hat 1.001, min ESS 2982
        web effect: sd = 1.58 on the logit scale; precision per web from 0.03 to 0.84
           class p_average p_new logit   sd webs predicted observed
           exact     0.441 0.458 -0.24 0.37   12       549      181
           genus     0.542 0.529  0.17 0.36   13      1725      806
       subfamily     0.317 0.367 -0.79 0.37    9      1020      134
          family     0.347 0.389 -0.65 0.36   13      4445     1410
           order     0.246 0.311 -1.15 0.36   13      5138     1564
           class     0.514 0.509  0.06 0.36    9       772      360
         generic     0.735 0.673  1.05 0.37    3       848      716
        studies left out together: 13 + 50
        classes without validation data: beyond

