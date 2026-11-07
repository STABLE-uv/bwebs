# print method is stable

    Code
      print(prepare_nodes(tiny_taxa(), rule = tiny_rule()))
    Output
      <bwebs_nodes> 9 nodes (6 consumers, 2 producers, 1 detritus)
        rules as consumer: R2 = 6, R3 = 0, R5 = 0, NT = 0, basal = 3, manual = 0
        rules as resource: R2 = 8, R3 = 0, R5 = 0, NT = 1
       node_id          rank     role resource_base consumer_base consumer_max
          cocc       species producer        family          none         none
          navi         genus producer        family          none         none
          cpom non_taxonomic detritus         exact          none         none
          baet       species consumer        family        family        class
          ecdy         genus consumer        family        family        class
          perl       species consumer        family        family        class
          simu         genus consumer        family        family        class
          gerr         genus consumer        family        family        class
          asel         genus consumer        family        family        class

