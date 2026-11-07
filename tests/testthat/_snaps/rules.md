# print method is stable

    Code
      print(matching_rule())
    Output
      <bwebs_rule>
        mode: climb; start at family, climb up to class; min_partners = 5
        predators climb as resources: no (predator if > 0.5 animal prey)
        registries: all, at every level
        downward inheritance: subfamily, family, order (roles: resource)
        generic resource groups: 5; exceptions: 1; manual-diet nodes: 0

---

    Code
      print(matching_rule(primary_registry = "gray2015"))
    Output
      <bwebs_rule>
        mode: climb; start at family, climb up to class; min_partners = 5
        predators climb as resources: no (predator if > 0.5 animal prey)
        primary registry: gray2015; other registries up to family
        downward inheritance: subfamily, family, order (roles: resource)
        generic resource groups: 5; exceptions: 1; manual-diet nodes: 0

