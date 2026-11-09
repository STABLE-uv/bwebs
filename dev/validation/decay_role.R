# Pairs not inferred: does the role of the resource, or the size of the web,
# predict the missed links better than the taxonomic distance alone?
# (version 0.3.0; not in the package)
#
# For every pair without an inferred link of the 12 validation webs (the
# default rule, each web rebuilt without the studies of its system), its
# evidence for the decay, the role of the resource (producer, detritus or
# animal) and whether the rule infers links of the consumer to resources of
# that kind (basal or animal; or none at all). Logistic models with a random
# effect of the web (lme4::glmer), scored as score_rule() does: each system
# left out in turn, probability in a new system (averaged over the web
# effect), squared error over its pairs not inferred. The size of the web is
# the log of its number of nodes (centred at 50). Gray registry, and Gray +
# GloBI (every registry at every level).
#
# Run from the root of the package (about 10 minutes; needs lme4):
#   Rscript dev/validation/decay_role.R <globi_registry.rds>
pkgload::load_all(".", quiet = TRUE)
g <- readRDS(commandArgs(trailingOnly = TRUE)[1])
out_dir <- file.path("dev", "validation")
gray <- gray2015[gray2015$study_id != "12", ]
both <- combine_registries(gray, g$globi_test)
rule <- matching_rule()
src <- setdiff(gray2015_sources$study_id[gray2015_sources$validation], "12")
key <- function(r, c) paste(r, c, sep = "\u001f")

pairs_of <- function(reg, excl) {
  out <- list()
  for (w in src) {
    vn <- registry_nodes(reg, w, rule)
    nodes <- vn$nodes
    left <- if (w %in% c("13", "50")) c("13", "50") else w
    extra <- if (is.null(excl)) character(0) else excl$study_id[excl$web == w]
    train <- reg[!reg$study_id %in% c(left, extra), ]
    mw <- metaweb_core(nodes, train, check_manual_links(NULL, nodes$node_id), prefilter = "both")
    studied <- nodes$node_id[vn$studied]
    e <- mw$edges[!mw$edges$cannibal & mw$edges$consumer %in% studied, ]
    role <- stats::setNames(nodes$role, nodes$node_id)
    e_basal <- role[e$resource] != "consumer"
    eats_basal <- tapply(e_basal, e$consumer, any)
    eats_animal <- tapply(!e_basal, e$consumer, any)
    p <- expand.grid(resource = nodes$node_id, consumer = studied, stringsAsFactors = FALSE)
    p <- p[p$resource != p$consumer &
             !key(p$resource, p$consumer) %in% key(e$resource, e$consumer), ]
    p$y <- key(p$resource, p$consumer) %in% key(vn$records$resource, vn$records$consumer)
    p$res_role <- unname(role[p$resource])
    eb <- unname(eats_basal[p$consumer])
    ea <- unname(eats_animal[p$consumer])
    eb[is.na(eb)] <- FALSE
    ea[is.na(ea)] <- FALSE
    basal <- p$res_role != "consumer"
    p$diet <- ifelse(!eb & !ea, "no_links",
                     ifelse(ifelse(basal, eb, ea), "kind_eaten", "kind_not_eaten"))
    p <- cbind(p, decay_evidence(decay_profiles(nodes, train, rule, consumers = studied),
                                 p$resource, p$consumer))
    p$web <- w
    out[[w]] <- p
  }
  do.call(rbind, out)
}

tx <- c("cl_order", "cl_class", "cl_phylum", "cl_any", "rg_family", "rg_order", "rg_class",
        "rg_phylum", "generic")
forms <- list(taxonomy = tx, taxonomy_web_size = c(tx, "log_nodes"),
              role = c(tx, "producer", "detritus"),
              role_diet = c(tx, "producer", "detritus", "animal_not_eaten", "basal_not_eaten",
                            "no_links"),
              role_diet_only = c("producer", "detritus", "animal_not_eaten", "basal_not_eaten",
                                 "no_links"))
tables <- list()
summaries <- list()
for (nm in c("gray", "both_all")) {
  p <- if (nm == "gray") pairs_of(gray, NULL) else pairs_of(both, g$exclusions)
  kind <- ifelse(p$res_role == "consumer", "animal", "basal")
  d <- data.frame(web = p$web, sys = ifelse(p$web %in% c("13", "50"), "13+50", p$web))
  for (lv in decay_consumer_levels[-1L]) d[[paste0("cl_", lv)]] <- as.numeric(p$consumer_level == lv)
  rg <- decay_resource_groups[decay_resource_group(match(p$resource_level, decay_levels))]
  for (lv in decay_resource_groups[-5L]) d[[paste0("rg_", lv)]] <- as.numeric(rg == lv)
  d$generic <- as.numeric(p$generic)
  d$producer <- as.numeric(p$res_role == "producer")
  d$detritus <- as.numeric(p$res_role == "detritus")
  d$animal_not_eaten <- as.numeric(kind == "animal" & p$diet == "kind_not_eaten")
  d$basal_not_eaten <- as.numeric(kind == "basal" & p$diet == "kind_not_eaten")
  d$no_links <- as.numeric(p$diet == "no_links")
  n_nodes <- tapply(p$resource, p$web, function(x) length(unique(x)))
  d$log_nodes <- log(as.numeric(n_nodes[p$web])) - log(50)
  d$y <- as.numeric(p$y)
  d$n <- 1
  a <- stats::aggregate(cbind(y, n) ~ ., data = d, FUN = sum)
  rates <- stats::aggregate(cbind(pairs = 1, missed = y) ~ res_role + diet, data = p, FUN = sum)
  rates$rate <- rates$missed / rates$pairs
  summaries[[nm]] <- rates
  tables[[nm]] <- do.call(rbind, lapply(names(forms), function(f) {
    rhs <- paste(forms[[f]], collapse = " + ")
    fo <- stats::as.formula(paste("cbind(y, n - y) ~", rhs, "+ (1 | web)"))
    fit <- suppressMessages(suppressWarnings(lme4::glmer(fo, data = a, family = stats::binomial)))
    sq <- 0
    for (s in unique(a$sys)) {
      m <- suppressMessages(suppressWarnings(lme4::glmer(fo, data = a[a$sys != s, ],
                                                          family = stats::binomial)))
      te <- a[a$sys == s, ]
      cf <- lme4::fixef(m)
      X <- stats::model.matrix(stats::as.formula(paste("~", rhs)), te)
      eta <- drop(X[, names(cf), drop = FALSE] %*% cf)
      ph <- p_new_system(eta, sqrt(unlist(lme4::VarCorr(m))))
      sq <- sq + sum(te$y * (1 - ph)^2 + (te$n - te$y) * ph^2)
    }
    data.frame(source = nm, model = f, aic = stats::AIC(fit),
               sigma = sqrt(unlist(lme4::VarCorr(fit))),
               log_nodes = if ("log_nodes" %in% names(lme4::fixef(fit))) {
                 lme4::fixef(fit)[["log_nodes"]]
               } else {
                 NA_real_
               },
               cv_squared_error = sq, squared_error_p0 = sum(a$y), stringsAsFactors = FALSE)
  }))
}
tables <- do.call(rbind, tables)
rownames(tables) <- NULL
print(summaries)
print(tables, digits = 5, row.names = FALSE)
saveRDS(list(rates = summaries, models = tables), file.path(out_dir, "decay_role.rds"))
