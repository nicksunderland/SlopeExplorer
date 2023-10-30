
# composite 1
foo = data.table::fread("/Users/xx20081/Documents/local_data/results/collider_bias/harmonised_effects.composite_1.autosomes.gz")
foo <- foo[PVAL.incidence<0.1, ]
cols <- c("SNP","BETA.incidence","PVAL.incidence","BETA.prognosis","PVAL.prognosis")
new_names <- c("SNP","i.BETA","i.P","p.BETA","p.P")
data.table::setnames(foo,cols, new_names)
foo[, names(foo)[!names(foo) %in% new_names] := NULL]
data.table::fwrite(foo, "/Users/xx20081/git/SlopeExplorer/inst/extdata/composite_1_autosomes.gz")


