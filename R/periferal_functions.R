
#The script includes the following non-base R functions:
#
#  fam_K_matrix_fun
# opticont




################################################################################
#function to run opticont depending on opticont_method
################################################################################

run_OC_max_EBV <- function(cand, kinship_constraint, ub, lb, opticont_method) {
  if(opticont_method == "max.EBV") {

    if(nrow(cand$classes[cand$classes$Group == "cohort"]) == 1 |
       nrow(cand$classes[cand$classes$Group != "cohort"]) == 2) { #overlapping = FALSE
      con <- list(ub.fPED = kinship_constraint,
                  ub = ub,
                  lb = lb)
    } else {

      con <- list(ub.fPED = kinship_constraint,
                  ub = ub,
                  lb = lb,
                  ub.oldest_age_r_vector = cand$classes[cand$classes$age == max(cand$classes$age)  ,"cont0"],
                  lb.oldest_age_r_vector_2 = cand$classes[cand$classes$age == max(cand$classes$age)  ,"cont0"]) #fix to r vector value of the oldest year class
    }

  }

  if(opticont_method == "min.fPED") {

    if(nrow(cand$classes[cand$classes$Group == "cohort"]) == 1 |
       nrow(cand$classes[cand$classes$Group != "cohort"]) == 2) { #overlapping = FALSE
      con <- list(ub = ub,
                  lb = lb)
    } else {

      con <- list(ub = ub,
                  lb = lb,
                  ub.oldest_age_r_vector = cand$classes[cand$classes$age == max(cand$classes$age)  ,"cont0"],
                  lb.oldest_age_r_vector_2 = cand$classes[cand$classes$age == max(cand$classes$age)  ,"cont0"]) #fix to the contribution of the oldest year class
    }

  }
  fit <- opticont(opticont_method, cand, con,
                  solver="cccp")
  return(fit)
}




################################################################################
#identify best individual from fams_to_retain
################################################################################

get_best_indiv <- function(fish, additional_fams_to_retain, candidate_parents) {
  tmp <- fish[order(fish$RANK, decreasing = TRUE),]
  tmp <- tmp[tmp$Fam %in% additional_fams_to_retain &
               tmp$Indiv %in% candidate_parents,]
  cand_parents_fixed <- tmp[!duplicated(tmp$Fam),]
  cand_parents_fixed <- cand_parents_fixed[order(cand_parents_fixed$Fam),"Indiv"]
  rm(tmp)
  return (cand_parents_fixed)
}






################################################################################
#get family K matrix
################################################################################

get_fam_K_matrix <- function(ped, cand_fams) {

  fam_K_matrix <-  unique(ped[,c("Fam", "Sire", "FAM_SIRE", "Dam", "FAM_DAM")])
  colnames(fam_K_matrix) <-  c("FAM", "SIRE", "FAM_SIRE", "DAM", "FAM_DAM")
  fam_K_matrix <- fam_K_matrix_fun(fam_K_matrix[!is.na(fam_K_matrix[,1]),])
  fam_K_matrix <- as.matrix(fam_K_matrix$K_matrix_families)
  fam_K_matrix <- fam_K_matrix[colnames(fam_K_matrix) %in% cand_fams, colnames(fam_K_matrix) %in% cand_fams]
  return(fam_K_matrix)
}
