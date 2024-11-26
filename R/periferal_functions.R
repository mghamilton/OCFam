
#The script includes the following non-base R functions:
#
# optiSel::opticont
# AGHmatrix::Amatrix



################################################################################
#function to run optiSel::opticont depending on opticont_method
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
  fit <- optiSel::opticont(opticont_method, cand, con,
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
  fam_K_matrix <- optiSelFam::fam_K_matrix_fun(fam_K_matrix[!is.na(fam_K_matrix[,1]),])
  fam_K_matrix <- as.matrix(fam_K_matrix$K_matrix_families)
  fam_K_matrix <- fam_K_matrix[colnames(fam_K_matrix) %in% cand_fams, colnames(fam_K_matrix) %in% cand_fams]
  return(fam_K_matrix)
}

################################################################################
#fam_K_matrix_fun
################################################################################

fam_K_matrix_fun <- function(family_dat) {
  #data
  # class: data.frame
  # fields:
  #   FAM
  #   SIRE
  #   FAM_SIRE
  #   DAM
  #   FAM_DAM

  #load required library - AGHmatrix
#  if("AGHmatrix" %in% installed.packages()[,"Package"] == FALSE) { install.packages("AGHmatrix" , repos='https://cran.csiro.au/')}  #Install nadiv package if not already installed
#  library(AGHmatrix)

  #sires
  sires <- matrix(unique(family_dat$SIRE), ncol = 1)
  colnames(sires)[1] <- "SIRE"
  sires <- merge(sires, family_dat[,c("SIRE","FAM_SIRE")], by = "SIRE", all.x = TRUE)
  colnames(sires) <- c("Indiv_id","FAM")
  sires <- merge(sires, family_dat[,c("FAM","SIRE","DAM")], by = "FAM", all.x = TRUE)
  sires <- subset(sires, select=-c(FAM)) #remove FAM column

  #dams
  dams <- matrix(unique(family_dat$DAM), ncol = 1)
  colnames(dams)[1] <- "DAM"
  dams <- merge(dams, family_dat[,c("DAM","FAM_DAM")], by = "DAM", all.x = TRUE)
  colnames(dams) <- c("Indiv_id","FAM")
  dams <- merge(dams, family_dat[,c("FAM","SIRE","DAM")], by = "FAM", all.x = TRUE)
  dams <- subset(dams, select=-c(FAM)) #remove FAM column

  parents <- rbind(sires,dams)
  parents <- unique(parents)
  parents <- parents[order(parents[,"Indiv_id"], decreasing = FALSE), ]

  families            <- unique(family_dat)
  tmp1 <- families
  tmp2 <- families
  tmp1$FAM <- paste0("1_", tmp1$FAM) #needs a unique id
  families <- rbind(tmp1, tmp2)

  fam_pedigree        <- rbind(as.matrix(parents),as.matrix(families[,c("FAM","SIRE","DAM")]))

  #Generate K Matrix
  # K_matrix_families      <- nadiv::makeA(fam_pedigree[,1:3])/2
  tmp <- fam_pedigree
  tmp[is.na(tmp)] <- 0
  K_matrix_families      <- AGHmatrix::Amatrix(tmp[,1:3], ploidy=2)/2
  rm(tmp)

  family_inbreeding <-  round((diag(as.matrix(K_matrix_families)) * 2) - 1, 10)
  family_inbreeding <-  family_inbreeding[(nrow(parents)+1):(nrow(parents)+nrow(families)/2)]

  K_matrix_families <- K_matrix_families[(nrow(parents)+1):(nrow(parents)+nrow(families)/2), (nrow(parents)+nrow(families)/2+1):nrow(K_matrix_families)]

  rownames(K_matrix_families) <- colnames(K_matrix_families)
  names(family_inbreeding) <- colnames(K_matrix_families)

  return(list(K_matrix_families = K_matrix_families,
              family_inbreeding = family_inbreeding))
}
