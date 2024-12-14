#' OCFam
#'
#' @description
#' This function generates a mating list for a set of parents.
#' The mating list can be generated i) to minimise the average inbreeding coefficient (F) of families generated or ii) according to assortative mating principles.
#' Inputs include a list of parents and a 3-column pedigree file specifying the ancestry of these candidates.
#' @param ped is a data frame with the following columns (class in parentheses):
#' \itemize{
#'  \item{'INDIV' (character).}
#'  \item{'SIRE'  (character).}
#'  \item{'DAM'  (character).}
#'  \item{'FAM'  (character).}
#'  \item{'BORN'  (numeric).}
#'  \item{'EBV'  (numeric).}
#'  \item{'N_AS_PARENT_CURRENT'  (numeric).} #number of families contributed to in the current spawn round
#'  \item{'AVAIL_BROOD'  (logical).}
#' }
#' @param indiv_contbn is the ...... (numeric between 0 and 1)
#' @param kinship_constraint is the ...... (numeric between 0 and 1)
#' @param step_interval is the ...... (numeric between 0 and 1)???
#' @param gene_flow_vector is a vector ......  as.numeric(NA) if overlapping_gens=FALSE (numeric). For example, gene_flow_vector = c(0.2, 0.8, 0, 0) - Year 4, 3, 2, 1, oldest to youngest
#' @param min_prop_fams is the min_prop_parent_fams_to_retain_in_youngest_age_class (numeric between 0 and 1)
#' @param max_parents_per_fam is the ...... (integer)
#' @return 'fam_contbn' is a data frame containing details of contributions by family:
#' \itemize{
#'  \item{ }
#' }
#' @return 'parent_contbn' is a data frame containing details of contributions by individual:
#' \itemize{
#'  \item{}
#' }
#' @return 'fit_out' is a vector containing details of the constraints applied:
#' \itemize{
#'  \item{}
#' }
#' @examples
#' #Retrieve example data
#' candidate_parents <- OCFam::candidate_parents
#' ped <- OCFam::ped
#'
#' output <- OCFam(
#' candidate_parents = candidate_parents,
#' ped = ped,
#' indiv_contbn = (1/180),
#' kinship_constraint = 0.09,
#' step_interval = 0.1,
#' overlapping_gens = FALSE,
#' gene_flow_vector = NA,
#' min_prop_fams = 1)
#' head(output$fam_contbn)
#' head(output$parent_contbn)
#' output$fit_out
#' @import optiSel
#' @import AGHmatrix
#' @import dplyr
#' @export

#Functions that are not base R functions
# dplyr::left_join
# optiSel::candes
# optiSel::pedIBD
# optiSel::prePed
# OCFam::run_OC_max_EBV
# OCFam::get_best_indiv
# OCFam::get_fam_K_matrix

#Required packages: optiSel, dplyr, AGHmatrix

#TEST#######################################################
#  load("C:/Users/MHamilton/OneDrive - CGIAR/Desktop/OCFam/R/Example_CC.RData")
#  ped <- ped[,c("INDIV", "SIRE", "DAM", "FAM", "EBV")]
#  candidate_parents <- candidate_parents_selection
#  indiv_contbn <- indiv_contbn_selection # 1/(N_fams_selection*2 + sum(ped[ped$LINE == "Selection" & !is.na(ped$N_AS_PARENT_CURRENT), "N_AS_PARENT_CURRENT"])) #parents of existing families
#  kinship_constraint <- 0.014
#  min_prop_fams <- min_prop_parent_fams_to_retain_selection



#Function to optimise contributions at family level
OCFam  <- function(#candidate_parents,
  ped,
  indiv_contbn,
  kinship_constraint,
  step_interval = 0.1,
  overlapping_gens,
  gene_flow_vector,
  min_prop_fams,
  max_parents_per_fam) {


  #candidate_parents - output of get_lb_ub function
  # INDIV                     N_AS_PARENT_CURRENT lb          ub                                Exclude_max_parents_per_fam
  # G202302_752E925           NA           0           0.005555556                       FALSE
  # G202302_7535BE7           NA           0           0.005555556                       FALSE
  # G202302_7536818           NA           0           0.005555556                       FALSE
  # G202302_753946A           NA           0           0.005555556                       FALSE
  # G202302_753DC36           NA           0           0.005555556                       FALSE
  # G202302_753E4AF           NA           0           0.005555556                       FALSE

  #indiv_contbn - scalar
  # e.g. 0.005555556 = 1/(N_fams_tilv*2 + sum(ped[ped$LINE == "TILV" & !is.na(ped$N_AS_PARENT_CURRENT), "N_AS_PARENT_CURRENT"])) #parents of existing families

  #kinship_constraint
  #e.g. 0.09

  #ped - pedigree
  # INDIV           SIRE            DAM             FAM             BORN     EBV      RANK FAM_SIRE
  # G202302_75FB917 G202202_7AE17C5 G202202_7ADC209 G202302_ON_0178 17.88489 -0.0565  135  G202202_ON_0081
  # G202302_75D0BE8 G202202_7AE17C5 G202202_7ADC209 G202302_ON_0178 17.88489 -0.0556  136  G202202_ON_0081
  # G202302_75D03F6 G202202_7AE17C5 G202202_7ADC209 G202302_ON_0178 17.88489 -0.0549  137  G202202_ON_0081
  # G202302_75D29F6 G202202_7AE17C5 G202202_7ADC209 G202302_ON_0178 17.88489 -0.0545  138  G202202_ON_0081
  # G202302_75FA452 G202202_7AE17C5 G202202_7ADC209 G202302_ON_0178 17.88489 -0.0545  139  G202202_ON_0081
  # G202302_75CF143 G202202_7AE17C5 G202202_7ADC209 G202302_ON_0178 17.88489 -0.0497  140  G202202_ON_0081
  # FAM_DAM         LINE COHORT        AVAIL_BROOD Breed        N_AS_PARENT_CURRENT     AVAIL_OR_PAST_BROOD
  # G202202_ON_0007 TiLV    <NA>       FALSE       ON           NA               FALSE
  # G202202_ON_0007 TiLV GROWOUT       FALSE       ON           NA               FALSE
  # G202202_ON_0007 TiLV GROWOUT       FALSE       ON           NA               FALSE
  # G202202_ON_0007 TiLV GROWOUT       FALSE       ON           NA               FALSE
  # G202202_ON_0007 TiLV GROWOUT       FALSE       ON           NA               FALSE
  # G202202_ON_0007 TiLV GROWOUT       FALSE       ON           NA               FALSE

  #step_interval
  #e.g. 0.1

  #overlapping_gens
  #e.g. FALSE

  #gene_flow_vector
  #e.g. NA - not required if overlapping_gens=FALSE

  #min_prop_fams
  #e.g. 1

  # cand_parents_fixed_initial <- candidate_parents[!is.na(candidate_parents$lb), "INDIV"]

  #  if(length(cand_parents_fixed_initial) >= (1 / indiv_contbn)) {break("Fixed contributions of indivuals equal to or greater than possible number based on indiv_contbn")}
  set.seed(12345)
  tmp <- ped[order(runif(nrow(ped))),c("INDIV", "EBV")]
  tmp <- tmp[order(tmp$EBV ), ]
  tmp$RANK <- nrow(tmp):1
  ped <- dplyr::left_join(ped, tmp, by = c("INDIV", "EBV"))
  rm(tmp)

  ped$AVAIL_OR_PAST_BROOD <- ped$AVAIL_BROOD | (ped$N_AS_PARENT_CURRENT > 0 & !is.na(ped$N_AS_PARENT_CURRENT))

  candidate_parents <- get_lb_ub(ped = ped,
                                 indiv_contbn = indiv_contbn,
                                 max_parents_per_fam = max_parents_per_fam)

  if(nrow(candidate_parents) < 3) {break("Not enough candidate parents.  Must be 3 or more.")}

  ################################################################################
  #Compute necessary columns
  ################################################################################

  #ped$Generation <- determine_generations(ped)

  fams <- ped[,c("INDIV", "FAM")]
  colnames(fams) <- c("SIRE", "FAM_SIRE")
  ped <- dplyr::left_join(ped, fams, by = "SIRE")
  colnames(fams) <- c("DAM", "FAM_DAM")
  ped <- dplyr::left_join(ped, fams, by = "DAM")
  rm(fams)

  ped$AVAIL_BROOD <- ped$INDIV %in% candidate_parents[!candidate_parents$Exclude_max_parents_per_fam, "INDIV"]

  ped$Breed <- "Ignored"

  #ped <- dplyr::left_join(ped, candidate_parents[,c("INDIV", "N_AS_PARENT_CURRENT")], by = "INDIV")

  ################################################################################
  #Inputs for optiSel
  ################################################################################


  if(overlapping_gens == TRUE) {
    #contributions
    gen_interval <- sum(gene_flow_vector * seq(length(gene_flow_vector),1))
    r_vector     <- matrix(NA, nrow = 1, ncol = length(gene_flow_vector)) #not separated by sex
    r_vector     <- as.vector(r_vector)
    for(i in 1:length(gene_flow_vector)) {
      r_vector[i] = sum(gene_flow_vector[1:i]/gen_interval) #oldest age class to youngest age class
    }

    cont <- data.frame(age = length(r_vector):1,                #oldest age class to youngest age class
                       male = r_vector / 2,
                       female = r_vector / 2)
    # cont <- cont[rowSums(cont[,c("male", "female")]) > 0,]
    cont <- cont[order(cont$age),]                              #youngest age class to oldest age class
    rownames(cont) <- cont$age

    #Kinships

    candidate_parents <- merge(candidate_parents, ped[,c("INDIV", "BORN")], by = "INDIV") #get BORN

    #keep one individual per family from immature generations
    tmp <- ped[ped$BORN > max(candidate_parents$BORN),]
    tmp <- tmp[!duplicated(tmp$FAM),]
    tmp$Mature <- FALSE
    tmp$INDIV <-  tmp$FAM
    tmp$EBV <- 0
    ped <- ped[ped$BORN <= max(candidate_parents$BORN),]
    ped$Mature <- TRUE
    ped <- rbind(ped, tmp)
    rm(tmp)

    ped$isCandidate <- FALSE
    ped$isCandidate <- ped$INDIV %in% candidate_parents$INDIV

    keep = ped[ped$isCandidate | !ped$Mature,"INDIV"] #keep if isCandidate or is immature
    Pedig <- optiSel::prePed(Pedig = ped, keep = keep) #keep if isCandidate or is immature
    fPED <- optiSel::pedIBD(Pedig, keep.only = keep)

    colnames(Pedig)[colnames(Pedig) %in% c("Indiv", "Sire", "Dam", "Sex", "Breed")] <- c("INDIV", "SIRE", "DAM", "SEX", "BREED")

    #replace immature year classes with family K matrix as per Hamilton paper
    fam_K <- OCFam::get_fam_K_matrix(Pedig, Pedig[!Pedig$Mature,"INDIV"])
    fPED[rownames(fam_K),colnames(fam_K)] <- fam_K

    # Pedig <- Pedig[Pedig$INDIV %in% colnames(fPED),]

    Sy <- summary(Pedig)
    Pedig <- merge(Pedig, Sy[, c("INDIV", "equiGen")], on="INDIV")

    #Population means

    #  Pedig[Pedig$EBV == 0,"EBV"] <- runif(nrow(Pedig[Pedig$EBV == 0,]))
    rownames(Pedig) <- Pedig$INDIV
    Pedig$oldest_age_r_vector <- 0
    Pedig[Pedig$BORN == min(Pedig[Pedig$AVAIL_OR_PAST_BROOD, "BORN"]), "oldest_age_r_vector"] <- 1 #oldest YC to 1
    Pedig$oldest_age_r_vector_2 <- Pedig$oldest_age_r_vector

    cand <- optiSel::candes(phen=Pedig[Pedig$INDIV %in% rownames(fPED),], #c("INDIV",	"SIRE",	"DAM",	"Sex",	"Breed",	"BORN",		"EBV",	"isCandidate", "oldest_age_r_vector", "oldest_age_r_vector_2")],
                            fPED=fPED,
                            cont = cont)

    #set initial ub and lb
    ub <- setNames(candidate_parents$ub, candidate_parents$INDIV)
    lb <- setNames(candidate_parents$lb, candidate_parents$INDIV)

  }

  if(overlapping_gens == FALSE) {

    Pedig <- optiSel::prePed(Pedig = ped, keep = candidate_parents$INDIV)

    fPED <- optiSel::pedIBD(Pedig, keep.only = candidate_parents$INDIV)

    Sy <- summary(Pedig)
    Pedig <- merge(Pedig, Sy[, c("INDIV", "equiGen")], on="INDIV")

    cont <- data.frame(age = 1,
                       male = 0.5,
                       female = 0.5)
    cand <- optiSel::candes(phen=Pedig[Pedig$INDIV %in% candidate_parents$INDIV,],
                            fPED=fPED,
                            cont = cont)

    #set initial ub and lb
    ub <- setNames(candidate_parents$ub, candidate_parents$INDIV)
    ub <- ub[cand$phen$INDIV]

    lb <- setNames(candidate_parents$lb, candidate_parents$INDIV)
    lb <- lb[cand$phen$INDIV]
  }


  if(is.na(kinship_constraint)) {
    opticont_method = "min.fPED"
  } else {
    opticont_method = "max.EBV"
  }

  ################################################################################
  #1. generate candidate list with one individual for Min_N_fams_selection_to_retain family (least related families)
  ################################################################################

  #identify candidate families



  #modified this code to limit only families in most youngest year class




  cand_fams  <- unlist(unique(ped[ped$INDIV %in% candidate_parents$INDIV, "FAM"]))
  cand_fams_in_youngest_age_class <- unlist(unique(ped[ped$INDIV %in%
                                                         candidate_parents[candidate_parents$BORN == max(candidate_parents$BORN),  "INDIV"],
                                                       "FAM"]))
  # cand_fams_in_oldest_age_class <- cand_fams[!cand_fams %in% cand_fams_in_youngest_age_class]

  fixed_fams <- unlist(unique(ped[ped$INDIV %in% names(lb[lb > 0.1*indiv_contbn]), "FAM"])) #contribution > 0
  fixed_fams_in_youngest_age_class <- fixed_fams[fixed_fams %in% cand_fams_in_youngest_age_class]

  # fixed_fams <- unique(ped[ped$INDIV %in% ped[!is.na(ped$N_AS_PARENT_CURRENT),"INDIV"], "FAM"])

  fam_K_matrix <- OCFam::get_fam_K_matrix(ped = ped,
                                               cand_fams = cand_fams)
  print(fam_K_matrix)
  mean(fam_K_matrix)
  write.csv(fam_K_matrix, "fam_K_matrix.csv")

  fam_K_matrix_in_youngest_age_class <- fam_K_matrix[rownames(fam_K_matrix) %in% cand_fams_in_youngest_age_class,
                                                     colnames(fam_K_matrix) %in% cand_fams_in_youngest_age_class]

  min_N_parent_fams_to_retain_in_youngest_age_class <- ceiling(length(cand_fams_in_youngest_age_class) * min_prop_fams)
  rm(cand_fams)

  additional_fams_to_retain <- NULL

  if(length(fixed_fams_in_youngest_age_class) < min_N_parent_fams_to_retain_in_youngest_age_class) {
    while(nrow(fam_K_matrix_in_youngest_age_class) > min_N_parent_fams_to_retain_in_youngest_age_class) {
      tmp <- rowSums(fam_K_matrix_in_youngest_age_class)
      tmp <- tmp[!names(tmp) %in% fixed_fams_in_youngest_age_class]
      tmp <- max(tmp)

      tmp <- rownames(fam_K_matrix_in_youngest_age_class)[rowSums(fam_K_matrix_in_youngest_age_class) == tmp]
      tmp <- tmp[!tmp %in% fixed_fams_in_youngest_age_class]

      fam_K_matrix_in_youngest_age_class <- fam_K_matrix_in_youngest_age_class[rownames(fam_K_matrix_in_youngest_age_class) != tmp[1], rownames(fam_K_matrix_in_youngest_age_class) != tmp[1]]
      rm(tmp)
    }

    additional_fams_to_retain <- rownames(fam_K_matrix_in_youngest_age_class)[!rownames(fam_K_matrix_in_youngest_age_class) %in% fixed_fams]
  }
  rm(fam_K_matrix, fam_K_matrix_in_youngest_age_class,
     fixed_fams, fixed_fams_in_youngest_age_class,
     min_N_parent_fams_to_retain_in_youngest_age_class,
     cand_fams_in_youngest_age_class)







  additional_parents_fixed <- OCFam::get_best_indiv(fish = ped[ped$AVAIL_BROOD ,], #exclude animals previously used as parents
                                                         additional_fams_to_retain = additional_fams_to_retain,
                                                         candidate_parents = candidate_parents[!candidate_parents$Exclude_max_parents_per_fam,"INDIV"])
  ub[names(ub) %in% additional_parents_fixed] <- indiv_contbn
  lb[names(lb) %in% additional_parents_fixed] <- indiv_contbn - 1e-10
  rm(additional_parents_fixed, additional_fams_to_retain)

  ub_orig <- ub
  lb_orig <- lb
  ################################################################################
  #2. identify individuals that contribute more than indiv_contbn and add to list of cand_parents_fixed
  ################################################################################

  fit <- OCFam::run_OC_max_EBV(cand = cand, kinship_constraint = kinship_constraint,
                                    ub = ub, lb = lb, opticont_method = opticont_method)





  ###############################################
  #NEED TO ENSURE IN OVERLAPPING GENERATIONS THAT ONE GENERATION IS NOT OVER-REPRESTENTED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





  ub <- setNames(floor(round(fit$parent$oc/indiv_contbn,4)) * indiv_contbn, fit$parent$INDIV)
  ub[ub == 0] <- indiv_contbn
  ub <- ub[cand$phen$INDIV]
  ub <- ub[names(ub) %in%  names(ub_orig)]

  lb <- setNames(floor(round(fit$parent$oc/indiv_contbn,4)) * indiv_contbn - 1e-10 , fit$parent$INDIV)
  lb[lb < 0] <- 0
  lb <- lb[cand$phen$INDIV]
  lb <- lb[names(lb) %in%  names(lb_orig)]

  ################################################################################
  #3.  add to list of cand_parents_fixed the best individual from families highly represented in OC
  #      - add individuals consecutively with families highly represented in OC fixed first
  ################################################################################







  cand_parents_fixed_current_iteration <- NA #lb[lb != 0]
  cand_parents_fixed_past_iteration <- NULL
  max_sum_oc <- 1

  for(adj in seq((1-step_interval),0,-step_interval)) {

    print(paste("Count parents previous iteration =",length(c(cand_parents_fixed_current_iteration, cand_parents_fixed_past_iteration )), "of", 1 / indiv_contbn))
    print(paste("indiv_contbn = ",indiv_contbn))
    print(paste("Iteration adjustment = ",adj))
    print(paste("max_sum_oc = ",max_sum_oc))
    print(paste("indiv_contbn * adj = ",indiv_contbn * adj))


    #  if(cand_parents_fixed_past_iteration < 1 / indiv_contbn) {


    if(max_sum_oc > (indiv_contbn * adj)) {

      if((sum(lb)) < 1) {


        cand_parents_fixed_past_iteration <- names(lb[lb!=0])

        cand_parents_not_fixed <- fit$parent[!fit$parent$INDIV %in% cand_parents_fixed_past_iteration,]

        #   cand_parents_fixed_past_iteration <- c(cand_parents_not_fixed[!cand_parents_not_fixed$AVAIL_BROOD, "INDIV"],
        #                                         cand_parents_fixed_past_iteration) #some families failed but parents in "PARENT" cohort (and pond)
        cand_parents_not_fixed <- cand_parents_not_fixed[cand_parents_not_fixed$AVAIL_BROOD, ]

        fam_contbn_not_fixed <- aggregate(cand_parents_not_fixed$oc, by = list(cand_parents_not_fixed$FAM), FUN = "sum")
        colnames(fam_contbn_not_fixed) <- c("FAM", "sum_oc")

        fams_to_retain <- fam_contbn_not_fixed[fam_contbn_not_fixed$sum_oc > (indiv_contbn * adj), "FAM"]

        if(length(fams_to_retain) > 0) { #probably should allow more than one individual to be fixed per iteration but any missed will be picked up in next iteration and so like that impact is minimal

          #get max_sum_oc
          fam_contbn_not_fixed[fam_contbn_not_fixed$FAM %in% fams_to_retain, "sum_oc"] <-
            fam_contbn_not_fixed[fam_contbn_not_fixed$FAM %in% fams_to_retain, "sum_oc"] - indiv_contbn # accounting for individuals fixed in this iteration
          max_sum_oc <- max(fam_contbn_not_fixed$sum_oc)

          cand_parents_fixed_current_iteration <- OCFam::get_best_indiv(fish = cand_parents_not_fixed,
                                                                             additional_fams_to_retain = fams_to_retain,
                                                                             candidate_parents = candidate_parents[is.na(candidate_parents$N_AS_PARENT_CURRENT),"INDIV"]) #exclude animals previously used as parents including those families that subsequently died and are now in the "PARENT" pond

          ub[names(ub) %in% cand_parents_fixed_current_iteration] <- indiv_contbn
          ub <- ub[cand$phen$INDIV]
          ub <- ub[names(ub) %in%  names(ub_orig)]

          lb[names(lb) %in% cand_parents_fixed_current_iteration] <- indiv_contbn -1e-10
          lb <- lb[cand$phen$INDIV]
          lb <- lb[names(lb) %in%  names(lb_orig)]

        }
        if(adj != 0) {
          try(fit <- OCFam::run_OC_max_EBV(cand = cand, kinship_constraint = kinship_constraint, ub = ub, lb = lb, opticont_method = opticont_method))
        }


        cand_parents_fixed_past_iteration <- cand_parents_fixed_past_iteration[!cand_parents_fixed_past_iteration %in% cand_parents_fixed_current_iteration] #this shouldn't make any difference but sometimes gets confused
        final_adj <- adj
      }
    }
  }
  #}
  #From last successful iteration with finer steps (step_interval^2)
  #  for(adj in  seq(final_adj - step_interval^2, final_adj - step_interval + step_interval^2, -step_interval^2)) {
  #    source("C:/Users/MHamilton/OneDrive - CGIAR/Current/Bangladesh Carp/Mate allocation and parent selection/CGIP optimal contributions/CGIP OC step 3 01.R")
  #  }

  print("FINAL ITERATION")
  print(paste("Count parents previous iteration =",length(c(cand_parents_fixed_current_iteration, cand_parents_fixed_past_iteration )), "of", 1 / indiv_contbn))
  print(paste("indiv_contbn = ",indiv_contbn))
  print(paste("Iteration adjustment = ",adj))
  print(paste("max_sum_oc = ",max_sum_oc))
  print(paste("indiv_contbn * adj = ",indiv_contbn * adj))

  fit$final.step <-  list(
    c("Count parents previous iteration =",length(c(cand_parents_fixed_current_iteration, cand_parents_fixed_past_iteration )), "of", 1 / indiv_contbn),
    c("indiv_contbn = ",indiv_contbn),
    c("Iteration adjustment = ",adj),
    c("max_sum_oc = ",max_sum_oc),
    c("indiv_contbn * adj = ",indiv_contbn * adj))

  fit_out <- fit[c("info", "summary", "mean", "bc", "obj.fun", "final.step")]

  ################################################################################
  #4. Make up the numbers by identifying the least related individuals
  ################################################################################

  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #Note - with overlapping generation this may need to be changed to ensure contributions from age classes are correct
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  cand_parents_fixed <- c(cand_parents_fixed_current_iteration, cand_parents_fixed_past_iteration )

  Pedig <- rbind(ped[ped$GENERATION !=  max(ped$BORN), ], #c("INDIV", "SIRE", "DAM")
                 ped[ped$INDIV %in% cand_parents_fixed,]) #c("INDIV", "SIRE", "DAM")

  Pedig <- optiSel::prePed(Pedig = Pedig, keep = cand_parents_fixed)

  A_mat <- optiSel::pedIBD(Pedig, keep.only = cand_parents_fixed)

  A_mat <-  A_mat * lb[colnames(A_mat)] / indiv_contbn #weight according to parent contributions

  tmp <- ped[!is.na(ped$N_AS_PARENT_CURRENT) & ped$N_AS_PARENT_CURRENT != 0 & ped$INDIV %in% rownames(A_mat), "N_AS_PARENT_CURRENT"]
  N_retain <- 1 / indiv_contbn +  length(tmp) - sum(tmp) #accounts for greater contribution of some past parents
  rm(tmp)

  while(nrow(A_mat) > N_retain) {
    tmp <- colSums(A_mat) + seq.int(0.00000001, 1e-10, length.out = nrow(A_mat)) #add small amount in case some are the same
    tmp <- tmp[names(tmp) %in% cand_parents_fixed_current_iteration]
    tmp <- names(tmp[tmp == max(tmp)])
    A_mat <- A_mat[rownames(A_mat) != tmp[1], rownames(A_mat) != tmp[1]]
    rm(tmp)
  }
  rm(N_retain)

  indivs_to_retain <- colnames(A_mat)
  rm(A_mat)

  #need to add the total contributions of all individuals to parents
  tmp <- data.frame(INDIV = names(lb)[names(lb) %in% indivs_to_retain],
                    N_TOTAL_AS_PARENT = lb[names(lb) %in% indivs_to_retain] / indiv_contbn)

  if(length(cand_parents_fixed_current_iteration) > 0) {
    #  tmp[tmp$INDIV  %in% cand_parents_fixed_current_iteration,"N_TOTAL_AS_PARENT"] <- 1 #lb not updated
    tmp[tmp$INDIV  %in% cand_parents_fixed_current_iteration,"ADDED_IN_LAST_ITERATION"] <- TRUE
  } else {
    tmp$ADDED_IN_LAST_ITERATION <- FALSE
  }
  tmp[is.na(tmp$ADDED_IN_LAST_ITERATION), "ADDED_IN_LAST_ITERATION"] <- FALSE

  parent_contbn <- ped[ped$INDIV %in% indivs_to_retain,]
  parent_contbn[is.na(parent_contbn$N_AS_PARENT_CURRENT),"N_AS_PARENT_CURRENT"] <- 0
  parent_contbn <- dplyr::left_join(parent_contbn, tmp, by = "INDIV")

  #if insufficient individuals then what???????????????????????????????????????????????

  fam_contbn_parents <- aggregate(parent_contbn$N_AS_PARENT_CURRENT, by = list(parent_contbn$FAM), FUN = "sum")
  colnames(fam_contbn_parents) <- c("FAM", "N_AS_PAST_PARENT")

  fam_contbn         <- aggregate(parent_contbn$N_TOTAL_AS_PARENT, by = list(parent_contbn$FAM), FUN = "sum")
  colnames(fam_contbn) <- c("FAM", "N_INDIV_TOTAL")
  fam_contbn$N_INDIV_TOTAL <- round(fam_contbn$N_INDIV_TOTAL,6)
  fam_contbn <- dplyr::left_join(fam_contbn, fam_contbn_parents, by = "FAM")
  fam_contbn$N_INDIV <- fam_contbn$N_INDIV_TOTAL - fam_contbn$N_AS_PAST_PARENT
  fam_contbn

  return(list(fam_contbn   = fam_contbn,
              parent_contbn = parent_contbn,
              candidate_parents = candidate_parents,
              fit_out      = fit_out))#, summary.opticont = summary.opticont
}





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
  tmp <- tmp[tmp$FAM %in% additional_fams_to_retain &
               tmp$INDIV %in% candidate_parents,]
  cand_parents_fixed <- tmp[!duplicated(tmp$FAM),]
  cand_parents_fixed <- cand_parents_fixed[order(cand_parents_fixed$FAM),"INDIV"]
  rm(tmp)
  return (cand_parents_fixed)
}

################################################################################
#get family K matrix
################################################################################

get_fam_K_matrix <- function(ped, cand_fams) {

  fam_K_matrix <-  unique(ped[,c("FAM", "SIRE", "FAM_SIRE", "DAM", "FAM_DAM")])
  colnames(fam_K_matrix) <-  c("FAM", "SIRE", "FAM_SIRE", "DAM", "FAM_DAM")
  fam_K_matrix <- OCFam::fam_K_matrix_fun(fam_K_matrix[!is.na(fam_K_matrix[,1]),])
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

determine_generations <- function(pedigree) {

  pedigree[is.na(pedigree[,2]) , 2] <- 0
  pedigree[is.na(pedigree[,3]) , 3] <- 0

  tmp <- pedigree[!duplicated(pedigree[,c(2, 3)]) |
                    pedigree[,2] == 0 |
                    pedigree[,3] == 0  , 1]
  trim_ped <- AGHmatrix::filterpedigree(tmp, pedigree[1:3])
  trim_ped <- pedigree[pedigree[,1] %in% trim_ped[,1], 1:3]
  rm(tmp)

  # Extract the INDIV, SIRE, and DAM columns by position
  INDIV <- trim_ped[,1]
  SIRE <- trim_ped[,2]
  DAM <- trim_ped[,3]


  # Initialize generation for each individual
  generation <- rep(NA, nrow(trim_ped))

  # Founders (those with no parents) are generation 0
  generation[SIRE == 0 & DAM == 0] <- 0

  # Iteratively calculate generations for non-founders
  # while (any(is.na(generation))) {
  for (i in seq_len(nrow(trim_ped))) {
    if (is.na(generation[i])) {
      # Fetch SIRE and DAM generations
      sire_gen <- ifelse(SIRE[i] %in% INDIV,
                         generation[INDIV == SIRE[i]], NA)
      dam_gen <- ifelse(DAM[i] %in% INDIV,
                        generation[INDIV == DAM[i]], NA)

      # Calculate generation as average of parents' generations + 1, if available
      if (!is.na(sire_gen) | !is.na(dam_gen)) {
        generation[i] <- mean(c(sire_gen, dam_gen), na.rm = TRUE) + 1
      }
    }
  }
  # }

  trim_ped[,"generation"] <- generation
  trim_ped     <-  trim_ped[!(trim_ped[,2] == 0 & trim_ped[,3] == 0),] #remove founders
  trim_ped <- unique(trim_ped[,2:4]) #retain unique families

  pedigree <- dplyr::left_join(pedigree,
                               trim_ped, by = colnames(pedigree[,2:3]))

  pedigree[pedigree[,2] == 0 & pedigree[,3] == 0,"generation"] <- 0


  return(pedigree$generation)
}

#get lower and upper bounds
get_lb_ub <- function(ped, indiv_contbn, max_parents_per_fam) {

  #exclude based on max_parents_per_fam
  candidate_parents <- ped[ ped$AVAIL_OR_PAST_BROOD,]

  include <- candidate_parents[candidate_parents$AVAIL_BROOD,]
  avail_brood <- include$INDIV
  include <- include[order(include$RANK, decreasing = FALSE),]
  past_brood <- candidate_parents[!candidate_parents$AVAIL_BROOD & candidate_parents$AVAIL_OR_PAST_BROOD,] #Past brood
  include <- rbind(past_brood, include) #past brood at top - can't change past parents
  include <- by(include, include["FAM"], head, n=max_parents_per_fam)
  include <- Reduce(rbind, include)
  include <- c(include[,"INDIV"])

  exlcude <- avail_brood[!avail_brood %in% include] #exclude from avail_brood only - can't change past parents
  rm(include, avail_brood, past_brood)

  candidate_parents <- candidate_parents[, c("INDIV", "N_AS_PARENT_CURRENT")]

  max_parents_per_fam

  candidate_parents$lb <- 0
  candidate_parents[!is.na(candidate_parents$N_AS_PARENT_CURRENT),"lb"] <-
    candidate_parents[!is.na(candidate_parents$N_AS_PARENT_CURRENT),"N_AS_PARENT_CURRENT"] * indiv_contbn - 1e-10

  candidate_parents$ub <- indiv_contbn
  candidate_parents[!is.na(candidate_parents$N_AS_PARENT_CURRENT),"ub"] <-
    candidate_parents[!is.na(candidate_parents$N_AS_PARENT_CURRENT),"N_AS_PARENT_CURRENT"] * indiv_contbn

  candidate_parents[candidate_parents$lb < 0 &
                      !is.na(candidate_parents$N_AS_PARENT_CURRENT), "ub"] <-
    candidate_parents[candidate_parents$lb < 0 &
                        !is.na(candidate_parents$N_AS_PARENT_CURRENT), "ub"] + 1e-10  # ub must be different to lb
  candidate_parents[candidate_parents$lb < 0 &
                      !is.na(candidate_parents$N_AS_PARENT_CURRENT), "lb"] <-  0 # can't be less than zero

  candidate_parents[candidate_parents$INDIV %in% exlcude,"lb"] <- 0
  candidate_parents[candidate_parents$INDIV %in% exlcude,"ub"] <- 1e-10
  candidate_parents[candidate_parents$INDIV %in% exlcude,"N_AS_PARENT_CURRENT"] <- 0

  candidate_parents$Exclude_max_parents_per_fam <- FALSE
  candidate_parents[candidate_parents$INDIV %in% exlcude,"Exclude_max_parents_per_fam"] <- TRUE

  return(candidate_parents)
}
