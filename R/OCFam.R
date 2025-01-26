#' OCFam
#'
#' @description
#' This function implements the optimal contributions method of Hamilton (2020) https://doi.org/10.1093/jhered/esaa051
#' It also addresses the rounding issue associated with standard optimal contributions, particularly in highly fecund species in which relatively few families are generated (or parents are used).
#'
#' @param ped is a data frame with the following columns (class in parentheses):
#' \itemize{
#'  \item{'INDIV' is the individual identifier (character).}
#'  \item{'SIRE' is the male parent identifier (character).}
#'  \item{'DAM' is the female parent identifier (character).}
#'  \item{'FAM' is a full-sibling family identifier (character).}
#'  \item{'SEX' is the sex of the individual. "male" (or 1 = male), "female" (or 2 = female), NA if unknown.  If one animal is of unknown sex, sex will be ignored (integer or character).}
#'  \item{'BORN' integer indicating age class.  May be the year of birth if one age class per year or an integer indicating the sequence of age classes (integer).}
#'  \item{'EBV' is the estimated breeding value (numeric).}
#'  \item{'N_AS_PARENT_PREV' is the number of families in the next age class previously contributed to (integer).}
#'  \item{'AVAIL_BROOD' is TRUE if the individual is candidate parent of next age class (logical).}
#'  }
#' @param N_fams is the number of families to be generated in the next age class including those already produced (see N_AS_PARENT_PREV in ped) (integer)
#' @param kinship_constraint is the maximum value of the mean kinship between families in the next age class, while mean EBV is maximised (max.EBV).  If NA there is no constraint placed on kinship and the average kinship is minimised while EBV is not considered (min.fPED) (numeric between 0 and 1)
#' @param step_interval is a parameter that controls rounding error in individual contributions.  The larger the value the greater the rounding error but the quicker it runs (numeric between 0 and 1)
#' @param gene_flow_vector is a applicable to overlapping generations (if NA discrete generations is assumed).  It is vector representing parental contributions by age class to the next age class. For example, gene_flow_vector = c(0.2, 0.8, 0, 0) - oldest age class to youngest age class.
#' @param min_prop_fams is the proportion of families to be retained (i.e. to contribute at least one parent) from the oldest age class with parental candidates (numeric between 0 and 1)
#' @param max_parents_per_fam is the maximum number of parents to contribute from each family (integer)
#' @return 'fam_contbn' is a data frame containing details of contributions by family:
#' \itemize{
#'  \item{FAM}
#'  \item{N_INDIV_TOTAL}
#'  \item{N_AS_PAST_PARENT}
#'  \item{N_INDIV}
#' }
#' @return 'parent_contbn' is a data frame containing details of contributions by individual:
#' \itemize{
#'  \item{'INDIV' is the individual identifier (character).}
#'  \item{'SIRE' is the male parent identifier (character).}
#'  \item{'DAM' is the female parent identifier (character).}
#'  \item{'FAM' is a full-sibling family identifier (character).}
#'  \item{'BORN' integer indicating age class.  May be the year of birth if one age class per year or an integer indicating the sequence of age classes (integer).}
#'  \item{'EBV' is the estimated breeding value (numeric).}
#'  \item{'N_AS_PARENT_PREV' is the number of families in the next age class previously contributed to (integer).}
#'  \item{'AVAIL_BROOD' is TRUE if the individual is candidate parent of next age class (logical).}
#'  \item{'RANK' individual rank based on EBV}
#'  \item{'AVAIL_OR_PAST_BROOD' is TRUE if AVAIL_BROOD is TRUE or the individual has previously been used as a parent in the next age class.}
#'  \item{'FAM_SIRE' is a full-sibling family identifier of the SIRE (character).}
#'  \item{'FAM_DAM' is a full-sibling family identifier of the DAM (character).}
#'  \item{'BREED' is the species.}
#'  \item{'MATURE' is sexually mature}
#'  \item{'IS_CANDIDATE' was considered a candidate parent in OCFam}
#'  \item{'N_TOTAL_AS_PARENT' is the total count as a parent in the next year class}
#'  \item{'ADDED_IN_LAST_ITERATION' individual was addes in the last iteration of OCFam.  Would be the first individuals to be removed in necessary.}
#' }
#' @return 'candidate_parents' is a data frame containing details of candidate parents:
#' \itemize{
#'  \item{'INDIV' is the individual identifier (character).}
#'  \item{'N_AS_PARENT_PREV' is the number of families contributed to in the next age class prior to running optimal contributions (integer).}
#'  \item{'LB' is the low bound of possilbe contribution to the next age class}
#'  \item{'UB' is the upper bound of possilbe contribution to the next age class}
#'  \item{EXCLUDE_MAX_PARENTS_PER_FAM' If TRUE individual is excluded on the basis that there are better parents available from the family to meet the max_parents_per_fam constraint.}
#' }
#' @return 'fit_out' is a vector containing details of the constraints applied in the implementation of optiSel:
#' @examples
#' #Retrieve example data
#' ped <- OCFam::ped
#' tail(ped)
#'
#' #Run OCFam function
#' OCFam_output <- OCFam::OCFam(ped = ped,
#'                              N_fams = 60,
#'                              kinship_constraint = 0.012,
#'                              step_interval = 0.1,
#'                              gene_flow_vector = NA,
#'                              min_prop_fams = 0.9,
#'                              max_parents_per_fam = 4
#' )
#'
#' OCFam_output$fit_out$summary
#' OCFam_output$fit_out$mean
#' head(OCFam_output$fam_contbn)
#' head(OCFam_output$parent_contbn)
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

#Function to optimise contributions at family level
OCFam  <- function(ped,
                   N_fams, #total families including those already produced
                   kinship_constraint, #if NA opticont_method = "min.fPED" else opticont_method = "max.EBV"
                   step_interval = 0.1,
                   gene_flow_vector, #if NA then discrete generations
                   min_prop_fams,
                   max_parents_per_fam) {

  indiv_contbn <- 1/(N_fams*2)

  overlapping_gens <-  !is.na(gene_flow_vector)[1]

  #change names for optiSel
  colnames(ped)[colnames(ped) == "INDIV"] <- "Indiv"
  colnames(ped)[colnames(ped) == "SIRE"] <- "Sire"
  colnames(ped)[colnames(ped) == "DAM"] <- "Dam"
  colnames(ped)[colnames(ped) == "SEX"] <- "Sex"
  colnames(ped)[colnames(ped) == "BORN"] <- "Born"

  if(!("Sex" %in% colnames(ped))) {
    ped$Sex <- NA
  }

  ped$Sex <- as.character(ped$Sex)
  ped[ped$Sex == "1" & !is.na(ped$Sex),"Sex"] <- "male"
  ped[ped$Sex == "2" & !is.na(ped$Sex),"Sex"] <- "female"
  ped[ped$Sex == "0" & !is.na(ped$Sex),"Sex"] <- NA

  #Data checks###############################################################################
  if(!(sum(ped$Sex %in% c("male","female")) == nrow(ped) |
       sum(is.na(ped$Sex)) == nrow(ped))) {
    break("SEX must be defined for all individuals or not defined for all individuals (i.e. = NA)")
  }







  #candidate_parents - output of get_lb_ub function
  # Indiv                     N_AS_PARENT_PREV lb          ub                                EXCLUDE_MAX_PARENTS_PER_FAM
  # G202302_752E925           NA           0           0.005555556                       FALSE
  # G202302_7535BE7           NA           0           0.005555556                       FALSE
  # G202302_7536818           NA           0           0.005555556                       FALSE
  # G202302_753946A           NA           0           0.005555556                       FALSE
  # G202302_753DC36           NA           0           0.005555556                       FALSE
  # G202302_753E4AF           NA           0           0.005555556                       FALSE

  #indiv_contbn - scalar
  # e.g. 0.005555556 = 1/(N_fams_tilv*2 + sum(ped[ped$LINE == "TILV" & !is.na(ped$N_AS_PARENT_PREV), "N_AS_PARENT_PREV"])) #parents of existing families

  #kinship_constraint
  #e.g. 0.09

  #ped - pedigree
  # Indiv           Sire            Dam             FAM             Born     EBV      RANK FAM_Sire
  # G202302_75FB917 G202202_7AE17C5 G202202_7ADC209 G202302_ON_0178 17.88489 -0.0565  135  G202202_ON_0081
  # G202302_75D0BE8 G202202_7AE17C5 G202202_7ADC209 G202302_ON_0178 17.88489 -0.0556  136  G202202_ON_0081
  # G202302_75D03F6 G202202_7AE17C5 G202202_7ADC209 G202302_ON_0178 17.88489 -0.0549  137  G202202_ON_0081
  # G202302_75D29F6 G202202_7AE17C5 G202202_7ADC209 G202302_ON_0178 17.88489 -0.0545  138  G202202_ON_0081
  # G202302_75FA452 G202202_7AE17C5 G202202_7ADC209 G202302_ON_0178 17.88489 -0.0545  139  G202202_ON_0081
  # G202302_75CF143 G202202_7AE17C5 G202202_7ADC209 G202302_ON_0178 17.88489 -0.0497  140  G202202_ON_0081
  # FAM_Dam         LINE COHORT        AVAIL_BROOD Breed        N_AS_PARENT_PREV     AVAIL_OR_PAST_BROOD
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

  # cand_parents_fixed_initial <- candidate_parents[!is.na(candidate_parents$lb), "Indiv"]

  #  if(length(cand_parents_fixed_initial) >= (1 / indiv_contbn)) {break("Fixed contributions of indivuals equal to or greater than possible number based on indiv_contbn")}
  set.seed(12345)
  tmp <- ped[order(runif(nrow(ped))),c("Indiv", "EBV")]
  tmp <- tmp[order(tmp$EBV ), ]
  tmp$RANK <- nrow(tmp):1
  ped <- dplyr::left_join(ped, tmp, by = c("Indiv", "EBV"))
  rm(tmp)

  ped$AVAIL_OR_PAST_BROOD <- ped$AVAIL_BROOD | (ped$N_AS_PARENT_PREV > 0 & !is.na(ped$N_AS_PARENT_PREV))

  candidate_parents <- OCFam::get_lb_ub(ped = ped,
                                        indiv_contbn = indiv_contbn,
                                        max_parents_per_fam = max_parents_per_fam)

  candidate_parents_orig <- candidate_parents #for output
  colnames(candidate_parents_orig) <- c("INDIV", "N_AS_PARENT_PREV", "LB", "UB", "EXCLUDE_MAX_PARENTS_PER_FAM")

  candidate_parents <- candidate_parents[!(candidate_parents$EXCLUDE_MAX_PARENTS_PER_FAM),]

  if(nrow(candidate_parents) < 3) {break("Not enough candidate parents.  Must be 3 or more.")}

  ################################################################################
  #Compute necessary columns
  ################################################################################

  #ped$Generation <- determine_generations(ped)

  fams <- ped[,c("Indiv", "FAM")]
  colnames(fams) <- c("Sire", "FAM_Sire")
  ped <- dplyr::left_join(ped, fams, by = "Sire")
  colnames(fams) <- c("Dam", "FAM_Dam")
  ped <- dplyr::left_join(ped, fams, by = "Dam")
  rm(fams)

  ped$AVAIL_BROOD <- ped$Indiv %in% candidate_parents[!candidate_parents$EXCLUDE_MAX_PARENTS_PER_FAM, "Indiv"]

  ped$Breed <- "Breed_ignored"

  #ped <- dplyr::left_join(ped, candidate_parents[,c("Indiv", "N_AS_PARENT_PREV")], by = "Indiv")

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

    candidate_parents <- merge(candidate_parents, ped[,c("Indiv", "Born")], by = "Indiv") #get Born

    #keep one individual per family from immature generations
    tmp <- ped[ped$Born > max(candidate_parents$Born),]
    tmp <- tmp[!duplicated(tmp$FAM),]
    if(nrow(tmp)>0) {
      tmp$Mature <- FALSE
      tmp$EBV <- 0
    } else {
      tmp$Mature <- NA[-1]
      tmp$EBV <- NA[-1]
    }
    tmp$Indiv <-  tmp$FAM
    ped <- ped[ped$Born <= max(candidate_parents$Born),]
    ped$Mature <- TRUE
    ped <- rbind(ped, tmp)
    rm(tmp)


    ped$isCandidate <- FALSE
    ped$isCandidate <- ped$Indiv %in% candidate_parents$Indiv

    keep = ped[ped$isCandidate | !ped$Mature,"Indiv"] #keep if isCandidate or is immature
    Pedig <- optiSel::prePed(Pedig = ped, keep = keep) #keep if isCandidate or is immature
    fPED <- optiSel::pedIBD(Pedig, keep.only = keep)


    #replace immature year classes with family K matrix as per Hamilton paper
    fam_K <- OCFam::get_fam_K_matrix(Pedig, Pedig[!Pedig$Mature,"Indiv"])
    fPED[rownames(fam_K),colnames(fam_K)] <- fam_K

    # Pedig <- Pedig[Pedig$Indiv %in% colnames(fPED),]

    Sy <- summary(Pedig)
    Pedig <- merge(Pedig, Sy[, c("Indiv", "equiGen")], on="Indiv")

    #Population means

    #  Pedig[Pedig$EBV == 0,"EBV"] <- runif(nrow(Pedig[Pedig$EBV == 0,]))
    rownames(Pedig) <- Pedig$Indiv
    Pedig$oldest_age_r_vector <- 0
    Pedig[Pedig$Born == min(Pedig[Pedig$AVAIL_OR_PAST_BROOD, "Born"]), "oldest_age_r_vector"] <- 1 #oldest YC to 1
    Pedig$oldest_age_r_vector_2 <- Pedig$oldest_age_r_vector

    cand <- optiSel::candes(phen=Pedig[Pedig$Indiv %in% rownames(fPED),], #c("Indiv",	"Sire",	"Dam",	"Sex",	"Breed",	"Born",		"EBV",	"isCandidate", "oldest_age_r_vector", "oldest_age_r_vector_2")],
                            fPED=fPED,
                            cont = cont)

    #set initial ub and lb
    ub <- setNames(candidate_parents$ub, candidate_parents$Indiv)
    lb <- setNames(candidate_parents$lb, candidate_parents$Indiv)

  }

  if(overlapping_gens == FALSE) {

    Pedig <- optiSel::prePed(Pedig = ped, keep = candidate_parents$Indiv)

    fPED <- optiSel::pedIBD(Pedig, keep.only = candidate_parents$Indiv)

    Sy <- summary(Pedig)
    Pedig <- merge(Pedig, Sy[, c("Indiv", "equiGen")], on="Indiv")

    cont <- data.frame(age = 1,
                       male = 0.5,
                       female = 0.5)
    cand <- optiSel::candes(phen=Pedig[Pedig$Indiv %in% candidate_parents$Indiv,],
                            fPED=fPED,
                            cont = cont)

    #set initial ub and lb
    ub <- setNames(candidate_parents$ub, candidate_parents$Indiv)
    ub <- ub[cand$phen$Indiv]

    lb <- setNames(candidate_parents$lb, candidate_parents$Indiv)
    lb <- lb[cand$phen$Indiv]
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



  #modified this code to limit only families in the youngest year class




  cand_fams  <- unlist(unique(ped[ped$Indiv %in% candidate_parents$Indiv, "FAM"]))
  cand_fams_in_youngest_age_class  <- unlist(unique(ped[ped$Indiv %in% candidate_parents$Indiv &
                                                          ped$Born == max(ped$Born), "FAM"]))

  #cand_fams_in_youngest_age_class <- unlist(unique(ped[ped$Indiv %in%
  #                                                       candidate_parents[candidate_parents$Born == max(candidate_parents$Born),  "Indiv"],
  #                                                     "FAM"]))
  # cand_fams_in_oldest_age_class <- cand_fams[!cand_fams %in% cand_fams_in_youngest_age_class]

  fixed_fams <- unlist(unique(ped[ped$Indiv %in% names(lb[lb > 0.1*indiv_contbn]), "FAM"])) #contribution > 0
  fixed_fams_in_youngest_age_class <- fixed_fams[fixed_fams %in% cand_fams_in_youngest_age_class]

  # fixed_fams <- unique(ped[ped$Indiv %in% ped[!is.na(ped$N_AS_PARENT_PREV),"Indiv"], "FAM"])

  fam_K_matrix <- OCFam::get_fam_K_matrix(ped = ped,
                                          cand_fams = cand_fams)
  print(fam_K_matrix)
  mean(fam_K_matrix)
  # write.csv(fam_K_matrix, "fam_K_matrix.csv")

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
                                                    candidate_parents = candidate_parents[!candidate_parents$EXCLUDE_MAX_PARENTS_PER_FAM,"Indiv"])
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





  ub <- setNames(floor(round(fit$parent$oc/indiv_contbn,4)) * indiv_contbn, fit$parent$Indiv)
  ub[ub == 0] <- indiv_contbn
  ub <- ub[cand$phen$Indiv]
  ub <- ub[names(ub) %in%  names(ub_orig)]

  lb <- setNames(floor(round(fit$parent$oc/indiv_contbn,4)) * indiv_contbn - 1e-10 , fit$parent$Indiv)
  lb[lb < 0] <- 0
  lb <- lb[cand$phen$Indiv]
  lb <- lb[names(lb) %in%  names(lb_orig)]

  ################################################################################
  #3.  add to list of cand_parents_fixed the best individual from families highly represented in OC
  #      - add individuals consecutively with families highly represented in OC fixed first
  ################################################################################

  cand_parents_fixed_current_iteration <- NA #lb[lb != 0]
  cand_parents_fixed_past_iteration <- NULL
  max_sum_oc <- 1

  for(adj in seq((1-step_interval),0,-step_interval)) {

    prev_count <- length(c(cand_parents_fixed_current_iteration, cand_parents_fixed_past_iteration ))

    print(paste("Count parents previous iteration =",length(c(cand_parents_fixed_current_iteration, cand_parents_fixed_past_iteration )), "of", 1 / indiv_contbn))
    print(paste("indiv_contbn = ",indiv_contbn))
    print(paste("Iteration adjustment = ",adj))
    print(paste("max_sum_oc = ",max_sum_oc))
    print(paste("indiv_contbn * adj = ",indiv_contbn * adj))


    #  if(cand_parents_fixed_past_iteration < 1 / indiv_contbn) {


    if(max_sum_oc > (indiv_contbn * adj)) {

      if((sum(lb)) < 1) {


        cand_parents_fixed_past_iteration <- names(lb[lb!=0])

        cand_parents_not_fixed <- fit$parent[!fit$parent$Indiv %in% cand_parents_fixed_past_iteration,]

        #   cand_parents_fixed_past_iteration <- c(cand_parents_not_fixed[!cand_parents_not_fixed$AVAIL_BROOD, "Indiv"],
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
                                                                        candidate_parents = candidate_parents[is.na(candidate_parents$N_AS_PARENT_PREV),"Indiv"]) #exclude animals previously used as parents including those families that subsequently died and are now in the "PARENT" pond

          ub[names(ub) %in% cand_parents_fixed_current_iteration] <- indiv_contbn
          ub <- ub[cand$phen$Indiv]
          ub <- ub[names(ub) %in%  names(ub_orig)]

          lb[names(lb) %in% cand_parents_fixed_current_iteration] <- indiv_contbn -1e-10
          lb <- lb[cand$phen$Indiv]
          lb <- lb[names(lb) %in%  names(lb_orig)]

        }
        if(adj != 0) {
          try(fit <- OCFam::run_OC_max_EBV(cand = cand, kinship_constraint = kinship_constraint, ub = ub, lb = lb, opticont_method = opticont_method))
        }


        cand_parents_fixed_past_iteration <- cand_parents_fixed_past_iteration[!cand_parents_fixed_past_iteration %in% cand_parents_fixed_current_iteration] #this shouldn't make any difference but sometimes gets confused
        final_adj <- adj
      }
    }

    #if already have enough parents the break after one addition iteration
    if(length(c(cand_parents_fixed_current_iteration, cand_parents_fixed_past_iteration )) >= (1 / indiv_contbn) &
       length(c(cand_parents_fixed_current_iteration, cand_parents_fixed_past_iteration )) == prev_count) {break}

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

  Pedig <- rbind(ped[ped$Born !=  max(ped$Born), ], #c("Indiv", "Sire", "Dam")          #Is this correct for overlapping generations???????
                 ped[ped$Indiv %in% cand_parents_fixed,]) #c("Indiv", "Sire", "Dam")

  Pedig <- optiSel::prePed(Pedig = Pedig, keep = cand_parents_fixed)

  A_mat <- optiSel::pedIBD(Pedig, keep.only = cand_parents_fixed)

  A_mat * lb[colnames(A_mat)] / indiv_contbn #weight according to parent contributions


  tmp <- ped[!is.na(ped$N_AS_PARENT_PREV) & ped$N_AS_PARENT_PREV != 0 & ped$Indiv %in% rownames(A_mat), "N_AS_PARENT_PREV"]
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


  ################################################################################
  # THE FOLLOWING CODE WORKS BUT CAN RESULT IN INSUFFICIENT INDIVIDUALS SELECTED.
  # CURRENTLY JUST ACCEPT THAT ADJUSTMENT FOR SEX CONTRIBUTIONS MAY BE REQUIRED.!!!!!!
  ################################################################################

  #  if(sum(is.na(ped$Sex)) > 0) { #if sex is not defined
  #  tmp <- ped[!is.na(ped$N_AS_PARENT_PREV) & ped$N_AS_PARENT_PREV != 0 & ped$Indiv %in% rownames(A_mat), "N_AS_PARENT_PREV"]
  #  N_retain <- 1 / indiv_contbn +  length(tmp) - sum(tmp) #accounts for greater contribution of some past parents
  #  rm(tmp)

  #  while(nrow(A_mat) > N_retain) {
  #    tmp <- colSums(A_mat) + seq.int(0.00000001, 1e-10, length.out = nrow(A_mat)) #add small amount in case some are the same
  #    tmp <- tmp[names(tmp) %in% cand_parents_fixed_current_iteration]
  #    tmp <- names(tmp[tmp == max(tmp)])
  #    A_mat <- A_mat[rownames(A_mat) != tmp[1], rownames(A_mat) != tmp[1]]
  #    rm(tmp)
  #  }
  #  rm(N_retain)

  #  indivs_to_retain <- colnames(A_mat)
  #  rm(A_mat)

  #  } else {

  #   #males##############################################
  #    males <- ped[ped$Sex == "male", "Indiv"]

  #    tmp <- ped[!is.na(ped$N_AS_PARENT_PREV) & ped$N_AS_PARENT_PREV != 0 & ped$Indiv %in% rownames(A_mat), "N_AS_PARENT_PREV"]
  #    N_retain_males <- 0.5 / indiv_contbn +  length(tmp) - sum(tmp) #accounts for greater contribution of some past parents
  #    rm(tmp)

  #    A_mat <- A_mat[rownames(A_mat) %in% males,colnames(A_mat) %in% males]

  #    while(nrow(A_mat) > N_retain_males) {
  #      tmp <- colSums(A_mat) + seq.int(0.00000001, 1e-10, length.out = nrow(A_mat)) #add small amount in case some are the same
  #      tmp <- tmp[names(tmp) %in% cand_parents_fixed_current_iteration]
  #      tmp <- tmp[names(tmp) %in% males]
  #      tmp <- names(tmp[tmp == max(tmp)])
  #      A_mat <- A_mat[rownames(A_mat) != tmp[1], rownames(A_mat) != tmp[1]]
  #      rm(tmp)
  #    }
  #    rm(N_retain_males)

  #    indivs_to_retain <- colnames(A_mat)
  #    rm(A_mat)


  #    #females##############################################

  #    males_removed <- cand_parents_fixed_current_iteration[cand_parents_fixed_current_iteration %in% males][!cand_parents_fixed_current_iteration[cand_parents_fixed_current_iteration %in% males] %in% indivs_to_retain]

  #    A_mat <- optiSel::pedIBD(Pedig, keep.only = cand_parents_fixed[!cand_parents_fixed %in% males_removed])

  #    A_mat * lb[colnames(A_mat)] / indiv_contbn #weight according to parent contributions

  #    females <- ped[ped$Sex == "female", "Indiv"]

  #    tmp <- ped[!is.na(ped$N_AS_PARENT_PREV) & ped$N_AS_PARENT_PREV != 0 & ped$Indiv %in% rownames(A_mat), "N_AS_PARENT_PREV"]
  #    N_retain_females <- 0.5 / indiv_contbn +  length(tmp) - sum(tmp) #accounts for greater contribution of some past parents
  #    rm(tmp)

  #    A_mat <- A_mat[rownames(A_mat) %in% females,colnames(A_mat) %in% females]
  #    while(nrow(A_mat) > N_retain_females) {
  #      tmp <- colSums(A_mat) + seq.int(0.00000001, 1e-10, length.out = nrow(A_mat)) #add small amount in case some are the same
  #      tmp <- tmp[names(tmp) %in% cand_parents_fixed_current_iteration]
  #      tmp <- tmp[names(tmp) %in% females]
  #      tmp <- names(tmp[tmp == max(tmp)])
  #      A_mat <- A_mat[rownames(A_mat) != tmp[1], rownames(A_mat) != tmp[1]]
  #      rm(tmp)
  #    }
  #    rm(N_retain_females)

  #    indivs_to_retain <- c(indivs_to_retain, colnames(A_mat))
  #    rm(A_mat)
  #  }







  #need to add the total contributions of all individuals to parents
  tmp <- data.frame(Indiv = names(lb)[names(lb) %in% indivs_to_retain],
                    N_TOTAL_AS_PARENT = lb[names(lb) %in% indivs_to_retain] / indiv_contbn)

  if(length(cand_parents_fixed_current_iteration) > 0) {
    #  tmp[tmp$Indiv  %in% cand_parents_fixed_current_iteration,"N_TOTAL_AS_PARENT"] <- 1 #lb not updated
    tmp[tmp$Indiv  %in% cand_parents_fixed_current_iteration,"ADDED_IN_LAST_ITERATION"] <- TRUE
  } else {
    tmp$ADDED_IN_LAST_ITERATION <- FALSE
  }
  tmp[is.na(tmp$ADDED_IN_LAST_ITERATION), "ADDED_IN_LAST_ITERATION"] <- FALSE

  parent_contbn <- ped[ped$Indiv %in% indivs_to_retain,]
  parent_contbn[is.na(parent_contbn$N_AS_PARENT_PREV),"N_AS_PARENT_PREV"] <- 0
  parent_contbn <- dplyr::left_join(parent_contbn, tmp, by = "Indiv")

  #if insufficient individuals then what???????????????????????????????????????????????

  fam_contbn_parents <- aggregate(parent_contbn$N_AS_PARENT_PREV, by = list(parent_contbn$FAM), FUN = "sum")
  colnames(fam_contbn_parents) <- c("FAM", "N_AS_PAST_PARENT")

  fam_contbn         <- aggregate(parent_contbn$N_TOTAL_AS_PARENT, by = list(parent_contbn$FAM), FUN = "sum")
  colnames(fam_contbn) <- c("FAM", "N_INDIV_TOTAL")
  fam_contbn$N_INDIV_TOTAL <- round(fam_contbn$N_INDIV_TOTAL,6)
  fam_contbn <- dplyr::left_join(fam_contbn, fam_contbn_parents, by = "FAM")
  fam_contbn$N_INDIV <- fam_contbn$N_INDIV_TOTAL - fam_contbn$N_AS_PAST_PARENT
  fam_contbn

  #change names from optiSel
  colnames(parent_contbn)[colnames(parent_contbn) == "Indiv"] <- "INDIV"
  colnames(parent_contbn)[colnames(parent_contbn) == "Sire"] <- "SIRE"
  colnames(parent_contbn)[colnames(parent_contbn) == "Dam"] <- "DAM"
  colnames(parent_contbn)[colnames(parent_contbn) == "Sex"] <- "SEX"
  colnames(parent_contbn)[colnames(parent_contbn) == "Born"] <- "BORN"
  colnames(parent_contbn)[colnames(parent_contbn) == "FAM_Sire"] <- "FAM_SIRE"
  colnames(parent_contbn)[colnames(parent_contbn) == "FAM_Dam"] <- "FAM_DAM"
  colnames(parent_contbn)[colnames(parent_contbn) == "Breed"] <- "BREED"
  colnames(parent_contbn)[colnames(parent_contbn) == "Mature"] <- "MATURE"
  colnames(parent_contbn)[colnames(parent_contbn) == "isCandidate"] <- "IS_CANDIDATE"

  colnames(candidate_parents)[colnames(candidate_parents) == "Indiv"] <- "INDIV"
  colnames(candidate_parents)[colnames(candidate_parents) == "lb"] <- "LB"
  colnames(candidate_parents)[colnames(candidate_parents) == "ub"] <- "UB"
  colnames(candidate_parents)[colnames(candidate_parents) == "Born"] <- "BORN"

  fam_contbn <- fam_contbn[order(fam_contbn$FAM),]
  parent_contbn <- parent_contbn[order(parent_contbn$INDIV),]
  candidate_parents <- candidate_parents[order(candidate_parents$INDIV),]

  return(list(fam_contbn   = fam_contbn,
              parent_contbn = parent_contbn,
              candidate_parents = candidate_parents_orig,
              fit_out      = fit_out))#, summary.opticont = summary.opticont
}





#The script includes the following non-base R functions:
#
# optiSel::opticont
# AGHmatrix::Amatrix



################################################################################
#function to run optiSel::opticont depending on opticont_method
################################################################################
#' @export
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

#' @export
get_best_indiv <- function(fish, additional_fams_to_retain, candidate_parents) {
  tmp <- fish[order(fish$RANK, decreasing = TRUE),]
  tmp <- tmp[tmp$FAM %in% additional_fams_to_retain &
               tmp$Indiv %in% candidate_parents,]
  cand_parents_fixed <- tmp[!duplicated(tmp$FAM),]
  cand_parents_fixed <- cand_parents_fixed[order(cand_parents_fixed$FAM),"Indiv"]
  rm(tmp)
  return (cand_parents_fixed)
}

################################################################################
#get family K matrix
################################################################################
#' @export
get_fam_K_matrix <- function(ped, cand_fams) {

  fam_K_matrix <-  unique(ped[,c("FAM", "Sire", "FAM_Sire", "Dam", "FAM_Dam")])
  colnames(fam_K_matrix) <-  c("FAM", "Sire", "FAM_Sire", "Dam", "FAM_Dam")
  fam_K_matrix <- fam_K_matrix[!is.na(fam_K_matrix$FAM),]
  fam_K_matrix <- fam_K_matrix[fam_K_matrix$FAM != "",]
  fam_K_matrix <- OCFam::fam_K_matrix_fun(fam_K_matrix[!is.na(fam_K_matrix[,1]),])
  fam_K_matrix <- as.matrix(fam_K_matrix$K_matrix_families)
  fam_K_matrix <- fam_K_matrix[colnames(fam_K_matrix) %in% cand_fams, colnames(fam_K_matrix) %in% cand_fams]
  return(fam_K_matrix)
}

################################################################################
#fam_K_matrix_fun
################################################################################
#' @export
fam_K_matrix_fun <- function(family_dat) {
  #data
  # class: data.frame
  # fields:
  #   FAM
  #   Sire
  #   FAM_Sire
  #   Dam
  #   FAM_Dam

  #load required library - AGHmatrix
  #  if("AGHmatrix" %in% installed.packages()[,"Package"] == FALSE) { install.packages("AGHmatrix" , repos='https://cran.csiro.au/')}  #Install nadiv package if not already installed
  #  library(AGHmatrix)

  #sires
  sires <- matrix(unique(family_dat$Sire), ncol = 1)
  colnames(sires)[1] <- "Sire"
  sires <- merge(sires, family_dat[,c("Sire","FAM_Sire")], by = "Sire", all.x = TRUE)
  colnames(sires) <- c("Indiv_id","FAM")
  sires <- merge(sires, family_dat[,c("FAM","Sire","Dam")], by = "FAM", all.x = TRUE)
  sires <- subset(sires, select=-c(FAM)) #remove FAM column

  #dams
  dams <- matrix(unique(family_dat$Dam), ncol = 1)
  colnames(dams)[1] <- "Dam"
  dams <- merge(dams, family_dat[,c("Dam","FAM_Dam")], by = "Dam", all.x = TRUE)
  colnames(dams) <- c("Indiv_id","FAM")
  dams <- merge(dams, family_dat[,c("FAM","Sire","Dam")], by = "FAM", all.x = TRUE)
  dams <- subset(dams, select=-c(FAM)) #remove FAM column

  parents <- rbind(sires,dams)
  parents <- unique(parents)
  parents <- parents[order(parents[,"Indiv_id"], decreasing = FALSE), ]

  families            <- unique(family_dat)
  tmp1 <- families
  tmp2 <- families
  tmp1$FAM <- paste0("1_", tmp1$FAM) #needs a unique id
  families <- rbind(tmp1, tmp2)

  fam_pedigree        <- rbind(as.matrix(parents),as.matrix(families[,c("FAM","Sire","Dam")]))

  #Generate K Matrix
  tmp <- fam_pedigree
  tmp <- tmp[tmp[,"Indiv_id"] != 0,]
  tmp[is.na(tmp[,"Sire"]),"Sire"] <- 0
  tmp[is.na(tmp[,"Dam"]),"Dam"] <- 0
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

################################################################################
#determine_generations
################################################################################

#' @export
determine_generations <- function(pedigree) {

  pedigree[is.na(pedigree[,2]) , 2] <- 0
  pedigree[is.na(pedigree[,3]) , 3] <- 0

  tmp <- pedigree[!duplicated(pedigree[,c(2, 3)]) |
                    pedigree[,2] == 0 |
                    pedigree[,3] == 0  , 1]
  trim_ped <- AGHmatrix::filterpedigree(tmp, pedigree[1:3])
  trim_ped <- pedigree[pedigree[,1] %in% trim_ped[,1], 1:3]
  rm(tmp)

  # Extract the Indiv, Sire, and Dam columns by position
  Indiv <- trim_ped[,1]
  Sire <- trim_ped[,2]
  Dam <- trim_ped[,3]


  # Initialize generation for each individual
  generation <- rep(NA, nrow(trim_ped))

  # Founders (those with no parents) are generation 0
  generation[Sire == 0 & Dam == 0] <- 0

  # Iteratively calculate generations for non-founders
  # while (any(is.na(generation))) {
  for (i in seq_len(nrow(trim_ped))) {
    if (is.na(generation[i])) {
      # Fetch Sire and Dam generations
      sire_gen <- ifelse(Sire[i] %in% Indiv,
                         generation[Indiv == Sire[i]], NA)
      dam_gen <- ifelse(Dam[i] %in% Indiv,
                        generation[Indiv == Dam[i]], NA)

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
#' @export
get_lb_ub <- function(ped, indiv_contbn, max_parents_per_fam) {

  #exclude based on max_parents_per_fam
  candidate_parents <- ped[ ped$AVAIL_OR_PAST_BROOD,]

  include <- candidate_parents[candidate_parents$AVAIL_BROOD,]
  avail_brood <- include$Indiv
  include <- include[order(include$RANK, decreasing = FALSE),]
  past_brood <- candidate_parents[!candidate_parents$AVAIL_BROOD & candidate_parents$AVAIL_OR_PAST_BROOD,] #Past brood
  include <- rbind(past_brood, include) #past brood at top - can't change past parents
  include <- by(include, include["FAM"], head, n=max_parents_per_fam)
  include <- Reduce(rbind, include)
  include <- c(include[,"Indiv"])

  exlcude <- avail_brood[!avail_brood %in% include] #exclude from avail_brood only - can't change past parents
  rm(include, avail_brood, past_brood)

  candidate_parents <- candidate_parents[, c("Indiv", "N_AS_PARENT_PREV")]

  candidate_parents[candidate_parents$N_AS_PARENT_PREV == 0,"N_AS_PARENT_PREV"] <- NA #if 0 make NA for following lines

  #max_parents_per_fam

  candidate_parents$lb <- 0
  candidate_parents[!is.na(candidate_parents$N_AS_PARENT_PREV),"lb"] <-
    candidate_parents[!is.na(candidate_parents$N_AS_PARENT_PREV),"N_AS_PARENT_PREV"] * indiv_contbn - 1e-10

  candidate_parents$ub <- indiv_contbn
  candidate_parents[!is.na(candidate_parents$N_AS_PARENT_PREV),"ub"] <-
    candidate_parents[!is.na(candidate_parents$N_AS_PARENT_PREV),"N_AS_PARENT_PREV"] * indiv_contbn

  candidate_parents[candidate_parents$lb < 0 &
                      !is.na(candidate_parents$N_AS_PARENT_PREV), "ub"] <-
    candidate_parents[candidate_parents$lb < 0 &
                        !is.na(candidate_parents$N_AS_PARENT_PREV), "ub"] + 1e-10  # ub must be different to lb
  candidate_parents[candidate_parents$lb < 0 &
                      !is.na(candidate_parents$N_AS_PARENT_PREV), "lb"] <-  0 # can't be less than zero

  candidate_parents[candidate_parents$Indiv %in% exlcude,"lb"] <- 0
  candidate_parents[candidate_parents$Indiv %in% exlcude,"ub"] <- 1e-10
  candidate_parents[candidate_parents$Indiv %in% exlcude,"N_AS_PARENT_PREV"] <- 0

  candidate_parents$EXCLUDE_MAX_PARENTS_PER_FAM <- FALSE
  candidate_parents[candidate_parents$Indiv %in% exlcude,"EXCLUDE_MAX_PARENTS_PER_FAM"] <- TRUE

  return(candidate_parents)
}
