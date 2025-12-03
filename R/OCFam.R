#' OCFam
#'
#' @description
#' This function addresses the rounding issue associated with standard optimal contributions, particularly in highly fecund species in which relatively few families are generated (or parents are used).
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
#'  \item{'AVAIL_PARENT' is TRUE if the individual is candidate parent of next age class (logical).}
#'  }
#' @param sire_count_by_age is scalar (discrete generations) or a vector (overlapping generations) representing required sires by age class. For example (overlapping generations), sire_count_by_age = c(20, 80, 0, 0) - oldest parental age class to youngest age class.
#' @param dam_count_by_age is scalar (discrete generations) or a vector (overlapping generations) representing required dams by age class. For example (overlapping generations), dam_count_by_age = c(20, 80, 0, 0) - oldest parental age class to youngest age class.
#' @param kinship_constraint is the maximum value of the mean kinship among families in the next age class, while mean EBV is maximised (max.EBV).  If NA there is no constraint placed on kinship and the average kinship is minimised while EBV is not considered (min.fPED) (numeric between 0 and 1; default = NA)
#' @param iterations specifies the number of iterations. The greater the number of iterations specified, the closer the solution will be to the optimum but the longer it will take to run (positive integer; default = 10)
#' @param min_prop_fams is the proportion of families to be retained (i.e. to contribute at least one parent) from the oldest age class with parental candidates (numeric between 0 and 1; default = 0)
#' @param max_parents_per_fam is the maximum number of parents to contribute from each family (integer; default = 10000)
#' @param fast when TRUE, only candidate parents with the highest estimated breeding values (EBVs) within each family are retained as candidate parents, up to the limit defined by max_parents_per_fam (logical; default = TRUE)

#' @return 'fam_contbn' is a data frame containing details of contributions by family:
#' \itemize{
#'  \item{'FAM' is a full-sibling family identifier (character).}
#'  \item{'N_INDIV_TOTAL' is the total count of parental contributions to the next age class (integer)}
#'  \item{N_AS_PAST_PARENT is the count of past parental contributions to the next age class  (integer).}
#'  \item{N_INDIV is the count of parental contributions to the next age class from the current OCFam run.  Equals N_INDIV_TOTAL minus N_AS_PAST_PARENT (integer)}
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
#'  \item{'AVAIL_PARENT' is TRUE if the individual is candidate parent of next age class (logical).}
#'  \item{'RANK' individual rank based on EBV}
#'  \item{'AVAIL_OR_PAST_PARENT' is TRUE if AVAIL_PARENT is TRUE or the individual has previously been used as a parent in the next age class.}
#'  \item{'FAM_SIRE' is a full-sibling family identifier of the SIRE (character).}
#'  \item{'FAM_DAM' is a full-sibling family identifier of the DAM (character).}
#'  \item{'BREED' is the species.}
#'  \item{'MATURE' is sexually mature}
#'  \item{'IS_CANDIDATE' was considered a candidate parent in OCFam}
#'  \item{'N_TOTAL_AS_PARENT' is the total count of parental contributions to the next age class}
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
#' @return 'fit_out' is a vector containing details of the constraints applied in the implementation of optiSel.
#' @examples
#' #Retrieve small example data
#' ped_small <- OCFam::ped_small
#' tail(ped_small)
#'
#' #Run OCFam function
#' OCFam_out_small <- OCFam::OCFam(ped = ped_small,
#'                              sire_count_by_age = 60,
#'                              dam_count_by_age = 60,
#'                              kinship_constraint = NA, #This is a small example.  Cannot constrain kinship_constraint.
#'                              iterations = 10,
#'                              min_prop_fams = 0.9,
#'                              max_parents_per_fam = 4,
#'                              fast = TRUE
#' )
#'
#' OCFam_out_small$fit_out$mean
#' head(OCFam_out_small$fam_contbn)
#' head(OCFam_out_small$parent_contbn)
#'
#' #Retrieve larger example (sex unknown)
#' ped_big <- OCFam::ped_big
#' tail(ped_big)
#'
#' #Run OCFam function
#' OCFam_out_big <- OCFam::OCFam(ped = ped_big,
#'                              sire_count_by_age = c(20, 100),
#'                              dam_count_by_age = c(10, 110),
#'                              kinship_constraint = 0.014,
#'                              iterations = 10,
#'                              min_prop_fams = 0,
#'                              max_parents_per_fam = 4,
#'                              fast = TRUE
#' )
#'
#' OCFam_out_big$fit_out$mean
#' head(OCFam_out_big$fam_contbn)
#' head(OCFam_out_big$parent_contbn)
#' @import optiSel
#' @import AGHmatrix
#' @import dplyr

#Functions that are not base R functions
# dplyr::left_join
# optiSel::candes
# optiSel::pedIBD
# optiSel::prePed
# OCFam::run_OC_max_EBV
# OCFam::get_best_indiv
# OCFam::get_fam_K_matrix

###############################################################################
#Currently OCFam doesn't allow integer contributions to deviate from 0 or 1
# - in mass spawning, may need to force minimum contributions to be > 1 for example???
###############################################################################

#' @export
#Function to optimise contributions at family level
OCFam  <- function(ped,
                   sire_count_by_age,
                   dam_count_by_age,
                   kinship_constraint = NA, #if NA opticont_method = "min.fPED" else opticont_method = "max.EBV"
                   iterations = 10,
                   min_prop_fams = 0,
                   max_parents_per_fam = 10000,
                   fast = TRUE) {

  ## ==== DATA CHECKS ==========================================================

  OCFam::check_OCFam_ped(ped)
  OCFam::check_OCFam_args(kinship_constraint = kinship_constraint,
                          iterations = iterations,
                          sire_count_by_age = sire_count_by_age,
                          dam_count_by_age = dam_count_by_age,
                          min_prop_fams = min_prop_fams,
                          max_parents_per_fam = max_parents_per_fam)

  ## ===========================================================================

  indiv_contbn_males <- 1/(sum(sire_count_by_age)*2)
  indiv_contbn_females <- 1/(sum(dam_count_by_age)*2)
  step_interval <- 1/iterations

  overlapping_gens <-  length(sire_count_by_age) > 1

  #Date check
  if(overlapping_gens == FALSE) {
    if(sum(ped[ped$BORN != max(ped$BORN),"AVAIL_PARENT"]) > 0) {
      stop("If generations are not overlapping, AVAIL_PARENT in 'ped' must be FALSE for all individuals not in the oldest age class")
    }
  }

  #  if(overlapping_gens == TRUE) {
  #    stop("OCFam does not yet cope with overlapping generations.  This will be fixed at some point.")
  #  }

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

  # Data checks on Sex vs AVAIL_PARENT ##########################################
  tmp <- ped[ped$AVAIL_PARENT, , drop = FALSE]

  if (nrow(tmp) > 0) {
    n_def <- sum(tmp$Sex %in% c("male", "female"), na.rm = TRUE)
    n_na  <- sum( is.na(tmp$Sex))

    ## Either all candidates have Sex defined, or none do
    if (!(n_def == nrow(tmp) || n_na == nrow(tmp))) {
      stop("SEX must either be defined as 'male'/'female' (or '1' for male and '2' for female) for all individuals ",
           "where AVAIL_PARENT is TRUE, or be NA for all such individuals.")
    }
  }

  set.seed(12345)
  tmp <- ped[order(runif(nrow(ped))),c("Indiv", "EBV")]
  tmp <- tmp[order(tmp$EBV ), ]
  tmp$RANK <- nrow(tmp):1
  ped <- dplyr::left_join(ped, tmp, by = c("Indiv", "EBV"))
  rm(tmp)

  ped$AVAIL_OR_PAST_PARENT <- ped$AVAIL_PARENT | (ped$N_AS_PARENT_PREV > 0 & !is.na(ped$N_AS_PARENT_PREV))


  if (all(is.na(ped[ped$AVAIL_PARENT, "Sex"]))) {
    sex_known <- FALSE
  } else if (any(is.na(ped[ped$AVAIL_PARENT, "Sex"]))) {
    stop("Some candidates have Sex = NA and others defined; please fix the input.")
  } else {
    sex_known <- TRUE
  }

  candidate_parents <- OCFam::get_lb_ub(ped = ped,
                                        indiv_contbn_males = indiv_contbn_males,
                                        indiv_contbn_females = indiv_contbn_females,
                                        max_parents_per_fam = max_parents_per_fam,  #max_parents_per_fam possibly shouldn't be specified in get_lb_ub if sex is known because it will potentially exclude individuals from a sex that might otherwise be included.
                                        sex_known = sex_known,
                                        fast = fast)

  candidate_parents_orig <- candidate_parents #for output
  colnames(candidate_parents_orig) <- c("INDIV", "N_AS_PARENT_PREV", "LB", "UB", "EXCLUDE_MAX_PARENTS_PER_FAM")

  candidate_parents <- candidate_parents[!(candidate_parents$EXCLUDE_MAX_PARENTS_PER_FAM),]

  if (nrow(candidate_parents) < 3) {
    stop("Not enough candidate parents: 'OCFam' requires at least 3 candidate ",
         "or past parents (rows with AVAIL_OR_PAST_PARENT == TRUE).")
  }


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

  ped$AVAIL_PARENT <- ped$Indiv %in% candidate_parents[!candidate_parents$EXCLUDE_MAX_PARENTS_PER_FAM, "Indiv"]

  ped$Breed <- "Breed_ignored"

  ################################################################################
  #Inputs for optiSel
  ################################################################################

  if(overlapping_gens == TRUE) {

    #contributions (averaged across sexes)
    #   avg_gene_flow_vector <- (sire_count_by_age + dam_count_by_age) / sum(sire_count_by_age, dam_count_by_age)
    #    avg_gen_interval <- sum(avg_gene_flow_vector * seq(length(avg_gene_flow_vector),1))
    #    avg_r_vector     <- matrix(NA, nrow = 1, ncol = length(avg_gene_flow_vector)) #not separated by sex
    #    avg_r_vector     <- as.vector(avg_r_vector)
    #    for(i in 1:length(avg_gene_flow_vector)) {
    #      avg_r_vector[i] = sum(avg_gene_flow_vector[1:i]/avg_gen_interval) #oldest age class to youngest age class
    #    }

    # contributions (sires)
    #   In optiSel:
    #    cont_age	Meaning
    #      1	        NEWBORN cohort (offspring produced by selected parents)
    #      2	        youngest breeding cohort
    #      3	        older breeding cohort

    sire_gene_flow_vector <- sire_count_by_age / sum(sire_count_by_age)
    sire_gen_interval <- sum(sire_gene_flow_vector * seq(length(sire_gene_flow_vector),1))
    sire_r_vector     <- matrix(NA, nrow = 1, ncol = length(sire_gene_flow_vector)) #not separated by sex
    sire_r_vector     <- as.vector(sire_r_vector)
    for(i in 1:length(sire_gene_flow_vector)) {
      sire_r_vector[i] = sum(sire_gene_flow_vector[1:i]/sire_gen_interval) #oldest age class to youngest age class
    }

    # contributions (dams)
    dam_gene_flow_vector <- dam_count_by_age / sum(dam_count_by_age)
    dam_gen_interval <- sum(dam_gene_flow_vector * seq(length(dam_gene_flow_vector),1))
    dam_r_vector     <- matrix(NA, nrow = 1, ncol = length(dam_gene_flow_vector)) #not separated by sex
    dam_r_vector     <- as.vector(dam_r_vector)
    for(i in 1:length(dam_gene_flow_vector)) {
      dam_r_vector[i] = sum(dam_gene_flow_vector[1:i]/dam_gen_interval) #oldest age class to youngest age class
    }

    #I don't know what should be done here really
    sire_r_vector <- c(sire_r_vector, 1/min(sire_gen_interval, dam_gen_interval))
    sire_r_vector <- sire_r_vector/sum(sire_r_vector)
    dam_r_vector <- c(dam_r_vector, 1/min(sire_gen_interval, dam_gen_interval))
    dam_r_vector <- dam_r_vector/sum(dam_r_vector)

    cont <- data.frame(age = length(sire_r_vector):1,                #oldest age class to youngest age class
                       male = sire_r_vector / 2,
                       female = dam_r_vector / 2)
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

    #replace immature age classes with family K matrix as per Hamilton paper
    fam_K <- OCFam::get_fam_K_matrix(Pedig, Pedig[!Pedig$Mature,"Indiv"])
    fPED[rownames(fam_K),colnames(fam_K)] <- fam_K

    # Pedig <- Pedig[Pedig$Indiv %in% colnames(fPED),]

    Sy <- summary(Pedig)
    Pedig <- merge(Pedig, Sy[, c("Indiv", "equiGen")], on="Indiv")

    #Population means

    #  Pedig[Pedig$EBV == 0,"EBV"] <- runif(nrow(Pedig[Pedig$EBV == 0,]))
    rownames(Pedig) <- Pedig$Indiv
  #  Pedig$oldest_age_r_vector <- 0
  #  Pedig[Pedig$Born == min(Pedig[Pedig$AVAIL_OR_PAST_PARENT, "Born"]), "oldest_age_r_vector"] <- 1 #oldest YC to 1
  #  Pedig$oldest_age_r_vector_2 <- Pedig$oldest_age_r_vector

    cand <- optiSel::candes(phen=Pedig[Pedig$Indiv %in% rownames(fPED),], #c("Indiv",	"Sire",	"Dam",	"Sex",	"Breed",	"Born",		"EBV",	"isCandidate", "oldest_age_r_vector", "oldest_age_r_vector_2")],
                            fPED=fPED,
                            cont = cont)

    #set initial ub and lb
    ub <- setNames(candidate_parents$ub, candidate_parents$Indiv)
    ub <- ub[cand$phen$Indiv]

    lb <- setNames(candidate_parents$lb, candidate_parents$Indiv)
    lb <- lb[cand$phen$Indiv]

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
  #generate candidate list with one individual for Min_N_fams_selection_to_retain family (least related families)
  ################################################################################

  #identify candidate families

  cand_fams  <- unlist(unique(ped[ped$Indiv %in% candidate_parents$Indiv, "FAM"]))
  cand_fams_in_youngest_age_class  <- unlist(unique(ped[ped$Indiv %in% candidate_parents$Indiv &
                                                          ped$Born == max(ped$Born), "FAM"]))

  #previously fixed
  fixed_fams <- unlist(unique(ped[ped$Indiv %in% names(lb[lb > 0.1*min(indiv_contbn_males,indiv_contbn_females)]), "FAM"])) #contribution > 0

  fixed_fams_in_youngest_age_class <- fixed_fams[fixed_fams %in% cand_fams_in_youngest_age_class]

  fam_K_matrix <- OCFam::get_fam_K_matrix(ped = ped,
                                          cand_fams = cand_fams)
  # print(fam_K_matrix)
  # mean(fam_K_matrix)
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

  tmp <- ped
  tmp$oc <- tmp$EBV
  tmp$Age <- max(tmp$Born) - tmp$Born + 1

  additional_parents_fixed <- OCFam::get_best_indiv(fish = tmp[tmp$AVAIL_PARENT ,], #exclude animals previously used as parents
                                                    additional_fams_to_retain = additional_fams_to_retain,
                                                    candidates = candidate_parents[!candidate_parents$EXCLUDE_MAX_PARENTS_PER_FAM,"Indiv"],
                                                    count_male_req_by_age   = sire_count_by_age,
                                                    count_female_req_by_age = dam_count_by_age,
                                                    sex_known = sex_known)

  if(sex_known) {

    additional_sires_fixed <- additional_parents_fixed[additional_parents_fixed %in% ped[ped$Sex == "male", "Indiv"]]
    ub[names(ub) %in% additional_sires_fixed] <- indiv_contbn_males
    lb[names(lb) %in% additional_sires_fixed] <- indiv_contbn_males - 1e-10

    additional_dams_fixed <- additional_parents_fixed[additional_parents_fixed %in% ped[ped$Sex == "female", "Indiv"]]
    ub[names(ub) %in% additional_dams_fixed] <- indiv_contbn_females
    lb[names(lb) %in% additional_dams_fixed] <- indiv_contbn_females - 1e-10

    rm(additional_sires_fixed, additional_dams_fixed, additional_parents_fixed, additional_fams_to_retain)

  } else {

    ub[names(ub) %in% additional_parents_fixed] <- indiv_contbn_males
    lb[names(lb) %in% additional_parents_fixed] <- indiv_contbn_males - 1e-10
    rm(additional_parents_fixed, additional_fams_to_retain)

  }
  ub_orig <- ub
  lb_orig <- lb

  ################################################################################
  #add to list of cand_parents_fixed the best individual from families highly represented in OC
  #      - add individuals consecutively with families highly represented in OC fixed first
  ################################################################################

  male_indiv   <- Pedig$Indiv[Pedig$Sex == "male"]   # Subset the Pedig dataframe for male individuals
  female_indiv <- Pedig$Indiv[Pedig$Sex == "female"] # Subset the Pedig dataframe for female individuals

  cand_parents_fixed_current_iteration <- NA #lb[lb != 0]
  cand_parents_fixed_past_iteration <- NULL

  max_sum_oc <- 1

  candidate_parents <- merge(candidate_parents, ped[,c("Indiv", "Sex")], by = "Indiv") #get Sex

  for(adj in seq(1, (0-step_interval), -step_interval)) {

    prev_count <- length(c(cand_parents_fixed_current_iteration, cand_parents_fixed_past_iteration ))

    print(paste("Total count of parents fixed in previous iteration =",length(c(cand_parents_fixed_current_iteration, cand_parents_fixed_past_iteration )), "of", 0.5 / (indiv_contbn_males) + 0.5 / (indiv_contbn_females)))
    print(paste("indiv_contbn_males (0.5 / N_sires) = ",indiv_contbn_males))
    print(paste("indiv_contbn_females (0.5 / N_dams) = ",indiv_contbn_females))
    print(paste("Iteration = ", iterations - adj/step_interval + 1))
    print(paste("Threshold adjustment (adj) = ",adj))
    print(paste("Threshold sum of individual contributions for fixing an individual from a family = ", min(c(indiv_contbn_males,indiv_contbn_females)) * adj))

    if(max_sum_oc > (min(c(indiv_contbn_males,indiv_contbn_females)) * adj)) { #It is the sum of oc within families that determines if an individual from that family is fix in any given iteration

      if((sum(lb)) < 1) {

        if(adj != 1) {
          cand_parents_fixed_past_iteration <- names(lb[lb!=0])

          cand_parents_not_fixed <- fit$parent[!fit$parent$Indiv %in% cand_parents_fixed_past_iteration,]

          cand_parents_not_fixed <- cand_parents_not_fixed[cand_parents_not_fixed$AVAIL_PARENT, ]

          fam_contbn_not_fixed <- aggregate(cand_parents_not_fixed$oc, by = list(cand_parents_not_fixed$FAM), FUN = "sum")
          colnames(fam_contbn_not_fixed) <- c("FAM", "sum_oc")

          fams_to_retain <- fam_contbn_not_fixed[fam_contbn_not_fixed$sum_oc > (min(c(indiv_contbn_males,indiv_contbn_females)) * adj), "FAM"]

          if(length(fams_to_retain) > 0) { #probably should allow more than one individual to be fixed per iteration but any missed will be picked up in next iteration and so likely that impact is minimal

            #get max_sum_oc
            fam_contbn_not_fixed[fam_contbn_not_fixed$FAM %in% fams_to_retain, "sum_oc"] <-
              fam_contbn_not_fixed[fam_contbn_not_fixed$FAM %in% fams_to_retain, "sum_oc"] -
              min(c(indiv_contbn_males,indiv_contbn_females)) # accounting for individuals fixed in this iteration
            max_sum_oc <- max(fam_contbn_not_fixed$sum_oc)

            #  fix contributions from sex-age class once there are enough

            count_male_req_by_age <- NULL
            count_male_req_by_age[1:length(sire_count_by_age)] <- NA
            count_female_req_by_age <- NULL
            count_female_req_by_age[1:length(dam_count_by_age)] <- NA

            if(sex_known) {
              for(age in length(sire_count_by_age):1) {

                sum_male_lb <- sum(lb[names(lb) %in% male_indiv &
                                        names(lb) %in% fit$parent[ fit$parent$Age == age, "Indiv"]], na.rm = TRUE) # Match male-age identifiers in lb and sum the corresponding values
                count_male_fixed <- sum_male_lb / indiv_contbn_males
                count_male_req_by_age[length(count_male_req_by_age) - age + 1] <-
                  round(sire_count_by_age[length(sire_count_by_age) - age + 1] - count_male_fixed)

                sum_female_lb <- sum(lb[names(lb) %in% female_indiv &
                                          names(lb) %in% fit$parent[ fit$parent$Age == age, "Indiv"]], na.rm = TRUE) # Match female-age identifiers in lb and sum the corresponding values
                count_female_fixed <- sum_female_lb / indiv_contbn_females
                count_female_req_by_age[length(count_female_req_by_age) - age + 1] <-
                  round(dam_count_by_age[length(dam_count_by_age) - age + 1] - count_female_fixed)
              }
            } else { #only use one vector if sex unknown - count_male_req_by_age
              for(age in length(sire_count_by_age):1) {
                sum_lb <- sum(lb[names(lb) %in% fit$parent[ fit$parent$Age == age, "Indiv"]], na.rm = TRUE) # Match age identifiers in lb and sum the corresponding values
                count_fixed <- sum_lb / indiv_contbn_males
                count_male_req_by_age[length(count_male_req_by_age) - age + 1] <-
                  round(sire_count_by_age[length(sire_count_by_age) - age + 1] * 2 - count_fixed) # *2 because only use male or female vector in get_best_indiv if sex unknown
              }
              count_female_req_by_age <- NA #only use male  vector in get_best_indiv if sex unknown
            }

            cand_parents_fixed_current_iteration <- OCFam::get_best_indiv(fish = cand_parents_not_fixed,
                                                                          additional_fams_to_retain = fams_to_retain,
                                                                          candidates = candidate_parents[is.na(candidate_parents$N_AS_PARENT_PREV),"Indiv"],
                                                                          count_male_req_by_age = count_male_req_by_age,
                                                                          count_female_req_by_age = count_female_req_by_age,
                                                                          sex_known = sex_known)

            if(sex_known) {

              cand_sires_fixed_current_iteration <- cand_parents_fixed_current_iteration[cand_parents_fixed_current_iteration %in% ped[ped$Sex == "male", "Indiv"]]
              ub[names(ub) %in% cand_sires_fixed_current_iteration] <- indiv_contbn_males
              lb[names(lb) %in% cand_sires_fixed_current_iteration] <- indiv_contbn_males - 1e-10

              cand_dams_fixed_current_iteration <- cand_parents_fixed_current_iteration[cand_parents_fixed_current_iteration %in% ped[ped$Sex == "female", "Indiv"]]
              ub[names(ub) %in% cand_dams_fixed_current_iteration] <- indiv_contbn_females
              lb[names(lb) %in% cand_dams_fixed_current_iteration] <- indiv_contbn_females - 1e-10

            } else {

              ub[names(ub) %in% cand_parents_fixed_current_iteration] <- indiv_contbn_males
              lb[names(lb) %in% cand_parents_fixed_current_iteration] <- indiv_contbn_males - 1e-10

            }
              ub <- ub[cand$phen$Indiv] #reorder
              ub <- ub[names(ub) %in%  names(ub_orig)]

              lb <- lb[cand$phen$Indiv] #reorder
              lb <- lb[names(lb) %in%  names(lb_orig)]

            #check if contributions within sex-age classes are met

            if(sex_known) {
              for (age in length(sire_count_by_age):1) {
                tmp <- cand$phen[(max(cand$phen$Born) - cand$phen$Born + 1) == age & cand$phen$Sex == "male", "Indiv"]
                if(sum(tmp %in% c(cand_parents_fixed_current_iteration, cand_parents_fixed_past_iteration)) ==
                   sire_count_by_age[length(sire_count_by_age) - age + 1]) {# reached max in age-sex class
                  ub[names(ub) %in% tmp] <- lb[names(lb) %in% tmp] + 1e-10 #fix upper bound to lower bound value
                }

                tmp <- cand$phen[(max(cand$phen$Born) - cand$phen$Born + 1) == age & cand$phen$Sex == "female", "Indiv"]
                if(sum(tmp %in% c(cand_parents_fixed_current_iteration, cand_parents_fixed_past_iteration)) ==
                   dam_count_by_age[length(dam_count_by_age) - age + 1]) {# reached max in age-sex class
                  ub[names(ub) %in% tmp] <- lb[names(lb) %in% tmp] + 1e-10 #fix upper bound to lower bound value
                }
              }
            } else {
              for (age in length(sire_count_by_age):1) {
                tmp <- cand$phen[(max(cand$phen$Born) - cand$phen$Born + 1) == age, "Indiv"]
                if(sum(tmp %in% c(cand_parents_fixed_current_iteration, cand_parents_fixed_past_iteration)) ==
                   sire_count_by_age[length(sire_count_by_age) - age + 1]) {# reached max in age class
                  ub[names(ub) %in% tmp] <- lb[names(lb) %in% tmp] + 1e-10 #fix upper bound to lower bound value
                }
              }
            }

              #check if contributions within families are met
              for (fam in unique(fit$parent$FAM)) {
                tmp <- cand$phen[cand$phen$FAM == fam, "Indiv"]
                if(length(tmp %in% c(cand_parents_fixed_current_iteration, cand_parents_fixed_past_iteration)) >= fam ) {# reached max in family
                  ub[names(ub) %in% tmp] <- lb[names(lb) %in% tmp] + 1e-10 #fix upper bound to lower bound value
                }
              }

              ub <- ub[cand$phen$Indiv] #reorder
              ub <- ub[names(ub) %in%  names(ub_orig)]

              lb <- lb[cand$phen$Indiv] #reorder
              lb <- lb[names(lb) %in%  names(lb_orig)]

              #exclude families that have reached max_parents_per_fam
              fam_counts <- ped[ped$Indiv %in% cand_parents_fixed_past_iteration,c("Indiv","FAM")]
              if(nrow(fam_counts) > 0) {
                fam_counts <- aggregate(Indiv ~ FAM, data = fam_counts, FUN = length)
                families_below_threshold <- fam_counts$FAM[fam_counts$Indiv < max_parents_per_fam]
                fams_to_retain <- fams_to_retain[fams_to_retain %in% families_below_threshold]
              }
              rm(fam_counts, families_below_threshold)




          }
        }
        #if(adj != 0) {
        if(adj > (0-step_interval)) {
          try(fit <- OCFam::run_OC_max_EBV(cand = cand, kinship_constraint = kinship_constraint, ub = ub, lb = lb, opticont_method = opticont_method))
          final_adj <- adj
        }

        cand_parents_fixed_past_iteration <- cand_parents_fixed_past_iteration[!cand_parents_fixed_past_iteration %in% cand_parents_fixed_current_iteration] #this shouldn't make any difference but sometimes gets confused

      }
    }

    #if already have enough parents then break after one additional iteration
    if(length(c(cand_parents_fixed_current_iteration, cand_parents_fixed_past_iteration )) >= (0.5 / (indiv_contbn_males) + 0.5 / (indiv_contbn_females)) &
       length(c(cand_parents_fixed_current_iteration, cand_parents_fixed_past_iteration )) == prev_count) {break}
    #if not identifying additional parents then break
    #     if(length(c(cand_parents_fixed_current_iteration, cand_parents_fixed_past_iteration )) == prev_count) {break}

  }

  print("FINAL ITERATION")
  print(paste("Total count of parents fixed in previous iteration =",length(c(cand_parents_fixed_current_iteration, cand_parents_fixed_past_iteration )), "of", 0.5 / (indiv_contbn_males) + 0.5 / (indiv_contbn_females)))
  print(paste("indiv_contbn_males (0.5 / N_sires) = ",indiv_contbn_males))
  print(paste("indiv_contbn_females (0.5 / N_dams) = ",indiv_contbn_females))
  print(paste("Iteration = ", iterations - final_adj/step_interval + 1))
  print(paste("Threshold adjustment (final_adj) = ",final_adj))
  print(paste("Threshold sum of individual contributions for fixing an individual from a family = ", min(c(indiv_contbn_males,indiv_contbn_females)) * final_adj))

  fit$final.step <-  list(
    c("Total count of parents fixed in previous iteration =",length(c(cand_parents_fixed_current_iteration, cand_parents_fixed_past_iteration )), "of", 0.5 / (indiv_contbn_males) + 0.5 / (indiv_contbn_females)),
    c("indiv_contbn_males (0.5 / N_sires) = ",indiv_contbn_males),
    c("indiv_contbn_females (0.5 / N_dams) = ",indiv_contbn_females),
    c("Iteration = ", iterations - final_adj/step_interval + 1),
    c("Threshold adjustment (final_adj) = ",final_adj),
    c("Threshold sum of individual contributions for fixing an individual from a family = ", min(c(indiv_contbn_males,indiv_contbn_females)) * final_adj))

  fit_out <- fit[c("info", "summary", "mean", "bc", "obj.fun", "final.step")]

  ################################################################################
  ################################################################################

  cand_parents_fixed <- c(cand_parents_fixed_current_iteration, cand_parents_fixed_past_iteration )

  #need to add the total contributions of all individuals to parents
  if(sex_known) {

    cand_sires_fixed <- cand_parents_fixed[cand_parents_fixed %in% ped[ped$Sex == "male", "Indiv"]]
    tmp1 <- data.frame(Indiv = names(lb)[names(lb) %in% cand_sires_fixed],
                       N_TOTAL_AS_PARENT = lb[names(lb) %in% cand_sires_fixed] / indiv_contbn_males)

    cand_dams_fixed <- cand_parents_fixed[cand_parents_fixed %in% ped[ped$Sex == "female", "Indiv"]]
    tmp2 <- data.frame(Indiv = names(lb)[names(lb) %in% cand_dams_fixed],
                       N_TOTAL_AS_PARENT = lb[names(lb) %in% cand_dams_fixed] / indiv_contbn_females)
    tmp <- rbind(tmp1,tmp2)
    rm(tmp1,tmp2)

  } else {

    indiv_contbn <- indiv_contbn_males
    tmp <- data.frame(Indiv = names(lb)[names(lb) %in% cand_parents_fixed],
                      N_TOTAL_AS_PARENT = lb[names(lb) %in% cand_parents_fixed] / indiv_contbn)
    rm(indiv_contbn)

  }

  if(length(cand_parents_fixed_current_iteration) > 0) {
    #  tmp[tmp$Indiv  %in% cand_parents_fixed_current_iteration,"N_TOTAL_AS_PARENT"] <- 1 #lb not updated
    tmp[tmp$Indiv  %in% cand_parents_fixed_current_iteration,"ADDED_IN_LAST_ITERATION"] <- TRUE
  } else {
    tmp$ADDED_IN_LAST_ITERATION <- FALSE
  }
  tmp[is.na(tmp$ADDED_IN_LAST_ITERATION), "ADDED_IN_LAST_ITERATION"] <- FALSE

  parent_contbn <- ped[ped$Indiv %in% cand_parents_fixed,]
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

#    if(nrow(cand$classes[cand$classes$Group == "cohort"]) == 1 |
#       nrow(cand$classes[cand$classes$Group != "cohort"]) == 2) { #overlapping = FALSE
      con <- list(ub.fPED = kinship_constraint,
                  ub = ub,
                  lb = lb)
#    } else {
#
#      con <- list(ub.fPED = kinship_constraint,
#                  ub = ub,
#                  lb = lb,
#                  ub.oldest_age_r_vector = cand$classes[cand$classes$age == max(cand$classes$age)  ,"cont0"],
#                  lb.oldest_age_r_vector_2 = cand$classes[cand$classes$age == max(cand$classes$age)  ,"cont0"]) #fix to r vector value of the oldest year class
#    }

  }

  if(opticont_method == "min.fPED") {

#    if(nrow(cand$classes[cand$classes$Group == "cohort"]) == 1 |
#       nrow(cand$classes[cand$classes$Group != "cohort"]) == 2) { #overlapping = FALSE
      con <- list(ub = ub,
                  lb = lb)
#    } else {
#
#      con <- list(ub = ub,
#                  lb = lb,
#                  ub.oldest_age_r_vector = cand$classes[cand$classes$age == max(cand$classes$age)  ,"cont0"],
#                  lb.oldest_age_r_vector_2 = cand$classes[cand$classes$age == max(cand$classes$age)  ,"cont0"]) #fix to the contribution of the oldest year class
#    }

  }

  fit <- optiSel::opticont(method = opticont_method, cand = cand, con = con,
                           solver="cccp") #cccp2 tends to start iterations again after finding a solution
  #  }

  return(fit)
}


################################################################################
#identify best individual from fams_to_retain
################################################################################


#' @export
get_best_indiv <- function(fish,
                           additional_fams_to_retain,
                           candidates,
                           count_male_req_by_age,
                           count_female_req_by_age,
                           sex_known)  {

  tmp <- fish[order(fish$oc, decreasing = TRUE),]
  tmp <- tmp[tmp$FAM %in% additional_fams_to_retain &
               tmp$Indiv %in% candidates,]
  cand_parents_fixed <- tmp[!duplicated(tmp$FAM),]  #Best individual (highest oc) per family independent of sex!!!!!!!
  cand_parents_fixed <- cand_parents_fixed[,"Indiv"]

  cand_rows <- fish$Indiv %in% candidates
  candidates_sex_age   <- fish[cand_rows, c("Indiv", "Sex", "Age")]

  if(sex_known) {

    males_fixed <- NULL
    for(age in length(count_male_req_by_age):1) {
      males_by_age <- candidates_sex_age[candidates_sex_age[,"Sex"] == "male" &
                                           !is.na(candidates_sex_age[,"Sex"]) &
                                           candidates_sex_age[,"Age"] == age , "Indiv"]
      males_by_age <- cand_parents_fixed[cand_parents_fixed %in% males_by_age]
      if(min(count_male_req_by_age[length(count_male_req_by_age) - age + 1],
             length(males_by_age)) <= 0) {
        males_by_age <- NULL
      } else {
        males_by_age <- males_by_age[1:min(count_male_req_by_age[length(count_male_req_by_age) - age + 1],
                                           length(males_by_age))]
      }
      males_fixed <- c(males_fixed, males_by_age)
    }

    females_fixed <- NULL
    for(age in length(count_female_req_by_age):1) {
      females_by_age <- candidates_sex_age[candidates_sex_age[,"Sex"] == "female" &
                                             !is.na(candidates_sex_age[,"Sex"]) &
                                             candidates_sex_age[,"Age"] == age , "Indiv"]
      females_by_age <- cand_parents_fixed[cand_parents_fixed %in% females_by_age]
      if(min(count_female_req_by_age[length(count_female_req_by_age) - age + 1],
             length(females_by_age)) <= 0) {
        females_by_age <- NULL
      } else {
        females_by_age <- females_by_age[1:min(count_female_req_by_age[length(count_female_req_by_age) - age + 1],
                                               length(females_by_age))]
      }
      females_fixed <- c(females_fixed, females_by_age)
    }

    cand_parents_fixed <- c(males_fixed, females_fixed)

  } else {

    fixed <- NULL
    for(age in length(count_male_req_by_age):1) {
      by_age <- candidates_sex_age[candidates_sex_age[,"Age"] == age , "Indiv"]
      by_age <- cand_parents_fixed[cand_parents_fixed %in% by_age]
      if(round(min(count_male_req_by_age[length(count_male_req_by_age) - age + 1],
                   length(by_age))) <= 0) {
        by_age <- NULL
      } else {
        by_age <- by_age[1:round(min(count_male_req_by_age[length(count_male_req_by_age) - age + 1],
                                     length(by_age)))]
      }
      fixed <- c(fixed, by_age)
    }

    cand_parents_fixed <- fixed
  }

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
get_lb_ub <- function(ped, indiv_contbn_males, indiv_contbn_females, max_parents_per_fam, sex_known, fast = TRUE) {


  candidate_parents <- ped[ ped$AVAIL_OR_PAST_PARENT,]

  include <- candidate_parents[candidate_parents$AVAIL_PARENT,]
  AVAIL_PARENT <- include$Indiv
  include <- include[order(include$RANK, decreasing = FALSE),]
  past_parent <- candidate_parents[!candidate_parents$AVAIL_PARENT & candidate_parents$AVAIL_OR_PAST_PARENT,] #Past parent
  include <- rbind(past_parent, include) #past parent at top - can't change past parents
  if(fast) {
  include <- by(include, include["FAM"], head, n=max_parents_per_fam)  #exclude based on max_parents_per_fam
  include <- Reduce(rbind, include)
  } else {
      include <- by(include, include["FAM"], head, n=10000)  #large number instead of max_parents_per_fam
      include <- Reduce(rbind, include)
    }
  include <- include[,"Indiv"]

  exlcude <- AVAIL_PARENT[!AVAIL_PARENT %in% include] #exclude from AVAIL_PARENT only - can't change past parents
  rm(include, AVAIL_PARENT, past_parent)

  candidate_parents <- candidate_parents[, c("Indiv", "N_AS_PARENT_PREV")]

  candidate_parents[candidate_parents$N_AS_PARENT_PREV == 0,"N_AS_PARENT_PREV"] <- NA #if 0 make NA for following lines

  if(sex_known){
    #males
    candidate_sires <- candidate_parents[candidate_parents$Indiv %in% ped[ped$Sex == "male", "Indiv"],]
    candidate_sires$lb <- 0
    candidate_sires[!is.na(candidate_sires$N_AS_PARENT_PREV) & candidate_sires$SEX == "male","lb"] <-
      candidate_sires[!is.na(candidate_sires$N_AS_PARENT_PREV) & candidate_sires$SEX == "male","N_AS_PARENT_PREV"] * indiv_contbn_males - 1e-10

    candidate_sires$ub <- indiv_contbn_males
    candidate_sires[!is.na(candidate_sires$N_AS_PARENT_PREV) & candidate_sires$SEX == "male","ub"] <-
      candidate_sires[!is.na(candidate_sires$N_AS_PARENT_PREV) & candidate_sires$SEX == "male","N_AS_PARENT_PREV"] * indiv_contbn_males

    #females
    candidate_dams <- candidate_parents[candidate_parents$Indiv %in% ped[ped$Sex == "female", "Indiv"],]
    candidate_dams$lb <- 0
    candidate_dams[!is.na(candidate_dams$N_AS_PARENT_PREV) & candidate_dams$SEX == "female","lb"] <-
      candidate_dams[!is.na(candidate_dams$N_AS_PARENT_PREV) & candidate_dams$SEX == "female","N_AS_PARENT_PREV"] * indiv_contbn_females - 1e-10

    candidate_dams$ub <- indiv_contbn_females
    candidate_dams[!is.na(candidate_dams$N_AS_PARENT_PREV) & candidate_dams$SEX == "female","ub"] <-
      candidate_dams[!is.na(candidate_dams$N_AS_PARENT_PREV) & candidate_dams$SEX == "female","N_AS_PARENT_PREV"] * indiv_contbn_females

    candidate_parents <- rbind(candidate_sires, candidate_dams)
  } else {

    indiv_contbn <- indiv_contbn_males

    candidate_parents$lb <- 0
    candidate_parents[!is.na(candidate_parents$N_AS_PARENT_PREV),"lb"] <-
      candidate_parents[!is.na(candidate_parents$N_AS_PARENT_PREV),"N_AS_PARENT_PREV"] * indiv_contbn - 1e-10

    candidate_parents$ub <- indiv_contbn
    candidate_parents[!is.na(candidate_parents$N_AS_PARENT_PREV),"ub"] <-
      candidate_parents[!is.na(candidate_parents$N_AS_PARENT_PREV),"N_AS_PARENT_PREV"] * indiv_contbn

    rm(indiv_contbn)

  }

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

#' Internal: basic checks on 'ped' for OCFam()
#'
#' Ensures that the pedigree has all required columns and basic
#' structure/typing for OCFam to work safely.
#' @export
check_OCFam_ped <- function(ped) {

  # Must be a data.frame
  if (!is.data.frame(ped)) {
    stop("'ped' must be a data.frame.")
  }

  # Required columns for OCFam
  required_cols <- c(
    "INDIV", "SIRE", "DAM", "FAM",
    "SEX", "BORN", "EBV",
    "AVAIL_PARENT", "N_AS_PARENT_PREV"
  )

  missing_cols <- setdiff(required_cols, colnames(ped))
  if (length(missing_cols) > 0) {
    stop("The following required columns are missing from 'ped': ",
         paste(missing_cols, collapse = ", "))
  }

  ## ID / family columns -------------------------------------------------------
  # Coerce to character to avoid factor issues
  ped$INDIV <- as.character(ped$INDIV)
  ped$SIRE  <- as.character(ped$SIRE)
  ped$DAM   <- as.character(ped$DAM)
  ped$FAM   <- as.character(ped$FAM)
  ped$SEX   <- as.character(ped$SEX)

  # INDIV must be non-missing and non-empty
  if (any(is.na(ped$INDIV) | ped$INDIV == "")) {
    stop("'INDIV' must not contain NA or empty strings.")
  }

  # Prefer unique IDs (this is usually essential for a pedigree)
  if (any(duplicated(ped$INDIV))) {
    stop("'INDIV' contains duplicated identifiers.")
  }

  ## BORN ----------------------------------------------------------------------
  if (!is.numeric(ped$BORN)) {
    stop("'BORN' must be numeric (integer or numeric).")
  }

  # It is OK for some BORN to be NA for ancestors,
  # but not for candidate parents (checked later in OCFam).
  # Here we just forbid *all* being NA, which would be useless.
  if (all(is.na(ped$BORN))) {
    stop("All 'BORN' values in 'ped' are NA.")
  }

  ## EBV -----------------------------------------------------------------------
  if (!is.numeric(ped$EBV)) {
    stop("'EBV' must be numeric (NA allowed).")
  }

  ## AVAIL_PARENT ---------------------------------------------------------------
  # Accept logical; otherwise try to convert 0/1 to logical with a clear error.
  if (!is.logical(ped$AVAIL_PARENT)) {
    # Common case: 0/1 coding
    if (is.numeric(ped$AVAIL_PARENT) &&
        all(ped$AVAIL_PARENT %in% c(0, 1, NA))) {
      ped$AVAIL_PARENT <- ped$AVAIL_PARENT == 1
    } else {
      stop("'AVAIL_PARENT' must be logical or numeric 0/1.")
    }
  }

  ## N_AS_PARENT_PREV ----------------------------------------------------------
  if (!is.numeric(ped$N_AS_PARENT_PREV)) {
    stop("'N_AS_PARENT_PREV' must be numeric/integer.")
  }

  if (any(ped$N_AS_PARENT_PREV < 0, na.rm = TRUE)) {
    stop("'N_AS_PARENT_PREV' must be >= 0.")
  }

  ## Sanity checks on content --------------------------------------------------
  # At least one candidate parent available in this round
  n_cand <- sum(ped$AVAIL_PARENT, na.rm = TRUE)
  if (n_cand == 0) {
    stop("There are no candidate parents: all 'AVAIL_PARENT' are FALSE/NA.")
  }

  # At least some defined family IDs somewhere
  if (all(is.na(ped$FAM) | ped$FAM == "")) {
    stop("No valid family identifiers found in 'FAM'.")
  }

  ## -------------------------------------------------------------------
  ## CHECK: Sex consistency with historical parental roles (Sire/Dam)
  ## -------------------------------------------------------------------

  sex_vec <- ped$SEX

  ## Rule 1: Either all NA OR all defined
  if (all(is.na(sex_vec))) {
    ## Sex-unknown scenario (allowed) — no parent-role checks possible
    warning("All SEX values are NA: operating in sex-unknown mode. Cannot verify consistency of historical Sire/Dam roles.",
            call. = FALSE)

  } else if (any(is.na(sex_vec))) {
    ## Mixed case — not allowed
    stop("Some individuals have SEX = NA while others have defined SEX. ",
         "Either define SEX for all individuals, or set all SEX = NA.",
         call. = FALSE)

  } else {

    ## Sex is fully defined — perform Sire and Dam consistency checks

    ## Create lookup: Indiv ? Sex
    sex_by_id <- ped$SEX
    names(sex_by_id) <- ped$INDIV

    ## ---------- Check Sires ----------
    sire_ids <- unique(ped$SIRE)
    sire_ids <- sire_ids[!is.na(sire_ids)]
    sire_ids <- sire_ids[sire_ids %in% names(sex_by_id)]

    sire_sex <- sex_by_id[as.character(sire_ids)]
    bad_sires <- sire_ids[sire_sex != "male"]

    ## ---------- Check Dams ----------
    dam_ids <- unique(ped$DAM)
    dam_ids <- dam_ids[!is.na(dam_ids)]
    dam_ids <- dam_ids[dam_ids %in% names(sex_by_id)]

    dam_sex <- sex_by_id[as.character(dam_ids)]
    bad_dams <- dam_ids[dam_sex != "female"]

    ## ---------- Issue warnings (but do NOT stop) ----------
    if (length(bad_sires) > 0 || length(bad_dams) > 0) {

      msg <- "SEX conflicts detected with historical parental roles:"

      if (length(bad_sires) > 0) {
        msg <- paste0(msg,
                      "\n  - ", length(bad_sires),
                      " individual(s) used as Sire but SEX != 'male' ",
                      "(e.g. ", paste(head(bad_sires, 5), collapse = ", "),
                      if (length(bad_sires) > 5) ", ..." else "",
                      ")."
        )
      }

      if (length(bad_dams) > 0) {
        msg <- paste0(msg,
                      "\n  - ", length(bad_dams),
                      " individual(s) used as Dam but SEX != 'female' ",
                      "(e.g. ", paste(head(bad_dams, 5), collapse = ", "),
                      if (length(bad_dams) > 5) ", ..." else "",
                      ")."
        )
      }

      stop(msg, call. = FALSE)
    }
  }


  invisible(ped)
}

#' Internal: basic checks on scalar arguments for OCFam()
#' @export
check_OCFam_args <- function(kinship_constraint,
                             iterations,
                             sire_count_by_age,
                             dam_count_by_age,
                             min_prop_fams,
                             max_parents_per_fam) {

  ## kinship_constraint --------------------------------------------------------
  if (!is.null(kinship_constraint) && !is.na(kinship_constraint)) {
    if (length(kinship_constraint) != 1 ||
        !is.numeric(kinship_constraint) ||
        kinship_constraint < 0 || kinship_constraint > 1) {
      stop("'kinship_constraint' must be a single numeric in [0,1] or NA.")
    }
  }

  ## iterations -------------------------------------------------------------
  if (length(iterations) != 1 || !is.numeric(iterations) ||
      length(iterations) != 1 || !isTRUE(all.equal(iterations, round(iterations))) ||
      is.na(iterations) || iterations <= 0) {
    stop("'iterations' must be a postive integer")
  }

  ## sire_count_by_age and dam_count_by_age ----------------------------------------------------------
  if (sum(is.na(sire_count_by_age)) > 0 | sum(is.na(dam_count_by_age)) > 0 ) {
    stop("'sire_count_by_age' and 'dam_count_by_age' must not contain NA")
  }

  is_valid_count_vector <- function(x, tol = 1e-8) {
    is.numeric(x) &&
      !any(is.na(x)) &&
      all(abs(x - round(x)) < tol) &&   # integer-valued
      all(x >= 0)                       # non-negative
  }

  if (!is_valid_count_vector(sire_count_by_age) ||
      !is_valid_count_vector(dam_count_by_age)) {

    stop("'sire_count_by_age' and 'dam_count_by_age' must be numeric, integer-valued (scalar or vector), and all values must be >= 0.")
  }

# if (sum(sire_count_by_age) != sum(dam_count_by_age)) {
#    stop("'sire_count_by_age' must sum to the same value as 'dam_count_by_age'")
#  }

  if (length(sire_count_by_age) != length(dam_count_by_age)) {
    stop("'sire_count_by_age' and 'dam_count_by_age' must be the same length")
  }

  ## min_prop_fams -------------------------------------------------------------
  if (length(min_prop_fams) != 1 || !is.numeric(min_prop_fams) ||
      is.na(min_prop_fams) || min_prop_fams < 0 || min_prop_fams > 1) {
    stop("'min_prop_fams' must be a single numeric in (0,1].")
  }

  ## max_parents_per_fam -------------------------------------------------------
  if (length(max_parents_per_fam) != 1 ||
      !is.numeric(max_parents_per_fam) ||
      is.na(max_parents_per_fam) ||
      max_parents_per_fam < 1 ||
      max_parents_per_fam != as.integer(max_parents_per_fam)) {
    stop("'max_parents_per_fam' must be a single positive integer (>= 1).")
  }

  invisible(TRUE)
}

