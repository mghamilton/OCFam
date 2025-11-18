

##########################################################################################
#OCFamPrep
##########################################################################################
#' OCFamPrep
#'
#' @description
#' This functions provides tables and plot that allows analysis of historic trends in average coancestry (k) and inbreeding (Wrights inbreeding coefficient, F)
#'
#' @param ped is a data frame with the following columns (class in parentheses):
#' \itemize{
#'  \item{'INDIV' is the individual identifier (character).}
#'  \item{'SIRE' is the male parent identifier (character). Unknown sires should be entered as '0' or NA.}
#'  \item{'DAM' is the female parent identifier (character). Unknown sires should be entered as '0' or NA.}
#'  \item{'FAM' is a full-sibling family identifier.  May be NA if parents are unknown (character).}
#'  \item{'BORN' integer indicating age class.  May be the year of birth if one age class per year or an integer indicating the sequence of age classes (integer).}
#'  \item{'EBV' is the estimated breeding value. This field may contain missing values or NA (numeric).}
#'  }
#'
#' @param age_class_names is a data frame specifying names for age classes (class in parentheses):
#' \itemize{
#'  \item{'AGE_CLASS_NAME' is a name for the age class (character).}
#'  \item{'BORN' integer indicating age class.  May be the year of birth if one age class per year or an integer indicating the sequence of age classes (integer).}
#'  }
#'
#' @param gene_flow_vector is a applicable to overlapping generations (if NA discrete generations is assumed).  It is vector representing parental contributions by age class to the next age class. For example, gene_flow_vector = c(0.2, 0.8, 0, 0) - oldest age class to youngest age class.
#'
#' @return 'plot' is plot of trends in coancestry and inbreeding coefficients
#' @return 'age_class_means' is a dataframe listing inbreeding coefficients for individuals within families.
#' \itemize{
#'  \item{'BORN' is age class (integer).}
#'  \item{'INBREEDING' is the average Wright's inbreeding coefficient of families in the age class (numeric).}
#'  \item{'COANCESTRY' is the average coancestry between families in the age class (numeric).}
#'  \item{'OVERLAPPING_COANCESTRY' is the average coancestry between families weighted by the age class contributions specified in the gene_flow_vector (numeric).}
#' }
#' @return 'age_class_K_mat' is a matrix containing average coancestries within age classes on the diagonals and between age classes on the off-diagonals.  Row and column names are the FAM identifiers from 'ped'.
#' @return 'fam_K_matrix' is a data frame containing the between family coancestry matrix (half of the Between-Family Relationship Matrix described in Hamilton (2020) 'Optimal Contribution Selection in Highly Fecund Species With Overlapping Generations').  Row and column names are the FAM identifiers from 'ped'.
#' @return 'family_inbreeding' is a dataframe listing inbreeding coefficients for individuals within families.'
#' \itemize{
#'  \item{'FAM' is the family identifier (character).}
#'  \item{'F' is the Wright's inbreeding coefficient for individuals in the family (numeric).}
#' }
#' @return 'ped' is a data frame with the following columns (class in parentheses):
#' \itemize{
#'  \item{'INDIV' is the individual identifier (character).}
#'  \item{'SIRE' is the male parent identifier (character).}
#'  \item{'DAM' is the female parent identifier (character).}
#'  \item{'FAM' is a full-sibling family identifier (character).}
#'  \item{'BORN' integer indicating age class.  May be the year of birth if one age class per year or an integer indicating the sequence of age classes (integer).}
#'  \item{'EBV' is the estimated breeding value (numeric).}
#'  \item{'FAM_SIRE' is a full-sibling family identifier of the SIRE (character).}
#'  \item{'EBV_SIRE' is the estimated breeding value of the SIRE (numeric).}
#'  \item{'FAM_DAM' is a full-sibling family identifi of the DAM (character).}
#'  \item{'EBV_DAM'  is the estimated breeding value of the DAM (numeric).}
#'  \item{'EBV_FAM' is the mean EBV of the FAM's SIRE and DAM}
#'  }
#' @return 'family_ped' is a data frame with the following columns (class in parentheses):
#' \itemize{
#'  \item{'FAM' is a full-sibling family identifier (character).}
#'  \item{'SIRE' is the male parent identifier (character).}
#'  \item{'DAM' is the female parent identifier (character).}
#'  \item{'BORN' integer indicating age class.  May be the year of birth if one age class per year or an integer indicating the sequence of age classes (integer).}
#'  \item{'FAM_SIRE' is a full-sibling family identifier of the SIRE (character).}
#'  \item{'EBV_SIRE' is the estimated breeding value of the SIRE (numeric).}
#'  \item{'FAM_DAM' is a full-sibling family identifi of the DAM (character).}
#'  \item{'EBV_DAM'  is the estimated breeding value of the DAM (numeric).}
#'  \item{'EBV_FAM' is the mean EBV of the FAM's SIRE and DAM}
#'  }
#' @return 'gene_flow_vector' is a applicable to overlapping generations (if NA discrete generations is assumed).  It is vector representing parental contributions by age class to the next age class. For example, gene_flow_vector = c(0.2, 0.8, 0, 0) - oldest age class to youngest age class.
#' @return 'r' is a weight vector of age classes; r[i] denotes the long-term contribution of age class i (until lifetime is reached).  See Meuwissen, T. H. E. and A. K. Sonesson. 1998. Maximizing the response of selection with a predefined rate of inbreeding: overlapping generations.
#' @return 'gen_interval' is the generation interval.

#' @examples
#' #Retrieve example data
#' ped_small <- OCFam::ped_small
#' age_class_names <-  data.frame(AGE_CLASS_NAME = c("Founders", "Gen_0", "Gen_1", "Gen_2"), BORN = c(0, 1, 2, 3))
#' gene_flow_vector <- c(0.2,0.8)
#'
#' #Show example data
#' tail(ped_small)
#' age_class_names
#' gene_flow_vector
#'
#' #Run OCFam function
#' OCFamPrep_output <- OCFam::OCFamPrep(ped = ped_small,
#'                              age_class_names = age_class_names,
#'                              gene_flow_vector = gene_flow_vector
#' )
#'
#' OCFamPrep_output$age_class_means
#' OCFamPrep_output$age_class_K_mat
#' OCFamPrep_output$fam_K_matrix[1:5,1:5]
#' head(OCFamPrep_output$family_inbreeding)
#'
#' @import AGHmatrix
#' @import dplyr
#' @import ggplot2
#' @export


OCFamPrep <- function(ped, age_class_names = NULL, gene_flow_vector = NULL) {
  library(AGHmatrix)
  library(dplyr)
  library(ggplot2)

  ## ==== DATA CHECKS ==========================================================
  check_OCFamPrep_inputs(ped = ped,
                         age_class_names = age_class_names,
                         gene_flow_vector = gene_flow_vector)
  ## ===========================================================================

  if(length(gene_flow_vector) == 1) {gene_flow_vector <- NULL}
  if(length(age_class_names) == 1) {age_class_names <- NULL}

  #replace NA in sires and dams
  ped[is.na(ped$SIRE), "SIRE"] <- NA
  ped[is.na(ped$DAM), "DAM"] <- NA

  #create family_ped
  family_ped <- unique(ped[,c("FAM", "SIRE", "DAM", 'BORN')])
  sires <- ped[,c('INDIV', 'FAM', "EBV")]
  colnames(sires) <- c("SIRE", "FAM_SIRE", "EBV_SIRE")
  family_ped <- dplyr::left_join(family_ped, sires, by = "SIRE")
  rm(sires)
  dams <- ped[,c('INDIV', 'FAM', "EBV")]
  colnames(dams) <- c("DAM", "FAM_DAM", "EBV_DAM")
  family_ped <- dplyr::left_join(family_ped, dams, by = "DAM")
  rm(dams)

  family_ped$EBV_FAM <- (family_ped$EBV_SIRE + family_ped$EBV_DAM)/2

  family_ped[family_ped$FAM == "" | is.na(family_ped$FAM), "FAM"] <- NA
  family_ped <- family_ped[!is.na(family_ped$FAM),]

  ped <- dplyr::left_join(ped, family_ped, by = c("FAM", "SIRE", "DAM", 'BORN'))

  fam_K_matrix <- OCFam_fam_K_matrix(family_ped)

  family_inbreeding <- data.frame(FAM = names(fam_K_matrix$family_inbreeding),
                                  F    = fam_K_matrix$family_inbreeding)
  fam_K_matrix <- fam_K_matrix$K_matrix_families
  fam_K_matrix <- as.data.frame(as.matrix(fam_K_matrix))

  fam_K_matrix$FAM <- rownames(fam_K_matrix)

  ####################################################################################
  #Generate summary data
  ####################################################################################

  age_classs <- unique(family_ped$BORN)

  age_class_K_mat <- matrix(data = NA, nrow = length(age_classs), ncol = length(age_classs))
  rownames(age_class_K_mat) <- age_classs
  colnames(age_class_K_mat) <- age_classs

  age_class_N_fams <- family_ped %>%
    group_by(BORN) %>%
    summarise(N_FAMS = n_distinct(FAM)) %>%
    arrange(BORN)
  age_class_N_fams <- as.data.frame(age_class_N_fams)

  #  age_class_F <- data.frame(BORN = NA, VARIABLE = NA, MEAN = NA)[-1,]
  #  for(ac in age_classs) {
  #    fams <- family_ped[family_ped$BORN == as.character(ac),"FAM"]
  #    tmp_F <- mean(family_inbreeding[family_inbreeding$FAM %in% fams, "F"])
  #    tmp <- data.frame(BORN = ac, VARIABLE = "Inbreeding (F)", MEAN = tmp_F)
  #    age_class_F <- rbind(age_class_F, tmp)
  #  }

  age_class_F <- left_join(family_ped, family_inbreeding, by = "FAM") %>%
    group_by(BORN) %>%
    summarise(MEAN = mean(F, na.rm = TRUE))
  age_class_F$VARIABLE <- "Inbreeding (F)"
  age_class_F <- as.data.frame(age_class_F[,c("BORN", "VARIABLE", "MEAN")])

  #coancestry matrix
  for(ac1 in age_classs) {
    for(ac2 in age_classs) {

      if(ac1 >= ac2){
        fams1 <- family_ped[family_ped$BORN == as.character(ac1),"FAM"]
        fams2 <- family_ped[family_ped$BORN == as.character(ac2),"FAM"]

        tmp_mean <- mean(as.matrix(fam_K_matrix[fams1,fams2]))
        age_class_K_mat[as.character(ac1), as.character(ac2)] <-  tmp_mean
      }
    }
  }
  age_class_K_mat[upper.tri(age_class_K_mat)] <- age_class_K_mat[lower.tri(age_class_K_mat)]
  age_class_K <- data.frame(BORN = rownames(age_class_K_mat), VARIABLE = "Coancestry (k)", MEAN = diag(age_class_K_mat))

  if(is.null(gene_flow_vector)) {
    age_class_means <- rbind(age_class_F, age_class_K)
    r <- 1
    gen_interval <- 1
  }

  age_class_K_overlap <-  data.frame(BORN = age_classs, VARIABLE = "Coancestry (overlapping; k)", MEAN = NA)

  if(!is.null(gene_flow_vector)) {

    gen_interval  <- gene_flow_vector %*% c(length(gene_flow_vector):1)

    r <- rep(NA,length(gene_flow_vector))
    for (i in 1:length(gene_flow_vector)) {
      r[i] <- sum(gene_flow_vector[1:i])/gen_interval
    }

    for(ac in age_classs[length(gene_flow_vector):length(age_classs)]) {
      age_class_K_overlap[age_class_K_overlap$BORN == ac,"MEAN"] <-
        t(r) %*% age_class_K_mat[as.character((ac-length(r)+1):ac), as.character((ac-length(r)+1):ac)] %*% r #r' A_bar r
    }
    age_class_means <- rbind(age_class_F, age_class_K, age_class_K_overlap)
  }

  ####################################################################################
  #Generate plots
  ####################################################################################

  age_class_means$BORN <-  as.numeric(age_class_means$BORN)

  age_class_means <- dplyr::left_join(age_class_means, age_class_N_fams, by = "BORN")

  if(is.null(age_class_names)) {
    p <- ggplot(age_class_means, aes(x=BORN, y=MEAN, group=VARIABLE))
  } else {
    age_class_means <- dplyr::left_join(age_class_means, age_class_names, by = "BORN")
    ped <- dplyr::left_join(ped, age_class_names, by = "BORN")
    family_ped <- dplyr::left_join(family_ped, age_class_names, by = "BORN")

    # Ensure AGE_CLASS_NAME is ordered by BORN
    age_class_means <- age_class_means[order(age_class_means$BORN), ]

    # Reassign AGE_CLASS_NAME as a factor, now ordered by BORN
    age_class_means$AGE_CLASS_NAME <- factor(age_class_means$AGE_CLASS_NAME,
                                             levels = unique(age_class_means$AGE_CLASS_NAME))

    p <- ggplot(age_class_means, aes(x=AGE_CLASS_NAME, y=MEAN, group=VARIABLE))

  }

  # Filter the data for Coancestry
  coancestry_data <- subset(age_class_means, VARIABLE == "Coancestry (k)")

  p <- p +
    geom_line(aes(color = VARIABLE), na.rm = TRUE) +
    geom_point(aes(color = VARIABLE), na.rm = TRUE) +

    # Add N_FAMS text labels above the points for Coancestry
    geom_text(data = coancestry_data,
              aes(label = N_FAMS),
              vjust = -0.5,  # Move text above the points
              color = "black",  # Text color
              size = 3, # Text size
              na.rm = TRUE) +

    # Customizing the theme
    theme_classic() +
    theme( # remove the vertical grid lines
      panel.grid.major.x = element_blank() ,
      # explicitly set the horizontal lines (or they will disappear too)
      panel.grid.major.y = element_line(linewidth=.1, color="grey70" ),
      axis.text.x = element_text(angle = 90, vjust = 0.5)
    )

  print(p)

  age_class_means_out <- data.frame(BORN = age_classs,
                                    N_FAMS = age_class_N_fams$N_FAMS,
                                    INBREEDING = age_class_F$MEAN,
                                    COANCESTRY = age_class_K$MEAN,
                                    OVERLAPPING_COANCESTRY = age_class_K_overlap$MEAN)
  if(!is.null(age_class_names )) {
    age_class_means_out <- dplyr::left_join(age_class_names, age_class_means_out, by = "BORN")
  }

  return(list(plot = p,
              age_class_means = age_class_means_out,
              age_class_K_mat = age_class_K_mat,
              fam_K_matrix = fam_K_matrix,
              family_inbreeding = family_inbreeding,
              ped = ped,
              family_ped = family_ped,
              gene_flow_vector = gene_flow_vector,
              r = r,
              gen_interval = gen_interval))
}


OCFam_fam_K_matrix <- function(family_ped) {
  #data
  # class: data.frame
  # fields:
  #   FAM
  #   SIRE
  #   FAM_SIRE
  #   DAM
  #   FAM_DAM

  #sires
  sires <- matrix(unique(family_ped$SIRE), ncol = 1)
  colnames(sires)[1] <- "SIRE"
  sires <- merge(sires, family_ped[,c("SIRE","FAM_SIRE")], by = "SIRE", all.x = TRUE)
  colnames(sires) <- c("Indiv_id","FAM")
  sires <- merge(sires, family_ped[,c("FAM","SIRE","DAM")], by = "FAM", all.x = TRUE)
  sires <- subset(sires, select=-c(FAM)) #remove FAM column

  #dams
  dams <- matrix(unique(family_ped$DAM), ncol = 1)
  colnames(dams)[1] <- "DAM"
  dams <- merge(dams, family_ped[,c("DAM","FAM_DAM")], by = "DAM", all.x = TRUE)
  colnames(dams) <- c("Indiv_id","FAM")
  dams <- merge(dams, family_ped[,c("FAM","SIRE","DAM")], by = "FAM", all.x = TRUE)
  dams <- subset(dams, select=-c(FAM)) #remove FAM column

  parents <- rbind(sires,dams)
  parents <- unique(parents)
  parents <- parents[order(parents[,"Indiv_id"], decreasing = FALSE), ]

  families            <- unique(family_ped)
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

get.r <- function(gene_flow_vector) {

  #generation interval


  #r
  gene_flow_vector$r <- NA
  gene_flow_vector$age_t_plus_1 <- gene_flow_vector$age - 1
  for(age in max(gene_flow_vector$age):1) {
    gene_flow_vector[gene_flow_vector$sex == "male" & gene_flow_vector$age == age,"r"] <- sum(gene_flow_vector[gene_flow_vector$sex == "male" & gene_flow_vector$age %in% (max(gene_flow_vector$age):age),"s_1_plus_s_2"]) / Lbar_male
  }
  for(age in max(gene_flow_vector$age):1) {
    gene_flow_vector[gene_flow_vector$sex == "female" & gene_flow_vector$age == age,"r"] <- sum(gene_flow_vector[gene_flow_vector$sex == "female" & gene_flow_vector$age %in% (max(gene_flow_vector$age):age),"s_1_plus_s_2"]) / Lbar_female
  }

  gene_flow_vector <- gene_flow_vector

  return(list(n_fams = n_fams,
              Lbar_male = Lbar_male,
              Lbar_female = Lbar_female,
              gene_flow_vector = gene_flow_vector))
}

#' Internal: basic checks on inputs for OCFamPrep()
#'
#' Ensures that the input pedigree and optional age-class and gene-flow
#' information are structurally consistent with what OCFamPrep expects.
#' @keywords internal
check_OCFamPrep_inputs <- function(ped, age_class_names, gene_flow_vector) {

  ## ped -----------------------------------------------------------------------
  if (!is.data.frame(ped)) {
    stop("'ped' must be a data.frame.")
  }

  required_cols <- c("INDIV", "SIRE", "DAM", "FAM", "BORN", "EBV")
  missing_cols  <- setdiff(required_cols, colnames(ped))
  if (length(missing_cols) > 0) {
    stop("The following required columns are missing from 'ped': ",
         paste(missing_cols, collapse = ", "))
  }

  # Coerce to character where appropriate
  ped$INDIV <- as.character(ped$INDIV)
  ped$SIRE  <- as.character(ped$SIRE)
  ped$DAM   <- as.character(ped$DAM)
  ped$FAM   <- as.character(ped$FAM)

  if (any(is.na(ped$INDIV) | ped$INDIV == "")) {
    stop("'INDIV' must not contain NA or empty strings.")
  }

  if (!is.numeric(ped$BORN)) {
    stop("'BORN' in 'ped' must be numeric (integer or numeric).")
  }
  if (all(is.na(ped$BORN))) {
    stop("All 'BORN' values in 'ped' are NA.")
  }

  if (!is.numeric(ped$EBV)) {
    stop("'EBV' in 'ped' must be numeric (NA allowed).")
  }

  # At least some families present
  if (all(is.na(ped$FAM) | ped$FAM == "")) {
    stop("No valid family identifiers found in 'ped$FAM'.")
  }

  ## age_class_names -----------------------------------------------------------
  if (!is.null(age_class_names)) {

    if (!is.data.frame(age_class_names)) {
      stop("'age_class_names' must be a data.frame when supplied.")
    }

    required_ac_cols <- c("BORN", "AGE_CLASS_NAME")
    missing_ac <- setdiff(required_ac_cols, colnames(age_class_names))
    if (length(missing_ac) > 0) {
      stop("The following required columns are missing from 'age_class_names': ",
           paste(missing_ac, collapse = ", "))
    }

    if (!is.numeric(age_class_names$BORN)) {
      stop("'BORN' in 'age_class_names' must be numeric.")
    }
    if (!is.character(age_class_names$AGE_CLASS_NAME)) {
      stop("'AGE_CLASS_NAME' in 'age_class_names' must be character.")
    }

    # All BORN values in age_class_names should appear in ped$BORN
    if (any(!age_class_names$BORN %in% ped$BORN)) {
      stop("All 'BORN' values in 'age_class_names' must also appear in 'ped$BORN'.")
    }
  }

  ## gene_flow_vector ----------------------------------------------------------
  if (!is.null(gene_flow_vector)) {

    # Allow a single NA as "no overlapping gens" marker
    if (length(gene_flow_vector) == 1 && is.na(gene_flow_vector[1])) {
      return(invisible(TRUE))
    }

    if (!is.numeric(gene_flow_vector)) {
      stop("'gene_flow_vector' must be numeric, NULL, or a single NA.")
    }

    if (any(is.na(gene_flow_vector))) {
      stop("'gene_flow_vector' must not contain NA, except the single-NA case.")
    }

    if (any(gene_flow_vector < 0)) {
      stop("'gene_flow_vector' must have non-negative entries.")
    }

    if (sum(gene_flow_vector) <= 0) {
      stop("'gene_flow_vector' must have a positive sum.")
    }

    n_age_classes <- length(unique(ped$BORN))
    if (length(gene_flow_vector) > n_age_classes) {
      stop("'gene_flow_vector' length (", length(gene_flow_vector),
           ") is greater than the number of age classes in 'ped' (",
           n_age_classes, ").")
    }
  }

  invisible(TRUE)
}
