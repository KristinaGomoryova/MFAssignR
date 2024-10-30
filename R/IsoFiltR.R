c13_mass <- 1.0033548380
ch2_mass <- 14.01565 # replaced in a code
s34_mass <- 1.995797 # replaced in a code
KMDr_34S_int <- 12
KMDr_13C_int <- 21

#' Filters and arranges data based on abundance and mass
#'
#' This function renames columns and filters data based on a specified abundance threshold.
#' It arranges the data based on the experimental mass.
#'
#' @param peaks data frame:
#' The input data frame containing at least two columns, where the first column is the experimental mass,
#' and the second column is the abundance.
#' @param sn numeric:
#' The abundance threshold. Peaks with abundance below this threshold will be filtered out.
#' @return data frame:
#' A filtered and arranged data frame with columns "Exp_mass" and "Abundance".
#' @export
filter_and_arrange_data <- function(peaks, sn) {
  names(peaks)[1:2] <- c("Exp_mass", "Abundance")
  if (ncol(peaks) == 3) {
    names(peaks)[3] <- "RT"
  } else if (ncol(peaks) == 2) {
    peaks$RT <- 1:nrow(peaks)
  }

  data1 <- peaks |>
    dplyr::filter(Abundance >= sn) |>
    dplyr::arrange(Exp_mass)

  return(data1)
}

#' Creates data chunks
#'
#' This function takes a data frame and divides it into chunks.
#'
#' @param data1 data frame:
#' The input data frame to be divided into chunks.
#' @return list:
#' A list containing two elements: "Carblist" and "Sulflist," each representing a list of data chunks.
#' @export
create_data_chunks <- function(data1) {
  sect <- ceiling(ceiling(nrow(data1))/10)
  over <- ceiling(round(sect) * 0.15)

  data_chunks <- list()
  for (i in 1:10) {
    start_idx <- max((i-1)*sect, 1)
    end_idx <- min((i * sect + over), nrow(data1))
    data_chunks[[i]] <- data1[start_idx:end_idx, ]
  }

  return(list(Carblist = data_chunks, Sulflist = data_chunks))
}

#' Binds and filters two data frames
#'
#' This function combines two data frames and removes duplicate rows. Optionally, it can remove rows with missing
#' experimental masses or masses less than or equal to zero.
#'
#' @param data1 data frame:
#' The first data frame to be combined.
#' @param data2 data frame:
#' The second data frame to be combined.
#' @param remove_na logical:
#' If TRUE, remove rows with missing experimental masses or masses less than zero.
#' @return data frame:
#'   The combined and filtered data frame.
#' @export
bind_and_filter_data <- function(data1, data2, remove_na) {
  result <- rbind(data1, data2)

  if (remove_na) {
    result <- result[!is.na(result$Exp_mass) & result$Exp_mass > 0, ]
  }

  result <- unique(result)

  return(result)
}

#' Combines isotopes
#' This function combines isotopes based on experimental mass.
#' @param Isotopes data frame:
#' The input data frame containing isotopes.
#' @return data frame:
#'  The combined data frame.
combine_isotopes <- function(Isotopes) {
  Isotopes <- Isotopes[Isotopes$Exp_mass > 0, ]
  Isotopes$order <- 1:nrow(Isotopes)

  Isotopes <- Isotopes[order(Isotopes$order), ]
  Isotopes$Dups1 <- duplicated(Isotopes$Exp_mass)

  Isotopes <- Isotopes[order(-Isotopes$order), ]
  Isotopes$Dups2 <- duplicated(Isotopes$Exp_mass)

  Doubles <- Isotopes[(Isotopes$Dups1 == TRUE | Isotopes$Dups2 == TRUE), ]
  DS34 <- Doubles[Doubles$Tag == "S34", ]
  DC13 <- Doubles[Doubles$Tag != "S34", ]

  NewDub <- merge(DS34, DC13, by = "Exp_mass")
  NewDub <- NewDub[c(1, 2, 3, 4, 10)]
  NewDub$Tag <- paste(NewDub$Tag.y, NewDub$Tag.x, sep = "_")

  NewDub <- NewDub[c(1, 2, 3, 6)]

  names(NewDub)[2:3] <- c("Abundance", "RT")

  Singles <- Isotopes[(Isotopes$Dups1 == FALSE & Isotopes$Dups2 == FALSE), ]
  Singles <- Singles[c(1, 2, 3, 4)]

  Iso_Out <- rbind(Singles, NewDub)
  return(Iso_Out)
}

#' Identifies and separates likely isotopic masses from monoisotopic masses
#'
#' IsoFiltR separates likely isotopic masses from monoisotopic masses in a
#' mass list. It can identify likely 13C and 34S isotopic
#' masses and put them in a separate mass list from the monoisotopic
#' masses. This should be done prior to formula assignment in order
#' to lessen the chances of incorrectly assigned formulas, where, for
#' example a 13C containing CHO formula can be assigned to a
#' monoistopic CHNOS formula.
#'
#' The only necessary input is a two column data frame with the abundance in the
#' first column and the measured ion mass in the second column. This should be the
#' raw mass list output.
#'
#' The output of this function is a list of dataframes. Dataframe 1 contains the flagged
#' monoisotopic masses, and all masses that did not have a matching isotopic mass. The
#' second dataframe contains the masses flagged as isotopic. These two dataframes should next be run in the
#'  \code{\link{MFAssignCHO}} or \code{\link{MFAssign}} to assign formulas to the masses.
#'
#' Note that the classification of isotopic or monoisotopic from this function is not
#' definitive.
#'
#' @param peaks data frame:
#' The input 2 column data frame containing abundance and peak mass
#' @param SN numeric:
#' Sets the noise cut for the data, peaks below this value will not be evaluated
#' @param Carbrat numeric:
#' Sets the maximum 13C/12C ratio that is allowed for matching, default is 60
#' @param Sulfrat numeric:
#' Sets the maximum 34S/32S ratio that is allowed for matching, default is 30
#' @param Sulferr numeric:
#' Sets the maximum allowed error (ppm) for 34S mass matching, default is 5
#' @param Carberr numeric:
#' Sets the maximum allowed error (ppm) for 13C mass matching, default is 5
#' @return list(Monolist, Isolist):
#'   Monolist - monoistopic and non-matched masses,
#'   Isolist - isotopic masses
#' @export
IsoFiltR <- function(
    peaks,
    SN = 0,
    Carbrat = 60,
    Sulfrat = 30,
    Sulferr = 5,
    Carberr = 5) {
  cols <- ncol(peaks)

  data1 <- filter_and_arrange_data(peaks, SN)

  data_chunks <- create_data_chunks(data1)

  # Data frame set up
  End <- data1

  # Carbon Isotoping
  IsoOutC1_final <- data.frame(Exp_mass = -42, Abundance = -42, RT = -42)
  MonoOutC_final <- data.frame(Exp_mass = -42, Abundance = -42, RT = -42)
  IsoOutC2_final <- data.frame(Exp_mass = -42, Abundance = -42, RT = -42)

  for (i in 1:10) {
    data <- data_chunks$Carblist[[i]]

    result <- process_carbon_isotoping(data, Carberr, Carbrat, End)

    MonoOutC_final <- bind_and_filter_data(MonoOutC_final, result$MonoC, TRUE)
    IsoOutC1_final <- bind_and_filter_data(IsoOutC1_final, result$IsoC1, TRUE)
    IsoOutC2_final <- bind_and_filter_data(IsoOutC2_final, result$IsoC2, TRUE)
  }

  # Sulfur Section
  IsoOutS_final <- data.frame(Exp_mass = -42, Abundance = -42, RT = -42)
  MonoOutS_final <- data.frame(Exp_mass = -42, Abundance = -42, RT = -42)

  for (i in 1:10) {
    data <- data_chunks$Sulflist[[i]]
    data_end <- data

    result <- process_sulfur(data, Sulferr, data_end, Sulfrat)

    IsoOutS_final <- bind_and_filter_data(IsoOutS_final, result$IsooutS, FALSE)
    MonoOutS_final <- bind_and_filter_data(MonoOutS_final, result$MonooutS, FALSE)
  }

  # End of Sulfur Section
  IsoOutC1_final <- unique(IsoOutC1_final)
  IsoOutC2_final <- unique(IsoOutC2_final)
  MonoOutC_final <- unique(MonoOutC_final)

  IsoOutC1_final$Tag <- "C13"
  IsoOutC2_final$Tag <- "2C13"
  MonoOutC_final$Tag <- "C"
  MonoOutS_final$Tag <- "S"
  IsoOutS_final$Tag <- "S34"

  # Sulfur Check
  sulfur_check_result <- sulfur_check(
    MonoOutC_final,
    MonoOutS_final,
    IsoOutC1_final,
    IsoOutC2_final
  )

  # Set up the isotopes
  Isotopes <- rbind(
    sulfur_check_result$IsoOutC1_final,
    sulfur_check_result$IsoOutC2_final,
    IsoOutS_final
  )

  Iso_Out <- tidyr::drop_na(combine_isotopes(Isotopes))

  # This section organizes the data and helps to ensure there are no peaks being considered as
  # both monoisotopic and isotopic peaks.
  Mono_out <- rbind(MonoOutC_final, sulfur_check_result$MonoOutS_final)

  Dup_mass <- merge(Iso_Out, Mono_out, by = "Exp_mass")

  Dup_rem <- Dup_mass[c(1, 2, 3)]
  Iso_Sup <- Dup_mass[c(1, 2, 3, 4)]

  names(Iso_Sup)[2:4] <- c("Abundance", "RT", "Tag")

  Iso_Out <- merge(Iso_Out, Dup_rem, by = "Exp_mass", all = TRUE)
  Iso_Out <- Iso_Out[is.na(Iso_Out$Abundance.x), ]
  Iso_Out <- Iso_Out[c(1, 2, 3, 4)]

  # The Iso_mass in this section is the final output isotope list
  Iso_mass <- unique(Iso_Out)

  # This data frame is necessary to remove the flagged isotope peaks from
  # the overall data frame.
  Iso_align <- Iso_mass[c(1, 2, 3)]

  Aligned <- merge(End, Iso_align, by = "Exp_mass", all = TRUE)

  Mono_final <- Aligned[is.na(Aligned$Abundance.y), ]
  Mono_final <- Mono_final[c(1, 2, 3)]

  names(Mono_final) <- c("exp_mass", "abundance", "RT")
  Mono_final <- unique(Mono_final)
  Mono_final <- Mono_final[!is.na(Mono_final$abundance), ]

  Iso_final <- bind_and_filter_data(Iso_mass, Iso_Sup, FALSE)

  names(Iso_final) <- c("exp_mass", "abundance", "RT", "tag")
  Iso_final <- Iso_final[!is.na(Iso_final$abundance), ]

  # making sure all masses are present
  Iso_check <- unique(Iso_final)
  Mono_check <- unique(Mono_final)
  All_check <- dplyr::bind_rows(Mono_check, Iso_check)
  All_check <- All_check[!duplicated(All_check[1:3]), ]
  names(End)[1:2] <- c("exp_mass", "abundance")

  All_check <- merge(End, All_check,
    by = c("exp_mass", "abundance", "RT"), all = TRUE
  )

  Mono_final <- All_check[is.na(All_check$tag), ]
  Mono_final <- Mono_final[c(1:3)]

  Iso_final <- All_check[!is.na(All_check$tag), ]

  if (cols == 2) {
    Mono_final <- Mono_final[c(1, 2)]
    Iso_final <- Iso_final[c(1, 2, 4)]
  }

  return(list(Mono = Mono_final, Iso = Iso_final))
}

#' Extracts mono and iso pairs from a data frame
#'
#' This function takes a data frame of pairs and extracts columns for mono and iso masses.
#'
#' @param pairs data frame:
#' The input data frame
#' @return list:
#' A list with two elements: "mono_pair" containing columns "Mono_mass" and "Iso_mass1,"
#' and "iso_pair" containing columns "Iso_mass1" and "Iso_mass2."
#' @export
extract_mono_and_iso_pairs <- function(pairs) {
  mono_pair <- pairs[c(1, 2)]
  names(mono_pair) <- c("Mono_mass", "Iso_mass1")

  iso_pair <- pairs[c(1, 2)]
  names(iso_pair) <- c("Iso_mass1", "Iso_mass2")

  return(list(mono_pair = mono_pair, iso_pair = iso_pair))
}

#' Merges and filters mono and iso pairs
#'
#' This function merges mono and iso pairs based on the "Iso_mass1" column and filters the result.
#'
#' @param mono_pair data frame:
#' The data frame containing columns "Mono_mass" and "Iso_mass1."
#' @param iso_pair data frame:
#' The data frame containing columns "Iso_mass1" and "Iso_mass2."
#' @return data frame:
#' The merged and filtered data frame containing columns "Mono_mass," "Iso_mass1," and "Iso_mass2."
#' @export
merge_and_filter_pairs <- function(mono_pair, iso_pair) {
  merged_pair <- merge(mono_pair, iso_pair, by = "Iso_mass1", all = TRUE)

  filtered_data <- merged_pair[complete.cases(
    merged_pair$Mono_mass,
    merged_pair$Iso_mass2
  ), ]

  return(filtered_data)
}

#' Merges and filters abundances for mono and iso pairs
#'
#' This function merges mono and iso pairs with abundance data and filters the result.
#'
#' @param final_pair2 data frame:
#' The data frame containing columns "Mono_mass", "Iso_mass1", and "Iso_mass2".
#' @param end_data data frame:
#' The data frame containing columns "Exp_mass", "Abundance", "RT".
#' @return data frame:
#' The merged and filtered data frame containing columns "Mono_mass," "Mono_Abund," "Mono_RT," "Iso_mass1," "Iso_Abund1," "Iso_RT1," "Iso_mass2," "Iso_Abund2," and "Iso_RT2."
#' @export
merge_and_filter_abundances <- function(final_pair2, end_data) {
  mono_abund <- setNames(end_data, c("Mono_mass", "Mono_Abund", "Mono_RT"))
  iso_abund1 <- setNames(end_data, c("Iso_mass1", "Iso_Abund1", "Iso_RT1"))
  iso_abund2 <- setNames(end_data, c("Iso_mass2", "Iso_Abund2", "Iso_RT2"))

  aligned_data <- merge(final_pair2, mono_abund, by = "Mono_mass")
  aligned_data <- merge(aligned_data, iso_abund1, by = "Iso_mass1")
  aligned_data <- merge(aligned_data, iso_abund2, by = "Iso_mass2", all = TRUE)
  aligned_data <- aligned_data[!is.na(aligned_data$Mono_mass), ]

  return(aligned_data)
}

#' Return TRUE if any of the elements in the list has a length of 0, otherwise FALSE
#' @param filtered_data Filtered data with mono-iso links.
#' @return TRUE if any of the links is 0, otherwise FALSE.
#' @export
filtered_data_is_empty<- function(filtered_data) {
  return(any(lapply(filtered_data, length) == 0))
}

#' Process carbon isotoping data
#'
#' This function processes carbon isotoping data using specified parameters and filtering criteria.
#'
#' @param raw_data data frame:
#' Raw data containing information for carbon isotoping.
#' @param carb_error numeric:
#' Error for carbon isotoping.
#' @param carb_ratio numeric:
#' Carbon ratio for filtering.
#' @param end_data data frame:
#' End data containing information for carbon isotoping.
#' @return list:
#' A list containing three data frames: "MonoC" for monoisotopic data, "IsoC1" for first isotope data, and "IsoC2" for second isotope data.
#' @export
process_carbon_isotoping <- function(raw_data, carb_error, carb_ratio, end_data) {
  # The numbers chosen for filtering are based on the results of positive mode ESI for BB burning aerosol
  # Before change on 6/25/19, the second KMDrdiff was -0.4975

  # 0.00149 explanation can be found on page 211 [here](https://core.ac.uk/download/pdf/217038739.pdf)
  # The limits for filtering of isotope pairs for 13C:
  # 0.4975 < KMDrDiff < -0.494501 and 0.501501 < KMDrDiff < 0.5045
  # either 0.498001 is a typo or a value in the thesis (the upper limit)

  pairs <- compute_pairs(
    data = raw_data,
    err = carb_error,
    mass_difference = c13_mass,
    normalization_factor = 21,
    multiplier = 1,
    pair_thresold = c(-1.0015, -1.005),
    filter_criteria = c(0.00149, -0.494501, -0.498001, 0.5045, 0.501501)
  )

  extracted_pair <- extract_mono_and_iso_pairs(pairs)

  filtered_data <- merge_and_filter_pairs(
    extracted_pair$mono_pair,
    extracted_pair$iso_pair
  )

  if(filtered_data_is_empty(filtered_data)) {
    return(list(
      MonoC = raw_data, # treat all data as mono and return it.
      IsoC1 = data.frame(Exp_mass = c(), Abundance = c(), RT = c()),
      IsoC2 = data.frame(Exp_mass = c(), Abundance = c(), RT = c())))
  }

  iso1 <- data.frame(Mono_mass = filtered_data$Iso_mass1, Tag = "Iso1")
  iso1 <- rbind(iso1, data.frame(Mono_mass = -42, Tag = "Iso1"))

  final_set <- merge(extracted_pair$mono_pair, iso1, by = "Mono_mass", all = TRUE)
  final_pair <- final_set[is.na(final_set$Tag), ]

  final_set2 <- merge(final_pair, extracted_pair$iso_pair, by = "Iso_mass1", all = TRUE)
  final_pair2 <- final_set2[!is.na(final_set2$Mono_mass), ]
  final_pair2 <- final_pair2[-3]
  final_pair2 <- final_pair2[c(2, 1, 3)]

  aligned_data <- merge_and_filter_abundances(final_pair2, end_data)

  aligned_data <- aligned_data[aligned_data$Iso_Abund1 <
    (carb_ratio / 100) * aligned_data$Mono_Abund, ]

  iso_pair_data <- aligned_data[is.na(aligned_data$Iso_mass2), ]
  iso_tri_data <- aligned_data[!is.na(aligned_data$Iso_mass2) &
    (aligned_data$Iso_Abund2 < (carb_ratio / 100) * aligned_data$Iso_Abund1), ]

  final_aligned_data <- rbind(iso_pair_data, iso_tri_data)
  mono_columns <- final_aligned_data[c(3, 4, 5)]
  names(mono_columns) <- c("Exp_mass", "Abundance", "RT")

  iso_columns1 <- final_aligned_data[c(2, 6, 7)]
  names(iso_columns1) <- c("Exp_mass", "Abundance", "RT")

  iso_columns2 <- final_aligned_data[c(1, 8, 9)]
  names(iso_columns2) <- c("Exp_mass", "Abundance", "RT")

  return(list(MonoC = mono_columns, IsoC1 = iso_columns1, IsoC2 = iso_columns2))
}

#' Process sulfur data
#'
#' This function processes sulfur data using specified parameters and filtering criteria.
#'
#' @param raw_data data frame:
#' Raw data containing information for sulfur.
#' @param sulferr numeric:
#' Error for sulfur.
#' @param data_end data frame:
#' End data containing information for sulfur.
#' @param sulfrat numeric:
#' Sulfur ratio for filtering.
#' @return list:
#' A list containing two data frames: "IsooutS" for isotope data and "MonooutS" for monoisotopic data.
#' @export
process_sulfur <- function(raw_data, sulferr, data_end, sulfrat) {
  pairs <- compute_pairs(
    data = raw_data,
    err = sulferr,
    mass_difference = s34_mass,
    normalization_factor = 12,
    multiplier = 2,
    pair_thresold = c(-1.990, -2),
    filter_criteria = c(0.00249, -0.29051, -0.29349, 0.70949, 0.7075)
  )

  data_end2 <- data_end

  names(data_end2) <- c("Exp_mass1", "Iso_Abund", "Iso_RT")

  merged_pairs <- merge(pairs, data_end, by = "Exp_mass")
  merged_pairs <- unique(merged_pairs)
  merged_pairs <- merge(merged_pairs, data_end2, by = "Exp_mass1")

  MonoS <- merged_pairs[c(2, 14, 15)]
  IsoS <- merged_pairs[c(1, 16, 17)]
  names(IsoS)[1] <- "Iso_mass"

  abund <- cbind(MonoS, IsoS)
  abund_dummy <- data.frame(
    Exp_mass = -42,
    Abundance = -1,
    RT = -42,
    Iso_mass = -42,
    Iso_Abund = -1,
    Iso_RT = -42
  )
  abund <- rbind(abund, abund_dummy)
  abund$ratio <- abund$Iso_Abund / abund$Abundance * 100

  filtered_abund <- abund[abund$ratio <= sulfrat, ]
  filtered_abund <- unique(filtered_abund)

  MonooutS <- filtered_abund[c(1, 2, 3)]
  IsooutS <- filtered_abund[c(4, 5, 6)]
  names(IsooutS) <- c("Exp_mass", "Abundance", "RT")

  return(list(IsooutS = IsooutS, MonooutS = MonooutS))
}

#' Create pairs based on mass difference criteria
#'
#' This function generates pairs of mass values based on specified criteria.
#'
#' @param exp_mass numeric vector:
#' Vector of experimental mass values.
#' @param less_than numeric:
#' pair thresold.
#' @param greater_than numeric:
#' pair thresold.
#' @return data frame:
#' A data frame containing pairs of experimental mass values that satisfy the specified mass difference criteria.
#' @export
create_pairs <- function(exp_mass, less_than, greater_than) {
  pairs_data <- expand.grid(Exp_mass = exp_mass, Exp_mass1 = exp_mass)
  pairs_data$mdiff <- pairs_data$Exp_mass - pairs_data$Exp_mass1
  pairs_data <- pairs_data[pairs_data$mdiff < less_than & pairs_data$mdiff > greater_than, ]
  return(pairs_data)
}

#' Filter pairs based on mass difference error
#'
#' This function filters pairs based on mass difference error.
#'
#' @param pairs_data data frame:
#' Data frame containing pairs of mass values.
#' @param mass_difference numeric:
#' Mass difference value.
#' @param threshold numeric:
#' Error threshold for filtering pairs.
#' @return data frame:
#' A data frame containing pairs of mass values that satisfy the specified error threshold.
#' @export
filter_pairs <- function(pairs_data, mass_difference, threshold) {
  error_col <- abs(((pairs_data$Exp_mass + mass_difference) -
    pairs_data$Exp_mass1) / pairs_data$Exp_mass1 * 10^6)
  pairs_data <- pairs_data[error_col <= threshold, ]
  pairs_data <- pairs_data[c(1:3)]
  return(pairs_data)
}

#' Calculate variables based on mass values
#'
#' This function calculates variables based on specified mass values and parameters.
#'
#' @param exp_mass_col numeric:
#' experimental mass values.
#' @param mass_difference numeric:
#' Mass difference value.
#' @param normalization_factor numeric:
#' Normalization factor.
#' @param multiplier numeric:
#' Multiplier value.
#' @return matrix:
#' calculated variables: km, kmd, kmr, and kmdr.
#' @export
calculate_variables <- function(
    exp_mass_col,
    mass_difference,
    normalization_factor,
    multiplier) {
  km <- exp_mass_col * (multiplier / mass_difference)
  kmd <- round((round(exp_mass_col) - km), 3)
  kmr <- exp_mass_col * ((round((ch2_mass / normalization_factor)) /
    ((ch2_mass / normalization_factor))))
  kmdr <- round((round(kmr) - kmr), 3)

  return(cbind(km, kmd, kmr, kmdr))
}

#' Compute pairs based on specified criteria and filtering
#'
#' This function computes pairs of mass values based on specified criteria and filtering parameters.
#'
#' @param data data frame:
#' Data frame containing experimental mass values.
#' @param err numeric:
#' Error value for filtering.
#' @param mass_difference numeric:
#' Mass difference value.
#' @param normalization_factor numeric:
#' Normalization factor.
#' @param multiplier numeric:
#' Multiplier value.
#' @param pair_thresold numeric vector:
#' Vector containing pair mass difference criteria.
#' @param filter_criteria numeric vector:
#' Vector containing filtering criteria.
#' @return data frame:
#'   A data frame containing computed pairs of mass values that satisfy the specified criteria and filtering parameters.
#' @export
compute_pairs <- function(
    data,
    err,
    mass_difference,
    normalization_factor,
    multiplier,
    pair_thresold,
    filter_criteria) {
  exp_mass <- unlist(data[1])

  pairs_data <- create_pairs(exp_mass, pair_thresold[1], pair_thresold[2])
  pairs_data <- filter_pairs(pairs_data, mass_difference, err)

  pairs_data[, c("KM", "KMD", "KMr", "KMDr")] <- calculate_variables(
    pairs_data$Exp_mass,
    mass_difference,
    normalization_factor,
    multiplier
  )
  pairs_data[, c("KM1", "KMD1", "KMr1", "KMDr1")] <- calculate_variables(
    pairs_data$Exp_mass1,
    mass_difference,
    normalization_factor,
    multiplier
  )

  pairs_data$KMDrdiff <- pairs_data$KMDr - pairs_data$KMDr1
  pairs_data$KMDdiff <- pairs_data$KMD - pairs_data$KMD1

  filtering_criteria <- abs(pairs_data$KMDdiff) <= filter_criteria[1] & (
    (pairs_data$KMDrdiff < filter_criteria[2] & pairs_data$KMDrdiff > filter_criteria[3]) |
      (pairs_data$KMDrdiff < filter_criteria[4] & pairs_data$KMDrdiff > filter_criteria[5])
  )

  pairs <- pairs_data[filtering_criteria, ]

  return(pairs)
}

#' Compare and filter isotopic and bad data
#'
#' This function compares isotopic and bad data frames based on mass, abundance, and retention time.
#' It returns a list containing filtered isotopic and bad data frames.
#'
#' @param iso_df data frame:
#' Data frame containing isotopic data.
#' @param bad_df data frame:
#' Data frame containing bad data for comparison.
#' @return list:
#' A list containing filtered isotopic and bad data frames.
#' @export
compare_bad_and_filter <- function(iso_df, bad_df) {
  merged_df <- merge(iso_df, bad_df,
    by = c("Exp_mass", "Abundance", "RT"),
    all = TRUE
  )
  iso_df <- merged_df[(!is.na(merged_df$Tag) & is.na(merged_df$Tag.y)) |
    (!is.na(merged_df$Tag) & !is.na(merged_df$Tag.y)), ]
  iso_df <- iso_df[c(1:4)]
  bad_df <- merged_df[(is.na(merged_df$Tag) & !is.na(merged_df$Tag.y)), ]
  return(list(iso_df, bad_df))
}

#' Check sulfur data
#'
#' This function performs sulfur-related checks and processes the data accordingly.
#'
#' @param MonoOutC_final data frame:
#' Data frame containing final monoisotopic carbon data.
#' @param MonoOutS_final data frame:
#' Data frame containing final monoisotopic sulfur data.
#' @param IsoOutC1_final data frame:
#' Data frame containing final isotope 1 carbon data.
#' @param IsoOutC2_final data frame:
#' Data frame containing final isotope 2 carbon data.
#' @return list:
#' A list containing three processed data frames: monoisotopic sulfur, isotope 1 carbon, and isotope 2 carbon.
#' @export
sulfur_check <- function(MonoOutC_final, MonoOutS_final, IsoOutC1_final, IsoOutC2_final) {
  # Filtering the 34S masses that do not have a corresponding 13C
  Sulfs <- merge(
    MonoOutC_final,
    MonoOutS_final,
    by = c("Exp_mass", "Abundance", "RT"),
    all = TRUE
  )
  Sulfgood <- Sulfs[complete.cases(Sulfs), ]
  Sulfbad <- Sulfs[!complete.cases(Sulfs), ]

  # Comparison of "bad" 34S to C13
  iso_out_c1_list <- compare_bad_and_filter(IsoOutC1_final, Sulfbad)
  IsoOutC1_final <- iso_out_c1_list[[1]]
  Sulfbad <- iso_out_c1_list[[2]]

  # Comparison of "bad" 34S to 2C13
  Sulfbad <- Sulfbad[-4]
  iso_out_c2_list <- compare_bad_and_filter(IsoOutC2_final, Sulfbad)
  IsoOutC2_final <- iso_out_c2_list[[1]]
  Sulfbad <- iso_out_c2_list[[2]]

  # Process Sulfbad and Sulfgood
  Sulfbad <- Sulfbad[c(1, 2, 3, 6)]
  names(Sulfbad)[4] <- "Tag"
  Sulfgood <- Sulfgood[c(1, 2, 3, 5)]
  names(Sulfgood)[4] <- "Tag"

  MonoOutS_final <- rbind(Sulfgood, Sulfbad)

  return(list(
    MonoOutS_final = MonoOutS_final,
    IsoOutC1_final = IsoOutC1_final,
    IsoOutC2_final = IsoOutC2_final
  ))
}
