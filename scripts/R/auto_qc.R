#!/usr/bin/env Rscript
options(scipen=999);

# suppress warning on R build version #
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))

# Parameters -------------------------------------------------------------------------------------------------------------

# limits for biological plausibility of a solution
qclim_purity_low_wGII_cn_pass = 0.95
qclim_wGII_high_purity_cn_pass = 0.01
qclim_strict_loh_cn_pass = 0.01
qclim_biggest_homdel_Mb = 50
qclim_frac_gn_homdel = 0.025
qclim_LOH_gII_ratio = 0.01
qclim_frac_genome_altered_lim = 0.15
qclim_purity = 0.1
qclim_clonal_cn_change_size = 20000000

# Flags to consider reviewing a chosen solution
flaglim_allele_fit = 50
flaglim_n_segments = 350

# Functions--------------------------------------------------------------------------------------------------------------------

calc.ploidy <- function(segments, input = c('ascat', 'refphase', 'auto_qc')) {

  # refphase/ascat/auto_qc objects have slight different col names - adjust to input
  if( input == 'ascat'){
    start_col <- 'startpos' ; end_col <- 'endpos' ; Acn_raw_col <- 'nAraw'
    Bcn_raw_col <- 'nBraw' ; sample_col <- 'sample' ; totcn_col <- 'cnTotal'
    chr_col <- 'chr'
  }

  if( input == 'refphase'){
    start_col <- 'start' ; end_col <- 'end' ; Acn_raw_col <- 'cn_major_raw'
    Bcn_raw_col <- 'cn_minor_raw' ; sample_col <- 'sample' ; totcn_col <- 'total_cn'
    chr_col <- 'chrom'
  }

  if( input == 'auto_qc'){
    start_col <- 'start' ; end_col <- 'end' ; Acn_raw_col <- 'nAraw'
    Bcn_raw_col <- 'nBraw' ; sample_col <- 'sample' ; totcn_col <- 'cnTotal'
    chr_col <- 'chromosome'
  }

  # Calculates various ploidy metrics for a given sample
  # Input:
  #   ASCAT segments array (ascat.output$segments_raw)
  # Output:
  #   A datatable with the following columns:
  #       sample - sample
  #       g_ploidy  - The most frequent ploidy across the genome [Previously 'Major']
  #       wchr_ploidy - The most predominant ploidy weighted by chromosome
  #                 (i.e. most frequent ploidy with chromosomes weighted equally) [Perviously wMajor]
  #       mean_raw - the mean Araw and Braw scores across the genome.
  #       0..10 - The frequency of ploidy from 0 to 10. Note that ploidy higher than 10
  #               is filtered out.
  segments <- as.data.table(segments)
  segments[, segment_len := as.numeric(get(end_col) - get(start_col)) ]
  segments[, weighted_raw := (get(totcn_col)) * segment_len, by = get(sample_col)]

  # Calculate the proportion of the genome classified with a given ploidy
  segments_prop <- segments[, 
                            .(ploidy_len = sum(segment_len)), by = .(get(totcn_col), get(sample_col))][,
                                                                                           .(get, ploidy = ploidy_len / sum(ploidy_len)), by = get.1][
                                                                                             order(-ploidy),
                                                                                             ]
  setnames(segments_prop, c('get', 'get.1'), c(totcn_col, sample_col))

  # Calculate predominate ploidy by sample across entire genome
  g_ploidy <- segments_prop[, .(g_ploidy = get(totcn_col)[which.max(ploidy)]), by = eval(sample_col) ]

  # Calculate predominate ploidy by asking what the most common ploidy is across
  # all chromosomes.
  wchr_ploidy <- segments[, .(ploidy_len = sum(segment_len)), by = .(get(totcn_col), get(sample_col), get(chr_col))][,
                                         .(get, ploidy = ploidy_len / sum(ploidy_len)), by = .(get.1, get.2)][,
                                                                                                  .SD[which.max(ploidy)], by = .(get.1, get.2)][,
                                                                                                                          .(n = .N), by = .(get.1, get)][,
                                                                                                                                               .(wchr_ploidy = get[which.max(n)]), by = get.1]
  setnames(wchr_ploidy, 'get.1', sample_col)

  # Calculate the mean raw intensity score is
  mean_raw = segments[, .(mean_raw = sum(weighted_raw / sum(segment_len))), by = get(sample_col)]
  setnames(mean_raw, 'get', sample_col)

  out <- dcast(segments_prop, get(sample_col) ~ get(totcn_col), value.var = "ploidy")
  setnames(out, 'sample_col', sample_col)
  out <- g_ploidy[out, on = sample_col]
  out <- wchr_ploidy[out, on = sample_col]
  out <- mean_raw[out, on = sample_col]
  # Fill missing cols with NA
  for(i in 0:10) {
    if (!(i %in% names(out))) {
      out[[as.symbol(i)]] = NA
    }
  }

  # Only consider ploidy 1-10
  max_cn <- suppressWarnings(max(as.numeric(names(out)),na.rm=T))
  if( max_cn > 10){
     remove_cols <- 11:max_cn
     remove_cols <- as.character( remove_cols[ remove_cols %in% names(out) ] )
     out[, (remove_cols) := NULL ]
  }
  setcolorder(out, c(sample_col, 'mean_raw', 'wchr_ploidy', 'g_ploidy', as.character(0:10) ) )
  return(out)
}

calc.wgii <- function(segments, input = c('ascat', 'refphase', 'auto_qc'),
                      threshold = 0.6, check.names = FALSE, include.sexchrom = FALSE) {

  # refphase/ascat have slight different col names - adjust to input
  if( input == 'ascat'){
    start_col <- 'startpos' ; end_col <- 'endpos' ; Maj_cn_col <- 'nMajor'
    Min_cn_col <- 'nMinor' ; sample_col <- 'sample' ; totcn_col <- 'cnTotal'
    chr_col <- 'chr'
  }
  
  if( input == 'refphase'){
    start_col <- 'start' ; end_col <- 'end' ; Maj_cn_col <- 'cn_major_raw'
    Min_cn_col <- 'cn_minor_raw' ; sample_col <- 'sample' ; totcn_col <- 'total_cn'
    chr_col <- 'chrom'
  }

  if( input == 'auto_qc'){
    start_col <- 'start' ; end_col <- 'end' ; Maj_cn_col <- 'nMajor'
    Min_cn_col <- 'nMinor' ; sample_col <- 'sample' ; totcn_col <- 'cnTotal'
    chr_col <- 'chromosome'
  }

  # Edit 20140826: gii & wgii: exclude sex chromosomes
  # Edit 20150316: Use wMajor as ploidy
  # Edit 20200407: Rewritten for clarity & documentation
  # Input
  # ascat.output$segments_raw
  # Output:
  # A data.table with the following columns
  #    gii = "Genome integrity index"; Measures the proportion of the genome
  #          with an aberrant ploidy.
  #    wgii = "Weighted genome integrity index"; Measures the proportion of the genome
  #          with an aberrant ploidy, weighted by chromosome.
  #    floh = "Fraction Loss of Heterozygosity; The proportion of the genome with loss-of-heterozygosity.
  #    wfloh = "Weighted Fraction Loss of Heterozygosity "; The proportion of the genome, weighted by chromosome, with loss-of-heterozygosity.
  #            Does not include sex chromosomes.

  if(!include.sexchrom) segments <- segments[ !get(chr_col) %in% c(23, 24, 'chrX', 'chrY', 'X', 'Y') ]

  segments <- as.data.table(segments)

  # Get ploidy information
  ploidy <- calc.ploidy(segments, input)

  # Use the weighted genome ploidy (wmajor)
  # If ploidy is less than 2, set it to 2.
  ploidy[ wchr_ploidy < 2, wchr_ploidy := 2 ]
  
  segments <- ploidy[segments, on = sample_col]
  
  # Calc. segment length
  segments[, segment_len := as.numeric(get(end_col) - get(start_col)) ]
  
  # Label abberrant segments; wchr_ploidy = predominant ploidy genome wide, weighted by chromosome
  segments[, aberrant := (get(totcn_col) < (wchr_ploidy - threshold)) | (get(totcn_col) > (wchr_ploidy + threshold))]

  # Fold in "measured chromosome length"; Calculate measured genome length by sample
  chrom_len <- segments[,.(chrom_len = max(get(end_col)) - min(get(start_col))), by = .(get(chr_col), get(sample_col))]
  setnames(chrom_len, c('get', 'get.1'), c(chr_col, sample_col))

  # Join segment and chrom lengths
  segments <- chrom_len[segments, on = c(sample_col, chr_col)]

  # Summarize aberrant segments
  segment_len <- segments[, .(aberrant_len = sum(ifelse(aberrant, segment_len, 0)),
                              chrom_len = data.table::first(chrom_len)), by = .(get(sample_col), get(chr_col))]
  setnames(segment_len, c('get', 'get.1'), c(sample_col, chr_col))
  
  sex_chrom <- c("X", "Y")
  # Calculate segment metrics by whole genome
  calc_wg <- segment_len[, .(gii = sum(aberrant_len) / sum(chrom_len)), by = get(sample_col)]
  setnames(calc_wg, c('get'), c(sample_col))

  # Calculate segment metrics weighted by chromsoome, and without sex chromosomes.
  calc_chrom <- segment_len[!(get(chr_col) %in% sex_chrom), .(wgii = mean(aberrant_len / chrom_len)),
                             by = get(sample_col)]
  setnames(calc_chrom, c('get'), c(sample_col))

  out <- calc_chrom[ calc_wg, on = sample_col ]
  setcolorder(out, c(sample_col, "gii", "wgii"))
  return(out)
}


assess_solution <- function( segments, homdel_lim = 0.8, subcl_lim = 0.25, tight50CCF_lim = 0.4, tight_clonal_lim = 0.15, bin_size = 10000, 
                             loose_clonal_lim = 1.2, abs_tight_cn_lim = 3, strict_loh_lim = 1.2, imbalence_lim = 0.1){


  start_col <- 'start' ; end_col <- 'end' ; total_cn_col <- 'cnTotal'
  sample_col <- 'sample' ; chr_col <- 'chromosome'
  
  #remove sex chrs
  segments <- segments[ !get(chr_col) %in% c(23, 24, 'chrX', 'chrY', 'X', 'Y') ]

  # remove missing data
  segments <- segments[ !is.na(get(total_cn_col)) ]

  # Calc. segment length
  segments[, segment_len := as.numeric(get(end_col) - get(start_col)) ]
  segments[, segment_gn_frac := segment_len / sum(segment_len, na.rm = TRUE), by = get(sample_col) ]
  setnames( segments, total_cn_col, 'total_raw')
  segments[, is_subclonal := abs(total_raw) > subcl_lim  ]
  
  # want to get a sense for how nicely areas of imbalence fit onto integers (like goodness of fit)
  # Use this to calculate CCFs (even for solutions where the CN isn't sitting on integers)
  # Use these CCFs to find tight 50% CCF that may indicate solution needs doubling
  allelle_cns <- melt(segments, id.vars = c('segment_len', sample_col), measure.vars = 'total_raw',
                      variable.name= 'allele', value.name = 'raw_cn' )
  allelle_cns[ raw_cn > 0.5, int_diff := raw_cn / round(raw_cn) ]
  allelle_cns[, prob_clonal := abs(int_diff) < loose_clonal_lim & raw_cn < 6 ]

  # Devide into bins and take a median so we're more resistant to outliers (additional segments that switch from 'probably clonal' in solutions with slightly shifted pliody)
  allelle_cns[, n_bins := round(segment_len / bin_size) ]
  int_offsets <- allelle_cns[ (prob_clonal),
                              .(integer_offset = median( rep(int_diff, n_bins), na.rm=T ) ),
                              by = get(sample_col) ]
  setnames(int_offsets, c('get'), c(sample_col))
  segments <- merge( segments, int_offsets, all.x = TRUE )

  segments[, total_cor := total_raw / integer_offset ]
  segments[, is_tight_50CCF := (abs(total_cor - round(total_cor)) > tight50CCF_lim & total_cor < 5) ] 
  segments[, is_tight_clonal := (abs(total_cor - round(total_cor)) < tight_clonal_lim & total_cor < 4) ] 

  # Calculate max size of the biggest clonal CN change
  segments[, is_pliody := round(total_raw, 0) == round( mean(total_raw[is_tight_clonal]), 0 ), by = get(sample_col)] #?
  max_clonal_seg <- segments[ (is_tight_clonal & !is_pliody), .(max_abberant_clonal_segment = max(segment_len)), by = eval(sample_col)]

  # Look at whether odd CNs are presented excluding >=5 as high CNs are noiser and lower should always be present
  segments[, total_diff_from_odd := min(abs(total_raw - c(1, 3))), 1:nrow(segments) ]
  
  # How much Allelic imabalence is there? Should go up with GDs but also decreases with purity
  # Homozygous deletions (even subclonal) should be infrequent and small events when they occur
  # Generally if we can avoid solutions with lots of subclonal CN we do as this is more parsimonious 
  # Very high number of segemnts may indicate oversegmentation 
  out <- segments[, .(frac_at_odd = sum(segment_gn_frac[ total_diff_from_odd < 0.25 ]),
                      frac_at_2 = sum(segment_gn_frac[ abs(total_raw - 2) < 0.25 ]),
                      frac_homdel = sum( segment_gn_frac[total_raw < homdel_lim] ),
                      floh_strict = sum( segment_gn_frac[ total_raw < strict_loh_lim ] ),
                      frac_subclonal = sum(segment_gn_frac[ is_subclonal ]),
                      frac_tight_50_CCF_subclonal = sum(segment_gn_frac[ is_tight_50CCF ]),
                      frac_tight_clonal = sum(segment_gn_frac[ is_tight_clonal ]),
                      N_segments = .N,
                      integer_offset = unique(integer_offset), 
                      biggest_homdel_Mb = ifelse( any(total_raw < homdel_lim), max(segment_len[total_raw < homdel_lim]) / 1000000, 0)), 
                   by = get(sample_col) ]
  setnames(out, c('get'), c(sample_col))

  out <- merge( out, max_clonal_seg, on = sample_col)
  
  return( out )
}

# Function to select the best solution
select_best_solution <- function(seg_files){
  segs <- fread( seg_files )
  # match columns to expected by code
  setnames(segs, 'absolute_copy_number', 'cnTotal')
  
  # Calculate metrics per region-------------------------------------------------------------------------------------------------------------
  
  ## Calculate metrics based on CN data alone ##
  
  # To run per solution with the same functions temporarily modify the sample column
  segs[, sample_sol := paste(sample, ploidy, cellularity, sep = ':')]
  segs[, sample_id := sample ]
  segs[, sample := sample_sol ] # for code expected inputs (should really change the function)
  segs[, segment_len := end - start + 1 ]
  segs[, num_gds := ifelse(ploidy < 2.7, 0, ifelse(ploidy < 3.5, 1, 2))]
  
  # Caluclate metrics
  sample_qc_metrics <- assess_solution( segs )
  setnames(sample_qc_metrics, 'sample', 'sample_sol')
  cin_metrics <- calc.wgii( segs, input = 'auto_qc' )
  setnames(cin_metrics, 'sample', 'sample_sol')

  # make a goodness of fit score
  # average deviation from interger (limited to 0.2 max) - as per Vanloo et al ASCAT paper
  segs[, int_deviation := abs(cnTotal - round(cnTotal)) ]
  segs[, int_deviation_scaled := (1 - (ifelse(int_deviation > 0.2, 0.2, int_deviation)/0.2))*100 ]
  gof <- segs[, .(goodness_of_fit = sum(int_deviation_scaled * segment_len) / sum(segment_len)), by = sample_sol ]
  
  print(nrow(sample_qc_metrics))
  print(nrow(cin_metrics))
  print(nrow(gof))
  # Bind together
  cn_qc <- merge( sample_qc_metrics,  
                  cin_metrics, by = 'sample_sol' )
  cn_qc <- merge( cn_qc,  
                  gof, by = 'sample_sol' )
  cn_qc <- merge( cn_qc,  
                  unique(segs[,.(sample_sol, sample = sample_id, ploidy, purity = cellularity, num_gds)]), by = 'sample_sol' ) 
   
  # make a logical order
  cn_qc <- cn_qc[, c(16, 1:15, 17:19)] 
  
  # Choose a solution ufor each region -------------------------------------------------------------------------------------
  
  # Rule out solutions that are biologically implausible #
  cn_qc[, `:=`( homdel_size_fail = biggest_homdel_Mb > qclim_biggest_homdel_Mb,
                homdel_frac_fail = frac_homdel > qclim_frac_gn_homdel,
                loh_fail = floh_strict / gii < qclim_LOH_gII_ratio,
                purity_fail = purity < qclim_purity | purity > qclim_purity_low_wGII_cn_pass & wgii < qclim_wGII_high_purity_cn_pass,
                size_clonal_cn_change_fail = max_abberant_clonal_segment < qclim_clonal_cn_change_size ) ]
  
  cn_qc[, bio_plausible_sol := !( homdel_size_fail | homdel_frac_fail | loh_fail | purity_fail | size_clonal_cn_change_fail ) ]
  
  # Now lets make a solution score for the solutions
  
  # Lets adjust the fraction of gn with tight 50% CCF (indicating solution needs doubling) by how much noise there is (non tight subclonal CN) 
  # scale this score by number of GDs. If 0 GD called then this score shouldn't be used as its a check for GD samples. Also more concnered about this
  # score if we've got 2 rather than 1 GDs as should be conservative on number of GDs where possible (hence (1 / (num_gds + 1)))
  # squaring and X50 is abirtary to make the change in score of a reaonable scale compared to the other factors going into the score
  cn_qc[, frac_subclonal_tight50CCF := frac_tight_50_CCF_subclonal / (1-frac_tight_clonal) ]
  cn_qc[, frac_subclonal_tight50CCF := ifelse( frac_subclonal_tight50CCF > 1, 1, frac_subclonal_tight50CCF )]
  cn_qc[, extra_int_score := frac_subclonal_tight50CCF ^ 2 * 50 * (1 / (num_gds + 1)) ]
  
  # Prevents NA (can happen with very poor solutions but everything needs a score)
  cn_qc[ is.na(extra_int_score), extra_int_score := 1 ]
  
  # Add punishment for missing CN states
  # Do this on a log scale so that 0.001 is punished much more than 0.01
  # But don't allow it punish more than 0.6 (equivilent to 0.001 frac odd)
  cn_qc[, frac_at_odd_lim := ifelse(frac_at_odd < 0.001, 0.001, frac_at_odd) ]
  cn_qc[, missing_cn_score := -(log10(frac_at_odd_lim) - log10(0.1 * num_gds)) * 10  ]
  cn_qc[ missing_cn_score < 0 , missing_cn_score := 0 ]
  
  # Round missing/extra integer scores to nearest 5% so that choices between reasonable solutions will be mostly made on goodness of fit 
  cn_qc[, solution_fit_score := goodness_of_fit -
          round(extra_int_score * 5) / 5 -
          round(missing_cn_score * 5) / 5  ]
  
  # Then select solution based on score after ruling out all bad solutions
  cn_qc[ (bio_plausible_sol),
         selected_solution_auto := sample_sol == sample_sol[ solution_fit_score == max(solution_fit_score[ bio_plausible_sol ], na.rm=T) ][1],
         by = sample ]
  
  # Lets flag chosen solutions that look a bit dodgy for manual review -----------------------------------------------------------------------
  
  cn_qc[ (selected_solution_auto), flag_goodness_of_fit := goodness_of_fit < flaglim_allele_fit ]
  cn_qc[ (selected_solution_auto), flag_segmentation := N_segments > flaglim_n_segments ]
  
  cn_qc[ (selected_solution_auto), flag_poor_solution :=  flag_segmentation | flag_goodness_of_fit ]

  return(cn_qc)
}
