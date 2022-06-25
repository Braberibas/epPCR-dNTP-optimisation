sim_epPCR = function(nt, AAs, S_tot, N_cycles, seq, seed, base_probs, base_transition_freqs, mark_nt_freqs, mark_AA_freqs, dir_out, verbose = TRUE, plot_ind = FALSE, n.cores = "auto") {
  
  # Perform simulations of epPCR reactions.
  
  # Arguments:
  #   `nt`: A vector of character strings containing exactly the characters "A", "C", "G" and "T" in any order.
  #   `AAs`: A vector of character strings containing exactly the characters "K", "N", "T", "R", "S", "I", "M", "Q", "H", "P", "L", "E", "D", "A", "G", "V", "*", "Y", "C", "W", "F" in any order.
  #   `S_tot`: The number of synthetic sequences to generate and evaluate per cycle.
  #   `N_cycles`: The number of cycles to perform.
  #   `seq`: Character string with the template DNA sequence. Must be comprised solely of the characters in `nt`.
  #   `seed`: Random number seed.
  #   `base_probs`: Named numeric vector of length four, containing the probabilities that the nts in `nt` mutate. Names must exactly match `nt`.
  #   `base_transition_freqs`: Four-by-four matrix containing the probabilities that the base in the row changes to the base in the column, given that the base mutates. Rownames and column names must exactly match `nt`.
  #   `mark_nt_freqs`: Four-by-four matrix analogous to `base_transition_freqs`, but containing frequencies which will be marked in nt transition distributions. Specify `NULL` if no marking is required.
  #   `mark_AA_freqs`: Analogous to `mark_nt_freqs`, but for AAs. Row and column names must exactly match `AAs`.
  #   `dir_out`: Output directory, in which all plots and the computed transition frequencies and counts are stored.
  #   `verbose`: Should user information be printed? (Default: TRUE.)
  #   `plot_ind`: Should summary plots of each cycle be generated? This is the main time consumer for a moderate number of sequences per cycle (Default: FALSE.)
  #   `n.cores`: Integer or "auto". Specifies the number of cores to use for some data structure conversions (the major time consumer for a number of cycles larger than 500). "auto" uses the number of cores of the system minus one.
  
  # Returns: Nothing.
  
  ### User info
  f = list.files(dir_out, full.names = TRUE)
  if (length(f) != 0) {
    print("WARNING: The specified output directory is not empty. All files in the directory will be deleted if execution is not terminated. Execution will re-commence in ten seconds.")
    Sys.sleep(10)
  }
  unlink(f, recursive = TRUE)
  print("Initialising simulation.")
  
  ### Create directory, if necessary
  if (!dir.exists(dir_out)) {
    dir.create(dir_out, recursive = TRUE)
  }
  
  ### Useful operator
  `%notin%` = Negate(`%in%`)
  
  ### Check input
  ## `nt`
  nt_good = TRUE
  for (i in nt) {
    if (i %notin% c("A", "C", "G", "T")) {
      nt_good = FALSE
    }
  }
  if (length(nt) != 4) {
    nt_good = FALSE
  }
  if (!nt_good) {
    stop("`nt` must comprise exactly the symbols 'A', 'C', 'G' and 'T', in any order.")
  }
  
  ## `AAs`
  AA_good = TRUE
  for (i in AAs) {
    if (i %notin% c("K", "N", "T", "R", "S", "I", "M", "Q", "H", "P", "L", "E", "D", "A", "G", "V", "*", "Y", "C", "W", "F")) {
      AA_good = FALSE
    }
  }
  if (length(AAs) != 21) {
    AA_good = FALSE
  }
  if (!AA_good) {
    stop("`AAs` must comprise exactly the symbols 'K', 'N', 'T', 'R', 'S', 'I', 'M', 'Q', 'H', 'P', 'L', 'E', 'D', 'A', 'G', 'V', 'Y', 'C', 'W', 'F' and '*', in any order.")
  }
  
  ## `seq`
  for (i in 1:nchar(seq)) {
    symb = substr(seq, i, i)
    if (symb %notin% nt) {
      stop("`seq` contains forbidden symbol '", symb, "' at position ", i, ". The only allowed symbols are 'A', 'C', 'G' and 'T'.")
    }
  }
  
  ## `seed`
  if (!is.numeric(seed)) {
    stop("`seed` must be numeric.")
  }
  
  ## `base_probs`
  for (prob in base_probs) {
    if (prob < 0 | prob > 1) {
      stop("All base mutation probabilities in `base_probs` must lie between 0 and 1.")
    }
  }
  if (!identical(names(base_probs), nt)) {
    stop("`base_probs` must be a named vector with names exactly identical to `nt`.")
  }
  
  ## `base_transition_freqs`
  if (!is.matrix(base_transition_freqs)) {
    stop("`base_transition_freqs` must be a matrix.")
  }
  if ((nrow(base_transition_freqs) != length(nt)) | (ncol(base_transition_freqs) != length(nt))) {
    stop("`base_transition_freqs` must have the same number of rows and columns as there are entries in `nt`.")
  }
  if (!identical(rownames(base_transition_freqs), nt) | !identical(colnames(base_transition_freqs), nt)) {
    stop("Row and column names of `base_transition_freqs` must be exactly identical to `nt`.")
  }
  
  ## `mark_nt_freqs`
  if (!is.null(mark_nt_freqs)){
    if (!is.matrix(mark_nt_freqs)) {
      stop("`mark_nt_freqs` must be a matrix, or `NULL`.")
    }
    if (nrow(mark_nt_freqs) != length(nt) | (ncol(mark_nt_freqs) != length(nt))) {
      stop("`mark_nt_freqs` must have the same number of rows and columns as there are entries in `nt`.")
    }
    if (!identical(rownames(mark_nt_freqs), nt) | !identical(colnames(mark_nt_freqs), nt)) {
      stop("Row and column names of `mark_nt_freqs` must be exactly identical to `nt`.")
    }
  }
  
  ## `mark_AA_freqs`
  if (!is.null(mark_AA_freqs)) {
    if (!is.matrix(mark_AA_freqs)) {
      stop("`mark_nt_freqs` must be a matrix or `NULL`.")
    }
    if (nrow(mark_AA_freqs) != length(AAs) | (ncol(mark_AA_freqs) != length(AAs))) {
      stop("`mark_AA_freqs` must have the same number of rows and columns as there are entries in `AAs`.")
    }
    if (!identical(rownames(mark_AA_freqs), AAs) | !identical(colnames(mark_AA_freqs), AAs)) {
      stop("Row and column names of `mark_AA_freqs` must be exactly identical to `AAs`.")
    }
  }
  
  ## `n.cores`
  require(doParallel)
  if (n.cores != "auto" | n.cores < 0) {
    stop(paste0("Invalid `n.cores` argument. Must be either 'auto' or a positive integer specifying the number of workers to use for parallel processing."))
  }
  if (n.cores != "auto" & n.cores > detectCores()) {
    stop("`n.cores` must be smaller or equal to the number of logical processors of your system. We suggest leaving at least one core free.")
  }
  
  ### Load packages
  require(stringr)
  require(ggplot2)
  require(cowplot)
  require(Biostrings)
  require(Rsamtools)
  require(viridis)
  require(tidyverse)
  require(vcd)
  require(foreach)
  
  ### Plotting defaults
  theme_set(theme_classic(base_size = 12, base_line_size = 0.7))
  A4_full <- c(height = 11.7, width = 8.3)
  
  ### Set random number seed
  set.seed(seed)
  
  ### Define entitities to store all transition matrices
  global_nt_trans_count_mat = vector(mode = "list", length = N_cycles)
  global_nt_trans_freq_mat = vector(mode = "list", length = N_cycles)
  global_AA_trans_count_mat = vector(mode = "list", length = N_cycles)
  global_AA_trans_freq_mat = vector(mode = "list", length = N_cycles)
  
  ### Function for barplots
  draw_bar = function(df, title, legend_title, xlabel) {
    bar = ggplot(df, aes(x = feature, y = count)) +
      geom_col() +
      xlab(xlabel) +
      ylab("count") +
      ggtitle(title)
    
    return(bar)
  }
  ### Function for drawing heatmaps
  draw_heat = function(mat, title, legend_title) {
    plot_dat = as.vector(mat) # Operation is performed column-wise.
    plot_df = data.frame(to = rep(colnames(mat), each = nrow(mat)),
                         from = rep(rownames(mat), times = ncol(mat)),
                         count = plot_dat)
    plot_df$from = factor(plot_df$from, levels = rownames(mat))
    plot_df$to = factor(plot_df$to, levels = colnames(mat))
    heat = ggplot(plot_df, aes(to, from)) +
      geom_tile(aes(fill = count)) +
      scale_fill_viridis(legend_title) +
      ggtitle(title)
    
    return(heat)
  }
  
  ### User information
  if (verbose) {
    print("Simulation started.")
  }
  
  ### Iterate
  for (c in 1:N_cycles) {
    
    ### Define entities to track mutation events
    ## Define transition matrices
    # DNA-level
    nt_trans_count_mat = matrix(rep(0, length(nt)^2),
                                nrow = length(nt),
                                ncol = length(nt))
    colnames(nt_trans_count_mat) = nt
    rownames(nt_trans_count_mat) = nt
    # AA-level
    AA_trans_count_mat = matrix(rep(0, length(AAs)^2),
                                nrow = length(AAs),
                                ncol = length(AAs))
    colnames(AA_trans_count_mat) = AAs
    rownames(AA_trans_count_mat) = AAs
    
    ### Simulation
    ## Generate new DNA sequences
    synth_S = c()
    for (S in 1:S_tot) {
      synth_S[S] = seq
      # For each base
      for (N_int in 1:nchar(seq)) {
        N = substr(seq, N_int, N_int)
        # Mutate it with given probability
        mutate = runif(1)
        if (mutate <= base_probs[N]) {
          M = sample(colnames(base_transition_freqs), 
                     1, 
                     prob = base_transition_freqs[N,])
          substr(synth_S[S], N_int, N_int) = M
          # Add the change to mutation transition count matrix
          nt_trans_count_mat[N, M] = nt_trans_count_mat[N, M] + 1
        }
      }
    }
    ## Translate all sequences
    synth_S = DNAStringSet(synth_S)
    synth_AA = translate(synth_S)
    ref_AA = translate(DNAString(seq))
    ## Count and tabulate AA changes
    # Iterate over the sequences
    for (sequenceNr in 1:length(synth_AA)) {
      sequence = synth_AA[sequenceNr]
      # Catalogue mismatches
      for (aaNr in 1:width(synth_AA[sequenceNr])) {
        # Get aa of a position
        from = as.character(subseq(ref_AA, start = aaNr, end = aaNr))
        to = as.character(subseq(sequence, start = aaNr, end = aaNr))
        # Store change
        AA_trans_count_mat[from, to] = AA_trans_count_mat[from, to] + 1
      }
    }
    ## Compute transition frequencies (row-normalised counts)
    nt_count_from = rowSums(nt_trans_count_mat)
    nt_trans_freq_mat = nt_trans_freq_mat = sweep(nt_trans_count_mat,
                                                  1,
                                                  nt_count_from,
                                                  "/")
    AA_count_from = rowSums(AA_trans_count_mat)
    AA_trans_freq_mat = sweep(AA_trans_count_mat,
                              1,
                              AA_count_from,
                              "/")
    
    ### Store all the matrices
    global_nt_trans_count_mat[[c]] = nt_trans_count_mat
    global_nt_trans_freq_mat[[c]] = nt_trans_freq_mat
    global_AA_trans_count_mat[[c]] = AA_trans_count_mat
    global_AA_trans_freq_mat[[c]] = AA_trans_freq_mat
    
    ### Get initial nt and AA distributions
    DNAseq = DNAString(seq)
    nt_init = letterFrequency(DNAseq, nt)
    AA_init = letterFrequency(ref_AA, AAs)
    nt_init = data.frame(feature = names(nt_init), count = nt_init)
    nt_init$feature = factor(nt_init$feature, levels = nt)
    AA_init = data.frame(feature = names(AA_init), count = AA_init)
    AA_init$feature = factor(AA_init$feature, levels = AAs)
    
    ### Visualisation
    if (plot_ind) {
      ## Barplots of initial nt and AA distributions
      # Generate barplots
      bar_nt_init = draw_bar(nt_init, 
                             "Initial nt counts", 
                             "count", 
                             "nt")
      bar_AA_init = draw_bar(AA_init, 
                             "Initial AA counts", 
                             "count", 
                             "AA")
      
      ## Generate heatmaps
      heat_count_DNA = draw_heat(nt_trans_count_mat, 
                                 "nt transition counts", 
                                 "count")
      heat_count_AA = draw_heat(AA_trans_count_mat, 
                                "AA transition counts", 
                                "count")
      heat_freq_DNA = draw_heat(nt_trans_freq_mat, 
                                "nt transition frequencies",
                                "frequency")
      heat_freq_AA = draw_heat(AA_trans_freq_mat, 
                               "AA transition frequencies",
                               "frequency")
      
      ## Panel
      graph = plot_grid(bar_nt_init,
                        bar_AA_init,
                        heat_count_DNA,
                        heat_count_AA,
                        heat_freq_DNA,
                        heat_freq_AA,
                        nrow = 3,
                        align = "v",
                        byrow = TRUE,
                        axis = "lr")
      
      filename = paste0(dir_out, "/sim_cycle_", c, ".png")
      ggsave(filename, 
             width = A4_full["width"], 
             height = A4_full["height"])
    }
  }
  
  ### Compute distributions
  if (verbose) {
    if (plot_ind) {
      print(paste0("Simulations completed. Summary plots of individual cycles have been stored in '", dir_out, "'."))
    }
    else {
      print("Simulations completed.")
    }
    print("Computing empirical transition distributions.")
  }
  
  ## Convert matrices to dataframes in long format
  # Function definition
  matlist_to_df = function(matlist) {
    df_tot = foreach(l = 1:length(matlist), .combine = rbind) %dopar% (
      data.frame(to = factor(rep(colnames(matlist[[l]]), each = nrow(matlist[[l]])), levels = colnames(matlist[[l]])),
                 from = factor(rep(rownames(matlist[[l]]), times = ncol(matlist[[l]])), levels = rownames(matlist[[l]])),
                 value = as.vector(matlist[[l]]), # operation column-wise
                 cycle = l) 
    )
    return(df_tot)
  }
  # Conversions
  if (n.cores == "auto") {
    n.cores <- parallel::detectCores() - 1
  }
  my.cluster <- parallel::makeCluster(
    n.cores, 
    type = "PSOCK" #only available option for Windows OS
  )
  doParallel::registerDoParallel(cl = my.cluster)
  
  if (verbose) {
    print(paste0("Parallel execution with ", foreach::getDoParWorkers(), " workers."))
    print("Converting nt transition count data.")
  }
  global_nt_trans_count_df = matlist_to_df(global_nt_trans_count_mat)
  if (verbose) {
    print("Converting nt transition frequency data.")
  }
  global_nt_trans_freq_df = matlist_to_df(global_nt_trans_freq_mat)
  if (verbose) {
    print("Converting AA transition count data.")
  }
  global_AA_trans_count_df = matlist_to_df(global_AA_trans_count_mat)
  if (verbose) {
    print("Converting AA transtion freqency data.")
  }
  global_AA_trans_freq_df = matlist_to_df(global_AA_trans_freq_mat)
  
  ## Create faceted histograms
  # Plot
  fac_hist = function(data, data_mark, filename, width, height, facet_lab_size, x_ax_text_size, axis_title_size, rotate_x_text = FALSE, density = TRUE) {
    graph = ggplot(data, aes(x = value)) +
      geom_histogram(aes(y = ..density..),
                     colour = "black",
                     fill = "White",
                     binwidth = 0.05)
    if (density) {
      graph = graph + geom_density()
    }
    graph = graph +
      facet_grid(from ~ to, switch = "y", scales = "free_y") +
      scale_y_continuous(position = "right") +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
      xlab("frequency") +
      expand_limits(x = c(0, 1)) +
      theme(strip.text.x = element_text(size = facet_lab_size),
            strip.text.y = element_text(size = facet_lab_size),
            axis.text.x = element_text(size = x_ax_text_size),
            axis.title = element_text(size = axis_title_size))
    if (rotate_x_text) {
      graph = graph +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
    }
    if (!is.null(data_mark)) {
      data_mark = matlist_to_df(list(data_mark))
      graph = graph +
        geom_vline(data = data_mark,
                   aes(xintercept = value),
                   colour = "red")
    }
    
    ggsave(filename, width = width, height = height)
  }
  theme_set(theme_bw(base_size = 12, base_line_size = 0.7))
  print("Generating nt transition frequency histograms.")
  filename = paste0(dir_out, "/histograms_nt.png")
  fac_hist(global_nt_trans_freq_df,
           data_mark = mark_nt_freqs,
           filename = filename,
           width = A4_full["width"],
           height = A4_full["width"]/2,
           facet_lab_size = 10,
           x_ax_text_size = 6,
           axis_title_size = 12)
  if (verbose) {
    print(paste0("Summary plot of nt transition frequency distributions stored as '", filename, "'."))
  }
  theme_set(theme_bw(base_size = 8, base_line_size = 0.4))
  print("Generating AA transition frequency histograms.")
  filename = paste0(dir_out, "/histograms_AA.png")
  fac_hist(global_AA_trans_freq_df,
           data_mark = mark_AA_freqs,
           filename = filename,
           width = A4_full["width"]*5,
           height = (A4_full["width"]/2)*5,
           facet_lab_size = 20,
           x_ax_text_size = 15,
           axis_title_size = 25,
           rotate_x_text = TRUE,
           density = FALSE)
  if (verbose) {
    print(paste0("Summary plot of AA transition frequency distributions stored as '", filename, "'."))
  }
  theme_set(theme_bw(base_size = 12, base_line_size = 0.7))
  ### Convergence of distributions: Fit binomial distribution to data up to that point after each cycle. Can only use for nt transition counts!
  ## Replace each entry in the count matrices with the fitted probability for a binomial distribution for all the data up to that point
  print("Generating plot to check for convergence.")
  tmp = global_nt_trans_count_mat
  for (id_from in 1:nrow(tmp[[1]])) {
    from = rownames(tmp[[1]])[id_from]
    counts_base = nt_init$count[nt_init$feature == from]
    for (id_to in 1:ncol(tmp[[1]])) {
      vec = c()
      for (l in 1:length(tmp)) {
        vec = append(vec, tmp[[l]][id_from, id_to])
        tmp[[l]][id_from, id_to] = goodfit(vec, 
                                           type = "binomial",
                                           par = list(size = counts_base * S_tot))$par$prob
      }
    }
  }
  ## Transform into dataframe
  probs_df = matlist_to_df(tmp)
  
  ## Plot
  graph = ggplot(probs_df, aes(x = cycle, y = value)) +
    geom_line() +
    facet_grid(from ~ to, switch = "y") +
    scale_y_continuous(position = "right") +
    xlab("Cycle") +
    ylab("cumulative P(binomial)")
  
  filename = paste0(dir_out, "/convergence.png")
  ggsave(filename, 
         width = A4_full["width"],
         height = A4_full["width"]/2)
  
  if (verbose) {
    print(paste0("If desired, check convergence in '", filename, "'."))
  }
  
  ### Reset ggplot2 theme
  theme_set(theme_classic(base_size = 12, base_line_size = 0.7))
  
  ### Create summary panel (counts summed, frequencies averaged)
  print("Generating summary panel of the simulation run.")
  ## Get matrices of summed counts
  # Function
  sum_matlist = function(matlist) {
    for (l in 2:length(matlist)) {
      matlist[[l]] = matlist[[l-1]] + matlist[[l]]
    }
    return(matlist[[l]])
  }
  # Compute sums
  global_nt_trans_count_mat = sum_matlist(global_nt_trans_count_mat)
  global_AA_trans_count_mat = sum_matlist(global_AA_trans_count_mat)
  
  ## Get matrices with averaged frequencies
  # Function
  avg_matlist = function(matlist) {
    for (l in 2:length(matlist)) {
      matlist[[l]] = matlist[[l-1]] + matlist[[l]]
    }
    return(matlist[[l]]/l)
  }
  # Compute averages
  global_nt_trans_freq_mat = avg_matlist(global_nt_trans_freq_mat)
  global_AA_trans_freq_mat = avg_matlist(global_AA_trans_freq_mat)
  
  ## Plots
  # Generate barplots
  nt_init$count = nt_init$count * S_tot * N_cycles
  bar_nt_init = draw_bar(nt_init, 
                         "Initial nt counts", 
                         "count", 
                         "nt")
  AA_init$count = AA_init$count * S_tot * N_cycles
  bar_AA_init = draw_bar(AA_init, 
                         "Initial AA counts", 
                         "count", 
                         "AA")
  
  # Generate heatmaps
  heat_count_DNA = draw_heat(global_nt_trans_count_mat, 
                             "nt transition counts", 
                             "count")
  heat_count_AA = draw_heat(global_AA_trans_count_mat, 
                            "AA transition counts", 
                            "count")
  heat_freq_DNA = draw_heat(global_nt_trans_freq_mat, 
                            "nt transition frequencies",
                            "frequency")
  heat_freq_AA = draw_heat(global_AA_trans_freq_mat, 
                           "AA transition frequencies",
                           "frequency")
  
  # Panel
  graph = plot_grid(bar_nt_init,
                    bar_AA_init,
                    heat_count_DNA,
                    heat_count_AA,
                    heat_freq_DNA,
                    heat_freq_AA,
                    nrow = 3,
                    align = "v",
                    byrow = TRUE,
                    axis = "lr")
  
  filename = paste0(dir_out, "/sim_summary.png")
  ggsave(filename, 
         width = A4_full["width"], 
         height = A4_full["height"])
  
  if (verbose) {
    print(paste0("A simulation summary plot has been created '", filename, "'."))
  }
  
  ### Save dataframes for import for future work
  filename = paste(dir_out, "sim_results.rda", sep = "/")
  save(global_nt_trans_count_mat,
       global_nt_trans_freq_mat,
       global_AA_trans_count_mat,
       global_AA_trans_freq_mat,
       global_nt_trans_count_df,
       global_nt_trans_freq_df,
       global_AA_trans_count_df,
       global_AA_trans_freq_df,
       file = filename)
  if (verbose) {
    print(paste0("Simulation completed. Matrices and dataframes with the simulation results have been stored as '", filename, "'."))
  }
  
  ### Close parallelisation backend
  parallel::stopCluster(cl = my.cluster)
  if (verbose) {
    print("Parallel backend closed.")
  }
}