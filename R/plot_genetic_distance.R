#' Plot Observed Genetic Distance vs Panmictic Expectation
#'
#' Plots the observed pairwise genetic distance distribution against the
#' panmictic expectation. Supports genome-wide, per-chromosome, panel,
#' and batch plots from \code{distances_by_chr()} or \code{simulate_panmixia()}
#' results.
#'
#'
#' @param result      A list returned by \code{simulate_panmixia()} or
#' \code{distances_by_chr()}. 
#' @param bw          Numeric. Bandwidth for kernel density estimation.
#'                    Default: \code{0.6}.
#' @param binwidth    Numeric. Bin width for the observed distance histogram. 
#'                    Default: \code{0.7}.
#' @param chr         Character or integer. Name or index of a single
#'                    chromosome to plot. 
#'                    Requires a \code{distances_by_chr()} result. 
#'                    Use \code{"all"} to generate one plot per chromosome.
#'                    Use a chromosome name (\code {"Chromosome 01"}) or index
#'                    (\code{1}) to plot a single chromosome. 
#'                    Default: \code{NULL}.
#' @param title       If \code{TRUE}, adds an auto-generated title.
#'                    Default: \code{FALSE}.
#' @param distribution Character. 
#'                    Default: \code{"both"} plots the observed pairwise
#'                    genetic distances and the panmictic null distribution.
#'                    \code{"observed"} plots only the ovserved distances/
#'                    \code{"panmictic"} plots only the panmictic null
#'                    distribution.
#' @param xlim        Optional \code{c(min, max)} for x-axis limits (%). 
#'                    If \code{NULL}, inferred from data.
#' @param save        If \code{TRUE}, saves the plot to \code{filename}.
#'                    Default: \code{FALSE}.
#' @param filename    Output filename. 
#'                    Default: \code{"plot_distance"}
#'
#'
#' @return A \code{ggplot} object.
#'
#' @details
#' Three plotting modes are available, evaluated in the following order:
#' \enumerate{
#'   \item \strong{Default}: Genome plot.
#'   \item \strong{Single chromosome} (\code{chr}): Plot for one chromosome,
#'     selected by name or index.
#'   \item \strong{All} (\code{chr = "all"}): Generates one plot
#'     per entry (each chromosome), printing and optionally saving each
#'     individually.
#' }
#'
#'
#' @importFrom ggplot2 ggplot aes geom_histogram geom_line coord_cartesian geom_density
#' @importFrom ggplot2 scale_fill_manual scale_color_manual labs theme_classic theme
#' @importFrom ggplot2 element_text element_blank unit after_stat
#'
#'
#' @export


plot_genetic_distance <- function(
    result,
    bw               = 0.2,
    binwidth         = 0.7,
    chr              = NULL,
    title            = FALSE,
    distribution     = c("both", "observed", "panmictic"), 
    xlim             = NULL,
    save             = FALSE,
    filename         = "plot_distance"
) {
  
  distribution = match.arg(distribution)
  
# Prepare observed and panmictic data for plotting: 
# Converts raw distances to percentages and estimates the panmictic density curve
  .prepare_plot_data <- function(result, xlim = NULL) {
    observed_data  <- result[[1]] * 100        # Observed distances (%)
    panmictic_data <- result[[2]] * 100        # Panmictic distances (%)
    
    
    # Set axis limits from data range or user-supplied xlim
    x_lower <- if (is.null(xlim)) 0 else xlim[1]
    x_upper <- if (is.null(xlim)) max(c(observed_data, panmictic_data)) else xlim[2]
    
    # Clip observed values to axis range
    observed_data <- observed_data[observed_data >= x_lower & observed_data <= x_upper]
    
    # Kernel density estimate for the panmictic expectation curve
#    panmictic_density <- density(
#      panmictic_data,
#      bw   = bw,
#      adjust = 0.4,
#      from = x_lower,
#      to   = x_upper,
#      n    = 256
#    )
    
    list(
      observed_df  = data.frame(distance = observed_data),
      panmictic_df = data.frame(distance = panmictic_data),
      x_lower      = x_lower,
      x_upper      = x_upper
    )
  }

  
# Plot:
# Combines a histogram of observed distances with the panmictic density curve
  .build_single_plot <- function(result, chr_label = NULL,
                                 title       = FALSE) {
    plot_data <- .prepare_plot_data(result, xlim)
    
    p <- ggplot(plot_data$observed_df, aes(x = distance))
      
    # Histogram of observed pairwise distances
    if (distribution %in% c("both", "observed")) {
      p <- p + geom_histogram(
      aes(y = after_stat(density), fill = "Observed"),
      binwidth = binwidth, color = "black", alpha = 0.6
    )
    } 
        
    # Panmictic expectation density curve
    if (distribution %in% c("both", "panmictic")) {
      p <- p +
        geom_density(
        data = plot_data$panmictic_df,
        aes(x = distance, color = "Panmictic"),
        bw = bw,
        inherit.aes = FALSE,
        linewidth = 0.7
      )
    } 
    
    p <- p + 
      coord_cartesian(xlim = c(plot_data$x_lower, plot_data$x_upper)) +
      labs(
        title = .build_plot_title(chr_label),
        x     = "Relative genetic distance (%)",
        y     = "Density"
      ) +
      theme_classic(base_family = "Times", base_size = 14) +
      theme(
        plot.title           = element_text(hjust  = 0.5, 
                                            face   = "bold",
                                            family = "Times", 
                                            size   = 12),
        legend.position      = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.spacing.y     = unit(0, "pt"),
        legend.text          = element_text(size = 12),
        axis.line.y          = element_blank(),
        axis.text.y          = element_blank(), 
        axis.ticks.y         = element_blank(),
        axis.text            = element_text(size = 14),
        axis.title           = element_text(size = 14)
      )
      
    if (distribution %in% c("both", "observed")) {
      p <- p + scale_fill_manual(values  = c("Observed"  = "skyblue3"), name = NULL)
    }
    
    if (distribution %in% c("both", "panmictic")) {
      p <- p + scale_color_manual(values = c("Panmictic" = "#E67E22"),  name = NULL)
    }
    
    p
  }

  
# Build the plot title: 
# Returns NULL (no title), a custom string or an auto-generated label depending on the 'title' 
  .build_plot_title <- function(chr_label = NULL) {
    if (!is.null(chr_label)){            
    
    # Auto-generate a label based on the chromosome argument
    scope_label <- if (is.numeric(chr_label)) paste("Chromosome", chr_label) else chr_label
    return(bquote(bold(.(scope_label))
    ))
    }
    
    if (!isTRUE(title)) return(NULL)
    bquote(bold("Genome")
    )
  } 
  
  # Save a finished plot:
  # Builds the filename from 'filename' and an optional label suffix (e.g. chromosome name)
  .save_plot_to_file <- function(plot_object, label = NULL) {
    suffix    <- if (!is.null(label)) paste0("_", gsub(" ", "_", label)) else ""
    file_path <- paste0(filename, suffix, ".tiff")
    
    if (is.null(file_path) || nchar(file_path) == 0)
      stop("Invalid file path: ", file_path)
    
    dir_path <- dirname(file_path)
      if (dir_path != ".") { 
        dir.create(dir_path, showWarnings = FALSE, recursive = TRUE)
      }
    
    ggplot2::ggsave(
      file_path, plot_object,
      width  = 7, 
      height = 4, 
      units  = "in",
      dpi    = 600,
      scale  = 1.5, 
      compression = "lzw"
    )
    message("Saved: ", file_path)
  }   


  
# 'chr': 
# Single or all chromosomes 
  if (!is.null(chr)) {
    if (is.null(result$chromosomes))
      stop("Chromosomes are NULL, requires a result from genetic_distribution(vcf_file, chr = 'by_chr').")
    
    chr_labels <- names(result$chromosomes)
    
    
    # chr = "all": One plot per chromosome, saved individually
    if (chr == "all") {
      message("Generating plots...")
      
      chr_entries      <- Filter(Negate(is.null), result$chromosomes)
      plot_list        <- vector("list", length(chr_entries))
      names(plot_list) <- names(chr_entries)
      
      for (entry_label in names(chr_entries)) {
        p <- .build_single_plot(chr_entries[[entry_label]],
                                chr_label = entry_label,
                                title     = TRUE)
        
        plot_list[[entry_label]] <- p
        print(p)
        
        if(save) .save_plot_to_file(p, 
                                    label = entry_label)
      }
      
      return(invisible(plot_list))
    }
    
  # Single chromosome: 
    # Accept chromosome by index (integer) or by name (string)
    target_chr <- if (is.numeric(chr)) {
      if (chr < 1 || chr > length(chr_labels))
        stop("chr = ", chr, " out of the range. Available:\n",
             paste(chr_labels, collapse = ", "))
      chr_labels[chr]
    } else {
      if (!chr %in% chr_labels)
        stop("chr '", chr, "' not found. Available:\n",
             paste(chr_labels, collapse = ", "))
      chr
    }
    
    p <- .build_single_plot(result$chromosomes[[target_chr]], 
                            chr_label = target_chr,
                            title     = TRUE)
    
    if (save) .save_plot_to_file(p, label = target_chr)
    return(p)
  }

  
# Default:
# Genome plot
  genome_entry <- if (!is.null(result$chromosomes)) result$genome else result
  
  p <- .build_single_plot(genome_entry, 
                          chr_label = NULL)
  
  if (save) .save_plot_to_file(p)
  return(p)
}

