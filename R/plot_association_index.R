#' Plot the Standardised Index of Association (rd)
#'
#' Plots the null distribution of the standardised index of association
#' (rd) against the panmictic expectation (rd = 0). 
#' Supports genome-wide, per-chromosome, panel, and batch plots from 
#' \code{distances_by_chr()} or \code{simulate_panmixia()} results.
#'
#'
#' @param result      A list returned by \code{simulate_panmixia()} or
#'                    \code{distances_by_chr()}.
#' @param sample_size Integer. Number of SNPs sampled per iteration.
#'                    Auto-detected from \code{result} when available. 
#'                    Default: \code{1000L}.
#' @param statistics  If \code{TRUE}, adds median, SD and p-value to the plot.
#'                    Default: \code{FALSE}.
#' @param chr         Character or integer. Name or index of a single chromosome to plot. 
#'                    Requires a \code{distances_by_chr()} result. 
#'                    Default: \code{NULL}.
#' @param title       If \code{TRUE}, adds an auto-generated title.
#'                    Default: \code{FALSE}.
#' @param save        If \code{TRUE}, saves the plot to \code{filename}.
#'                    Default: \code{FALSE}.
#' @param filename    Output filename. 
#'                    Default: \code{"plot_sIA"}
#'
#'
#'@return A \code{ggplot} object.
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
#' @importFrom ggplot2 ggplot aes geom_line geom_vline annotate coord_cartesian
#' @importFrom ggplot2 scale_x_continuous labs theme_classic theme element_text element_blank facet_wrap
#' 
#'
#' @export


plot_association_index <- function(
    result,
    sample_size  = 1000L,
    statistics   = TRUE,
    chr          = NULL,
    title        = FALSE,
    save         = FALSE,
    filename     = "plot_rd"
) {
  
  
# Auto-detect sample_size from result object when available
  if (!is.null(result$sample_size)) {
    sample_size <- result$sample_size
  } else if (!is.null(result$genome$sample_size)) {
    sample_size <- result$genome$sample_size
  }
  
  
# Compute rd statistics from a result entry:
  # Returns mean, median, SD, p-value and x-axis range for plotting
  .calc_stats <- function(result) {
    
    rd_values <- result[[3]]
    rd_values <- rd_values[!is.na(rd_values)]
    
    rd_mean <- mean(rd_values)
    rd_sd   <- sd(rd_values)

    if (is.na(rd_sd)) 
      stop("Just one value to calculate, therefore rd is indertemined as the variances are 0.")
    
    # P-value: 
    rd_p <- pnorm(0, mean = rd_mean, sd = rd_sd, lower.tail = TRUE)
    
    list(
      values = rd_values,
      mean   = rd_mean,
      median = median(rd_values),
      sd     = rd_sd,
      p      = rd_p,
      x_upper = rd_mean + 4 * rd_sd,
      x_lower = min(0, rd_mean - 4 * rd_sd)
    #  x_upper = quantile(rd_values, 0.995, na.rm = TRUE),
    #  x_lower = quantile(rd_values, 0.005, na.rm = TRUE)
    )
  }
  
# Format p-value for display:
# Returns a list with a formatted label and a significance flag
  .format_pvalue <- function(p) {
    if (is.null(p) || length(p) == 0 || is.na(p)) return("NA")
    if (p < 0.001) format(p, scientific = TRUE, digits = 2) else as.character(round(p, 4))
  }
  
  
# rd plot:
# Gaussian curve fitted to rd values with a vertical panmixia reference line

  .build_single_plot <- function(stats, 
                                 chr_label   = NULL,
                                 show_x_axis = TRUE,
                                 title       = TRUE) {
    pf      <- .format_pvalue(stats$p)
    x_range <- seq(stats$x_lower, stats$x_upper, length.out = 1000)
    
    # Fitted normal distribution over the rd range
    gaussian_df <- data.frame(
      x = x_range,
      y = dnorm(x_range, mean = stats$mean, sd = stats$sd)
    )
    
    p <- ggplot(data.frame(values = stats$values), aes(x = values)) +
      # Fitted Gaussian curve
      geom_line(
        data = gaussian_df, aes(x = x, y = y),
        color = "#3A7EBF", linewidth = 0.7, inherit.aes = FALSE
      ) +
      # Vertical reference line at rd = 0 (panmixia)
      geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
      annotate(
        "text", x = 0, y = Inf, vjust = 1.5, hjust = -0.1,
        label = "Panmixia~(bar(r)[d] == 0)", parse = TRUE,
        color = "red", size = 3.8, family = "Times"
      ) +
      scale_x_continuous(expand = c(0.02, 0)) +
      coord_cartesian(xlim = c(0, stats$x_upper)) +
      labs(
        title = .build_plot_title(chr_label),
        x     = if (show_x_axis) expression(bar(r)[d]) else NULL,
        y     = "Density"
      ) +
      theme_classic(base_family = "Times") +
      theme(
        plot.title   = element_text(hjust = 0.5, face = "plain",
                                    family = "Times", size = 12),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y  = element_blank(),
        axis.title   = element_text(size = 14),
        axis.text    = element_text(size = 14)
      )
    
    # Optionally overlay summary statistics on the plot
    if (statistics)
      p <- p + annotate(
        "text",
        x = stats$x_upper * 0.85, y = Inf,
        hjust = 0, vjust = 1.4, size = 3.5,
        family = "Times", color = "grey20", lineheight = 1.3,
        label = paste0(
          "Median:  ", round(stats$median, 6), "\n",
          "SD:         ", round(stats$sd,     6), "\n",
          "p-value: ",   pf
        )
      )
    p
  }
  
# Build the plot title:
# Returns NULL (no title), a custom string or an auto-generated label depending on the 'title' and 'custom_title' arguments
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
  
  
  # Print summary statistics to the console
  .print_stats <- function(stats, label = NULL) {
    pf <- .format_pvalue(stats$p)
    cat(sprintf(
      "%-15s | Median: %10.6f | SD: %10.6f | p-value: %s\n",
      label, stats$median, stats$sd, pf
    ))
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
        stats                  <- .calc_stats(chr_entries[[entry_label]])
        .print_stats(stats, entry_label)
        
        p <- .build_single_plot(stats,
                                chr_label = entry_label,
                                title     = TRUE)
        
        plot_list[[entry_label]] <- p
        print(p)
        
        if(save) .save_plot_to_file(p, 
                                    label = entry_label)
      }
      
      return(invisible(plot_list))
    }
    
    # Accept chromosome by index (integer) or by name (string)
    target_chr <- if (is.numeric(chr)) {
      if (chr < 1 || chr > length(chr_labels))
        stop("chr = ", chr, " out of range. Available: 1 to ", length(chr_labels))
      chr_labels[chr]
    } else {
      if (!chr %in% chr_labels)
        stop("chr '", chr, "' not found. Available:\n",
             paste(chr_labels, collapse = ", "))
      chr
    }

    if (is.null(result$chromosomes[[target_chr]])) stop(target_chr, " is NULL.")

    stats <- .calc_stats(result$chromosomes[[target_chr]])
    .print_stats(stats, target_chr)
    p <- .build_single_plot(stats, 
                            chr_label = target_chr,
                            title     = TRUE)
    
    if (save) .save_plot_to_file(p, 
                                 label = target_chr)
    
    return(p)
  }
  
  
# Default:
# Genome plot
  if (is.null(result$genome)) {
     for(it in seq_along(result$chromosomes)) {
       if(!is.null(result$chromosomes[[it]])){
          plot_association_index(result, chr = it)
          return(invisible(NULL))
       }
     }
  }

  stats <- .calc_stats(result$genome)
  .print_stats(stats, "Genome")
  p <- .build_single_plot(stats, chr_label = NULL)
  if (save) .save_plot_to_file(p)
  return(p)
}
