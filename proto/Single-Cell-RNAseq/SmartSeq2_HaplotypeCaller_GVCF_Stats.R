

#!/usr/bin/env Rscript

library(Cairo)
library(setwidth)
grDevices::X11.options(type='cairo')
options(device='x11')
options(bitmapType='cairo')

args <- commandArgs(trailingOnly = TRUE)

input_file <- args[1]
output_file <- args[2]

dat <- read.delim(input_file, na.strings = c("", "NA", "NaN", "."),check.names = FALSE, stringsAsFactors = FALSE)
annotations     <- c("QD","FS","SOR","MQ","MQRankSum","ReadPosRankSum")
missing_columns <- setdiff(annotations, colnames(dat))

if (length(missing_columns) > 0) { stop( "Missing columns: ", paste(missing_columns, collapse = ", ") ) }

# GATK generic SNP hard-filter reference thresholds.
thresholds <- list(
    QD             = list(value = 2.0, direction = "low"),
    FS             = list(value = 60.0, direction = "high"),
    SOR            = list(value = 3.0, direction = "high"),
    MQ             = list(value = 40.0, direction = "low"),
    MQRankSum      = list(value = -12.5, direction = "low"),
    ReadPosRankSum = list(value = -8.0, direction = "low")
)

# Convert annotation columns safely to numeric.
for (field in annotations) {
    dat[[field]] <- suppressWarnings(as.numeric(dat[[field]]))
}

cat("Total biallelic SNPs:", nrow(dat), "\n\n")

summary_table <- do.call(
    rbind,
    lapply(annotations, function(field) {
        x <- dat[[field]]
        threshold <- thresholds[[field]]$value
        direction <- thresholds[[field]]$direction
        fail <- if (direction == "low") { x < threshold }else{ x > threshold }

        data.frame(
            Annotation      = field,
            Total_sites     = length(x),
            Non_missing     = sum(!is.na(x)),
            Missing         = sum(is.na(x)),
            Missing_percent = round(mean(is.na(x)) * 100, 2),
            Median = if (all(is.na(x))) NA else median(x, na.rm = TRUE),
            Q01 = if (all(is.na(x))) NA else unname(quantile(x, 0.01, na.rm = TRUE)),
            Q05 = if (all(is.na(x))) NA else unname(quantile(x, 0.05, na.rm = TRUE)),
            Q95 = if (all(is.na(x))) NA else unname(quantile(x, 0.95, na.rm = TRUE)),
            Q99 = if (all(is.na(x))) NA else unname(quantile(x, 0.99, na.rm = TRUE)),
            Reference_cutoff = threshold,
            Failed_sites     = sum(fail, na.rm = TRUE),
            Failed_percent_nonmissing = if (sum(!is.na(x)) == 0) { NA 
            }else{ 
                round(sum(fail, na.rm = TRUE) / sum(!is.na(x)) * 100, 2)
            },
            stringsAsFactors = FALSE
        )
    })
)

#summary_file <- sub( "\\.png$", ".summary.tsv", output_png )
write.table( summary_table, file=output_file, sep = "\t", quote = FALSE, row.names = FALSE )

output_png <- gsub("tsv$", "png", output_file)
png(output_png,  width=18, height=12, units='in', res=150)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
for (field in annotations) {
    x <- dat[[field]]
    x <- x[is.finite(x)]

    if (length(x) < 2) { 
        plot.new()
        title(main = paste(field, "(insufficient data)"))
        next
    }
    # Restrict only the plotted x-axis to the central 99.6%.
    # Full data remain in the summary table.
    x_limits <- quantile( x, probs = c(0.002, 0.998), na.rm = TRUE, names = FALSE )
    x_plot   <- x[x >= x_limits[1] & x <= x_limits[2]]

    #if (length(unique(x_plot)) >= 2) {
        density_result <- density(x_plot, na.rm = TRUE)
        plot( density_result, main = field, xlab = field, ylab = "Density", lwd = 2 )
    #}else{
    #    hist( x_plot, breaks = 30, main = field, xlab = field )
    #}
    abline( v = thresholds[[field]]$value, lty = 2, lwd = 2 )
    legend( "topright", legend = paste0("Reference cutoff: ", thresholds[[field]]$value ),
        lty = 2, lwd = 2, bty = "n", cex = 0.8
    )
}
dev.off()

#print(summary_table)
