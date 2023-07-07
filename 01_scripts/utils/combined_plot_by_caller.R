

MERGED <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/01_short_reads/02_SV/SVs_short_reads/08_merged_union/merged_delly_manta_smoove.table"
FILT <- "/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/01_short_reads/02_SV/SVs_short_reads/09_filtered_union/merged_delly_manta_smoove_SUPP2.table"

CALLER1 <- 'delly'
CALLER2 <- 'manta'
CALLER3 <- 'smoove'

# Import and assign required class to each variable
merged <- read.table(MERGED, header = FALSE, 
                     col.names = c('CHROM', 'POS', 'ID', 'SVTYPE', 'SVLEN',
                                   'END', 'SUPP', 'SUPP_VEC'),
                     colClasses = c('character')
)

merged[, c('POS', 'SVLEN', 'END', 'SUPP')] <- sapply(merged[, c('POS', 'SVLEN', 'END', 'SUPP')],
                                                     as.numeric)

# Convert SV lengths to num and bins
SVLEN_breaks <- c(-Inf, 50, 100, 250, 500, 1000, 2500, 5000, 10000, Inf)
SVLEN_names <- c('[0-50[',
                 '[50-100[',
                 '[100-250[',
                 '[250-500[',
                 '[500-1000[',
                 '[1000-2500[',
                 '[2500-5000[',
                 '[5000-10000[',
                 '[10000+')

merged$SVLEN_bin <-
  cut(abs(merged$SVLEN), breaks = SVLEN_breaks, labels = SVLEN_names, right = FALSE)


# Add explicit caller names
merged$callers <-
  sapply(X = merged$SUPP_VEC, FUN = function(x){
    switch(x,
           '100' = CALLER1,
           '010' = CALLER2,
           '001' = CALLER3,
           '110' = paste0(CALLER1, ' + ', CALLER2),
           '101' = paste0(CALLER1, ' + ', CALLER3),
           '011' = paste0(CALLER2, ' + ', CALLER3),
           '111' = paste0(CALLER1, ' + ', CALLER2, ' + ', CALLER3))
  }
  )

# Add a variable for each caller
for (i in c(CALLER1, CALLER2, CALLER3)){
  merged[i] <- ifelse(test = grepl(i, merged$callers),
                      yes = 1,
                      no = 0)
}

# Make one dataset per caller
library(tidyr)
merged_long <- 
pivot_longer(data = merged, 
             cols = c('delly', 'manta', 'smoove'),
             values_to = 'single_caller_yes',
             names_to = 'single_caller'
             )

merged_long_bycaller <- subset(merged_long, single_caller_yes == 1)



# 2.2 Choose fun color scheme for sv types --------------------------------
# Get SVTYPE values
svtypes <- sort(unique(merged$SVTYPE)) 
#### we sort so that INV falls at the end of vector and 
#### is assigned the most divergent color from DELs, 
#### as INVs are rare and hard to distinguish bar plots

# Get hex code for as many colors as svtypes for a given viridis palette
hex_svtypes <- viridisLite::viridis(n = length(svtypes), option = 'D')
show_col(hex_svtypes)

# Assign a color to each svtype in a named vector
cols_svtypes <- vector(mode = 'character', length = length(svtypes))
for (i in 1:length(svtypes)) {
  names(cols_svtypes)[i] <- svtypes[i]
  cols_svtypes[i] <- hex_svtypes[i]
}


ggplot(data = merged_long_bycaller) + 
  facet_grid(factor(single_caller, levels = c('delly', 'manta', 'smoove')) ~ ., 
             scales = 'free_y') +
  geom_bar(aes(x = SVLEN_bin, fill = SVTYPE)) + 
  theme(
    ## Plot title
    plot.title = element_text(size = 10, face = 'bold', hjust = 0.5),
    ## Axis
    axis.text.x = element_text(angle = 45, size = 6, hjust = 1),
    axis.text.y = element_text(size = 6, hjust = 1),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    ## Legend
    legend.title = element_text(size = 8, hjust = 0.5),
    legend.text = element_text(size = 7),
    legend.key.size = unit(5, 'mm')
  ) +
  labs(
    x = "SV size (bp)",
    y = "SV count",
    fill = "SV type"
  ) + 
  scale_fill_manual(values = cols_svtypes)
