# Summarize and plot merged SVs from short reads 
# This script is made to be executed by the 01_scripts/utils/summarize_plot.sh script
library(ggplot2)
library(ggpubr)
library(data.table)
library(scales)

# 1. Access files from command line ---------------------------------------
argv <- commandArgs(T)
MERGED <- argv[1]
FILT <- argv[2]
CALLER1 <- argv[3]
CALLER2 <- argv[4]
CALLER3 <- argv[5]


# 2. Explore raw merged SV calls ------------------------------------------

# 2.1 Import and format ---------------------------------------------------
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



# 2.3 Summarize SV counts in tables ---------------------------------------
## Number of merged SVs
nrow(merged)
## By SVTYPE
table(merged$SVTYPE)
## By callers
table(merged$callers)

## By number of supporting tools
table(merged$SUPP)
table(merged$SUPP, merged$SVTYPE)

for (i in c(CALLER1, CALLER2, CALLER3)){
  nb_calls <- nrow(merged[merged[i] == 1, ])
  print(paste(nb_calls, 'SVs called by', i))
  print(table(merged[merged[i] == 1, 'SVTYPE']))
}

# 2.4 Plot SV counts ------------------------------------------------------

# Plot SV counts by type and size for each caller separately
for (i in c(CALLER1, CALLER2, CALLER3)){
  caller_plot <- 
    ggplot(data = merged[merged[i] == 1, ]) + 
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
      fill = "SV type",
      title = paste(i)
    ) + 
    scale_fill_manual(values = cols_svtypes)
  
  ## Write to env
  assign(x = paste0(i, '_plot'), caller_plot)
  
  ### Export plot as rds object
  #saveRDS(get(paste0(i, '_plot')), file = paste0(strsplit(MERGED, '.table')[[1]], '_barplot_', i, '.rds'))
  
}

#### On superdome, loading the .rds files throws an error, 
#### but no issue if downloaded locally and opened on PC

# Combine the 3 single-caller plots in a single figure
single_callers_3 <- ggarrange(get(paste0(CALLER1, '_plot')), 
                              get(paste0(CALLER2, '_plot')),
                              get(paste0(CALLER3, '_plot')), 
                              ncol = 2, nrow = 2, common.legend = TRUE, legend = "right")
single_callers_3

## Save to rds object
saveRDS(single_callers_3, file = paste0(strsplit(MERGED, '.table')[[1]], '_3single_callers_barplot.rds'))


# Plot SV counts by caller(s) combination
## First reorder
combi <- unique(merged$callers)

reordered_callers <- c(CALLER1, CALLER2, CALLER3,
                       combi[! combi %in% c(CALLER1, CALLER2, CALLER3)]
)

merged$callers_reordered <- factor(merged$callers, 
                           levels = reordered_callers)
                          

merged_byCALLERS <- 
ggplot(data = merged) +
  geom_bar(aes(x = SVLEN_bin, fill = SVTYPE)) +
  facet_wrap(~ callers_reordered, scales = 'free_y') +
  #facet_grid(rows = vars(SVTYPE)) +
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
    fill = "SV type",
    title = "SV count by caller(s)"
  ) + 
  scale_fill_manual(values = cols_svtypes)

merged_byCALLERS

## Save to rds object
saveRDS(merged_byCALLERS, 
        file = paste0(strsplit(MERGED, '.table')[[1]], 'merged_callers_byCALLERS.rds'))

# Plot SV count by SVTYPE and callers
merged_bySVTYPE <- 
ggplot(data = merged) +
  geom_bar(aes(x = SVLEN_bin, fill = callers)) +
  facet_wrap(~SVTYPE, scales = 'free_y') +
  #facet_grid(rows = vars(SVTYPE)) +
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
    fill = "Caller(s)",
    title = "SV count by type"
  ) + 
  scale_fill_viridis_d(option = "C")

merged_bySVTYPE

## Save to rds object
saveRDS(merged_bySVTYPE, 
        file = paste0(strsplit(MERGED, '.table')[[1]], 'merged_callers_bySVTYPE_barplot.rds'))



# 2.5 Summarize and plot calls supported by at least 2 tools --------------
# Extract calls where SUPP > 1
merged_SUPP2 <- subset(merged, SUPP != 1)

# What fraction of raw SVs are shared by multiple tools ?
nrow(merged_SUPP2)/nrow(merged)
table(merged_SUPP2$SVTYPE)

# Plot 
merged_byCALLERS_SUPP2 <- 
  ggplot(data = merged_SUPP2) +
  geom_bar(aes(x = SVLEN_bin, fill = SVTYPE)) +
  facet_wrap(~ callers, scales = 'free_y') +
  #facet_grid(rows = vars(SVTYPE)) +
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
    fill = "SV type",
    title = "SV count by callers combination for calls supported by > 1 tool"
  ) + 
  scale_fill_manual(values = cols_svtypes)

merged_byCALLERS_SUPP2

## Save to rds object
saveRDS(merged_byCALLERS_SUPP2, 
        file = paste0(strsplit(MERGED, '.table')[[1]], 'merged_callers_byCALLERS_SUPP2_barplot.rds'))


merged_bySVTYPE_SUPP2 <- 
  ggplot(data = merged_SUPP2) +
  geom_bar(aes(x = SVLEN_bin, fill = callers)) +
  facet_wrap(~SVTYPE, scales = 'free_y') +
  #facet_grid(rows = vars(SVTYPE)) +
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
    fill = "Callers combination",
    title = "SV count by type for calls supported by > 1 tool"
  ) + 
  scale_fill_viridis_d(option = "C")

merged_bySVTYPE_SUPP2

## Save to rds object
saveRDS(merged_bySVTYPE_SUPP2, 
        file = paste0(strsplit(MERGED, '.table')[[1]], 'merged_callers_bySVTYPE_SUPP2_barplot.rds'))



# 3. Explore filtered SV calls --------------------------------------------

# 3.1 Import and format ---------------------------------------------------
# Import and assign required class to each variable
filtered <- read.table(FILT, header = FALSE, 
                     col.names = c('CHROM', 'POS', 'ID', 'SVTYPE', 'SVLEN',
                                   'END', 'SUPP', 'SUPP_VEC'),
                     colClasses = c('character')
)

filtered[, c('POS', 'SVLEN', 'END', 'SUPP')] <- sapply(filtered[, c('POS', 'SVLEN', 'END', 'SUPP')],
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

filtered$SVLEN_bin <-
  cut(abs(filtered$SVLEN), breaks = SVLEN_breaks, labels = SVLEN_names, right = FALSE)


# Add explicit caller names
filtered$callers <-
  sapply(X = filtered$SUPP_VEC, FUN = function(x){
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
  filtered[i] <- ifelse(test = grepl(i, filtered$callers),
                      yes = 1,
                      no = 0)
}


# 3.2 Summarize SV counts in tables ---------------------------------------
# By SVTYPE
table(filtered$SVTYPE)
## By callers
table(filtered$callers)
table(filtered$callers, filtered$SVTYPE)



# 3.3 Plot SV counts ------------------------------------------------------

# Plot SV counts by caller(s) combination
filtered_byCALLERS <- 
  ggplot(data = filtered) +
  geom_bar(aes(x = SVLEN_bin, fill = SVTYPE)) +
  facet_wrap(~ callers, scales = 'free_y') +
  #facet_grid(rows = vars(SVTYPE)) +
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
    fill = "SV type",
    title = "Filtered SV count by caller(s)"
  ) + 
  scale_fill_manual(values = cols_svtypes)

filtered_byCALLERS

## Save to rds object
saveRDS(filtered_byCALLERS, 
        file = paste0(strsplit(FILT, '.table')[[1]], 'filtered_callers_byCALLERS,', '.rds'))

# Plot SV count by SVTYPE and callers
filtered_bySVTYPE <- 
  ggplot(data = filtered) +
  geom_bar(aes(x = SVLEN_bin, fill = callers)) +
  facet_wrap(~SVTYPE, scales = 'free_y') +
  #facet_grid(rows = vars(SVTYPE)) +
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
    fill = "Caller(s)",
    title = "Filtered SR SV count by type"
  ) + 
  scale_fill_viridis_d(option = "C")

filtered_bySVTYPE

## Save to rds object
saveRDS(filtered_bySVTYPE, 
        file = paste0(strsplit(FILT, '.table')[[1]], 'filtered_callers_bySVTYPE_barplot.rds'))

