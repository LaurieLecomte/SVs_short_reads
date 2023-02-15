# Format sample names in merged VCF i.e. remove leading 0_, 1_, 2_, ...

argv <- commandArgs(T)
NAMES_ORI <- argv[1]
CALLER1 <- argv[2]
CALLER2 <- argv[3]
CALLER3 <- argv[4]
OUTPUT <- argv[5]

# 1. Import original names
names_ori <- read.table(NAMES_ORI, col.names = 'ori')

# 2. Split original names : extract first character representing caller 
names_ori$caller_num <- sapply(X = strsplit(x = as.character(names_ori$ori), split = '_'), FUN = '[', 1) 

# 3. Split original names : extract sample name
names_ori$ID <- sapply(X = strsplit(x = as.character(names_ori$ori), split = "_"), FUN = '[', 2) 

# 4. Replace 0,1,2 by caller's name
names_ori$caller_name <- 
  sapply(names_ori$caller_num, FUN = function(x) {
    switch(x, 
           "0" = as.character(CALLER1), 
           "1" = as.character(CALLER2),
           "2" = as.character(CALLER3))
  })

# 5. Paste sample ID and caller's name together
names_ori$final <- paste0(names_ori$ID, '_', names_ori$caller_name)

# 6. Export table
write.table(names_ori[, c('ori', 'final')], OUTPUT, sep = "\t", col.names = F, row.names = F, quote = F)