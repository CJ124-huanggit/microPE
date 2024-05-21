
## Load necessary R packages
library(lme4)
library(car)

### Read bacterial data
# Read the OTU table from the specified path
otutab <- read.csv("data/bac_otutab.txt", sep = "\t", row.names = 1)
# Read sample metadata from the specified path
treat <- read.csv("data/metadata.csv", row.names = 1)

### Filter high-frequency data
# Define a frequency cutoff value to filter OTUs
cutoff = .6

# Filter the OTU table to retain OTUs with a frequency of at least 'cutoff' in samples
bac_filter <- data.frame(
  otutab[which(
    apply(otutab, 1, function(x) {
      length(which(x != 0)) / length(x)
    }) >= cutoff),
  ]
)

# Transpose the dataframe so that rows are samples and columns are OTUs
comm <- t(bac_filter)

## Normalization
# Normalize the transposed bacterial data
divindex <- scale(comm)

divs2 <- sapply(1:ncol(divindex), function(j) {
  message("Now j=", j, " in ", ncol(divindex), ". ", date())
  
  # Check if the current column has fewer than 3 unique values
  if (length(unique(divindex[, j])) < 3) {
    result <- rep(NA, 38)
  } else {
    div <- data.frame(divtest = divindex[, j], treat)
    
    # Fit the mixed-effects model
    fm1 <- lmer(divtest ~ PE + (1 | Site), data = div)
    
    # Perform ANOVA on the mixed-effects model
    presult <- car::Anova(fm1, type = "II")
    
    # Extract coefficients
    coefs <- coef(summary(fm1))[ , "Estimate"]  ## four coefficients
    names(coefs) <- paste0(names(coefs), ".mean")
    
    # Extract standard errors
    SEvalues <- coef(summary(fm1))[ , "Std. Error"] ## standard errors
    names(SEvalues) <- paste0(names(SEvalues), ".se")
    
    # Extract t-values
    tvalues <- coef(summary(fm1))[ , "t value"] ## t values
    names(tvalues) <- paste0(names(tvalues), ".t")
    
    # Extract Chi-square and P-values
    chisqP <- c(presult[,1], presult[,3])
    names(chisqP) <- c(paste0(row.names(presult), ".chisq"), paste0(row.names(presult), ".P"))
    
    # Combine all results
    result <- c(coefs, tvalues, SEvalues, chisqP)
  }
  result
})



## Load necessary R packages
library(lme4)
library(car)

### Read bacterial data
# Read the OTU table from the specified path
otutab <- read.csv("data/fun_otutab.txt", sep = "\t", row.names = 1)
# Read sample metadata from the specified path
treat <- read.csv("data/metadata.csv", row.names = 1)

### Filter high-frequency data
# Define a frequency cutoff value to filter OTUs
cutoff = .2

# Filter the OTU table to retain OTUs with a frequency of at least 'cutoff' in samples
fun_filter <- data.frame(
  otutab[which(
    apply(otutab, 1, function(x) {
      length(which(x != 0)) / length(x)
    }) >= cutoff),
  ]
)

# Transpose the dataframe so that rows are samples and columns are OTUs
comm <- t(fun_filter)

## Normalization
# Normalize the transposed bacterial data
divindex <- scale(comm)
divs2 <- sapply(1:ncol(divindex), function(j) {
  message("Now j=", j, " in ", ncol(divindex), ". ", date())
  
  # Check if the current column has fewer than 3 unique values
  if (length(unique(divindex[, j])) < 3) {
    result <- rep(NA, 38)
  } else {
    div <- data.frame(divtest = divindex[, j], treat)
    
    # Fit the mixed-effects model
    fm1 <- lmer(divtest ~ PE + (1 | Site), data = div)
    
    # Perform ANOVA on the mixed-effects model
    presult <- car::Anova(fm1, type = "II")
    
    # Extract coefficients
    coefs <- coef(summary(fm1))[ , "Estimate"]  ## four coefficients
    names(coefs) <- paste0(names(coefs), ".mean")
    
    # Extract standard errors
    SEvalues <- coef(summary(fm1))[ , "Std. Error"] ## standard errors
    names(SEvalues) <- paste0(names(SEvalues), ".se")
    
    # Extract t-values
    tvalues <- coef(summary(fm1))[ , "t value"] ## t values
    names(tvalues) <- paste0(names(tvalues), ".t")
    
    # Extract Chi-square and P-values
    chisqP <- c(presult[,1], presult[,3])
    names(chisqP) <- c(paste0(row.names(presult), ".chisq"), paste0(row.names(presult), ".P"))
    
    # Combine all results
    result <- c(coefs, tvalues, SEvalues, chisqP)
  }
  result
})
