### Boutros R Training ############################################################################

# Load Library and Input 1&2 Files
rm(list = ls());
library(reshape);
setwd('D:/Documents/Life Now/Bio Stuff/Boutros/Tutorial/R/R_practice');
input1 <- read.delim('input1.txt');
input2 <- read.delim('input2.txt');

# Reorder and Combine Input 1&2
input1 <- input1[order(input1$GeneID), ];
input2 <- input2[order(input2$GeneID), ];
input = merge(input1,input2,by = "GeneID");

# Method 2
input <- cbind(input1, input2[ ,2:ncol(input2)]); 


### T test ########################################################################################

# Create storage variable "tvalues" and remove "GeneID" Column for Input1 & Input2
tvalues <- rep(x = 0, 
               times = nrow(input));
input1 <- input1[,-1];
input2 <- input2[,-1];

# Calculate T-value for each row and store into "tvalues" variable
for (i in 1:nrow(input)) {
  tvalue <- t.test( 
    x = as.numeric(input1[i,]),
    y = as.numeric(input2[i,])
    );
  tvalues[i] <- as.numeric(tvalue$statistic);
  } 

# Create Histogram and Plot
tvalues = data.frame(input$GeneID,tvalues);
hist( 
  x = tvalues$tvalues, 
  main = "Histogram for T-Values", 
  breaks = 20
  );
plot(
  tvalues$tvalues, 
  log = "y",
  type = 'h',
  las = 1,
  xlab = "Genes",
  ylab = "T Values"
  );

# Since the critical value for 10 degrees of freedom t-test is 1.812,
# and all the values here is below e^0.2 (1.22), we can conclude that
# there isn't a different in population mean and that patients 1-3 and 4-12
# are from the same population.

### Wilcoxon Test and Fold Changes ################################################################
group1dat <- melt(input1);
group2dat <- melt(input2);
wilcox.test(
  x = group1dat$value, 
  y = group2dat$value
  );

# The P value is around 0.5, so the neutral hypothesis is not rejected
# The two samples are from the same population


### Permutation Test ##############################################################################

# Below are some storage variables
p_vals <- rep(
  x = 0,
  times = nrow(input));
observed_medians <- p_vals;
expected_medians <- p_vals;

# Perform the Test in for loop
for (gene in 1:nrow(input)) {
  
  # Specify # of Iterations and random sample size.
  iterations <- 1000;
  sample_size <- as.numeric(ncol(input1));
  
  # Create counters for how many larger/smaller than observed medians there are.
  larger_count <- 0;
  
  # Calculate expected median.
  med_expected <- median(as.numeric(input1));
  
  # Randomly sample from all the input data and calculate observed median.
  # Then, check if expected median is larger/smaller than observed median.
  for (rep in 1:iterations) {
    med_random <- median(as.numeric(sample(
      x = input[gene,2:ncol(input)],
      size = sample_size
      )));
    if (med_random >= med_expected) { larger_count = larger_count + 1; }
  }
  
  # Caculate p-value and store the results into storage variables.
  p_val <- larger_count / iterations;
  p_vals[gene] <- p_val;
  observed_medians[gene] <- med_random;
  expected_medians[gene] <- med_expected;
}

# Calculate adjusted p-value and integrate into output file.
p_adjst <- p.adjust(p_vals);
output1 <- table(input$GeneID,observed_medians,expected_medians,p_vals,p_adjst);


### Permutation Test with Fold Changes ############################################################

# Below are some storage variables
p_vals <- rep(x = 0, 
              times = nrow(input));
observed_medians <- p_vals;
expected_medians <- p_vals;

# Perform the Test in for loop
for (gene in 1:nrow(input)) {
  
  # This is a counter for how many 2-fold up&down there are
  fold_changes <- 0;
  
  # Specify # of Iterations and random sample size.
  iterations <- 1000;
  sample_size <- as.numeric(ncol(input1));
  med_observed <- median(as.numeric(input1));
  for (rep in 1:iterations) {
    medrandom <- median(as.numeric(sample(
      x = input[gene,2:ncol(input)],
      size = sample_size
      )));
    fold_change <- log2(max(medrandom, med_observed) / min(medrandom, med_observed));
    # Calculate fold change for this specific iteration
    if (fold_change >= 2) {fold_changes <- fold_changes + 1;}
    # Counter +1 if fold changes reaches 2-folds
  }
  
  # Store the results into temporary variables
  observed_medians[gene] <- med_observed;
  expected_medians[gene] <- medrandom;
  p_val <- abs(fold_changes) / iterations;
  # Number of >= 2-fold divided by total iteration
  p_vals[gene] <- p_val;
}
# Adjust P-values
p_adjst <- (p.adjust(p_vals));

# Create output file
output2 <- table(input$GeneID,observed_medians,expected_medians,p_vals,p_adjst);

