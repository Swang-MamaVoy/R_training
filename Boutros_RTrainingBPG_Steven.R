
### Boutros Lab Package Implementation ############################################################

# Load Library and Input 1&2 Files
rm(list = ls());
library(reshape);
library(BoutrosLab.dist.overload);
library(BoutrosLab.plotting.general);
library(BoutrosLab.statistics.general);
.libPaths(c("~/BL", .libPaths()))
setwd("D:/Documents/Life Now/Bio Stuff/Boutros/Tutorial/R/R_practice");
input1 <- read.delim("input1.txt");
input2 <- read.delim("input2.txt");

# Reorder and Combine Input 1&2
input1 <- input1[order(input1$GeneID), ];
input2 <- input2[order(input2$GeneID), ];
input = merge(input1,input2,by = "GeneID");
input1 <- input1[,-1];
input2 <- input2[,-1];

### T test ########################################################################################

# Create storage variable "pvalues"
pvalues <- rep(x = 0, 
               times = nrow(input));

for (i in 1:nrow(input)) {
  group1 <- t(input1[i,]);
  group2 <- t(input2[i,]);
  group1 <- cbind(group1,rep("group1", nrow(group1)));
  group2 <- cbind(group2,rep("group2", nrow(group2)));
  t_data <- rbind(group1, group2);
  t_result <- ttest.analysis(
    values = as.numeric(t_data[,1]),
    groups = as.factor(t_data[,2])
    );
  
  # Store T-value into matrix
  pvalues[i] <- t_result$pvalue;
  } 
hist <- create.histogram(
  filename = "Histogram_Ttest.tiff",
  x = pvalues,
  xlab.label = "P Values",
  ylab.label = "Count",
  main = "Histogram of T-test p-values",
  main.cex = 1,
  xlab.cex = 1,
  ylab.cex = 1,
  width = 4,
  height= 4,
  top.padding = 0.7,
  description = "Histogram of T Test, R Training Course"
  );


### Wilcoxon Test and Fold Changes ################################################################


pvals_utest <- rep(x = 0, 
                   times = nrow(input));
fold_changes <- pvals_utest;
for (i in 1:nrow(input)) {
  group1dat <- melt(
    data = input1[i,], 
    id = NULL
    );
  group2dat <- melt(
    data = input2[i,], 
    id = NULL
    );
  group1dat$variable <- rep(
    x = "group1",
    times = nrow(group1dat)
    ); 
  group2dat$variable <- rep(
    x = "group2",
    times = nrow(group2dat));
  dat <- rbind(group1dat, group2dat);
  group1bol <- dat$variable == "group1";
  group2bol <- dat$variable == "group2";
  pval_utest <- get.utest.p(
    x = dat$value,
    group1 = group1bol,
    group2 = group2bol);
  pvals_utest[i] <- pval_utest;
  fold_change <- get.foldchange(
    x = dat$value,
    group1 = group1bol,
    group2 = group2bol,
    logged = FALSE
    );
  fold_changes[i] <- fold_change;
  }

hist <- create.histogram(
  filename = "Histogram_Wilcoxon.tiff",
  x = pvals_utest,
  xlab.label = "P Values",
  ylab.label = "Count",
  main = "Histogram of Wilcoxon Test p-values",
  main.cex = 1,
  xlab.cex = 1,
  ylab.cex = 1,
  top.padding = 0.7,
  height = 4,
  width = 4,
  description = "Histogram of T Test, R Training Course"
  );

scatterplot_dataframe <- data.frame(gene = input$GeneID, fold_changes);
scat <- create.scatterplot(
  filename = "ScatterPlot_FoldChanges.tiff",
  main = "Scatter Plot of Fold Changes",
  data = scatterplot_dataframe,
  formula = fold_changes ~ gene,
  xlab.label = "Gene ID",
  main.cex = 1,
  xlab.cex = 1,
  ylab.cex = 1,
  xaxis.cex = 0.25, 
  xaxis.rot = 90,
  ylab.label = "Fold Changes",
  cex = 0.1,
  top.padding = 0.7,
  height = 4,
  width = 10
  );

### Permutation Test ##############################################################################

# Some storage variables
p_vals <- rep(x = 0, 
              times = nrow(input));
observed_medians <- p_vals;
expected_medians <- p_vals;
sample_size <- as.numeric(ncol(input1));
for (gene in 1:nrow(input)) {
  larger_count <- 0;
  smaller_count <- 0;
  iterations <- 1000;
  med1 <- median(
    x = as.numeric(input1[gene,])
    );
  for (rep in 1:iterations) {
    medrandom <- median(as.numeric(sample(
      x = input[gene,2:ncol(input)],
      size = sample_size
      )));
    if (medrandom >= med1) {larger_count = larger_count + 1;}
    }
  observed_medians[gene] <- med1;
  expected_medians[gene] <- medrandom;
  
  # Store the results into temporary variables
  p_val <- larger_count / iterations;
  p_vals[gene] <- p_val;
  }
p_adjst <- p.adjust(p_vals);
hist <- create.histogram(
  filename = "Histogram_Bootstrap.tiff",
  x = p_vals,
  xlab.label = "P Values",
  ylab.label = "Count",
  main = "Histogram of Bootstrap Test P-values",
  main.cex = 1,
  xlab.cex = 1,
  ylab.cex = 1,
  top.padding = 0.7,
  height = 4,
  width = 4,
  description = "Histogram of T Test, R Training Course"
  );

#### Permutation Test with Fold Changes ###########################################################

p_vals <- rep(0, nrow(input));
fold_changes <- p_vals;
sample_size <- as.numeric(ncol(input1));

# This for loop executes the Permutation Test for each row in "input"
# Calculate the expected median using data from "input1"
# Repeat 1000 random samplings for each row
# Calculate fold changes between observed and expected medians
# Add 1 to counter if fold change exceeds the threshold of 2-folds up and down
# Calculate p-value by pval = counter / iterations

for (gene in 1:nrow(input)) {
  iterations <- 1000;
  group1dat <- melt(
    data = input1[gene,], 
    id = NULL
  );
  foldchanges = rep(0,iterations);
  group1dat$variable <- rep("group1",nrow(group1dat)); 
  for (iter in 1:iterations) {
    group2dat <- melt(
      data = as.numeric(sample(input[gene,2:ncol(input)], sample_size)),
      id = NULL
    );
    group2dat$variable <- rep("group2",nrow(group2dat));
    dat <- rbind(group1dat, group2dat);
    group1bol <- dat$variable == "group1";
    group2bol <- dat$variable == "group2";
    fold_change <- get.foldchange(
      x = dat$value,
      group1 = group1bol,
      group2 = group2bol,
      logged = FALSE
    );
    foldchanges[iter] <- fold_change;
  }
  fold_changes[gene] <- mean(foldchanges);
  }
scatterplot_dataframe <- data.frame(gene = input$GeneID, fold_changes);
scat <- create.scatterplot(
  filename = "ScatterPlot_Permutation_FoldChanges.tiff",
  main = "Scatter Plot of Permutation Fold Changes",
  data = scatterplot_dataframe,
  formula = fold_changes ~ gene,
  xlab.label = "Gene ID",
  main.cex = 1,
  xlab.cex = 1,
  ylab.cex = 1,
  xaxis.cex = 0.25, 
  xaxis.rot = 90,
  ylab.label = "Fold Changes",
  cex = 0.1,
  top.padding = 0.7,
  height = 4,
  width = 10
)
