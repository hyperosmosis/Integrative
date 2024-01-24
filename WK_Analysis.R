setwd("~/Documents/Wu/Weighted_Kernel/")

library(tidyverse)
library(ggpubr)
library(xtable)

n <- c(50, 100, 200, 300, 400)

readData <- function(directory, pathway){
  l <- list()
  n <- c(50, 100, 200, 300, 400)
  for (i in 1:length(n)){
    filename <- paste(directory, "WK_n", n[i], pathway, "_Bin.csv", sep = "")
    l[[i]] <- read.csv(filename)
  }
  return(l)
}

# Type I Error
n50 <- read.csv("Type1_Real/WK_n50map00564.csv")
n100 <- read.csv("Type1_Real/WK_n100map00564.csv")
n200 <- read.csv("Type1_Real/WK_n200map00564.csv")
n300 <- read.csv("Type1_Real/WK_n300map00564.csv")
n400 <- read.csv("Type1_Real/WK_n400map00564.csv")
n500 <- read.csv("Type1_Real/WK_n500map00564.csv")

w_map00564 <- list(n50, n100, n200, n300, n400, n500)

n50 <- read.csv("Type1_Real/WK_n50map02010.csv")
n100 <- read.csv("Type1_Real/WK_n100map02010.csv")
n200 <- read.csv("Type1_Real/WK_n200map02010.csv")
n300 <- read.csv("Type1_Real/WK_n300map02010.csv")
n400 <- read.csv("Type1_Real/WK_n400map02010.csv")
n500 <- read.csv("Type1_Real/WK_n500map02010.csv")

w_map02010 <- list(n50, n100, n200, n300, n400, n500)

names <- NULL
for(i in 1:5){
  names[i] <- paste("N: ", n[i])
}

names(w_map00564) <- names
names(w_map02010) <- names

makeQQPlot <- function(list, filename){
#  jpeg(filename, width = 8, height = 10, units = "in", res = 150)
  par(mfrow = c(3, 2))
  lapply(names(list), function(x) {qqplot(-log10(list[[x]][,2]), -log10(runif(5000)),
                                                main = x, xlab = "Observed", ylab = "Expected",xlim = c(0, 3),
                                                ylim = c(0, 3.5), pch = 20)
    abline(a = 0, b = 1, lwd = 3, lty = 2, col = "red")})
#  dev.off()
}

makeHistogram <- function(list, filename){
  pdf(filename, width = 8, height = 10)
  par(mfrow = c(3, 2))
  lapply(names(list), function(x) {hist(list[[x]][,2],
                                        main = x, xlab = "Observed")})
  dev.off()
}
jpeg("Figures/QQPlot_map00564.jpg", width = 8, height = 10, units = "in", res = 150)
par(mfrow = c(3, 2))
lapply(names(w_map02010), function(x) {qqplot(-log10(w_map02010[[x]][,2]), -log10(runif(5000)),
                                              main = x, xlab = "Observed", ylab = "Expected",xlim = c(0, 3),
                                              ylim = c(0, 3.5), pch = 20)
  abline(a = 0, b = 1, lwd = 3, lty = 2, col = "red")})
dev.off()

jpeg("Figures/QQPlot_map02010.jpg", width = 8, height = 10, units = "in", res = 150)
par(mfrow = c(3, 2))
lapply(names(w_map02010), function(x) {qqplot(-log10(w_map02010[[x]][,2]), -log10(runif(5000)),
                                              main = x, xlab = "Observed", ylab = "Expected",xlim = c(0, 3),
                                              ylim = c(0, 3.5), pch = 20)
  abline(a = 0, b = 1, lwd = 3, lty = 2, col = "red")})
dev.off()

findType1 <- function(data, level = 0.05){
  result <- NULL
  n <- dim(data[[1]])[1]
  for (i in 1:length(data)){
    if (level == 0.05){
      result[i] <- sum(data[[i]][2]< 0.05)/n
    }
    else{
      result[i] <- sum(data[[i]][3] < 0.01)/n
    }
  }
  return(result)
}

type1_map00564 <- data.frame(n = c(50, 100, 200, 300, 400, 500),
                             a05 = findType1(w_map00564, 0.05),
                             a01 = findType1(w_map00564, 0.01))

type1_map02010 <- data.frame(n = c(50, 100, 200, 300, 400, 500),
                             a05 = findType1(w_map02020, 0.05),
                             a01 = findType1(w_map02020, 0.01))

print(xtable(type1_map02010, digits = 4), include.rownames = F)


## Type I Error with Cauchy Combination Test

cauchy.map00564 <- readData("Type1_Cauchy/", "map00564")
cauchy.map02010 <- readData("Type1_Cauchy/", "map02010")
names(cauchy.map00564) <- names[1:5]
names(cauchy.map02010) <- names[1:5]

makeQQPlot(cauchy.map00564, "Figures/QQPlot_Cauchy_map00564.jpg")
makeQQPlot(cauchy.map02010, "Figures/QQPlot_Cauchy_map02010.jpg")
makeHistogram(cauchy.map00564, "Figures/Hist_Cauchy_map00564.pdf")
makeHistogram(cauchy.map02010, "Figures/Hist_Cauchy_map02010.pdf")


type1_map00564 <- data.frame(n = c(50, 100, 200, 300, 400),
                             a05 = findType1(cauchy.map00564, 0.05),
                             a01 = findType1(cauchy.map00564, 0.01))

type1_map02010 <- data.frame(n = c(50, 100, 200, 300, 400),
                             a05 = findType1(cauchy.map02010, 0.05),
                             a01 = findType1(cauchy.map02010, 0.01))

print(xtable(type1_map02010, digits = 4), include.rownames = F)

## Longitudinal

long <- readData("~/Documents/Wu/Weighted_Kernel/Longitudinal/Type1/", "")
names(long) <- names

