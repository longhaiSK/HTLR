library(corrplot)

diabetes <- read.csv("data-raw/diabetes.csv", header = TRUE)

summary(diabetes)
str(diabetes)

# removing those observation rows with 0 in any of the variables
for (i in 2:6) {
  diabetes <- diabetes[-which(diabetes[, i] == 0), ]
}

# scale the covariates for easier comparison of coefficient posteriors
for (i in 1:8) {
  diabetes[i] <- scale(diabetes[i])
}

# modify the data column names slightly for easier typing
names(diabetes)[7] <- "PedigreeFunction"
names(diabetes) <- tolower(names(diabetes))

str(diabetes)

corrplot(cor(diabetes[, c(9,1:8)]))

# preparing the dataset
diabetes392 <- list("X" = model.matrix(outcome ~ . - 1, data = diabetes),
                    "y" = factor(diabetes$outcome))

usethis::use_data(diabetes392)
