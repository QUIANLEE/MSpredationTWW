# Sample dataset(DAY28-CHX ecoli number)
# data <- data.frame(
#   group = rep(c("CS", "OS"), each = 5),
#   value = c(1000, 6000, 9000, 1000, 1000, 3000, 8000, 9000, 10000, 1000)  # Replace with real data
# )

# Sample dataset(-CHX: day0-day28 ecoli number)
rm(list = ls())
data <- data.frame(
  group = rep(c("CS", "OS"), each = 5),
  value = c(7279000, 5874000, 5471000, 7239000, 5079000, 4957000, 5752000, 6471000, 6510000, 4039000)  # Replace with real data
)

# Sample dataset(DAY28+CHX ecoli number)
# data <- data.frame(
#   group = rep(c("CS", "OS"), each = 5),
#   value = c(151000, 83000, 161000, 85000, 43000, 104000, 176000, 236000, 244000, 296000)  # Replace with real data
# )

# Sample dataset(+CHX: day0-day28 ecoli number)
rm(list = ls())
data <- data.frame(
  group = rep(c("CS", "OS"), each = 5),
  value = c(7009000, 7837000, 6439000, 8315000, 9677000, 6176000, 10864000, 6164000, 5876000, 6344000)  # Replace with real data
)


# Check normality
shapiro.test(data$value)  # If p < 0.05, data is non-normal distribution, then perform Wilcoxon test (Mann-Whitney U test)


# Perform an independent t-test----normal distribution
t_test_result <- t.test(value ~ group, data = data, var.equal = TRUE)  # Set var.equal = TRUE for equal variance
print(t_test_result)

# Perform Wilcoxon test (Mann-Whitney U test)----non-normal distribution
wilcox_test_result <- wilcox.test(value ~ group, data = data)
print(wilcox_test_result)