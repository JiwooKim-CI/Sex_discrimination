# ------------------------------------------------------------
# Sex discrimination evidence analysis
# ------------------------------------------------------------
# gender: 1 = male, 2 = female
#
# Dataset 1 job coding:
#   1 = non-production, promoted
#   2 = production, promoted
#   3 = non-production, not promoted
#   4 = production, not promoted
#
# Dataset 2 job coding:
#   1 = high school, regular
#   2 = college+, regular
#   3 = high school, non-regular
#   4 = college+, non-regular
# ------------------------------------------------------------

library(dplyr)
library(boot)
library(ggplot2)
library(tikzDevice)

# ------------------------------------------------------------
# 1. Data construction
# ------------------------------------------------------------

make_dataset <- function(job_counts, gender_counts) {
  stopifnot(length(job_counts) == length(gender_counts))
  
  df <- data.frame(
    job = unlist(Map(rep, seq_along(job_counts), job_counts)),
    gender = unlist(Map(rep, gender_counts$gender, gender_counts$n))
  )
  
  df$gender <- factor(df$gender, levels = c(1, 2), labels = c("male", "female"))
  df
}

# Dataset 1
job_1 <- data.frame(
  job = c(
    rep(1, 205), rep(2, 182), rep(3, 7), rep(4, 20),
    rep(1, 7),   rep(2, 19),  rep(4, 151)
  ),
  gender = c(rep(1, 414), rep(2, 177))
)
job_1$gender <- factor(job_1$gender, levels = c(1, 2), labels = c("male", "female"))

# Dataset 2
job_2 <- data.frame(
  job = c(
    rep(1, 5),   rep(2, 880), rep(3, 59),  rep(4, 195),
    rep(1, 55),  rep(2, 527), rep(3, 303), rep(4, 527)
  ),
  gender = c(rep(1, 1139), rep(2, 1412))
)
job_2$gender <- factor(job_2$gender, levels = c(1, 2), labels = c("male", "female"))

# ------------------------------------------------------------
# 2. Variable coding
# ------------------------------------------------------------

add_binary_variables <- function(df) {
  df %>%
    mutate(
      Z = gender,
      X = ifelse(job %in% c(2, 4), 1, 0),  # production / higher education indicator
      Y = ifelse(job %in% c(1, 2), 1, 0)   # promoted / regular indicator
    )
}

job_1 <- add_binary_variables(job_1)
job_2 <- add_binary_variables(job_2)

# ------------------------------------------------------------
# 3. Point estimate
# ------------------------------------------------------------

compute_point_estimate <- function(df) {
  n_male <- sum(df$gender == "male")
  n_female <- sum(df$gender == "female")
  
  max(
    max(sum(df$job == 1 & df$gender == "male") / n_male,
        sum(df$job == 1 & df$gender == "female") / n_female) +
      max(sum(df$job == 3 & df$gender == "male") / n_male,
          sum(df$job == 3 & df$gender == "female") / n_female),
    
    max(sum(df$job == 2 & df$gender == "male") / n_male,
        sum(df$job == 2 & df$gender == "female") / n_female) +
      max(sum(df$job == 4 & df$gender == "male") / n_male,
          sum(df$job == 4 & df$gender == "female") / n_female)
  )
}

point_est_1 <- compute_point_estimate(job_1)
point_est_2 <- compute_point_estimate(job_2)

print(point_est_1)
print(point_est_2)

# ------------------------------------------------------------
# 4. Q_dy construction and hypothesis tests
# ------------------------------------------------------------

make_Qdy <- function(df, d, y) {
  ifelse(
    df$Z == "male",
    as.integer(df$X == d & df$Y == y),
    1 - as.integer(df$X == d & df$Y == 1 - y)
  )
}

run_qdy_tests <- function(df) {
  results <- data.frame()
  
  for (d in 0:1) {
    for (y in 0:1) {
      tmp <- df
      tmp$Qdy <- make_Qdy(tmp, d, y)
      
      tbl <- table(tmp$Z, tmp$Qdy)
      
      fisher_p <- NA_real_
      if (all(dim(tbl) == c(2, 2))) {
        fisher_p <- fisher.test(tbl, alternative = "greater")$p.value
      }
      
      n1 <- sum(tbl[1, ])
      n0 <- sum(tbl[2, ])
      p1 <- tbl[1, 2] / n1
      p0 <- tbl[2, 2] / n0
      se <- sqrt(p1 * (1 - p1) / n1 + p0 * (1 - p0) / n0)
      wald_z <- (p1 - p0) / se
      wald_p <- 1 - pnorm(wald_z)
      
      results <- rbind(
        results,
        data.frame(
          d = d,
          y = y,
          Z1_Q1 = tbl[1, 2],
          Z1_n = n1,
          Z0_Q1 = tbl[2, 2],
          Z0_n = n0,
          p_fisher = fisher_p,
          p_wald = wald_p
        )
      )
    }
  }
  
  results
}

test_results_1 <- run_qdy_tests(job_1)
test_results_2 <- run_qdy_tests(job_2)

print(test_results_1)
print(test_results_2)

# ------------------------------------------------------------
# 5. Bootstrap statistic
# ------------------------------------------------------------

boot_statistic <- function(data, indices) {
  d <- data[indices, ]
  n_male <- sum(d$gender == "male")
  n_female <- sum(d$gender == "female")
  
  max(
    max(sum(d$job == 1 & d$gender == "male") / n_male,
        sum(d$job == 1 & d$gender == "female") / n_female) +
      max(sum(d$job == 3 & d$gender == "male") / n_male,
          sum(d$job == 3 & d$gender == "female") / n_female),
    
    max(sum(d$job == 2 & d$gender == "male") / n_male,
        sum(d$job == 2 & d$gender == "female") / n_female) +
      max(sum(d$job == 4 & d$gender == "male") / n_male,
          sum(d$job == 4 & d$gender == "female") / n_female)
  )
}

run_bootstrap <- function(df, R = 10000) {
  boot(
    data = df,
    statistic = boot_statistic,
    R = R,
    strata = df$gender
  )
}

boot_result_1 <- run_bootstrap(job_1, R = 10000)
boot_result_2 <- run_bootstrap(job_2, R = 10000)

print(summary(boot_result_1))
print(summary(boot_result_2))

print(boot.ci(boot_result_1))
print(boot.ci(boot_result_2))

# ------------------------------------------------------------
# 6. Plot bootstrap distribution
# ------------------------------------------------------------

plot_boot_histogram <- function(boot_result,
                                output_tex = NULL,
                                dashed_lines = NULL,
                                red_line = 1,
                                xlab = "Result of Statistical Evidence for Sex Discrimination",
                                ylab = "Frequency") {
  t_df <- data.frame(stat = boot_result$t[, 1])
  
  g <- ggplot(t_df, aes(x = stat)) +
    geom_histogram(color = "black", fill = "gray", bins = 30) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold")
    ) +
    xlab(xlab) +
    ylab(ylab)
  
  if (!is.null(dashed_lines)) {
    g <- g + geom_vline(
      xintercept = dashed_lines,
      color = "black",
      linetype = "dashed",
      linewidth = 1
    )
  }
  
  if (!is.null(red_line)) {
    g <- g + geom_vline(
      xintercept = red_line,
      color = "red",
      linewidth = 1
    )
  }
  
  if (!is.null(output_tex)) {
    tikz(output_tex, standAlone = TRUE, width = 5, height = 4)
    print(g)
    dev.off()
  }
  
  g
}

# Example plot for dataset 2
g2 <- plot_boot_histogram(
  boot_result = boot_result_2,
  output_tex = "figure_plot_2.tex",
  dashed_lines = c(1.111, 1.181),
  red_line = 1
)

print(g2)

# Base boot plot if needed
plot(boot_result_2)
plot.boot(boot_result_2)