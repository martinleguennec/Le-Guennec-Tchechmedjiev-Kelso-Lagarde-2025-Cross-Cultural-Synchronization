## ──── GLOBAL PLOT THEME ──────────────────────────────────────────────────────────────────────────

theme_set(
  theme_light() +
    theme(
      text = element_text(family = "serif"),
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.title = element_text(size = 22),
      axis.text = element_text(size = 18),
      legend.key.size = unit(1, "cm"),
      legend.text = element_text(size = 20),
      legend.title = element_blank(),
      panel.background = element_rect(fill = "white"),
      panel.border = element_rect(color = "black"),
      panel.grid = element_blank(),
      panel.spacing = unit(1.5, "cm"),
      plot.subtitle = element_text(size = 20, hjust = 0.5),
      plot.title = element_text(face = "bold", size = 22, hjust = 0.5),
      plot.tag.position = "top",
      plot.tag = element_text(
        face = "bold",
        margin = margin(0, 0, 10, 0)
      ),
      strip.background = element_rect(color = "black", fill = "white"),
      strip.text = element_text(color = "black", size = 22)
    )
)

# Color palette
FR <- "#2596be"
FR_dark <- "#046281"
IN <- "#d55f32"
IN_dark <- "#892702"


## ──── PREDICTION UTILITIES ───────────────────────────────────────────────────────────────────────

#' Calculate standard error of regression predictions
#' 
#' This function calculates the standard error of predictions from a linear regression model.
#' It allows to plot the regression line with confidence intervals. The 95% confidence interval is 
#' calculated as the predicted value plus or minus 1.96 times the standard error of the prediction 
#' (based on the assumption that the residual are normally distributed).
#' 
#' #' @param beta0 Intercept of the regression model.
#' #' @param se_beta0 Standard error of the intercept.
#' #' @param beta1 Slope of the regression model.
#' #' @param se_beta1 Standard error of the slope.
#' #' @param x_range Range of x values for which to calculate predictions.
#' #' @param n_est Number of points to estimate within the x_range.
#' #' @return A data frame containing the predicted values, standard errors, and confidence intervals.
calculate_se_regression <- function(
    beta0, se_beta0, 
    beta1, se_beta1,
    x_range, n_est = 100) {
  
  x_values <- seq(x_range[1], x_range[2], length.out = n_est)
  
  pred_data <- data.frame(
    x = x_values,
    pred = beta0 + beta1 * x_values
  ) %>%
    mutate(
      se_pred = sqrt(se_beta0^2 + (x_values^2) * se_beta1^2),
      lower = pred - 1.96 * se_pred,
      upper = pred + 1.96 * se_pred
    )
  
  return(pred_data)
}


## ──── LINEAR MIXED MODEL UTILITIES ───────────────────────────────────────────────────────────────

#' Calculate Intraclass Correlation Coefficient (ICC)
#' 
#' This function calculates the Intraclass Correlation Coefficient (ICC) from a linear mixed model, 
#' only when the random effect is a random intercept alone. It is calculated as the ratio of the
#' variance of the random intercept to the total variance (random effects variance + residual variance).
#' 
#' #' @param lmm A fitted linear mixed model object (from `lme4` or similar).
#' #' @return The ICC value.
calculate_icc <- function(lmm) {
  random_effects <- as.data.frame(VarCorr(lmm))
  n_effects <- length(random_effects$vcov)
  icc <- random_effects$vcov[1] / sum(random_effects$vcov)
  return(icc)
}


#' Diagnostic plots for linear models
#' 
#' This function generates diagnostic plots for a linear model to check the assumptions of 
#' linearity, normality and homoscedasticity of residuals. If the model is a linear mixed model, the 
#' distribution of random effects are also plotted.
#' 
#' @param lm A fitted linear model object.
#' @param n_breaks A vector of two integers specifying the number of breaks for the histogram of residuals
lm_assumptions_check <- function(lm, n_breaks = c(30, 10)) {
  
  if (class(lm) == "lmerModLmerTest") {
    # Extract random effects as a data frame and identify them
    ranefs <- as.data.frame(ranef(lm))
    ranefs_names <- levels(ranefs$term)
    nrows <- length(ranefs_names) + 1
  } else {
    # If the model is not a linear mixed model, we only have one row for the residuals
    nrows <- 1
  }
  
  # Set up the plotting area
  par(mfrow = c(nrows, 3), cex.main = 2, cex.lab = 2, cex.axis = 2)
  
  # A. Residuals vs fitted values plot
  plot(fitted(lm), residuals(lm),
       xlab = "Fitted values", ylab = "Residuals",
       main = "A. Residuals vs. fitted values"
  )
  abline(h = 0, col = "grey40", lty = 2)
  
  
  # B. Normal Q-Q plot for residuals
  car::qqPlot(
    residuals(lm), 
    pch = 19, col.lines = "grey40", grid = FALSE,
    xlab = "Theoretical Quantiles", ylab = "Residuals", 
    main = "B. Normal Q-Q Plot for residuals"
  )
  
  # C. Histogram of residuals
  hist(residuals(lm),
       breaks = n_breaks[1],
       xlab = "Residuals",
       main = "C. Histogram of the residuals"
  )
  abline(v = 0, lwd = 3, lty = "dashed", col = "red")
  
  
  if (class(lm) == "lmerModLmerTest") {
    # Create a "group" column based on rownames, subject's ID is contained in it
    # If the intercept varies by frequency, the grp column is ID:freq
    # Identify subject's group based on the last 2 characters of their ID: FR = French; IN = Indian
    # Identifying groups is necessary to plot the random effect afterward
    if (grepl(":", ranefs$grp[1], fixed = TRUE)) {
      ranefs <- ranefs %>%
        separate(grp, into = c("subject", "frequency"), sep = ":") %>%
        mutate(group = ifelse(str_sub(subject, -2) == "Fr", "French", "Indian"))
    } else if ("grp" %in% colnames(ranefs)) {
      ranefs$group <- ifelse(str_sub(ranefs$grp, -2) == "Fr", "French", "Indian")
    } else {
      ranefs$group <- ifelse(
        str_sub(rownames(ranefs), -2) == "Fr", "French", "Indian"
      )
    }
    
    # Compute 95% confidence intervals for random effects 
    # Computed as 1.96 × SD, based on the assumption that the random effects are normally distributed
    ranefs <- ranefs %>%
      mutate(
        CI_low = condval - 1.96 * condsd,
        CI_high = condval + 1.96 * condsd
      ) %>%
      arrange(term, group, condval) %>%
      mutate(ID = row_number())
    
    
    letter_num <- 4
    
    for (the_ranef in ranefs_names) {
      ranef_data <- ranefs[ranefs$term == the_ranef, ]
      color_map <- c("French" = "#2596be", "Indian" = "#d55f32")
      colors <- color_map[ranef_data[, "group"]]
      
      # Compute y-limits based on the range of CI_low, CI_high, and condval
      ylim_vals <- range(ranef_data[, "CI_low"], ranef_data[, "CI_high"], ranef_data[, "condval"], na.rm = TRUE)
      
      # Random effect with confidence interval by group
      plot(
        ranef_data[, "ID"], ranef_data[, "condval"],
        col = colors,
        pch = 16, ylim = ylim_vals,
        xlab = "", ylab = "Random effect", xaxt = "n",
        main = sprintf(
          "%s. Individual coefficients for '%s'",
          LETTERS[letter_num], as.character(the_ranef)
        )
      )
      segments(ranef_data[, "ID"], ranef_data[, "CI_low"],
               ranef_data[, "ID"], ranef_data[, "CI_high"],
               col = colors
      )
      abline(h = 0, lty = "dashed")
      legend("topleft",
             legend = c("French", "Indian"), col = c("#2596be", "#d55f32"),
             lty = 1, pch = 16
      )
      
      
      # Normal Q-Q plot for residuals for random effect
      car::qqPlot(
        ranef_data[, "condval"],
        pch = 19, col.lines = "grey40", grid = FALSE,
        xlab = "Theoretical Quantiles", ylab = "Random Effects", 
        main = sprintf(
          "%s. Normal Q-Q Plot for '%s'",
          LETTERS[letter_num+1], as.character(the_ranef)
        )
      )
      
      # Histogram for random effect
      hist(ranef_data[, "condval"],
           breaks = n_breaks[2],
           xlab = "Random effects",
           main = sprintf(
             "%s. Histogram of '%s' ",
             LETTERS[letter_num + 2], as.character(the_ranef)
           )
      )
      abline(v = 0, lwd = 3, lty = "dashed", col = "red")
      
      
      letter_num <- letter_num + 3
    }
    
  }
  
  par(mfrow = c(1, 1))
}


## ──── MODEL SUMMARY AND DATAFRAME FORMATTING ─────────────────────────────────────────────────────

create_df_ranef_intercept <- function(
    model,
    n_decimals = c(2, 2)) {
  format_variance <- paste0("%0.", n_decimals[1], "f")
  format_sd <- paste0("%0.", n_decimals[2], "f")
  
  # Get the random effects estimates
  # First row contains variance about random intercept
  # Second row contains residual variance
  as.data.frame(VarCorr(model)) %>%
    mutate(
      # Name of the component
      grp = c("Random Effects: Intercept", "Residual"),
      # Variance of the component
      vcov = sprintf(format_variance, vcov),
      # Standard deviation of the component
      sdcor = sprintf(format_sd, sdcor),
      # Intraclass correlation (only for random effect)
      ICC = c(sprintf("%0.2f", calculate_icc(model)), ""),
      # Marginal and condition R-squared, set to nothing for the moment
      mR2 = "",
      cR2 = ""
    ) %>%
    dplyr::select(-c(var1, var2)) %>%
    # Add a row for the model fit: it only displays marginal and conditional R-squared
    bind_rows(
      data.frame(
        grp = "Model Fit",
        vcov = "",
        sdcor = "",
        ICC = "",
        mR2 = sprintf("%0.2f", r.squaredGLMM(model)[1]),
        cR2 = sprintf("%0.2f", r.squaredGLMM(model)[2])
      )
    ) %>%
    setNames(c(
      "Component", "Variance", "SD", "ICC",
      "$R^2$ (marg.)", "$R^2$ (cond.)"
    ))
}


format_model_comparison <- function(
    comparison_results,
    ic_digits = 1,
    r2_digits = 2,
    error_digits = 2) {
  as.data.frame(comparison_results) %>%
    mutate(
      AICc = sprintf(
        paste0("%0.", ic_digits, "f (%s)"),
        AICc,
        ifelse(AICc_wt < 0.001, "<.001",
               ifelse(AICc_wt > 0.999, ">.999",
                      sprintf("%0.3f", AICc_wt)
               )
        )
      ),
      R2_conditional = sprintf(
        paste0("%0.", r2_digits, "f"),
        R2_conditional
      ),
      R2_marginal = sprintf(
        paste0("%0.", r2_digits, "f"),
        R2_marginal
      )
    ) %>%
    select(Name, AICc, R2_conditional, R2_marginal, ) %>%
    setNames(c(
      "Model", "AICc (weights)", "$R^2$ (cond.)", "$R^2$ (marg.)"
    ))
}


format_fixed_coefficients <- function(
    model,
    n_digits = 2) {
  # Get the fixed effects coefficients
  fixef <- as.data.frame((summary(model))$coefficients)
  
  fixef_formatted <- fixef %>%
    # Modify columns names so that they are easier to modify
    setNames(c("estimate", "SE", "DF", "t_value", "p_value")) %>%
    # Modify formatting of the values
    mutate(
      estimate = sprintf(paste0("%0.", n_digits, "f"), estimate),
      SE = sprintf(paste0("%0.", n_digits, "f"), SE),
      DF = sprintf("%0.1f", DF),
      t_value = sprintf("%0.2f", t_value),
      p_value = ifelse(
        p_value >= 0.001,
        sprintf("%0.3f", p_value),
        "<.001"
      )
    ) %>%
    # Format column names
    setNames(c("Estimate", "SE", "DF", "$t$-value", "$p$-value"))
  
  return(fixef_formatted)
}


format_df_emmeans <- function(
    model_emmeans,
    condition = "condition",
    grouping = "group",
    n_digits = 2) {
  
  model_emmeans %>%
    mutate(
      df = sprintf("%0.1f", df),
      across(
        where(is.numeric),
        ~ sprintf(paste0("%0.", n_digits, "f"), .)
      ),
      CI = paste0("[", lower.CL, ", ", upper.CL, "]")
    ) %>%
    select(all_of(condition), all_of(grouping), emmean, CI, SE) %>%
    pivot_wider(
      names_from = all_of(grouping),
      values_from = c(emmean, CI, SE),
      names_vary = "slowest"
    )
}

format_df_contrast <- function(df_emmeans, lmm, condition, n_digits = 2){
  
  df_pairs <- as.data.frame(pairs(
    df_emmeans, 
    by = all_of(condition), 
    adjust="bonferroni"
  ))
  
  df_effsize <- 
    df_formatted <- merge(
      as.data.frame(
        pairs(
          df_emmeans,
          by = all_of(condition),
          adjust="bonferroni"
        )
      ),
      as.data.frame(
        eff_size(
          df_emmeans,
          sigma = sigma(lmm),
          edf = 26
        )
      ) %>%
        select(
          contrast, all_of(condition), effect.size,
          lower.CL, upper.CL
        )
    ) %>%
    mutate(
      p.value = ifelse(p.value < 0.001, "<.001",
                       sprintf("%0.3f", p.value)
      ),
      estimate = sprintf(paste0("%0.", n_digits, "f"), estimate),
      SE = sprintf(paste0("%0.", n_digits, "f"), SE),
      across(where(is.numeric), ~ sprintf(paste0("%0.", n_digits, "f"), .)),
      CI = sprintf("[%s, %s]", lower.CL, upper.CL)
    ) %>%
    dplyr::select(-c(contrast, lower.CL, upper.CL))
  
  return(df_formatted)
}

format_df_means_and_comparison <- function(
    emmeans_grid, 
    lmm) {
  
  # Dataframe with the estimated marginal means for each group
  df_emmeans <- as.data.frame(emmeans_grid) %>%
    mutate(
      across(
        .cols = where(is.numeric) & !matches("^df$"),
        ~ sprintf("%0.2f", .)
      ),
      CI = sprintf("[%s, %s]", lower.CL, upper.CL)
    ) %>%
    select(group, emmean, CI, SE) %>%
    setNames(c("Component", "Estimate", "95\\% CI (Estimate)", "SE"))
  
  # Dataframe with the results of the comparison
  df_comparison <- as.data.frame(pairs(emmeans_grid, adjust = "bonferroni")) %>%
    mutate(
      p.value = ifelse(p.value < 0.001, "<.001",
                       sprintf("%0.3f", p.value)
      ),
      df = sprintf("%g", df),
      across(where(is.numeric), ~ sprintf("%0.2f", .))
    )
  
  # Dataframe with the effect size of the comparison
  df_effsize <- as.data.frame(
    eff_size(
      emmeans_grid,
      sigma = sigma(lmm),
      edf = 26
    )
  ) %>%
    mutate(
      across(where(is.numeric), ~sprintf("%0.2f", .)),
      CI = sprintf("[%s, %s]", lower.CL, upper.CL)
    ) %>%
    select(contrast, effect.size, CI)
  
  # Combine the two dataframes with values about the comparison
  df_comparison_effsize <- merge(df_comparison, df_effsize) %>%
    setNames(c(
      "Component", "Estimate", "SE", "DF", 
      "$t$-value", "$p$-value", "$d$", "95\\% CI ($d$)"
    ))
  
  # Merge dataframe with the mean and the dataframe with the comparison
  df_formatted <- bind_rows(df_emmeans, df_comparison_effsize)
  
  return(df_formatted)
}


format_df_lm_results <- function(model, n_decimals = 2) {
  model_summary <- summary(model)
  
  # Extract coefficients as a data frame
  df <- as.data.frame(model_summary$coefficients) %>%
    setNames(c("Estimate", "SE", "t_value", "p_value")) %>%
    mutate(
      DF = model_summary$df[2], # Residual degrees of freedom
      term = rownames(model_summary$coefficients),
      r_squared = NA,
      sigma = NA
    )
  
  # Move 'term' to the first column
  df <- df %>%
    select(term, Estimate, SE, DF, t_value, p_value, r_squared, sigma)
  
  # Create model fit row
  model_fit <- tibble::tibble(
    term = "Model fit",
    Estimate = NA_real_,
    SE = NA_real_,
    DF = NA,
    t_value = NA_real_,
    p_value = NA_real_,
    r_squared = if (!is.null(model_summary$r.squared)) {
      if (model_summary$r.squared >= 0.01) {
        sprintf("%0.2f", model_summary$r.squared)
      } else {
        "<.01"
      }
    } else {
      NA
    },
    sigma = if (!is.null(model_summary$sigma)) {
      sprintf("%0.2f", model_summary$sigma)
    } else {
      NA
    }
  )
  
  # Bind and rename columns for final output
  df_out <- bind_rows(df, model_fit) %>%
    setNames(c(" ", "Estimate", "SE", "DF", "$t$-value", "$p$-value", "$R^2$", "$\\sigma_r$")) %>%
    select(` `, Estimate, SE, DF, `$t$-value`, `$p$-value`, `$R^2$`, `$\\sigma_r$`) %>%
    mutate(
      Estimate = ifelse(is.na(Estimate), NA, sprintf(paste0("%0.", n_decimals, "f"), Estimate)),
      SE = ifelse(is.na(SE), NA, sprintf(paste0("%0.", n_decimals, "f"), SE)),
      `$t$-value` = ifelse(is.na(`$t$-value`), NA, sprintf("%0.2f", as.numeric(`$t$-value`))),
      `$p$-value` = ifelse(
        is.na(`$p$-value`), NA,
        ifelse(`$p$-value` < 0.001, "<.001", sprintf("%0.3f", `$p$-value`))
      )
    )
  
  return(df_out)
}
