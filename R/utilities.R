#' Normal Distribution Data Fit
#'
#' This function fits a normal distribution to the given data and provides various statistics and graphs.
#'
#' @param data A numeric vector of data to be fitted.
#' @param pertinent_graphs A logical value indicating whether to generate pertinent graphs.
#' @param desired_quantiles A numeric vector of quantiles to be computed.
#' @param value_of_alpha_for_CI A numeric value for the level of significance for confidence intervals.
#' @param value_of_alpha_for_CR A numeric value for the level of significance for confidence regions.
#'
#' @return A list containing various summary statistics, test results, and graphs (if pertinent_graphs = TRUE).
#'
#' @examples
#' \dontrun{
#' data <- rnorm(100)
#' Normal_Distribution_Data_Fit(data, TRUE, c(0.1, 0.5, 0.9), 0.05, 0.05)
#' }
#'
#' @export
Normal_Distribution_Data_Fit <- function(data, pertinent_graphs, desired_quantiles, value_of_alpha_for_CI, value_of_alpha_for_CR) {
  
  
  
  #       ----------------   Actuarial Coding Portfolio  ----------------
  
  #                         Ioannis (Yanni) Papadopoulos
  
  #       =================================================================
  
  
  #   Before we begin, we will load all the necessary libraries.
  
  
  
  
  #   We will now begin by loading the data.
  #   We will proceed with the following method to load the data for reproducibility reasons:
  
  
  
  
  #       ======    Part I:    Model Selection    ======
  
  
  # *           -----   Descriptive Statistics  -----
  
  
  #   Preliminary summary statistics:
  
  preliminary_summary_statistics = base::summary(data)
  
  
  #   Observing the sample size:
  
  sample_size = base::length(data)
  
  
  #   Observing the minimum, maximum, and range:
  
  minimum = base::min(data)
  
  maximum = base::max(data)
  
  
  #   The minimum and maximum of the data can also be observed using the following function:
  
  min_and_max_of_data = base::range(data)
  
  #   To find the range, we will take the difference of the maximum and the minimum:
  
  range_difference = maximum - minimum
  
  
  
  #   Observing the mean of the data:
  
  mean_of_data = base::mean(data)
  
  
  
  #   Observing the median of the data:
  
  median_of_data = stats::median(data)
  
  
  
  #   Observing the variance and standard deviation of the data:
  
  variance_of_data = actuar::var(data)
  
  
  standard_deviation_of_data = actuar::sd(data)
  
  
  
  #   Observing the coefficient of variation for the data:
  
  coefficient_of_variation = standard_deviation_of_data/mean_of_data
  
  
  
  #   Calculation of Quantiles:
  
  #       The following are the required quantiles:
  
  Required_Quantiles = desired_quantiles
  
  
  #      We now observe the required quantile values:
  
  quantile_values = stats::quantile(data, Required_Quantiles)
  
  
  #   Observation of skewness and kurtosis.
  #     For this, the 'moments' package is used.
  
  Skewness_of_data = EnvStats::skewness(data)
  
  
  Kurtosis_of_data = EnvStats::kurtosis(data, excess = FALSE)
  
  
  
  
  
  # *           -----   Supporting Graphs  -----
  
  
  #   In this section, the required exploratory graphs are produced.
  
  #   We will also sort the data in increasing order for when such an operation would be appropriate.
  
  
  data_increasing_order = base::sort(data)
  
  
  #   Boxplot:
  
  
  graph_of_boxplot = function(data) {
    
    graphics::boxplot(data)
    graphics::title(main = "Boxplot", ylab = "Data")
    
  }
  
  
  
  #   Histogram:
  
  
  graph_of_histogram = function(data) {
    
    graphics::hist(data, main="Histogram for Data", xlab="Data", col="green")
    
  }
  
  
  
  
  #   Empirical CDF:
  
  eCDF_of_data = stats::ecdf(data)
  
  graph_of_eCDF = function(data) {
    
    #    eCDF_of_data = ecdf(data)
    graphics::plot(eCDF_of_data, main="Empirical CDF for Data", xlab="Data", ylab="F_n")
    
  }
  
  
  
  
  #   In preparation of the Empirical MRL (Mean Residual Life) function, a CDF function is defined.
  
  CDF = function(x){
    
    a = x
    
    for(i in 1:base::length(x) ) {
      
      a[i] = base::sum(x<=x[i]) / base::length(x)
      
    }
    return(a)
  }
  
  
  #   The Empirical MRL function is as follows:
  
  eMRL = function(x){
    
    a = x
    
    x.Fn = CDF(x)
    
    for(i in 1:base::length(x)){
      
      if(x[i] < base::max(x)){
        
        a[i] = base::sum(x[x - x[i] > 0] - x[i]) / (base::length(x)) / (1 - x.Fn[i])
        
      } else {
        
        a[i] = 0
        
      }
      
    }
    
    return(a)
    
  }
  
  
  
  
  #   The graph of the Empirical MRL is as follows:
  
  data_eMRL = eMRL(data_increasing_order)
  
  graph_of_eMRL = function(data) {
    
    #    data_eMRL = eMRL(data_increasing_order)
    graphics::plot(data_increasing_order, data_eMRL, type="l", xlab="Data", ylab="e_n", main="Empirical MRL")
    
  }
  
  
  
  
  
  #   The graph of the Empirical LEV is as follows:
  
  data_eLEV = base::mean(data_increasing_order) - data_eMRL * (1 - eCDF_of_data(data_increasing_order))
  
  graph_of_eLEV = function(data) {
    
    #    data_eLEV = mean(data_increasing_order) - data_eMRL * (1 - eCDF_of_data(data_increasing_order))
    graphics::plot(data_increasing_order, data_eLEV, type="l", xlab="Data", ylab="Limited Expected Value", main="Empirical LEV", las = 1)
    
  }
  
  
  
  
  
  #   To verify the validity of the assumptions, an IID test is performed:
  
  
  graph_of_IID_Test = function(data) {
    
    IID_Test = stats::acf(x = data, lag.max = sample_size, type = "correlation", main="ACF of Data")
    
  }
  
  
  
  
  
  #       ======    Part II:    Estimation    ======
  
  
  
  # *           -----   Method of Moments (MME) Parameter Estimation  -----
  
  
  #   We will proceed with fitting the model distribution candidates with respect to the Method of Moments (MME) approach:
  
  
  Normal_Fit_MME = fitdistrplus::fitdist(data, "norm", method = "mme")
  
  
  
  
  
  # *           -----   Maximum Likelihood Estimator (MLE) Parameter Estimation  -----
  
  
  #   Now, we will proceed with fitting the model distribution candidates with respect to the Maximum Likelihood Estimator (MLE) approach:
  
  
  Normal_Fit_MLE = fitdistrplus::fitdist(data, "norm")
  
  
  
  
  
  
  # *           -----   Estimation of Asymptotic Covariance Matrices  -----
  
  
  #   The Information Matrices with respect to each candidate distribution are observed as follows:
  
  
  #     Normal:
  
  Normal_Fit_MLE_Estimates = Normal_Fit_MLE$estimate
  
  Normal_Fit_MLE_Estimated_mean = base::as.numeric(Normal_Fit_MLE_Estimates[1])
  Normal_Fit_MLE_Estimated_sd = base::as.numeric(Normal_Fit_MLE_Estimates[2])
  
  Normal_vcov_Matrix = Normal_Fit_MLE$vcov
  
  Normal_Information_Matrix = matlib::Inverse(Normal_vcov_Matrix)
  
  
  
  
  
  
  
  
  
  
  # *           -----   Confidence Regions  -----
  
  
  #   Note:   The Confidence Regions (CR) will be constructed after the Formal Goodness of Fit test have been conducted, which is at the end of the script.
  #               This section will be denoted as "Confidence Regions, revisited"
  
  #   Construction of a Confidence Interval (CI) with respect to each candidate distribution's estimators is observed as follows:
  
  #   The chosen level of confidence is (value_of_alpha_for_CI)% with respect to a two-tailed test.
  
  alpha_for_CI = value_of_alpha_for_CI
  Level_of_Confidence.2_T = 1-(alpha_for_CI/2)
  
  Z_Value_at_LoC.2_T = stats::qnorm(Level_of_Confidence.2_T)
  
  
  #   We now proceed.
  
  
  #     Normal:
  
  Lower_bound_of_CI_for_MLE_Normal_mean = Normal_Fit_MLE_Estimated_mean - Z_Value_at_LoC.2_T*base::sqrt(Normal_vcov_Matrix[1,1])
  Upper_bound_of_CI_for_MLE_Normal_mean = Normal_Fit_MLE_Estimated_mean + Z_Value_at_LoC.2_T*base::sqrt(Normal_vcov_Matrix[1,1])
  
  
  CI_of_MLE_Normal_mean = base::as.numeric(
    base::c(Lower_bound_of_CI_for_MLE_Normal_mean, Upper_bound_of_CI_for_MLE_Normal_mean)
  )
  
  
  
  
  Lower_bound_of_CI_for_MLE_Normal_sd = Normal_Fit_MLE_Estimated_sd - Z_Value_at_LoC.2_T*base::sqrt(Normal_vcov_Matrix[2,2])
  Upper_bound_of_CI_for_MLE_Normal_sd = Normal_Fit_MLE_Estimated_sd + Z_Value_at_LoC.2_T*base::sqrt(Normal_vcov_Matrix[2,2])
  
  
  CI_of_MLE_Normal_sd = base::as.numeric(
    base::c(Lower_bound_of_CI_for_MLE_Normal_sd, Upper_bound_of_CI_for_MLE_Normal_sd)
  )
  
  
  
  
  
  
  #       ======    Part III:    Goodness-of-fit Tests    ======
  
  
  
  #     Creating the required functions for each distribution before proceeding.
  
  
  #      Normal
  
  #   CDF:
  
  CDF_Normal_Fit = function(x){
    
    stats::pnorm(x, mean = Normal_Fit_MLE_Estimated_mean, sd = Normal_Fit_MLE_Estimated_sd)
    
  }
  
  #   Survival Function:
  
  Survival_Function_Normal_Fit = function(x){
    
    1 - CDF_Normal_Fit(x)
    
  }
  
  
  
  #   Probability Density Function:
  
  PDF_Normal_Fit = function(x){
    
    stats::dnorm(x, mean = Normal_Fit_MLE_Estimated_mean, sd = Normal_Fit_MLE_Estimated_sd)
    
  }
  
  
  
  #   Limited Expected Value:
  
  LEV_Normal_Fit = function(t){
    
    stats::integrate(f = Survival_Function_Normal_Fit, lower = 0, upper = t)$value
    
  }
  
  
  
  
  #   Mean Residual Life:
  
  MRL_Normal_Fit = function(t){
    
    ( LEV_Normal_Fit(Inf) - LEV_Normal_Fit(t) ) / Survival_Function_Normal_Fit(t)
    
  }
  
  
  
  
  
  
  # *           -----   Graphical Comparisons  -----
  
  
  #   Now, we will create a graph comparing the Empirical CDF with the Fitted Distributions' CDFs.
  
  #     Firstly, we will create a sequence to aid in the creation of the graph.
  
  Sequence_for_Graphs = base::seq(minimum, maximum, length = sample_size)
  
  #     Also, we will introduce a multiplier for "par("cex")" as well as a line width parameter for flexibility:
  
  Multiplier_for_par.cex = 3
  
  Line_Width = Multiplier_for_par.cex*graphics::par("cex")
  
  
  Empirical_Distribution.text = "Empirical Distribution"
  
  
  Collection_of_Fitted_Distributions = base::c(
    "Fitted Normal Distribution"
    #    ,"Fitted Gamma Distribution"
    #    ,"Fitted Weibull Distribution"
    #    ,"Fitted Inverse Gamma Distribution"
    #    ,"Fitted Inverse Weibull Distribution"
  )
  
  
  Collection_of_Colours = base::c(
    "blue"
    #    ,"orange"
    #    ,"purple"
    #    ,"red"
    #    ,"turquoise"
  )
  
  
  
  #   Graphing the Empirical CDF with the Normal Fitted CDF overlay:
  
  
  graph_of_eCDF_with_Normal_Fitted_CDF_overlay = function() {
    
    
    #   The Empirical CDF is graphed as follows:
    
    EnvStats::ecdfPlot(data_increasing_order, main = "Graphical Comparison of the eCDF with the Fitted CDFs", xlab = "Data")
    
    
    #   We will now overlay the fitted distributions' CDFs.
    
    Normal_Fitted_CDF = base::sapply(X = Sequence_for_Graphs, CDF_Normal_Fit)
    graphics::lines(x = Sequence_for_Graphs, y = Normal_Fitted_CDF, lwd = Line_Width, col = Collection_of_Colours[1])
    
    
    #   Allowing the legend to be placed outside the plot area:
    graphics::par(xpd = TRUE)
    
    
    # Call locator() in the console and click on the plot to choose the legend position
    legend_position = graphics::locator(1)
    
    # Use the returned coordinates to place the legend
    graphics::legend(x = legend_position$x, y = legend_position$y, base::c(Empirical_Distribution.text, Collection_of_Fitted_Distributions), fill = base::c("black", Collection_of_Colours))
    
    
    # Revoking the allowance of the legend to be placed outside the plot area:
    graphics::par(xpd = FALSE)
    
  }
  
  
  
  
  
  
  
  
  #   Graphing the Empirical LEV with the Normal Fitted LEV overlay:
  
  graph_of_eLEV_with_Normal_Fitted_LEV_overlay = function() {
    
    #   The Empirical LEV is graphed as follows:
    
    graphics::plot(actuar::elev(data), type="l", xlab = "Data", ylab = "Limited Expected Value", main="Graphical Comparison of the eLEV with the Fitted LEVs", lwd = Line_Width, las = 1)
    
    
    #   We will now overlay the fitted distributions' LEVs.
    
    Normal_Fitted_LEV = base::sapply(X = Sequence_for_Graphs, LEV_Normal_Fit)
    graphics::lines(x = Sequence_for_Graphs, y = Normal_Fitted_LEV, lwd = Line_Width, col = Collection_of_Colours[1])
    
    
    #   Allowing the legend to be placed outside the plot area:
    graphics::par(xpd = TRUE)
    
    
    # Call locator() in the console and click on the plot to choose the legend position
    legend_position = graphics::locator(1)
    
    # Use the returned coordinates to place the legend
    graphics::legend(x = legend_position$x, y = legend_position$y, base::c(Empirical_Distribution.text, Collection_of_Fitted_Distributions), fill = base::c("black", Collection_of_Colours))
    
    
    # Revoking the allowance of the legend to be placed outside the plot area:
    graphics::par(xpd = FALSE)
    
  }
  
  
  
  
  
  
  
  
  #   Graphing the Empirical MRL with the Normal Fitted MRL overlay:
  
  graph_of_eMRL_with_Normal_Fitted_MRL_overlay = function() {
    
    #   The Empirical MRL is graphed as follows:
    
    graphics::plot(data_increasing_order, data_eMRL, type="l", xlab="Data", ylab="e_n", main="Graphical Comparison of the eMRL with the Fitted MRLs", lwd = Line_Width)
    
    
    #   We will now overlay the fitted distributions' MRLs.
    
    Normal_Fitted_MRL = base::sapply(X = Sequence_for_Graphs, MRL_Normal_Fit)
    graphics::lines(x = Sequence_for_Graphs, y = Normal_Fitted_MRL, lwd = Line_Width, col = Collection_of_Colours[1])
    
    
    #   Allowing the legend to be placed outside the plot area:
    graphics::par(xpd = TRUE)
    
    
    # Call locator() in the console and click on the plot to choose the legend position
    legend_position = graphics::locator(1)
    
    # Use the returned coordinates to place the legend
    graphics::legend(x = legend_position$x, y = legend_position$y, base::c(Empirical_Distribution.text, Collection_of_Fitted_Distributions), fill = base::c("black", Collection_of_Colours))
    
    
    # Revoking the allowance of the legend to be placed outside the plot area:
    graphics::par(xpd = FALSE)
    
  }
  
  
  
  
  
  
  
  
  
  #   Graphing the Empirical Histogram with the Normal Fitted PDF overlay:
  
  graph_of_histogram_with_Normal_Fitted_PDF_overlay = function() {
    
    #   The histogram of the empirical data is graphed as follows:
    
    graphics::hist(data, prob = T, main="Empirical Histogram compared to the Fitted Empirical Densities", xlab="Data", col="green")
    
    #   We will now overlay the fitted distributions' PDFs altogether.
    
    graphics::lines(data_increasing_order, stats::dnorm(data_increasing_order, Normal_Fit_MLE_Estimated_mean, Normal_Fit_MLE_Estimated_sd), xpd = T, lwd = Line_Width, col = Collection_of_Colours[1])
    
    
    #   Allowing the legend to be placed outside the plot area:
    graphics::par(xpd = TRUE)
    
    
    # Call locator() in the console and click on the plot to choose the legend position
    legend_position = graphics::locator(1)
    
    # Use the returned coordinates to place the legend
    graphics::legend(x = legend_position$x, y = legend_position$y, base::c(Empirical_Distribution.text, Collection_of_Fitted_Distributions), fill = base::c("black", Collection_of_Colours))
    
    
    # Revoking the allowance of the legend to be placed outside the plot area:
    graphics::par(xpd = FALSE)
    
  }
  
  
  
  
  
  
  
  #   Graphing the P-P Plot with the X = Y line overlay:
  
  graph_of_P.P_Plot_with_X.Y_line_overlay = function() {
    
    #     P-P Plots
    
    Legend_Content_for_P.P_Plots = base::c("(Fn, F_Fitted)","X = Y")
    
    Colour_of_X.Y_Line_for_P.P_Plots = "blue"
    Legend_Content_Fill_for_P.P_Plots = base::c("black", Colour_of_X.Y_Line_for_P.P_Plots)
    
    
    
    #   Normal P-P Plot:
    
    ppoints_vector = stats::ppoints(sample_size)
    pp_plot = stats::qnorm(ppoints_vector, Normal_Fit_MLE_Estimated_mean, Normal_Fit_MLE_Estimated_sd)
    graphics::plot(data_increasing_order, pp_plot, main = "Normal P-P Plot for Data", xlab = "Empirical CDF", ylab = "Fitted Normal CDF")
    graphics::abline(0, 1, lwd = Line_Width, col = Colour_of_X.Y_Line_for_P.P_Plots)
    
    
    #   Allowing the legend to be placed outside the plot area:
    graphics::par(xpd = TRUE)
    
    
    # Call locator() in the console and click on the plot to choose the legend position
    legend_position = graphics::locator(1)
    
    # Use the returned coordinates to place the legend
    graphics::legend(x = legend_position$x, y = legend_position$y, Legend_Content_for_P.P_Plots, fill = Legend_Content_Fill_for_P.P_Plots)
    
    
    # Revoking the allowance of the legend to be placed outside the plot area:
    graphics::par(xpd = FALSE)
    
  }
  
  
  
  
  
  
  
  #     Q-Q Plots
  
  
  
  
  #   The following Q-Q Plots compare the graphs with respect to a line of least squares, the 0-1 (y-intercept = 0, slope = 1) line, and the robust line (a line that is fit between the first and third quartiles of the data).
  
  
  Q.Q_Plot.L_S.function = function() {
    
    EnvStats::qqPlot(data_increasing_order, distribution = "norm", param.list = base::list(mean = Normal_Fit_MLE_Estimated_mean, sd = Normal_Fit_MLE_Estimated_sd), add.line = TRUE, qq.line.type = "least squares")
    
  }
  
  
  
  
  
  
  Q.Q_Plot.0_1.function = function() {
    
    EnvStats::qqPlot(data_increasing_order, distribution = "norm", param.list = base::list(mean = Normal_Fit_MLE_Estimated_mean, sd = Normal_Fit_MLE_Estimated_sd), add.line = TRUE, qq.line.type = "0-1")
    
  }
  
  
  
  
  
  
  Q.Q_Plot.Robust.function = function() {
    
    EnvStats::qqPlot(data_increasing_order, distribution = "norm", param.list = base::list(mean = Normal_Fit_MLE_Estimated_mean, sd = Normal_Fit_MLE_Estimated_sd), add.line = TRUE, qq.line.type = "robust")
    
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  # Define a function to calculate summary statistics
  Fitted_Distribution_Summary_Statistics = function(dist, par, quantiles, stats_fun) {
    q = base::sapply(quantiles, function(x) base::do.call(dist, base::c(base::list(p = x), par)))
    stats = base::do.call(stats_fun, par)
    base::c(quantiles = q, stats)
  }
  
  # Define parameters
  
  fitted_normal_par = base::list(mean = Normal_Fit_MLE_Estimated_mean, sd = Normal_Fit_MLE_Estimated_sd)
  
  
  # Define functions to calculate additional statistics
  
  fitted_normal_stats_fun = function(mean, sd) {
    base::list(
      mean = mean,
      sd = sd,
      cv = sd / mean,
      skew = 0,
      kurt = 3
    )
  }
  
  
  
  
  
  # Calculate summary statistics
  
  fitted_normal_stats = Fitted_Distribution_Summary_Statistics(stats::qnorm, fitted_normal_par, Required_Quantiles, fitted_normal_stats_fun)
  
  
  
  
  
  
  
  
  # *           -----   Formal Goodness-of-Fit Tests  -----
  
  
  
  
  #   Conducting the Kolmogorov-Smirnov Test on each fitted distribution:
  
  KS_Test = stats::ks.test(data, stats::pnorm, Normal_Fit_MLE_Estimated_mean, Normal_Fit_MLE_Estimated_sd)
  
  
  
  
  #   Conducting the Anderson-Darling Test on each fitted distribution:
  
  AD_Test = ADGofTest::ad.test(data, stats::pnorm, Normal_Fit_MLE_Estimated_mean, Normal_Fit_MLE_Estimated_sd)
  
  
  
  
  
  #   Conducting the AIC / BIC Tests on each fitted distribution:
  
  #   Preparing to conduct the tests:
  
  
  Normal_Fit_Preliminary_Summary = fitdistrplus::fitdist(data, "norm")
  
  
  
  
  #   The test results are observable as follows:
  
  Normal_Fit_Summary = base::summary(Normal_Fit_Preliminary_Summary)
  
  
  
  
  
  
  
  
  # *           -----   Confidence Regions, revisited  -----
  
  
  #   Construction of a Confidence Region (CR) with respect to each candidate distribution's estimators is observed as follows:
  
  #   The chosen level of confidence is (value_of_alpha_for_CR)% with respect to the Chi-Square test.
  
  
  alpha_for_CR = value_of_alpha_for_CR
  Level_of_Confidence.1_T = 1-(alpha_for_CR)
  
  Z_Value_at_LoC.1_T = stats::qnorm(Level_of_Confidence.1_T)
  
  
  #   We now proceed.
  
  
  
  graph_of_confidence_region_for_Normal_Fit_parameters = function() {
    
    
    #       Normal:
    
    
    #   Generating the confidence ellipse:
    
    Normal_Fit_MLE_Estimates_CR = ellipse::ellipse(x = Normal_vcov_Matrix, center = Normal_Fit_MLE_Estimates, shape = Normal_vcov_Matrix, radius = Z_Value_at_LoC.1_T)
    
    #   Plotting the confidence ellipse:
    
    graphics::plot(Normal_Fit_MLE_Estimates_CR, type = "l", xlab = "Mean", ylab = "Standard Deviation")
    graphics::polygon(Normal_Fit_MLE_Estimates_CR, col = grDevices::rgb(0, 0, 1, 0.5))  # Adding the shaded region
    
    graphics::points(Normal_Fit_MLE_Estimated_mean, Normal_Fit_MLE_Estimated_sd, pch = 19)
    
    graphics::title(main = "CR for the ML Estimators of the Fitted Normal Distribution")
    
    
  }
  
  
  #   Number of parameters in the Normal distribution is observed as follows:
  
  Number_of_Parameters_Normal_Fit_MLE_Estimates = base::length(stats::coef(Normal_Fit_Summary))
  
  
  #   Quantile of the Chi-Square distribution:
  
  Quantile_of_Chi.Square_Dist_Normal_Fit = stats::qchisq(1 - alpha_for_CR, df = Number_of_Parameters_Normal_Fit_MLE_Estimates)
  
  # Calculating the threshold on the log-likelihood for the likelihood ratio test at the alpha = 5% significance level:
  
  Log.Likelihood_Threshold_Normal_Fit = Normal_Fit_Summary$loglik - 0.5 * Quantile_of_Chi.Square_Dist_Normal_Fit
  
  
  
  
  
  Useful_List = base::list(
    preliminary_summary_statistics = base::summary(data)
    , sample_size = base::length(data)
    , minimum = base::min(data)
    , maximum = base::max(data)
    , min_and_max_of_data = base::range(data)
    , range_difference = maximum - minimum
    , mean_of_data = base::mean(data)
    , median_of_data = stats::median(data)
    , variance_of_data = actuar::var(data)
    , standard_deviation_of_data = actuar::sd(data)
    , coefficient_of_variation = standard_deviation_of_data/mean_of_data
    , Required_Quantiles = desired_quantiles
    , quantile_values = stats::quantile(data, Required_Quantiles)
    , Skewness_of_data = EnvStats::skewness(data)
    , Kurtosis_of_data = EnvStats::kurtosis(data, excess = FALSE)
    , Normal_Fit_MME = fitdistrplus::fitdist(data, "norm", method = "mme")
    , Normal_Fit_MLE = fitdistrplus::fitdist(data, "norm")
    , Normal_Fit_MLE_Estimates = Normal_Fit_MLE$estimate
    , Normal_Fit_MLE_Estimated_mean = base::as.numeric(Normal_Fit_MLE_Estimates[1])
    , Normal_Fit_MLE_Estimated_sd = base::as.numeric(Normal_Fit_MLE_Estimates[2])
    , Normal_vcov_Matrix = Normal_Fit_MLE$vcov
    , Normal_Information_Matrix = matlib::Inverse(Normal_vcov_Matrix)
    , alpha_for_CI = value_of_alpha_for_CI
    , Level_of_Confidence.2_T = 1-(alpha_for_CI/2)
    , Z_Value_at_LoC.2_T = stats::qnorm(Level_of_Confidence.2_T)
    , Lower_bound_of_CI_for_MLE_Normal_mean = Normal_Fit_MLE_Estimated_mean - Z_Value_at_LoC.2_T*base::sqrt(Normal_vcov_Matrix[1,1])
    , Upper_bound_of_CI_for_MLE_Normal_mean = Normal_Fit_MLE_Estimated_mean + Z_Value_at_LoC.2_T*base::sqrt(Normal_vcov_Matrix[1,1])
    , CI_of_MLE_Normal_mean = base::as.numeric(
      base::c(Lower_bound_of_CI_for_MLE_Normal_mean, Upper_bound_of_CI_for_MLE_Normal_mean)
    )
    , Lower_bound_of_CI_for_MLE_Normal_sd = Normal_Fit_MLE_Estimated_sd - Z_Value_at_LoC.2_T*base::sqrt(Normal_vcov_Matrix[2,2])
    , Upper_bound_of_CI_for_MLE_Normal_sd = Normal_Fit_MLE_Estimated_sd + Z_Value_at_LoC.2_T*base::sqrt(Normal_vcov_Matrix[2,2])
    , CI_of_MLE_Normal_sd = base::as.numeric(
      base::c(Lower_bound_of_CI_for_MLE_Normal_sd, Upper_bound_of_CI_for_MLE_Normal_sd)
    )
    , fitted_normal_par = base::list(mean = Normal_Fit_MLE_Estimated_mean, sd = Normal_Fit_MLE_Estimated_sd)
    , fitted_normal_stats = Fitted_Distribution_Summary_Statistics(stats::qnorm, fitted_normal_par, Required_Quantiles, fitted_normal_stats_fun)
    , KS_Test = stats::ks.test(data, stats::pnorm, Normal_Fit_MLE_Estimated_mean, Normal_Fit_MLE_Estimated_sd)
    , AD_Test = ADGofTest::ad.test(data, stats::pnorm, Normal_Fit_MLE_Estimated_mean, Normal_Fit_MLE_Estimated_sd)
    , Normal_Fit_Preliminary_Summary = fitdistrplus::fitdist(data, "norm")
    , Normal_Fit_Summary = base::summary(Normal_Fit_Preliminary_Summary)
    , alpha_for_CR = value_of_alpha_for_CR
    , Level_of_Confidence.1_T = 1-(alpha_for_CR)
    , Z_Value_at_LoC.1_T = stats::qnorm(Level_of_Confidence.1_T)
    , Number_of_Parameters_Normal_Fit_MLE_Estimates = base::length(stats::coef(Normal_Fit_Summary))
    , Quantile_of_Chi.Square_Dist_Normal_Fit = stats::qchisq(1 - alpha_for_CR, df = Number_of_Parameters_Normal_Fit_MLE_Estimates)
    , Log.Likelihood_Threshold_Normal_Fit = Normal_Fit_Summary$loglik - 0.5 * Quantile_of_Chi.Square_Dist_Normal_Fit
    
    , if(pertinent_graphs == TRUE) {
      
      base::c(
        Graph_of_eCDF_with_Normal_Fitted_CDF_Overlay = graph_of_eCDF_with_Normal_Fitted_CDF_overlay()
        , Graph_of_eLEV_with_Normal_Fitted_LEV_Overlay = graph_of_eLEV_with_Normal_Fitted_LEV_overlay()
        , Graph_of_eMRL_with_Normal_Fitted_MRL_Overlay = graph_of_eMRL_with_Normal_Fitted_MRL_overlay()
        , Graph_of_Histogram_with_Normal_Fitted_PDF_Overlay = graph_of_histogram_with_Normal_Fitted_PDF_overlay()
        , Graph_of_P.P_Plot_with_X.Y_Line_Overlay = graph_of_P.P_Plot_with_X.Y_line_overlay()
        , Q.Q_Plot.0_1 = Q.Q_Plot.0_1.function()
        , Q.Q_Plot.L_S = Q.Q_Plot.L_S.function()
        , Q.Q_Plot.Robust = Q.Q_Plot.Robust.function()
        , Graph_of_Confidence_Region_for_Normal_Fit_Parameters = graph_of_confidence_region_for_Normal_Fit_parameters()
      )
      
      
    } else {}
    
    
  )
  
  
  Useful_List
  
  
}

#' Gamma Distribution Data Fit
#'
#' This function fits a gamma distribution to the given data and provides various statistics and graphs.
#'
#' @param data A numeric vector of data to be fitted.
#' @param pertinent_graphs A logical value indicating whether to generate pertinent graphs.
#' @param desired_quantiles A numeric vector of quantiles to be computed.
#' @param value_of_alpha_for_CI A numeric value for the level of significance for confidence intervals.
#' @param value_of_alpha_for_CR A numeric value for the level of significance for confidence regions.
#'
#' @return A list containing various summary statistics, test results, and graphs (if pertinent_graphs = TRUE).
#'
#' @examples
#' \dontrun{
#' data <- rgamma(100)
#' Gamma_Distribution_Data_Fit(data, TRUE, c(0.1, 0.5, 0.9), 0.05, 0.05)
#' }
#'
#' @export
Gamma_Distribution_Data_Fit <- function(data, pertinent_graphs, desired_quantiles, value_of_alpha_for_CI, value_of_alpha_for_CR) {
  
  
  
  #       ----------------   Actuarial Coding Portfolio  ----------------
  
  #                         Ioannis (Yanni) Papadopoulos
  
  #       =================================================================
  
  
  #   Before we begin, we will load all the necessary libraries.
  
  
  
  
  #   We will now begin by loading the data.
  #   We will proceed with the following method to load the data for reproducibility reasons:
  
  
  
  
  #       ======    Part I:    Model Selection    ======
  
  
  # *           -----   Descriptive Statistics  -----
  
  
  #   Preliminary summary statistics:
  
  preliminary_summary_statistics = base::summary(data)
  
  
  #   Observing the sample size:
  
  sample_size = base::length(data)
  
  
  #   Observing the minimum, maximum, and range:
  
  minimum = base::min(data)
  
  maximum = base::max(data)
  
  
  #   The minimum and maximum of the data can also be observed using the following function:
  
  min_and_max_of_data = base::range(data)
  
  #   To find the range, we will take the difference of the maximum and the minimum:
  
  range_difference = maximum - minimum
  
  
  
  #   Observing the mean of the data:
  
  mean_of_data = base::mean(data)
  
  
  
  #   Observing the median of the data:
  
  median_of_data = stats::median(data)
  
  
  
  #   Observing the variance and standard deviation of the data:
  
  variance_of_data = actuar::var(data)
  
  
  standard_deviation_of_data = actuar::sd(data)
  
  
  
  #   Observing the coefficient of variation for the data:
  
  coefficient_of_variation = standard_deviation_of_data/mean_of_data
  
  
  
  #   Calculation of Quantiles:
  
  #       The following are the required quantiles:
  
  Required_Quantiles = desired_quantiles
  
  
  #      We now observe the required quantile values:
  
  quantile_values = stats::quantile(data, Required_Quantiles)
  
  
  #   Observation of skewness and kurtosis.
  #     For this, the 'moments' package is used.
  
  Skewness_of_data = EnvStats::skewness(data)
  
  
  Kurtosis_of_data = EnvStats::kurtosis(data, excess = FALSE)
  
  
  
  
  
  # *           -----   Supporting Graphs  -----
  
  
  #   In this section, the required exploratory graphs are produced.
  
  #   We will also sort the data in increasing order for when such an operation would be appropriate.
  
  
  data_increasing_order = base::sort(data)
  
  
  #   Boxplot:
  
  
  graph_of_boxplot = function(data) {
    
    graphics::boxplot(data)
    graphics::title(main = "Boxplot", ylab = "Data")
    
  }
  
  
  
  #   Histogram:
  
  
  graph_of_histogram = function(data) {
    
    graphics::hist(data, main="Histogram for Data", xlab="Data", col="green")
    
  }
  
  
  
  
  #   Empirical CDF:
  
  eCDF_of_data = stats::ecdf(data)
  
  graph_of_eCDF = function(data) {
    
    #    eCDF_of_data = ecdf(data)
    graphics::plot(eCDF_of_data, main="Empirical CDF for Data", xlab="Data", ylab="F_n")
    
  }
  
  
  
  
  #   In preparation of the Empirical MRL (Mean Residual Life) function, a CDF function is defined.
  
  CDF = function(x){
    
    a = x
    
    for(i in 1:base::length(x) ) {
      
      a[i] = base::sum(x<=x[i]) / base::length(x)
      
    }
    return(a)
  }
  
  
  #   The Empirical MRL function is as follows:
  
  eMRL = function(x){
    
    a = x
    
    x.Fn = CDF(x)
    
    for(i in 1:base::length(x)){
      
      if(x[i] < base::max(x)){
        
        a[i] = base::sum(x[x - x[i] > 0] - x[i]) / (base::length(x)) / (1 - x.Fn[i])
        
      } else {
        
        a[i] = 0
        
      }
      
    }
    
    return(a)
    
  }
  
  
  
  
  #   The graph of the Empirical MRL is as follows:
  
  data_eMRL = eMRL(data_increasing_order)
  
  graph_of_eMRL = function(data) {
    
    #    data_eMRL = eMRL(data_increasing_order)
    graphics::plot(data_increasing_order, data_eMRL, type="l", xlab="Data", ylab="e_n", main="Empirical MRL")
    
  }
  
  
  
  
  
  #   The graph of the Empirical LEV is as follows:
  
  data_eLEV = base::mean(data_increasing_order) - data_eMRL * (1 - eCDF_of_data(data_increasing_order))
  
  graph_of_eLEV = function(data) {
    
    #    data_eLEV = mean(data_increasing_order) - data_eMRL * (1 - eCDF_of_data(data_increasing_order))
    graphics::plot(data_increasing_order, data_eLEV, type="l", xlab="Data", ylab="Limited Expected Value", main="Empirical LEV", las = 1)
    
  }
  
  
  
  
  
  #   To verify the validity of the assumptions, an IID test is performed:
  
  
  graph_of_IID_Test = function(data) {
    
    IID_Test = stats::acf(x = data, lag.max = sample_size, type = "correlation", main="ACF of Data")
    
  }
  
  
  
  
  
  #       ======    Part II:    Estimation    ======
  
  
  
  # *           -----   Method of Moments (MME) Parameter Estimation  -----
  
  
  #   We will proceed with fitting the model distribution candidates with respect to the Method of Moments (MME) approach:
  
  
  Gamma_Fit_MME = fitdistrplus::fitdist(data, "gamma", method = "mme")
  
  
  
  
  
  # *           -----   Maximum Likelihood Estimator (MLE) Parameter Estimation  -----
  
  
  #   Now, we will proceed with fitting the model distribution candidates with respect to the Maximum Likelihood Estimator (MLE) approach:
  
  
  Gamma_Fit_MLE = fitdistrplus::fitdist(data, "gamma")
  
  
  
  
  
  
  # *           -----   Estimation of Asymptotic Covariance Matrices  -----
  
  
  #   The Information Matrices with respect to each candidate distribution are observed as follows:
  
  
  #     Gamma:
  
  Gamma_Fit_MLE_Estimates = Gamma_Fit_MLE$estimate
  
  Gamma_Fit_MLE_Estimated_shape = base::as.numeric(Gamma_Fit_MLE_Estimates[1])
  Gamma_Fit_MLE_Estimated_rate = base::as.numeric(Gamma_Fit_MLE_Estimates[2])
  
  Gamma_vcov_Matrix = Gamma_Fit_MLE$vcov
  
  Gamma_Information_Matrix = matlib::Inverse(Gamma_vcov_Matrix)
  
  
  
  
  
  
  
  
  
  
  # *           -----   Confidence Regions  -----
  
  
  #   Note:   The Confidence Regions (CR) will be constructed after the Formal Goodness of Fit test have been conducted, which is at the end of the script.
  #               This section will be denoted as "Confidence Regions, revisited"
  
  #   Construction of a Confidence Interval (CI) with respect to each candidate distribution's estimators is observed as follows:
  
  #   The chosen level of confidence is (value_of_alpha_for_CI)% with respect to a two-tailed test.
  
  alpha_for_CI = value_of_alpha_for_CI
  Level_of_Confidence.2_T = 1-(alpha_for_CI/2)
  
  Z_Value_at_LoC.2_T = stats::qnorm(Level_of_Confidence.2_T)
  
  
  #   We now proceed.
  
  
  #     Gamma:
  
  Lower_bound_of_CI_for_MLE_Gamma_shape = Gamma_Fit_MLE_Estimated_shape - Z_Value_at_LoC.2_T*base::sqrt(Gamma_vcov_Matrix[1,1])
  Upper_bound_of_CI_for_MLE_Gamma_shape = Gamma_Fit_MLE_Estimated_shape + Z_Value_at_LoC.2_T*base::sqrt(Gamma_vcov_Matrix[1,1])
  
  
  CI_of_MLE_Gamma_shape = base::as.numeric(
    base::c(Lower_bound_of_CI_for_MLE_Gamma_shape, Upper_bound_of_CI_for_MLE_Gamma_shape)
  )
  
  
  
  
  Lower_bound_of_CI_for_MLE_Gamma_rate = Gamma_Fit_MLE_Estimated_rate - Z_Value_at_LoC.2_T*base::sqrt(Gamma_vcov_Matrix[2,2])
  Upper_bound_of_CI_for_MLE_Gamma_rate = Gamma_Fit_MLE_Estimated_rate + Z_Value_at_LoC.2_T*base::sqrt(Gamma_vcov_Matrix[2,2])
  
  
  CI_of_MLE_Gamma_rate = base::as.numeric(
    base::c(Lower_bound_of_CI_for_MLE_Gamma_rate, Upper_bound_of_CI_for_MLE_Gamma_rate)
  )
  
  
  
  
  
  
  #       ======    Part III:    Goodness-of-fit Tests    ======
  
  
  
  #     Creating the required functions for each distribution before proceeding.
  
  
  #      Gamma
  
  #   CDF:
  
  CDF_Gamma_Fit = function(x){
    
    stats::pgamma(x, shape = Gamma_Fit_MLE_Estimated_shape, rate = Gamma_Fit_MLE_Estimated_rate)
    
  }
  
  #   Survival Function:
  
  Survival_Function_Gamma_Fit = function(x){
    
    1 - CDF_Gamma_Fit(x)
    
  }
  
  
  
  #   Probability Density Function:
  
  PDF_Gamma_Fit = function(x){
    
    stats::dgamma(x, shape = Gamma_Fit_MLE_Estimated_shape, rate = Gamma_Fit_MLE_Estimated_rate)
    
  }
  
  
  
  #   Limited Expected Value:
  
  LEV_Gamma_Fit = function(t){
    
    stats::integrate(f = Survival_Function_Gamma_Fit, lower = 0, upper = t)$value
    
  }
  
  
  
  
  #   Mean Residual Life:
  
  MRL_Gamma_Fit = function(t){
    
    ( LEV_Gamma_Fit(Inf) - LEV_Gamma_Fit(t) ) / Survival_Function_Gamma_Fit(t)
    
  }
  
  
  
  
  
  
  # *           -----   Graphical Comparisons  -----
  
  
  #   Now, we will create a graph comparing the Empirical CDF with the Fitted Distributions' CDFs.
  
  #     Firstly, we will create a sequence to aid in the creation of the graph.
  
  Sequence_for_Graphs = base::seq(minimum, maximum, length = sample_size)
  
  #     Also, we will introduce a multiplier for "par("cex")" as well as a line width parameter for flexibility:
  
  Multiplier_for_par.cex = 3
  
  Line_Width = Multiplier_for_par.cex*graphics::par("cex")
  
  
  Empirical_Distribution.text = "Empirical Distribution"
  
  
  Collection_of_Fitted_Distributions = base::c(
    #     "Fitted Normal Distribution"
    "Fitted Gamma Distribution"
    #    ,"Fitted Weibull Distribution"
    #    ,"Fitted Inverse Gamma Distribution"
    #    ,"Fitted Inverse Weibull Distribution"
  )
  
  
  Collection_of_Colours = base::c(
    #     "blue"
    "orange"
    #    ,"purple"
    #    ,"red"
    #    ,"turquoise"
  )
  
  
  
  #   Graphing the Empirical CDF with the Gamma Fitted CDF overlay:
  
  
  graph_of_eCDF_with_Gamma_Fitted_CDF_overlay = function() {
    
    
    #   The Empirical CDF is graphed as follows:
    
    EnvStats::ecdfPlot(data_increasing_order, main = "Graphical Comparison of the eCDF with the Fitted CDFs", xlab = "Data")
    
    
    #   We will now overlay the fitted distributions' CDFs.
    
    Gamma_Fitted_CDF = base::sapply(X = Sequence_for_Graphs, CDF_Gamma_Fit)
    graphics::lines(x = Sequence_for_Graphs, y = Gamma_Fitted_CDF, lwd = Line_Width, col = Collection_of_Colours[1])
    
    
    #   Allowing the legend to be placed outside the plot area:
    graphics::par(xpd = TRUE)
    
    
    # Call locator() in the console and click on the plot to choose the legend position
    legend_position = graphics::locator(1)
    
    # Use the returned coordinates to place the legend
    graphics::legend(x = legend_position$x, y = legend_position$y, base::c(Empirical_Distribution.text, Collection_of_Fitted_Distributions), fill = base::c("black", Collection_of_Colours))
    
    
    # Revoking the allowance of the legend to be placed outside the plot area:
    graphics::par(xpd = FALSE)
    
  }
  
  
  
  
  
  
  
  
  #   Graphing the Empirical LEV with the Gamma Fitted LEV overlay:
  
  graph_of_eLEV_with_Gamma_Fitted_LEV_overlay = function() {
    
    #   The Empirical LEV is graphed as follows:
    
    graphics::plot(actuar::elev(data), type="l", xlab = "Data", ylab = "Limited Expected Value", main="Graphical Comparison of the eLEV with the Fitted LEVs", lwd = Line_Width, las = 1)
    
    
    #   We will now overlay the fitted distributions' LEVs.
    
    Gamma_Fitted_LEV = base::sapply(X = Sequence_for_Graphs, LEV_Gamma_Fit)
    graphics::lines(x = Sequence_for_Graphs, y = Gamma_Fitted_LEV, lwd = Line_Width, col = Collection_of_Colours[1])
    
    
    #   Allowing the legend to be placed outside the plot area:
    graphics::par(xpd = TRUE)
    
    
    # Call locator() in the console and click on the plot to choose the legend position
    legend_position = graphics::locator(1)
    
    # Use the returned coordinates to place the legend
    graphics::legend(x = legend_position$x, y = legend_position$y, base::c(Empirical_Distribution.text, Collection_of_Fitted_Distributions), fill = base::c("black", Collection_of_Colours))
    
    
    # Revoking the allowance of the legend to be placed outside the plot area:
    graphics::par(xpd = FALSE)
    
  }
  
  
  
  
  
  
  
  
  #   Graphing the Empirical MRL with the Gamma Fitted MRL overlay:
  
  graph_of_eMRL_with_Gamma_Fitted_MRL_overlay = function() {
    
    #   The Empirical MRL is graphed as follows:
    
    graphics::plot(data_increasing_order, data_eMRL, type="l", xlab="Data", ylab="e_n", main="Graphical Comparison of the eMRL with the Fitted MRLs", lwd = Line_Width)
    
    
    #   We will now overlay the fitted distributions' MRLs.
    
    Gamma_Fitted_MRL = base::sapply(X = Sequence_for_Graphs, MRL_Gamma_Fit)
    graphics::lines(x = Sequence_for_Graphs, y = Gamma_Fitted_MRL, lwd = Line_Width, col = Collection_of_Colours[1])
    
    
    #   Allowing the legend to be placed outside the plot area:
    graphics::par(xpd = TRUE)
    
    
    # Call locator() in the console and click on the plot to choose the legend position
    legend_position = graphics::locator(1)
    
    # Use the returned coordinates to place the legend
    graphics::legend(x = legend_position$x, y = legend_position$y, base::c(Empirical_Distribution.text, Collection_of_Fitted_Distributions), fill = base::c("black", Collection_of_Colours))
    
    
    # Revoking the allowance of the legend to be placed outside the plot area:
    graphics::par(xpd = FALSE)
    
  }
  
  
  
  
  
  
  
  
  
  #   Graphing the Empirical Histogram with the Gamma Fitted PDF overlay:
  
  graph_of_histogram_with_Gamma_Fitted_PDF_overlay = function() {
    
    #   The histogram of the empirical data is graphed as follows:
    
    graphics::hist(data, prob = T, main="Empirical Histogram compared to the Fitted Empirical Densities", xlab="Data", col="green")
    
    #   We will now overlay the fitted distributions' PDFs altogether.
    
    graphics::lines(data_increasing_order, stats::dgamma(data_increasing_order, Gamma_Fit_MLE_Estimated_shape, Gamma_Fit_MLE_Estimated_rate), xpd = T, lwd = Line_Width, col = Collection_of_Colours[1])
    
    
    #   Allowing the legend to be placed outside the plot area:
    graphics::par(xpd = TRUE)
    
    
    # Call locator() in the console and click on the plot to choose the legend position
    legend_position = graphics::locator(1)
    
    # Use the returned coordinates to place the legend
    graphics::legend(x = legend_position$x, y = legend_position$y, base::c(Empirical_Distribution.text, Collection_of_Fitted_Distributions), fill = base::c("black", Collection_of_Colours))
    
    
    # Revoking the allowance of the legend to be placed outside the plot area:
    graphics::par(xpd = FALSE)
    
  }
  
  
  
  
  
  
  
  #   Graphing the P-P Plot with the X = Y line overlay:
  
  graph_of_P.P_Plot_with_X.Y_line_overlay = function() {
    
    #     P-P Plots
    
    Legend_Content_for_P.P_Plots = base::c("(Fn, F_Fitted)","X = Y")
    
    Colour_of_X.Y_Line_for_P.P_Plots = "blue"
    Legend_Content_Fill_for_P.P_Plots = base::c("black", Colour_of_X.Y_Line_for_P.P_Plots)
    
    
    
    #   Gamma P-P Plot:
    
    ppoints_vector = stats::ppoints(sample_size)
    pp_plot = stats::qgamma(ppoints_vector, Gamma_Fit_MLE_Estimated_shape, Gamma_Fit_MLE_Estimated_rate)
    graphics::plot(data_increasing_order, pp_plot, main = "Gamma P-P Plot for Data", xlab = "Empirical CDF", ylab = "Fitted Gamma CDF")
    graphics::abline(0, 1, lwd = Line_Width, col = Colour_of_X.Y_Line_for_P.P_Plots)
    
    
    #   Allowing the legend to be placed outside the plot area:
    graphics::par(xpd = TRUE)
    
    
    # Call locator() in the console and click on the plot to choose the legend position
    legend_position = graphics::locator(1)
    
    # Use the returned coordinates to place the legend
    graphics::legend(x = legend_position$x, y = legend_position$y, Legend_Content_for_P.P_Plots, fill = Legend_Content_Fill_for_P.P_Plots)
    
    
    # Revoking the allowance of the legend to be placed outside the plot area:
    graphics::par(xpd = FALSE)
    
  }
  
  
  
  
  
  
  
  #     Q-Q Plots
  
  
  
  
  #   The following Q-Q Plots compare the graphs with respect to a line of least squares, the 0-1 (y-intercept = 0, slope = 1) line, and the robust line (a line that is fit between the first and third quartiles of the data).
  
  
  Q.Q_Plot.L_S.function = function() {
    
    EnvStats::qqPlot(data_increasing_order, distribution = "gamma", param.list = base::list(shape = Gamma_Fit_MLE_Estimated_shape, scale = 1/Gamma_Fit_MLE_Estimated_rate), add.line = TRUE, qq.line.type = "least squares")
    
  }
  
  
  
  
  
  
  Q.Q_Plot.0_1.function = function() {
    
    EnvStats::qqPlot(data_increasing_order, distribution = "gamma", param.list = base::list(shape = Gamma_Fit_MLE_Estimated_shape, scale = 1/Gamma_Fit_MLE_Estimated_rate), add.line = TRUE, qq.line.type = "0-1")
    
  }
  
  
  
  
  
  
  Q.Q_Plot.Robust.function = function() {
    
    EnvStats::qqPlot(data_increasing_order, distribution = "gamma", param.list = base::list(shape = Gamma_Fit_MLE_Estimated_shape, scale = 1/Gamma_Fit_MLE_Estimated_rate), add.line = TRUE, qq.line.type = "robust")
    
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  # Define a function to calculate summary statistics
  Fitted_Distribution_Summary_Statistics = function(dist, par, quantiles, stats_fun) {
    q = base::sapply(quantiles, function(x) base::do.call(dist, base::c(base::list(p = x), par)))
    stats = base::do.call(stats_fun, par)
    base::c(quantiles = q, stats)
  }
  
  # Define parameters
  
  fitted_gamma_par = base::list(shape = Gamma_Fit_MLE_Estimated_shape, rate = Gamma_Fit_MLE_Estimated_rate)
  
  
  # Define functions to calculate additional statistics
  
  fitted_gamma_stats_fun = function(shape, rate) {
    base::list(
      mean = shape / rate,
      sd = base::sqrt(shape / rate^2),
      cv = base::sqrt(shape) / shape,
      skew = 2 / base::sqrt(shape),
      kurt = 3 + 6 / shape
    )
  }
  
  
  
  
  
  # Calculate summary statistics
  
  fitted_gamma_stats = Fitted_Distribution_Summary_Statistics(stats::qgamma, fitted_gamma_par, Required_Quantiles, fitted_gamma_stats_fun)
  
  
  
  
  
  
  
  
  # *           -----   Formal Goodness-of-Fit Tests  -----
  
  
  
  
  #   Conducting the Kolmogorov-Smirnov Test on each fitted distribution:
  
  KS_Test = stats::ks.test(data, stats::pgamma, Gamma_Fit_MLE_Estimated_shape, Gamma_Fit_MLE_Estimated_rate)
  
  
  
  
  #   Conducting the Anderson-Darling Test on each fitted distribution:
  
  AD_Test = ADGofTest::ad.test(data, stats::pgamma, Gamma_Fit_MLE_Estimated_shape, Gamma_Fit_MLE_Estimated_rate)
  
  
  
  
  
  #   Conducting the AIC / BIC Tests on each fitted distribution:
  
  #   Preparing to conduct the tests:
  
  
  Gamma_Fit_Preliminary_Summary = fitdistrplus::fitdist(data, "gamma")
  
  
  
  
  #   The test results are observable as follows:
  
  Gamma_Fit_Summary = base::summary(Gamma_Fit_Preliminary_Summary)
  
  
  
  
  
  
  
  
  # *           -----   Confidence Regions, revisited  -----
  
  
  #   Construction of a Confidence Region (CR) with respect to each candidate distribution's estimators is observed as follows:
  
  #   The chosen level of confidence is (value_of_alpha_for_CR)% with respect to the Chi-Square test.
  
  
  alpha_for_CR = value_of_alpha_for_CR
  Level_of_Confidence.1_T = 1-(alpha_for_CR)
  
  Z_Value_at_LoC.1_T = stats::qnorm(Level_of_Confidence.1_T)
  
  
  #   We now proceed.
  
  
  
  graph_of_confidence_region_for_Gamma_Fit_parameters = function() {
    
    
    #       Gamma:
    
    
    #   Generating the confidence ellipse:
    
    Gamma_Fit_MLE_Estimates_CR = ellipse::ellipse(x = Gamma_vcov_Matrix, center = Gamma_Fit_MLE_Estimates, shape = Gamma_vcov_Matrix, radius = Z_Value_at_LoC.1_T)
    
    #   Plotting the confidence ellipse:
    
    graphics::plot(Gamma_Fit_MLE_Estimates_CR, type = "l", xlab = "Shape", ylab = "Rate")
    graphics::polygon(Gamma_Fit_MLE_Estimates_CR, col = grDevices::rgb(0, 0, 1, 0.5))  # Adding the shaded region
    
    graphics::points(Gamma_Fit_MLE_Estimated_shape, Gamma_Fit_MLE_Estimated_rate, pch = 19)
    
    graphics::title(main = "CR for the ML Estimators of the Fitted Gamma Distribution")
    
    
  }
  
  
  
  #   Number of parameters in the Gamma distribution is observed as follows:
  
  Number_of_Parameters_Gamma_Fit_MLE_Estimates = base::length(stats::coef(Gamma_Fit_Summary))
  
  
  #   Quantile of the Chi-Square distribution:
  
  Quantile_of_Chi.Square_Dist_Gamma_Fit = stats::qchisq(1 - alpha_for_CR, df = Number_of_Parameters_Gamma_Fit_MLE_Estimates)
  
  # Calculating the threshold on the log-likelihood for the likelihood ratio test at the alpha = 5% significance level:
  
  Log.Likelihood_Threshold_Gamma_Fit = Gamma_Fit_Summary$loglik - 0.5 * Quantile_of_Chi.Square_Dist_Gamma_Fit
  
  
  
  
  
  Useful_List = base::list(
    preliminary_summary_statistics = base::summary(data)
    , sample_size = base::length(data)
    , minimum = base::min(data)
    , maximum = base::max(data)
    , min_and_max_of_data = base::range(data)
    , range_difference = maximum - minimum
    , mean_of_data = base::mean(data)
    , median_of_data = stats::median(data)
    , variance_of_data = actuar::var(data)
    , standard_deviation_of_data = actuar::sd(data)
    , coefficient_of_variation = standard_deviation_of_data/mean_of_data
    , Required_Quantiles = desired_quantiles
    , quantile_values = stats::quantile(data, Required_Quantiles)
    , Skewness_of_data = EnvStats::skewness(data)
    , Kurtosis_of_data = EnvStats::kurtosis(data, excess = FALSE)
    , Gamma_Fit_MME = fitdistrplus::fitdist(data, "gamma", method = "mme")
    , Gamma_Fit_MLE = fitdistrplus::fitdist(data, "gamma")
    , Gamma_Fit_MLE_Estimates = Gamma_Fit_MLE$estimate
    , Gamma_Fit_MLE_Estimated_shape = base::as.numeric(Gamma_Fit_MLE_Estimates[1])
    , Gamma_Fit_MLE_Estimated_rate = base::as.numeric(Gamma_Fit_MLE_Estimates[2])
    , Gamma_vcov_Matrix = Gamma_Fit_MLE$vcov
    , Gamma_Information_Matrix = matlib::Inverse(Gamma_vcov_Matrix)
    , alpha_for_CI = value_of_alpha_for_CI
    , Level_of_Confidence.2_T = 1-(alpha_for_CI/2)
    , Z_Value_at_LoC.2_T = stats::qnorm(Level_of_Confidence.2_T)
    , Lower_bound_of_CI_for_MLE_Gamma_shape = Gamma_Fit_MLE_Estimated_shape - Z_Value_at_LoC.2_T*base::sqrt(Gamma_vcov_Matrix[1,1])
    , Upper_bound_of_CI_for_MLE_Gamma_shape = Gamma_Fit_MLE_Estimated_shape + Z_Value_at_LoC.2_T*base::sqrt(Gamma_vcov_Matrix[1,1])
    , CI_of_MLE_Gamma_shape = base::as.numeric(
      base::c(Lower_bound_of_CI_for_MLE_Gamma_shape, Upper_bound_of_CI_for_MLE_Gamma_shape)
    )
    , Lower_bound_of_CI_for_MLE_Gamma_rate = Gamma_Fit_MLE_Estimated_rate - Z_Value_at_LoC.2_T*base::sqrt(Gamma_vcov_Matrix[2,2])
    , Upper_bound_of_CI_for_MLE_Gamma_rate = Gamma_Fit_MLE_Estimated_rate + Z_Value_at_LoC.2_T*base::sqrt(Gamma_vcov_Matrix[2,2])
    , CI_of_MLE_Gamma_rate = base::as.numeric(
      base::c(Lower_bound_of_CI_for_MLE_Gamma_rate, Upper_bound_of_CI_for_MLE_Gamma_rate)
    )
    , fitted_gamma_par = base::list(shape = Gamma_Fit_MLE_Estimated_shape, rate = Gamma_Fit_MLE_Estimated_rate)
    , fitted_gamma_stats = Fitted_Distribution_Summary_Statistics(stats::qgamma, fitted_gamma_par, Required_Quantiles, fitted_gamma_stats_fun)
    , KS_Test = stats::ks.test(data, stats::pgamma, Gamma_Fit_MLE_Estimated_shape, Gamma_Fit_MLE_Estimated_rate)
    , AD_Test = ADGofTest::ad.test(data, stats::pgamma, Gamma_Fit_MLE_Estimated_shape, Gamma_Fit_MLE_Estimated_rate)
    , Gamma_Fit_Preliminary_Summary = fitdistrplus::fitdist(data, "gamma")
    , Gamma_Fit_Summary = base::summary(Gamma_Fit_Preliminary_Summary)
    , alpha_for_CR = value_of_alpha_for_CR
    , Level_of_Confidence.1_T = 1-(alpha_for_CR)
    , Z_Value_at_LoC.1_T = stats::qnorm(Level_of_Confidence.1_T)
    , Number_of_Parameters_Gamma_Fit_MLE_Estimates = base::length(stats::coef(Gamma_Fit_Summary))
    , Quantile_of_Chi.Square_Dist_Gamma_Fit = stats::qchisq(1 - alpha_for_CR, df = Number_of_Parameters_Gamma_Fit_MLE_Estimates)
    , Log.Likelihood_Threshold_Gamma_Fit = Gamma_Fit_Summary$loglik - 0.5 * Quantile_of_Chi.Square_Dist_Gamma_Fit
    
    , if(pertinent_graphs == TRUE) {
      
      base::c(
        Graph_of_eCDF_with_Gamma_Fitted_CDF_Overlay = graph_of_eCDF_with_Gamma_Fitted_CDF_overlay()
        , Graph_of_eLEV_with_Gamma_Fitted_LEV_Overlay = graph_of_eLEV_with_Gamma_Fitted_LEV_overlay()
        , Graph_of_eMRL_with_Gamma_Fitted_MRL_Overlay = graph_of_eMRL_with_Gamma_Fitted_MRL_overlay()
        , Graph_of_Histogram_with_Gamma_Fitted_PDF_Overlay = graph_of_histogram_with_Gamma_Fitted_PDF_overlay()
        , Graph_of_P.P_Plot_with_X.Y_Line_Overlay = graph_of_P.P_Plot_with_X.Y_line_overlay()
        , Q.Q_Plot.0_1 = Q.Q_Plot.0_1.function()
        , Q.Q_Plot.L_S = Q.Q_Plot.L_S.function()
        , Q.Q_Plot.Robust = Q.Q_Plot.Robust.function()
        , Graph_of_Confidence_Region_for_Gamma_Fit_Parameters = graph_of_confidence_region_for_Gamma_Fit_parameters()
      )
      
      
    } else {}
    
    
  )
  
  
  Useful_List
  
  
}

#' Weibull Distribution Data Fit
#'
#' This function fits a Weibull distribution to the given data and provides various statistics and graphs.
#'
#' @param data A numeric vector of data to be fitted.
#' @param pertinent_graphs A logical value indicating whether to generate pertinent graphs.
#' @param desired_quantiles A numeric vector of quantiles to be computed.
#' @param value_of_alpha_for_CI A numeric value for the level of significance for confidence intervals.
#' @param value_of_alpha_for_CR A numeric value for the level of significance for confidence regions.
#'
#' @return A list containing various summary statistics, test results, and graphs (if pertinent_graphs = TRUE).
#'
#' @examples
#' \dontrun{
#' data <- rweibull(100)
#' Weibull_Distribution_Data_Fit(data, TRUE, c(0.1, 0.5, 0.9), 0.05, 0.05)
#' }
#'
#' @export
Weibull_Distribution_Data_Fit <- function(data, pertinent_graphs, desired_quantiles, value_of_alpha_for_CI, value_of_alpha_for_CR) {
  
  
  
  #       ----------------   Actuarial Coding Portfolio  ----------------
  
  #                         Ioannis (Yanni) Papadopoulos
  
  #       =================================================================
  
  
  #   Before we begin, we will load all the necessary libraries.
  
  
  
  
  #   We will now begin by loading the data.
  #   We will proceed with the following method to load the data for reproducibility reasons:
  
  
  
  
  #       ======    Part I:    Model Selection    ======
  
  
  # *           -----   Descriptive Statistics  -----
  
  
  #   Preliminary summary statistics:
  
  preliminary_summary_statistics = base::summary(data)
  
  
  #   Observing the sample size:
  
  sample_size = base::length(data)
  
  
  #   Observing the minimum, maximum, and range:
  
  minimum = base::min(data)
  
  maximum = base::max(data)
  
  
  #   The minimum and maximum of the data can also be observed using the following function:
  
  min_and_max_of_data = base::range(data)
  
  #   To find the range, we will take the difference of the maximum and the minimum:
  
  range_difference = maximum - minimum
  
  
  
  #   Observing the mean of the data:
  
  mean_of_data = base::mean(data)
  
  
  
  #   Observing the median of the data:
  
  median_of_data = stats::median(data)
  
  
  
  #   Observing the variance and standard deviation of the data:
  
  variance_of_data = actuar::var(data)
  
  
  standard_deviation_of_data = actuar::sd(data)
  
  
  
  #   Observing the coefficient of variation for the data:
  
  coefficient_of_variation = standard_deviation_of_data/mean_of_data
  
  
  
  #   Calculation of Quantiles:
  
  #       The following are the required quantiles:
  
  Required_Quantiles = desired_quantiles
  
  
  #      We now observe the required quantile values:
  
  quantile_values = stats::quantile(data, Required_Quantiles)
  
  
  #   Observation of skewness and kurtosis.
  #     For this, the 'moments' package is used.
  
  Skewness_of_data = EnvStats::skewness(data)
  
  
  Kurtosis_of_data = EnvStats::kurtosis(data, excess = FALSE)
  
  
  
  
  
  # *           -----   Supporting Graphs  -----
  
  
  #   In this section, the required exploratory graphs are produced.
  
  #   We will also sort the data in increasing order for when such an operation would be appropriate.
  
  
  data_increasing_order = base::sort(data)
  
  
  #   Boxplot:
  
  
  graph_of_boxplot = function(data) {
    
    graphics::boxplot(data)
    graphics::title(main = "Boxplot", ylab = "Data")
    
  }
  
  
  
  #   Histogram:
  
  
  graph_of_histogram = function(data) {
    
    graphics::hist(data, main="Histogram for Data", xlab="Data", col="green")
    
  }
  
  
  
  
  #   Empirical CDF:
  
  eCDF_of_data = stats::ecdf(data)
  
  graph_of_eCDF = function(data) {
    
    #    eCDF_of_data = ecdf(data)
    graphics::plot(eCDF_of_data, main="Empirical CDF for Data", xlab="Data", ylab="F_n")
    
  }
  
  
  
  
  #   In preparation of the Empirical MRL (Mean Residual Life) function, a CDF function is defined.
  
  CDF = function(x){
    
    a = x
    
    for(i in 1:base::length(x) ) {
      
      a[i] = base::sum(x<=x[i]) / base::length(x)
      
    }
    return(a)
  }
  
  
  #   The Empirical MRL function is as follows:
  
  eMRL = function(x){
    
    a = x
    
    x.Fn = CDF(x)
    
    for(i in 1:base::length(x)){
      
      if(x[i] < base::max(x)){
        
        a[i] = base::sum(x[x - x[i] > 0] - x[i]) / (base::length(x)) / (1 - x.Fn[i])
        
      } else {
        
        a[i] = 0
        
      }
      
    }
    
    return(a)
    
  }
  
  
  
  
  #   The graph of the Empirical MRL is as follows:
  
  data_eMRL = eMRL(data_increasing_order)
  
  graph_of_eMRL = function(data) {
    
    #    data_eMRL = eMRL(data_increasing_order)
    graphics::plot(data_increasing_order, data_eMRL, type="l", xlab="Data", ylab="e_n", main="Empirical MRL")
    
  }
  
  
  
  
  
  #   The graph of the Empirical LEV is as follows:
  
  data_eLEV = base::mean(data_increasing_order) - data_eMRL * (1 - eCDF_of_data(data_increasing_order))
  
  graph_of_eLEV = function(data) {
    
    #    data_eLEV = mean(data_increasing_order) - data_eMRL * (1 - eCDF_of_data(data_increasing_order))
    graphics::plot(data_increasing_order, data_eLEV, type="l", xlab="Data", ylab="Limited Expected Value", main="Empirical LEV", las = 1)
    
  }
  
  
  
  
  
  #   To verify the validity of the assumptions, an IID test is performed:
  
  
  graph_of_IID_Test = function(data) {
    
    IID_Test = stats::acf(x = data, lag.max = sample_size, type = "correlation", main="ACF of Data")
    
  }
  
  
  
  
  
  #       ======    Part II:    Estimation    ======
  
  
  
  # *           -----   Method of Moments (MME) Parameter Estimation  -----
  
  
  #   We will proceed with fitting the model distribution candidates with respect to the Method of Moments (MME) approach:
  
  
  Weibull_Fit_MME = EnvStats::eweibull(data, method = "mme")
  
  
  
  
  
  # *           -----   Maximum Likelihood Estimator (MLE) Parameter Estimation  -----
  
  
  #   Now, we will proceed with fitting the model distribution candidates with respect to the Maximum Likelihood Estimator (MLE) approach:
  
  
  Weibull_Fit_MLE = fitdistrplus::fitdist(data, "weibull")
  
  
  
  
  
  
  # *           -----   Estimation of Asymptotic Covariance Matrices  -----
  
  
  #   The Information Matrices with respect to each candidate distribution are observed as follows:
  
  
  #     Weibull:
  
  Weibull_Fit_MLE_Estimates = Weibull_Fit_MLE$estimate
  
  Weibull_Fit_MLE_Estimated_shape = base::as.numeric(Weibull_Fit_MLE_Estimates[1])
  Weibull_Fit_MLE_Estimated_scale = base::as.numeric(Weibull_Fit_MLE_Estimates[2])
  
  Weibull_vcov_Matrix = Weibull_Fit_MLE$vcov
  
  Weibull_Information_Matrix = matlib::Inverse(Weibull_vcov_Matrix)
  
  
  
  
  
  
  
  
  
  
  # *           -----   Confidence Regions  -----
  
  
  #   Note:   The Confidence Regions (CR) will be constructed after the Formal Goodness of Fit test have been conducted, which is at the end of the script.
  #               This section will be denoted as "Confidence Regions, revisited"
  
  #   Construction of a Confidence Interval (CI) with respect to each candidate distribution's estimators is observed as follows:
  
  #   The chosen level of confidence is (value_of_alpha_for_CI)% with respect to a two-tailed test.
  
  alpha_for_CI = value_of_alpha_for_CI
  Level_of_Confidence.2_T = 1-(alpha_for_CI/2)
  
  Z_Value_at_LoC.2_T = stats::qnorm(Level_of_Confidence.2_T)
  
  
  #   We now proceed.
  
  
  #     Weibull:
  
  Lower_bound_of_CI_for_MLE_Weibull_shape = Weibull_Fit_MLE_Estimated_shape - Z_Value_at_LoC.2_T*base::sqrt(Weibull_vcov_Matrix[1,1])
  Upper_bound_of_CI_for_MLE_Weibull_shape = Weibull_Fit_MLE_Estimated_shape + Z_Value_at_LoC.2_T*base::sqrt(Weibull_vcov_Matrix[1,1])
  
  
  CI_of_MLE_Weibull_shape = base::as.numeric(
    base::c(Lower_bound_of_CI_for_MLE_Weibull_shape, Upper_bound_of_CI_for_MLE_Weibull_shape)
  )
  
  
  
  
  Lower_bound_of_CI_for_MLE_Weibull_scale = Weibull_Fit_MLE_Estimated_scale - Z_Value_at_LoC.2_T*base::sqrt(Weibull_vcov_Matrix[2,2])
  Upper_bound_of_CI_for_MLE_Weibull_scale = Weibull_Fit_MLE_Estimated_scale + Z_Value_at_LoC.2_T*base::sqrt(Weibull_vcov_Matrix[2,2])
  
  
  CI_of_MLE_Weibull_scale = base::as.numeric(
    base::c(Lower_bound_of_CI_for_MLE_Weibull_scale, Upper_bound_of_CI_for_MLE_Weibull_scale)
  )
  
  
  
  
  
  
  #       ======    Part III:    Goodness-of-fit Tests    ======
  
  
  
  #     Creating the required functions for each distribution before proceeding.
  
  
  #      Weibull
  
  #   CDF:
  
  CDF_Weibull_Fit = function(x){
    
    stats::pweibull(x, shape = Weibull_Fit_MLE_Estimated_shape, scale = Weibull_Fit_MLE_Estimated_scale)
    
  }
  
  #   Survival Function:
  
  Survival_Function_Weibull_Fit = function(x){
    
    1 - CDF_Weibull_Fit(x)
    
  }
  
  
  
  #   Probability Density Function:
  
  PDF_Weibull_Fit = function(x){
    
    stats::dweibull(x, shape = Weibull_Fit_MLE_Estimated_shape, scale = Weibull_Fit_MLE_Estimated_scale)
    
  }
  
  
  
  #   Limited Expected Value:
  
  LEV_Weibull_Fit = function(t){
    
    stats::integrate(f = Survival_Function_Weibull_Fit, lower = 0, upper = t)$value
    
  }
  
  
  
  
  #   Mean Residual Life:
  
  MRL_Weibull_Fit = function(t){
    
    ( LEV_Weibull_Fit(Inf) - LEV_Weibull_Fit(t) ) / Survival_Function_Weibull_Fit(t)
    
  }
  
  
  
  
  
  
  # *           -----   Graphical Comparisons  -----
  
  
  #   Now, we will create a graph comparing the Empirical CDF with the Fitted Distributions' CDFs.
  
  #     Firstly, we will create a sequence to aid in the creation of the graph.
  
  Sequence_for_Graphs = base::seq(minimum, maximum, length = sample_size)
  
  #     Also, we will introduce a multiplier for "par("cex")" as well as a line width parameter for flexibility:
  
  Multiplier_for_par.cex = 3
  
  Line_Width = Multiplier_for_par.cex*graphics::par("cex")
  
  
  Empirical_Distribution.text = "Empirical Distribution"
  
  
  Collection_of_Fitted_Distributions = base::c(
    #     "Fitted Normal Distribution"
    #    ,"Fitted Gamma Distribution"
    "Fitted Weibull Distribution"
    #    ,"Fitted Inverse Gamma Distribution"
    #    ,"Fitted Inverse Weibull Distribution"
  )
  
  
  Collection_of_Colours = base::c(
    #     "blue"
    #    ,"orange"
    "purple"
    #    ,"red"
    #    ,"turquoise"
  )
  
  
  
  #   Graphing the Empirical CDF with the Weibull Fitted CDF overlay:
  
  
  graph_of_eCDF_with_Weibull_Fitted_CDF_overlay = function() {
    
    
    #   The Empirical CDF is graphed as follows:
    
    EnvStats::ecdfPlot(data_increasing_order, main = "Graphical Comparison of the eCDF with the Fitted CDFs", xlab = "Data")
    
    
    #   We will now overlay the fitted distributions' CDFs.
    
    Weibull_Fitted_CDF = base::sapply(X = Sequence_for_Graphs, CDF_Weibull_Fit)
    graphics::lines(x = Sequence_for_Graphs, y = Weibull_Fitted_CDF, lwd = Line_Width, col = Collection_of_Colours[1])
    
    
    #   Allowing the legend to be placed outside the plot area:
    graphics::par(xpd = TRUE)
    
    
    # Call locator() in the console and click on the plot to choose the legend position
    legend_position = graphics::locator(1)
    
    # Use the returned coordinates to place the legend
    graphics::legend(x = legend_position$x, y = legend_position$y, base::c(Empirical_Distribution.text, Collection_of_Fitted_Distributions), fill = base::c("black", Collection_of_Colours))
    
    
    # Revoking the allowance of the legend to be placed outside the plot area:
    graphics::par(xpd = FALSE)
    
  }
  
  
  
  
  
  
  
  
  #   Graphing the Empirical LEV with the Weibull Fitted LEV overlay:
  
  graph_of_eLEV_with_Weibull_Fitted_LEV_overlay = function() {
    
    #   The Empirical LEV is graphed as follows:
    
    graphics::plot(actuar::elev(data), type="l", xlab = "Data", ylab = "Limited Expected Value", main="Graphical Comparison of the eLEV with the Fitted LEVs", lwd = Line_Width, las = 1)
    
    
    #   We will now overlay the fitted distributions' LEVs.
    
    Weibull_Fitted_LEV = base::sapply(X = Sequence_for_Graphs, LEV_Weibull_Fit)
    graphics::lines(x = Sequence_for_Graphs, y = Weibull_Fitted_LEV, lwd = Line_Width, col = Collection_of_Colours[1])
    
    
    #   Allowing the legend to be placed outside the plot area:
    graphics::par(xpd = TRUE)
    
    
    # Call locator() in the console and click on the plot to choose the legend position
    legend_position = graphics::locator(1)
    
    # Use the returned coordinates to place the legend
    graphics::legend(x = legend_position$x, y = legend_position$y, base::c(Empirical_Distribution.text, Collection_of_Fitted_Distributions), fill = base::c("black", Collection_of_Colours))
    
    
    # Revoking the allowance of the legend to be placed outside the plot area:
    graphics::par(xpd = FALSE)
    
  }
  
  
  
  
  
  
  
  
  #   Graphing the Empirical MRL with the Weibull Fitted MRL overlay:
  
  graph_of_eMRL_with_Weibull_Fitted_MRL_overlay = function() {
    
    #   The Empirical MRL is graphed as follows:
    
    graphics::plot(data_increasing_order, data_eMRL, type="l", xlab="Data", ylab="e_n", main="Graphical Comparison of the eMRL with the Fitted MRLs", lwd = Line_Width)
    
    
    #   We will now overlay the fitted distributions' MRLs.
    
    Weibull_Fitted_MRL = base::sapply(X = Sequence_for_Graphs, MRL_Weibull_Fit)
    graphics::lines(x = Sequence_for_Graphs, y = Weibull_Fitted_MRL, lwd = Line_Width, col = Collection_of_Colours[1])
    
    
    #   Allowing the legend to be placed outside the plot area:
    graphics::par(xpd = TRUE)
    
    
    # Call locator() in the console and click on the plot to choose the legend position
    legend_position = graphics::locator(1)
    
    # Use the returned coordinates to place the legend
    graphics::legend(x = legend_position$x, y = legend_position$y, base::c(Empirical_Distribution.text, Collection_of_Fitted_Distributions), fill = base::c("black", Collection_of_Colours))
    
    
    # Revoking the allowance of the legend to be placed outside the plot area:
    graphics::par(xpd = FALSE)
    
  }
  
  
  
  
  
  
  
  
  
  #   Graphing the Empirical Histogram with the Weibull Fitted PDF overlay:
  
  graph_of_histogram_with_Weibull_Fitted_PDF_overlay = function() {
    
    #   The histogram of the empirical data is graphed as follows:
    
    graphics::hist(data, prob = T, main="Empirical Histogram compared to the Fitted Empirical Densities", xlab="Data", col="green")
    
    #   We will now overlay the fitted distributions' PDFs altogether.
    
    graphics::lines(data_increasing_order, stats::dweibull(data_increasing_order, Weibull_Fit_MLE_Estimated_shape, Weibull_Fit_MLE_Estimated_scale), xpd = T, lwd = Line_Width, col = Collection_of_Colours[1])
    
    
    #   Allowing the legend to be placed outside the plot area:
    graphics::par(xpd = TRUE)
    
    
    # Call locator() in the console and click on the plot to choose the legend position
    legend_position = graphics::locator(1)
    
    # Use the returned coordinates to place the legend
    graphics::legend(x = legend_position$x, y = legend_position$y, base::c(Empirical_Distribution.text, Collection_of_Fitted_Distributions), fill = base::c("black", Collection_of_Colours))
    
    
    # Revoking the allowance of the legend to be placed outside the plot area:
    graphics::par(xpd = FALSE)
    
  }
  
  
  
  
  
  
  
  #   Graphing the P-P Plot with the X = Y line overlay:
  
  graph_of_P.P_Plot_with_X.Y_line_overlay = function() {
    
    #     P-P Plots
    
    Legend_Content_for_P.P_Plots = base::c("(Fn, F_Fitted)","X = Y")
    
    Colour_of_X.Y_Line_for_P.P_Plots = "blue"
    Legend_Content_Fill_for_P.P_Plots = base::c("black", Colour_of_X.Y_Line_for_P.P_Plots)
    
    
    
    #   Weibull P-P Plot:
    
    ppoints_vector = stats::ppoints(sample_size)
    pp_plot = stats::qweibull(ppoints_vector, Weibull_Fit_MLE_Estimated_shape, Weibull_Fit_MLE_Estimated_scale)
    graphics::plot(data_increasing_order, pp_plot, main = "Weibull P-P Plot for Data", xlab = "Empirical CDF", ylab = "Fitted Weibull CDF")
    graphics::abline(0, 1, lwd = Line_Width, col = Colour_of_X.Y_Line_for_P.P_Plots)
    
    
    #   Allowing the legend to be placed outside the plot area:
    graphics::par(xpd = TRUE)
    
    
    # Call locator() in the console and click on the plot to choose the legend position
    legend_position = graphics::locator(1)
    
    # Use the returned coordinates to place the legend
    graphics::legend(x = legend_position$x, y = legend_position$y, Legend_Content_for_P.P_Plots, fill = Legend_Content_Fill_for_P.P_Plots)
    
    
    # Revoking the allowance of the legend to be placed outside the plot area:
    graphics::par(xpd = FALSE)
    
  }
  
  
  
  
  
  
  
  #     Q-Q Plots
  
  
  
  
  #   The following Q-Q Plots compare the graphs with respect to a line of least squares, the 0-1 (y-intercept = 0, slope = 1) line, and the robust line (a line that is fit between the first and third quartiles of the data).
  
  
  Q.Q_Plot.L_S.function = function() {
    
    EnvStats::qqPlot(data_increasing_order, distribution = "weibull", param.list = base::list(shape = Weibull_Fit_MLE_Estimated_shape, scale = Weibull_Fit_MLE_Estimated_scale), add.line = TRUE, qq.line.type = "least squares")
    
  }
  
  
  
  
  
  
  Q.Q_Plot.0_1.function = function() {
    
    EnvStats::qqPlot(data_increasing_order, distribution = "weibull", param.list = base::list(shape = Weibull_Fit_MLE_Estimated_shape, scale = Weibull_Fit_MLE_Estimated_scale), add.line = TRUE, qq.line.type = "0-1")
    
  }
  
  
  
  
  
  
  Q.Q_Plot.Robust.function = function() {
    
    EnvStats::qqPlot(data_increasing_order, distribution = "weibull", param.list = base::list(shape = Weibull_Fit_MLE_Estimated_shape, scale = Weibull_Fit_MLE_Estimated_scale), add.line = TRUE, qq.line.type = "robust")
    
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  # Define a function to calculate summary statistics
  Fitted_Distribution_Summary_Statistics = function(dist, par, quantiles, stats_fun) {
    q = base::sapply(quantiles, function(x) base::do.call(dist, base::c(base::list(p = x), par)))
    stats = base::do.call(stats_fun, par)
    base::c(quantiles = q, stats)
  }
  
  # Define parameters
  
  fitted_weibull_par = base::list(shape = Weibull_Fit_MLE_Estimated_shape, scale = Weibull_Fit_MLE_Estimated_scale)
  
  
  # Define functions to calculate additional statistics
  
  fitted_weibull_stats_fun = function(shape, scale) {
    list(
      mean = (scale * base::gamma(1 + 1 / shape)),
      sd = (base::sqrt((scale^2) * (base::gamma(1 + 2 / shape) - (base::gamma(1 + 1 / shape))^2))),
      cv = base::sqrt(base::gamma(1 + 2 / shape) - (base::gamma(1 + 1 / shape))^2) / base::gamma(1 + 1 / shape),
      skew = (((base::gamma(1 + 3 / shape) * scale^3) - (3 * (scale * base::gamma(1 + 1 / shape)) * (base::sqrt((scale^2) * (base::gamma(1 + 2 / shape) - (base::gamma(1 + 1 / shape))^2)))^2) - ((scale * base::gamma(1 + 1 / shape))^3)) / ((base::sqrt((scale^2) * (base::gamma(1 + 2 / shape) - (base::gamma(1 + 1 / shape))^2)))^3)),
      kurt = (
        (
          ((base::gamma(1 + 4 / shape) * scale^4) - (4 * (((base::gamma(1 + 3 / shape) * scale^3) - (3 * (scale * base::gamma(1 + 1 / shape)) * (base::sqrt((scale^2) * (base::gamma(1 + 2 / shape) - (base::gamma(1 + 1 / shape))^2)))^2) - ((scale * base::gamma(1 + 1 / shape))^3)) / ((base::sqrt((scale^2) * (base::gamma(1 + 2 / shape) - (base::gamma(1 + 1 / shape))^2)))^3)) * (base::sqrt((scale^2) * (base::gamma(1 + 2 / shape) - (base::gamma(1 + 1 / shape))^2)))^3 * (scale * base::gamma(1 + 1 / shape))) - (6 * (scale * base::gamma(1 + 1 / shape))^2 * (base::sqrt((scale^2) * (base::gamma(1 + 2 / shape) - (base::gamma(1 + 1 / shape))^2)))^2) - (scale * base::gamma(1 + 1 / shape))^4) / ((base::sqrt((scale^2) * (base::gamma(1 + 2 / shape) - (base::gamma(1 + 1 / shape))^2)))^4))
      )
    )
  }
  
  
  
  
  
  # Calculate summary statistics
  
  fitted_weibull_stats = Fitted_Distribution_Summary_Statistics(stats::qweibull, fitted_weibull_par, Required_Quantiles, fitted_weibull_stats_fun)
  
  
  
  
  
  
  
  
  # *           -----   Formal Goodness-of-Fit Tests  -----
  
  
  
  
  #   Conducting the Kolmogorov-Smirnov Test on each fitted distribution:
  
  KS_Test = stats::ks.test(data, stats::pweibull, Weibull_Fit_MLE_Estimated_shape, Weibull_Fit_MLE_Estimated_scale)
  
  
  
  
  #   Conducting the Anderson-Darling Test on each fitted distribution:
  
  AD_Test = ADGofTest::ad.test(data, stats::pweibull, Weibull_Fit_MLE_Estimated_shape, Weibull_Fit_MLE_Estimated_scale)
  
  
  
  
  
  #   Conducting the AIC / BIC Tests on each fitted distribution:
  
  #   Preparing to conduct the tests:
  
  
  Weibull_Fit_Preliminary_Summary = fitdistrplus::fitdist(data, "weibull")
  
  
  
  
  #   The test results are observable as follows:
  
  Weibull_Fit_Summary = base::summary(Weibull_Fit_Preliminary_Summary)
  
  
  
  
  
  
  
  
  # *           -----   Confidence Regions, revisited  -----
  
  
  #   Construction of a Confidence Region (CR) with respect to each candidate distribution's estimators is observed as follows:
  
  #   The chosen level of confidence is (value_of_alpha_for_CR)% with respect to the Chi-Square test.
  
  
  alpha_for_CR = value_of_alpha_for_CR
  Level_of_Confidence.1_T = 1-(alpha_for_CR)
  
  Z_Value_at_LoC.1_T = stats::qnorm(Level_of_Confidence.1_T)
  
  
  #   We now proceed.
  
  
  
  graph_of_confidence_region_for_Weibull_Fit_parameters = function() {
    
    
    #       Weibull:
    
    
    #   Generating the confidence ellipse:
    
    Weibull_Fit_MLE_Estimates_CR = ellipse::ellipse(x = Weibull_vcov_Matrix, center = Weibull_Fit_MLE_Estimates, shape = Weibull_vcov_Matrix, radius = Z_Value_at_LoC.1_T)
    
    #   Plotting the confidence ellipse:
    
    graphics::plot(Weibull_Fit_MLE_Estimates_CR, type = "l", xlab = "Shape", ylab = "Scale")
    graphics::polygon(Weibull_Fit_MLE_Estimates_CR, col = grDevices::rgb(0, 0, 1, 0.5))  # Adding the shaded region
    
    graphics::points(Weibull_Fit_MLE_Estimated_shape, Weibull_Fit_MLE_Estimated_scale, pch = 19)
    
    graphics::title(main = "CR for the ML Estimators of the Fitted Weibull Distribution")
    
    
  }
  
  
  
  #   Number of parameters in the Weibull distribution is observed as follows:
  
  Number_of_Parameters_Weibull_Fit_MLE_Estimates = base::length(stats::coef(Weibull_Fit_Summary))
  
  
  #   Quantile of the Chi-Square distribution:
  
  Quantile_of_Chi.Square_Dist_Weibull_Fit = stats::qchisq(1 - alpha_for_CR, df = Number_of_Parameters_Weibull_Fit_MLE_Estimates)
  
  # Calculating the threshold on the log-likelihood for the likelihood ratio test at the alpha = 5% significance level:
  
  Log.Likelihood_Threshold_Weibull_Fit = Weibull_Fit_Summary$loglik - 0.5 * Quantile_of_Chi.Square_Dist_Weibull_Fit
  
  
  
  
  
  
  Useful_List = base::list(
    preliminary_summary_statistics = base::summary(data)
    , sample_size = base::length(data)
    , minimum = base::min(data)
    , maximum = base::max(data)
    , min_and_max_of_data = base::range(data)
    , range_difference = maximum - minimum
    , mean_of_data = base::mean(data)
    , median_of_data = stats::median(data)
    , variance_of_data = actuar::var(data)
    , standard_deviation_of_data = actuar::sd(data)
    , coefficient_of_variation = standard_deviation_of_data/mean_of_data
    , Required_Quantiles = desired_quantiles
    , quantile_values = stats::quantile(data, Required_Quantiles)
    , Skewness_of_data = EnvStats::skewness(data)
    , Kurtosis_of_data = EnvStats::kurtosis(data, excess = FALSE)
    , Weibull_Fit_MME = EnvStats::eweibull(data, method = "mme")
    , Weibull_Fit_MLE = fitdistrplus::fitdist(data, "weibull")
    , Weibull_Fit_MLE_Estimates = Weibull_Fit_MLE$estimate
    , Weibull_Fit_MLE_Estimated_shape = base::as.numeric(Weibull_Fit_MLE_Estimates[1])
    , Weibull_Fit_MLE_Estimated_scale = base::as.numeric(Weibull_Fit_MLE_Estimates[2])
    , Weibull_vcov_Matrix = Weibull_Fit_MLE$vcov
    , Weibull_Information_Matrix = matlib::Inverse(Weibull_vcov_Matrix)
    , alpha_for_CI = value_of_alpha_for_CI
    , Level_of_Confidence.2_T = 1-(alpha_for_CI/2)
    , Z_Value_at_LoC.2_T = stats::qnorm(Level_of_Confidence.2_T)
    , Lower_bound_of_CI_for_MLE_Weibull_shape = Weibull_Fit_MLE_Estimated_shape - Z_Value_at_LoC.2_T*base::sqrt(Weibull_vcov_Matrix[1,1])
    , Upper_bound_of_CI_for_MLE_Weibull_shape = Weibull_Fit_MLE_Estimated_shape + Z_Value_at_LoC.2_T*base::sqrt(Weibull_vcov_Matrix[1,1])
    , CI_of_MLE_Weibull_shape = base::as.numeric(
      base::c(Lower_bound_of_CI_for_MLE_Weibull_shape, Upper_bound_of_CI_for_MLE_Weibull_shape)
    )
    , Lower_bound_of_CI_for_MLE_Weibull_scale = Weibull_Fit_MLE_Estimated_scale - Z_Value_at_LoC.2_T*base::sqrt(Weibull_vcov_Matrix[2,2])
    , Upper_bound_of_CI_for_MLE_Weibull_scale = Weibull_Fit_MLE_Estimated_scale + Z_Value_at_LoC.2_T*base::sqrt(Weibull_vcov_Matrix[2,2])
    , CI_of_MLE_Weibull_scale = base::as.numeric(
      base::c(Lower_bound_of_CI_for_MLE_Weibull_scale, Upper_bound_of_CI_for_MLE_Weibull_scale)
    )
    , fitted_weibull_par = base::list(shape = Weibull_Fit_MLE_Estimated_shape, scale = Weibull_Fit_MLE_Estimated_scale)
    , fitted_weibull_stats = Fitted_Distribution_Summary_Statistics(stats::qweibull, fitted_weibull_par, Required_Quantiles, fitted_weibull_stats_fun)
    , KS_Test = stats::ks.test(data, stats::pweibull, Weibull_Fit_MLE_Estimated_shape, Weibull_Fit_MLE_Estimated_scale)
    , AD_Test = ADGofTest::ad.test(data, stats::pweibull, Weibull_Fit_MLE_Estimated_shape, Weibull_Fit_MLE_Estimated_scale)
    , Weibull_Fit_Preliminary_Summary = fitdistrplus::fitdist(data, "weibull")
    , Weibull_Fit_Summary = base::summary(Weibull_Fit_Preliminary_Summary)
    , alpha_for_CR = value_of_alpha_for_CR
    , Level_of_Confidence.1_T = 1-(alpha_for_CR)
    , Z_Value_at_LoC.1_T = stats::qnorm(Level_of_Confidence.1_T)
    , Number_of_Parameters_Weibull_Fit_MLE_Estimates = base::length(stats::coef(Weibull_Fit_Summary))
    , Quantile_of_Chi.Square_Dist_Weibull_Fit = stats::qchisq(1 - alpha_for_CR, df = Number_of_Parameters_Weibull_Fit_MLE_Estimates)
    , Log.Likelihood_Threshold_Weibull_Fit = Weibull_Fit_Summary$loglik - 0.5 * Quantile_of_Chi.Square_Dist_Weibull_Fit
    
    , if(pertinent_graphs == TRUE) {
      
      base::c(
        Graph_of_eCDF_with_Weibull_Fitted_CDF_Overlay = graph_of_eCDF_with_Weibull_Fitted_CDF_overlay()
        , Graph_of_eLEV_with_Weibull_Fitted_LEV_Overlay = graph_of_eLEV_with_Weibull_Fitted_LEV_overlay()
        , Graph_of_eMRL_with_Weibull_Fitted_MRL_Overlay = graph_of_eMRL_with_Weibull_Fitted_MRL_overlay()
        , Graph_of_Histogram_with_Weibull_Fitted_PDF_Overlay = graph_of_histogram_with_Weibull_Fitted_PDF_overlay()
        , Graph_of_P.P_Plot_with_X.Y_Line_Overlay = graph_of_P.P_Plot_with_X.Y_line_overlay()
        , Q.Q_Plot.0_1 = Q.Q_Plot.0_1.function()
        , Q.Q_Plot.L_S = Q.Q_Plot.L_S.function()
        , Q.Q_Plot.Robust = Q.Q_Plot.Robust.function()
        , Graph_of_Confidence_Region_for_Weibull_Fit_Parameters = graph_of_confidence_region_for_Weibull_Fit_parameters()
      )
      
      
    } else {}
    
    
  )
  
  
  Useful_List
  
  
}

#' Inverse Gamma Distribution Data Fit
#'
#' This function fits an inverse gamma distribution to the given data and provides various statistics and graphs.
#'
#' @param data A numeric vector of data to be fitted.
#' @param pertinent_graphs A logical value indicating whether to generate pertinent graphs.
#' @param desired_quantiles A numeric vector of quantiles to be computed.
#' @param value_of_alpha_for_CI A numeric value for the level of significance for confidence intervals.
#' @param value_of_alpha_for_CR A numeric value for the level of significance for confidence regions.
#'
#' @return A list containing various summary statistics, test results, and graphs (if pertinent_graphs = TRUE).
#'
#' @examples
#' \dontrun{
#' install.packages("actuar")    # See below.
#' library(actuar)               # This is done to ensure the required functions are defined.
#' data <- rinvgamma(100)
#' Inverse_Gamma_Distribution_Data_Fit(data, TRUE, c(0.1, 0.5, 0.9), 0.05, 0.05)
#' }
#'
#' @export
Inverse_Gamma_Distribution_Data_Fit <- function(data, pertinent_graphs, desired_quantiles, value_of_alpha_for_CI, value_of_alpha_for_CR) {
  
  
  
  #       ----------------   Actuarial Coding Portfolio  ----------------
  
  #                         Ioannis (Yanni) Papadopoulos
  
  #       =================================================================
  
  
  #   Before we begin, we will load all the necessary libraries.
  
  
  
  
  #   We will now begin by loading the data.
  #   We will proceed with the following method to load the data for reproducibility reasons:
  
  
  
  
  #       ======    Part I:    Model Selection    ======
  
  
  # *           -----   Descriptive Statistics  -----
  
  
  #   Preliminary summary statistics:
  
  preliminary_summary_statistics = base::summary(data)
  
  
  #   Observing the sample size:
  
  sample_size = base::length(data)
  
  
  #   Observing the minimum, maximum, and range:
  
  minimum = base::min(data)
  
  maximum = base::max(data)
  
  
  #   The minimum and maximum of the data can also be observed using the following function:
  
  min_and_max_of_data = base::range(data)
  
  #   To find the range, we will take the difference of the maximum and the minimum:
  
  range_difference = maximum - minimum
  
  
  
  #   Observing the mean of the data:
  
  mean_of_data = base::mean(data)
  
  
  
  #   Observing the median of the data:
  
  median_of_data = stats::median(data)
  
  
  
  #   Observing the variance and standard deviation of the data:
  
  variance_of_data = actuar::var(data)
  
  
  standard_deviation_of_data = actuar::sd(data)
  
  
  
  #   Observing the coefficient of variation for the data:
  
  coefficient_of_variation = standard_deviation_of_data/mean_of_data
  
  
  
  #   Calculation of Quantiles:
  
  #       The following are the required quantiles:
  
  Required_Quantiles = desired_quantiles
  
  
  #      We now observe the required quantile values:
  
  quantile_values = stats::quantile(data, Required_Quantiles)
  
  
  #   Observation of skewness and kurtosis.
  #     For this, the 'moments' package is used.
  
  Skewness_of_data = EnvStats::skewness(data)
  
  
  Kurtosis_of_data = EnvStats::kurtosis(data, excess = FALSE)
  
  
  
  
  
  # *           -----   Supporting Graphs  -----
  
  
  #   In this section, the required exploratory graphs are produced.
  
  #   We will also sort the data in increasing order for when such an operation would be appropriate.
  
  
  data_increasing_order = base::sort(data)
  
  
  #   Boxplot:
  
  
  graph_of_boxplot = function(data) {
    
    graphics::boxplot(data)
    graphics::title(main = "Boxplot", ylab = "Data")
    
  }
  
  
  
  #   Histogram:
  
  
  graph_of_histogram = function(data) {
    
    graphics::hist(data, main="Histogram for Data", xlab="Data", col="green")
    
  }
  
  
  
  
  #   Empirical CDF:
  
  eCDF_of_data = stats::ecdf(data)
  
  graph_of_eCDF = function(data) {
    
    #    eCDF_of_data = ecdf(data)
    graphics::plot(eCDF_of_data, main="Empirical CDF for Data", xlab="Data", ylab="F_n")
    
  }
  
  
  
  
  #   In preparation of the Empirical MRL (Mean Residual Life) function, a CDF function is defined.
  
  CDF = function(x){
    
    a = x
    
    for(i in 1:base::length(x) ) {
      
      a[i] = base::sum(x<=x[i]) / base::length(x)
      
    }
    return(a)
  }
  
  
  #   The Empirical MRL function is as follows:
  
  eMRL = function(x){
    
    a = x
    
    x.Fn = CDF(x)
    
    for(i in 1:base::length(x)){
      
      if(x[i] < base::max(x)){
        
        a[i] = base::sum(x[x - x[i] > 0] - x[i]) / (base::length(x)) / (1 - x.Fn[i])
        
      } else {
        
        a[i] = 0
        
      }
      
    }
    
    return(a)
    
  }
  
  
  
  
  #   The graph of the Empirical MRL is as follows:
  
  data_eMRL = eMRL(data_increasing_order)
  
  graph_of_eMRL = function(data) {
    
    #    data_eMRL = eMRL(data_increasing_order)
    graphics::plot(data_increasing_order, data_eMRL, type="l", xlab="Data", ylab="e_n", main="Empirical MRL")
    
  }
  
  
  
  
  
  #   The graph of the Empirical LEV is as follows:
  
  data_eLEV = base::mean(data_increasing_order) - data_eMRL * (1 - eCDF_of_data(data_increasing_order))
  
  graph_of_eLEV = function(data) {
    
    #    data_eLEV = mean(data_increasing_order) - data_eMRL * (1 - eCDF_of_data(data_increasing_order))
    graphics::plot(data_increasing_order, data_eLEV, type="l", xlab="Data", ylab="Limited Expected Value", main="Empirical LEV", las = 1)
    
  }
  
  
  
  
  
  #   To verify the validity of the assumptions, an IID test is performed:
  
  
  graph_of_IID_Test = function(data) {
    
    IID_Test = stats::acf(x = data, lag.max = sample_size, type = "correlation", main="ACF of Data")
    
  }
  
  
  
  
  
  #       ======    Part II:    Estimation    ======
  
  
  
  # *           -----   Method of Moments (MME) Parameter Estimation  -----
  
  
  #   We will proceed with fitting the model distribution candidates with respect to the Method of Moments (MME) approach:
  
  
  #   Inverse Gamma Distribution
  beta_Inverse_Gamma_Fit_MME = (((mean_of_data)^3) / variance_of_data) + mean_of_data
  alpha_Inverse_Gamma_Fit_MME = (beta_Inverse_Gamma_Fit_MME / mean_of_data) + 1
  
  #   Printing the parameters:
  #cat("Inverse Gamma Parameters:\n")
  #cat("Shape (alpha):", alpha_Inverse_Gamma_Fit_MME, "\n")
  #cat("Scale (beta):", beta_Inverse_Gamma_Fit_MME, "\n")
  
  Inverse_Gamma_Fit_MME = base::c(alpha_Inverse_Gamma_Fit_MME, beta_Inverse_Gamma_Fit_MME)
  
  
  
  
  
  # *           -----   Maximum Likelihood Estimator (MLE) Parameter Estimation  -----
  
  
  #   Now, we will proceed with fitting the model distribution candidates with respect to the Maximum Likelihood Estimator (MLE) approach:
  
  
  Inverse_Gamma_Fit_MLE = fitdistrplus::fitdist(data, "invgamma")
  
  
  
  
  
  
  # *           -----   Estimation of Asymptotic Covariance Matrices  -----
  
  
  #   The Information Matrices with respect to each candidate distribution are observed as follows:
  
  
  #     Inverse Gamma:
  
  Inverse_Gamma_Fit_MLE_Estimates = Inverse_Gamma_Fit_MLE$estimate
  
  Inverse_Gamma_Fit_MLE_Estimated_shape = base::as.numeric(Inverse_Gamma_Fit_MLE_Estimates[1])
  Inverse_Gamma_Fit_MLE_Estimated_scale = base::as.numeric(Inverse_Gamma_Fit_MLE_Estimates[2])
  
  Inverse_Gamma_vcov_Matrix = Inverse_Gamma_Fit_MLE$vcov
  
  Inverse_Gamma_Information_Matrix = matlib::Inverse(Inverse_Gamma_vcov_Matrix)
  
  
  
  
  
  
  
  
  
  
  # *           -----   Confidence Regions  -----
  
  
  #   Note:   The Confidence Regions (CR) will be constructed after the Formal Goodness of Fit test have been conducted, which is at the end of the script.
  #               This section will be denoted as "Confidence Regions, revisited"
  
  #   Construction of a Confidence Interval (CI) with respect to each candidate distribution's estimators is observed as follows:
  
  #   The chosen level of confidence is (value_of_alpha_for_CI)% with respect to a two-tailed test.
  
  alpha_for_CI = value_of_alpha_for_CI
  Level_of_Confidence.2_T = 1-(alpha_for_CI/2)
  
  Z_Value_at_LoC.2_T = stats::qnorm(Level_of_Confidence.2_T)
  
  
  #   We now proceed.
  
  
  #     Inverse Gamma:
  
  Lower_bound_of_CI_for_MLE_Inverse_Gamma_shape = Inverse_Gamma_Fit_MLE_Estimated_shape - Z_Value_at_LoC.2_T * base::sqrt(Inverse_Gamma_vcov_Matrix[1, 1])
  Upper_bound_of_CI_for_MLE_Inverse_Gamma_shape = Inverse_Gamma_Fit_MLE_Estimated_shape + Z_Value_at_LoC.2_T * base::sqrt(Inverse_Gamma_vcov_Matrix[1, 1])
  
  CI_of_MLE_Inverse_Gamma_shape = base::as.numeric(base::c(Lower_bound_of_CI_for_MLE_Inverse_Gamma_shape, Upper_bound_of_CI_for_MLE_Inverse_Gamma_shape))
  
  Lower_bound_of_CI_for_MLE_Inverse_Gamma_scale = Inverse_Gamma_Fit_MLE_Estimated_scale - Z_Value_at_LoC.2_T * base::sqrt(Inverse_Gamma_vcov_Matrix[2, 2])
  Upper_bound_of_CI_for_MLE_Inverse_Gamma_scale = Inverse_Gamma_Fit_MLE_Estimated_scale + Z_Value_at_LoC.2_T * base::sqrt(Inverse_Gamma_vcov_Matrix[2, 2])
  
  CI_of_MLE_Inverse_Gamma_scale = base::as.numeric(base::c(Lower_bound_of_CI_for_MLE_Inverse_Gamma_scale, Upper_bound_of_CI_for_MLE_Inverse_Gamma_scale))
  
  
  
  
  
  
  #       ======    Part III:    Goodness-of-fit Tests    ======
  
  
  
  #     Creating the required functions for each distribution before proceeding.
  
  
  #      Inverse Gamma
  
  #   CDF:
  
  CDF_Inverse_Gamma_Fit = function(x){
    
    actuar::pinvgamma(x, shape = Inverse_Gamma_Fit_MLE_Estimated_shape, scale = Inverse_Gamma_Fit_MLE_Estimated_scale)
    
  }
  
  #   Survival Function:
  
  Survival_Function_Inverse_Gamma_Fit = function(x){
    
    1 - CDF_Inverse_Gamma_Fit(x)
    
  }
  
  
  
  #   Probability Density Function:
  
  PDF_Inverse_Gamma_Fit = function(x){
    
    actuar::dinvgamma(x, shape = Inverse_Gamma_Fit_MLE_Estimated_shape, scale = Inverse_Gamma_Fit_MLE_Estimated_scale)
    
  }
  
  
  
  #   Limited Expected Value:
  
  LEV_Inverse_Gamma_Fit = function(t){
    
    stats::integrate(f = Survival_Function_Inverse_Gamma_Fit, lower = 0, upper = t)$value
    
  }
  
  
  
  
  #   Mean Residual Life:
  
  MRL_Inverse_Gamma_Fit = function(t){
    
    ( LEV_Inverse_Gamma_Fit(Inf) - LEV_Inverse_Gamma_Fit(t) ) / Survival_Function_Inverse_Gamma_Fit(t)
    
  }
  
  
  
  
  
  
  # *           -----   Graphical Comparisons  -----
  
  
  #   Now, we will create a graph comparing the Empirical CDF with the Fitted Distributions' CDFs.
  
  #     Firstly, we will create a sequence to aid in the creation of the graph.
  
  Sequence_for_Graphs = base::seq(minimum, maximum, length = sample_size)
  
  #     Also, we will introduce a multiplier for "par("cex")" as well as a line width parameter for flexibility:
  
  Multiplier_for_par.cex = 3
  
  Line_Width = Multiplier_for_par.cex*graphics::par("cex")
  
  
  Empirical_Distribution.text = "Empirical Distribution"
  
  
  Collection_of_Fitted_Distributions = base::c(
    #     "Fitted Normal Distribution"
    #    ,"Fitted Gamma Distribution"
    #    ,"Fitted Weibull Distribution"
    "Fitted Inverse Gamma Distribution"
    #    ,"Fitted Inverse Weibull Distribution"
  )
  
  
  Collection_of_Colours = base::c(
    #     "blue"
    #    ,"orange"
    #    ,"purple"
    "red"
    #    ,"turquoise"
  )
  
  
  
  #   Graphing the Empirical CDF with the Inverse Gamma Fitted CDF overlay:
  
  
  graph_of_eCDF_with_Inverse_Gamma_Fitted_CDF_overlay = function() {
    
    
    #   The Empirical CDF is graphed as follows:
    
    EnvStats::ecdfPlot(data_increasing_order, main = "Graphical Comparison of the eCDF with the Fitted CDFs", xlab = "Data")
    
    
    #   We will now overlay the fitted distributions' CDFs.
    
    Inverse_Gamma_Fitted_CDF = base::sapply(X = Sequence_for_Graphs, CDF_Inverse_Gamma_Fit)
    graphics::lines(x = Sequence_for_Graphs, y = Inverse_Gamma_Fitted_CDF, lwd = Line_Width, col = Collection_of_Colours[1])
    
    
    #   Allowing the legend to be placed outside the plot area:
    graphics::par(xpd = TRUE)
    
    
    # Call locator() in the console and click on the plot to choose the legend position
    legend_position = graphics::locator(1)
    
    # Use the returned coordinates to place the legend
    graphics::legend(x = legend_position$x, y = legend_position$y, base::c(Empirical_Distribution.text, Collection_of_Fitted_Distributions), fill = base::c("black", Collection_of_Colours))
    
    
    # Revoking the allowance of the legend to be placed outside the plot area:
    graphics::par(xpd = FALSE)
    
  }
  
  
  
  
  
  
  
  
  #   Graphing the Empirical LEV with the Inverse Gamma Fitted LEV overlay:
  
  graph_of_eLEV_with_Inverse_Gamma_Fitted_LEV_overlay = function() {
    
    #   The Empirical LEV is graphed as follows:
    
    graphics::plot(actuar::elev(data), type="l", xlab = "Data", ylab = "Limited Expected Value", main="Graphical Comparison of the eLEV with the Fitted LEVs", lwd = Line_Width, las = 1)
    
    
    #   We will now overlay the fitted distributions' LEVs.
    
    Inverse_Gamma_Fitted_LEV = base::sapply(X = Sequence_for_Graphs, LEV_Inverse_Gamma_Fit)
    graphics::lines(x = Sequence_for_Graphs, y = Inverse_Gamma_Fitted_LEV, lwd = Line_Width, col = Collection_of_Colours[1])
    
    
    #   Allowing the legend to be placed outside the plot area:
    graphics::par(xpd = TRUE)
    
    
    # Call locator() in the console and click on the plot to choose the legend position
    legend_position = graphics::locator(1)
    
    # Use the returned coordinates to place the legend
    graphics::legend(x = legend_position$x, y = legend_position$y, base::c(Empirical_Distribution.text, Collection_of_Fitted_Distributions), fill = base::c("black", Collection_of_Colours))
    
    
    # Revoking the allowance of the legend to be placed outside the plot area:
    graphics::par(xpd = FALSE)
    
  }
  
  
  
  
  
  
  
  
  #   Graphing the Empirical MRL with the Inverse Gamma Fitted MRL overlay:
  
  graph_of_eMRL_with_Inverse_Gamma_Fitted_MRL_overlay = function() {
    
    #   The Empirical MRL is graphed as follows:
    
    graphics::plot(data_increasing_order, data_eMRL, type="l", xlab="Data", ylab="e_n", main="Graphical Comparison of the eMRL with the Fitted MRLs", lwd = Line_Width)
    
    
    #   We will now overlay the fitted distributions' MRLs.
    
    Inverse_Gamma_Fitted_MRL = base::sapply(X = Sequence_for_Graphs, MRL_Inverse_Gamma_Fit)
    graphics::lines(x = Sequence_for_Graphs, y = Inverse_Gamma_Fitted_MRL, lwd = Line_Width, col = Collection_of_Colours[1])
    
    
    #   Allowing the legend to be placed outside the plot area:
    graphics::par(xpd = TRUE)
    
    
    # Call locator() in the console and click on the plot to choose the legend position
    legend_position = graphics::locator(1)
    
    # Use the returned coordinates to place the legend
    graphics::legend(x = legend_position$x, y = legend_position$y, base::c(Empirical_Distribution.text, Collection_of_Fitted_Distributions), fill = base::c("black", Collection_of_Colours))
    
    
    # Revoking the allowance of the legend to be placed outside the plot area:
    graphics::par(xpd = FALSE)
    
  }
  
  
  
  
  
  
  
  
  
  #   Graphing the Empirical Histogram with the Inverse Gamma Fitted PDF overlay:
  
  graph_of_histogram_with_Inverse_Gamma_Fitted_PDF_overlay = function() {
    
    #   The histogram of the empirical data is graphed as follows:
    
    graphics::hist(data, prob = T, main="Empirical Histogram compared to the Fitted Empirical Densities", xlab="Data", col="green")
    
    #   We will now overlay the fitted distributions' PDFs altogether.
    
    graphics::lines(data_increasing_order, actuar::dinvgamma(data_increasing_order, shape = Inverse_Gamma_Fit_MLE_Estimated_shape, scale = Inverse_Gamma_Fit_MLE_Estimated_scale), xpd = T, lwd = Line_Width, col = Collection_of_Colours[1])
    
    
    #   Allowing the legend to be placed outside the plot area:
    graphics::par(xpd = TRUE)
    
    
    # Call locator() in the console and click on the plot to choose the legend position
    legend_position = graphics::locator(1)
    
    # Use the returned coordinates to place the legend
    graphics::legend(x = legend_position$x, y = legend_position$y, base::c(Empirical_Distribution.text, Collection_of_Fitted_Distributions), fill = base::c("black", Collection_of_Colours))
    
    
    # Revoking the allowance of the legend to be placed outside the plot area:
    graphics::par(xpd = FALSE)
    
  }
  
  
  
  
  
  
  
  #   Graphing the P-P Plot with the X = Y line overlay:
  
  graph_of_P.P_Plot_with_X.Y_line_overlay = function() {
    
    #     P-P Plots
    
    Legend_Content_for_P.P_Plots = base::c("(Fn, F_Fitted)","X = Y")
    
    Colour_of_X.Y_Line_for_P.P_Plots = "blue"
    Legend_Content_Fill_for_P.P_Plots = base::c("black", Colour_of_X.Y_Line_for_P.P_Plots)
    
    
    
    #   Inverse Gamma P-P Plot:
    
    ppoints_vector = stats::ppoints(sample_size)
    pp_plot = actuar::qinvgamma(ppoints_vector, shape = Inverse_Gamma_Fit_MLE_Estimated_shape, scale = Inverse_Gamma_Fit_MLE_Estimated_scale)
    graphics::plot(data_increasing_order, pp_plot, main = "Inverse Gamma P-P Plot for Data", xlab = "Empirical CDF", ylab = "Fitted Inverse Gamma CDF")
    graphics::abline(0, 1, lwd = Line_Width, col = Colour_of_X.Y_Line_for_P.P_Plots)
    
    
    #   Allowing the legend to be placed outside the plot area:
    graphics::par(xpd = TRUE)
    
    
    # Call locator() in the console and click on the plot to choose the legend position
    legend_position = graphics::locator(1)
    
    # Use the returned coordinates to place the legend
    graphics::legend(x = legend_position$x, y = legend_position$y, Legend_Content_for_P.P_Plots, fill = Legend_Content_Fill_for_P.P_Plots)
    
    
    # Revoking the allowance of the legend to be placed outside the plot area:
    graphics::par(xpd = FALSE)
    
  }
  
  
  
  
  
  
  
  #     Q-Q Plots
  
  
  
  
  #   The following Q-Q Plots compare the graphs with respect to a line of least squares, the 0-1 (y-intercept = 0, slope = 1) line, and the robust line (a line that is fit between the first and third quartiles of the data).
  
  
  Q.Q_Plot.function = function() {
    
    #   The Q-Q Plot for the Fitted Inverse Gamma distribution is as follows:
    inv_gamma_quant <- base::vector()
    for (i in 1:sample_size){
      inv_gamma_quant[i] = actuar::qinvgamma(i/(sample_size+1), shape = Inverse_Gamma_Fit_MLE_Estimated_shape, scale = Inverse_Gamma_Fit_MLE_Estimated_scale)
    }
    graphics::plot(x = data_increasing_order, y = inv_gamma_quant, main="Q-Q Plot", xlab = "Empirical Percentile", ylab = "Inverse Gamma Percentile")
    graphics::lines(x = data_increasing_order, y = data_increasing_order, lwd = Line_Width, col = "blue")
    graphics::legend("bottomright", base::c("(Empirical, Inverse Gamma)","X = Y"), cex = 1, fill=base::c("black","blue"))
    
  }
  
  
  
  
  
  
  
  
  # Define a function to calculate summary statistics
  Fitted_Distribution_Summary_Statistics = function(dist, par, quantiles, stats_fun) {
    q = base::sapply(quantiles, function(x) base::do.call(dist, base::c(base::list(p = x), par)))
    stats = base::do.call(stats_fun, par)
    base::c(quantiles = q, stats)
  }
  
  # Define parameters
  
  fitted_inv_gamma_par = base::list(shape = Inverse_Gamma_Fit_MLE_Estimated_shape, scale = Inverse_Gamma_Fit_MLE_Estimated_scale)
  
  
  # Define functions to calculate additional statistics
  
  fitted_inv_gamma_stats_fun = function(shape, scale) {
    base::list(
      mean = scale / (shape - 1),
      sd = base::sqrt(scale^2 / ((shape - 1)^2 * (shape - 2))),
      cv = base::sqrt((shape - 2) / (shape - 1)),
      skew = 4 * base::sqrt(2) / base::sqrt(shape - 3),
      kurt = 3 + 12 / (shape - 3)
    )
  }
  
  
  
  
  
  # Calculate summary statistics
  
  fitted_inv_gamma_stats = Fitted_Distribution_Summary_Statistics(actuar::qinvgamma, fitted_inv_gamma_par, Required_Quantiles, fitted_inv_gamma_stats_fun)
  
  
  
  
  
  
  
  
  # *           -----   Formal Goodness-of-Fit Tests  -----
  
  
  
  
  #   Conducting the Kolmogorov-Smirnov Test on each fitted distribution:
  
  KS_Test = stats::ks.test(data, actuar::pinvgamma, shape = Inverse_Gamma_Fit_MLE_Estimated_shape, scale = Inverse_Gamma_Fit_MLE_Estimated_scale)
  
  
  
  
  #   Conducting the Anderson-Darling Test on each fitted distribution:
  
  AD_Test = ADGofTest::ad.test(data, actuar::pinvgamma, shape = Inverse_Gamma_Fit_MLE_Estimated_shape, scale = Inverse_Gamma_Fit_MLE_Estimated_scale)
  
  
  
  
  
  #   Conducting the AIC / BIC Tests on each fitted distribution:
  
  #   Preparing to conduct the tests:
  
  
  Inverse_Gamma_Fit_Preliminary_Summary = fitdistrplus::fitdist(data, "invgamma")
  
  
  
  
  #   The test results are observable as follows:
  
  Inverse_Gamma_Fit_Summary = base::summary(Inverse_Gamma_Fit_Preliminary_Summary)
  
  
  
  
  
  
  
  
  # *           -----   Confidence Regions, revisited  -----
  
  
  #   Construction of a Confidence Region (CR) with respect to each candidate distribution's estimators is observed as follows:
  
  #   The chosen level of confidence is (value_of_alpha_for_CR)% with respect to the Chi-Square test.
  
  
  alpha_for_CR = value_of_alpha_for_CR
  Level_of_Confidence.1_T = 1-(alpha_for_CR)
  
  Z_Value_at_LoC.1_T = stats::qnorm(Level_of_Confidence.1_T)
  
  
  #   We now proceed.
  
  
  
  graph_of_confidence_region_for_Inverse_Gamma_Fit_parameters = function() {
    
    
    #       Inverse Gamma:
    
    
    #   Generating the confidence ellipse:
    
    Inverse_Gamma_Fit_MLE_Estimates_CR = ellipse::ellipse(x = Inverse_Gamma_vcov_Matrix, center = Inverse_Gamma_Fit_MLE_Estimates, shape = Inverse_Gamma_vcov_Matrix, radius = Z_Value_at_LoC.1_T)
    
    #   Plotting the confidence ellipse:
    
    graphics::plot(Inverse_Gamma_Fit_MLE_Estimates_CR, type = "l", xlab = "Shape", ylab = "Scale")
    graphics::polygon(Inverse_Gamma_Fit_MLE_Estimates_CR, col = grDevices::rgb(0, 0, 1, 0.5))  # Adding the shaded region
    
    graphics::points(Inverse_Gamma_Fit_MLE_Estimated_shape, Inverse_Gamma_Fit_MLE_Estimated_scale, pch = 19)
    
    graphics::title(main = "CR for the ML Estimators of the Fitted Inverse Gamma Distribution")
    
    
  }
  
  
  
  #   Number of parameters in the Inverse Gamma distribution is observed as follows:
  
  Number_of_Parameters_Inverse_Gamma_Fit_MLE_Estimates = base::length(stats::coef(Inverse_Gamma_Fit_Summary))
  
  
  #   Quantile of the Chi-Square distribution:
  
  Quantile_of_Chi.Square_Dist_Inverse_Gamma_Fit = stats::qchisq(1 - alpha_for_CR, df = Number_of_Parameters_Inverse_Gamma_Fit_MLE_Estimates)
  
  # Calculating the threshold on the log-likelihood for the likelihood ratio test at the alpha = 5% significance level:
  
  Log.Likelihood_Threshold_Inverse_Gamma_Fit = Inverse_Gamma_Fit_Summary$loglik - 0.5 * Quantile_of_Chi.Square_Dist_Inverse_Gamma_Fit
  
  
  
  
  
  
  Useful_List = base::list(
    preliminary_summary_statistics = base::summary(data)
    , sample_size = base::length(data)
    , minimum = base::min(data)
    , maximum = base::max(data)
    , min_and_max_of_data = base::range(data)
    , range_difference = maximum - minimum
    , mean_of_data = base::mean(data)
    , median_of_data = stats::median(data)
    , variance_of_data = actuar::var(data)
    , standard_deviation_of_data = actuar::sd(data)
    , coefficient_of_variation = standard_deviation_of_data/mean_of_data
    , Required_Quantiles = desired_quantiles
    , quantile_values = stats::quantile(data, Required_Quantiles)
    , Skewness_of_data = EnvStats::skewness(data)
    , Kurtosis_of_data = EnvStats::kurtosis(data, excess = FALSE)
    , Inverse_Gamma_Fit_MME = base::c(alpha_Inverse_Gamma_Fit_MME, beta_Inverse_Gamma_Fit_MME)
    , Inverse_Gamma_Fit_MLE = fitdistrplus::fitdist(data, "invgamma")
    , Inverse_Gamma_Fit_MLE_Estimates = Inverse_Gamma_Fit_MLE$estimate
    , Inverse_Gamma_Fit_MLE_Estimated_shape = base::as.numeric(Inverse_Gamma_Fit_MLE_Estimates[1])
    , Inverse_Gamma_Fit_MLE_Estimated_scale = base::as.numeric(Inverse_Gamma_Fit_MLE_Estimates[2])
    , Inverse_Gamma_vcov_Matrix = Inverse_Gamma_Fit_MLE$vcov
    , Inverse_Gamma_Information_Matrix = matlib::Inverse(Inverse_Gamma_vcov_Matrix)
    , alpha_for_CI = value_of_alpha_for_CI
    , Level_of_Confidence.2_T = 1-(alpha_for_CI/2)
    , Z_Value_at_LoC.2_T = stats::qnorm(Level_of_Confidence.2_T)
    , Lower_bound_of_CI_for_MLE_Inverse_Gamma_shape = Inverse_Gamma_Fit_MLE_Estimated_shape - Z_Value_at_LoC.2_T * base::sqrt(Inverse_Gamma_vcov_Matrix[1, 1])
    , Upper_bound_of_CI_for_MLE_Inverse_Gamma_shape = Inverse_Gamma_Fit_MLE_Estimated_shape + Z_Value_at_LoC.2_T * base::sqrt(Inverse_Gamma_vcov_Matrix[1, 1])
    , CI_of_MLE_Inverse_Gamma_shape = base::as.numeric(base::c(Lower_bound_of_CI_for_MLE_Inverse_Gamma_shape, Upper_bound_of_CI_for_MLE_Inverse_Gamma_shape))
    , Lower_bound_of_CI_for_MLE_Inverse_Gamma_scale = Inverse_Gamma_Fit_MLE_Estimated_scale - Z_Value_at_LoC.2_T * base::sqrt(Inverse_Gamma_vcov_Matrix[2, 2])
    , Upper_bound_of_CI_for_MLE_Inverse_Gamma_scale = Inverse_Gamma_Fit_MLE_Estimated_scale + Z_Value_at_LoC.2_T * base::sqrt(Inverse_Gamma_vcov_Matrix[2, 2])
    , CI_of_MLE_Inverse_Gamma_scale = base::as.numeric(base::c(Lower_bound_of_CI_for_MLE_Inverse_Gamma_scale, Upper_bound_of_CI_for_MLE_Inverse_Gamma_scale))
    , fitted_inv_gamma_par = base::list(shape = Inverse_Gamma_Fit_MLE_Estimated_shape, scale = Inverse_Gamma_Fit_MLE_Estimated_scale)
    , fitted_inv_gamma_stats = Fitted_Distribution_Summary_Statistics(actuar::qinvgamma, fitted_inv_gamma_par, Required_Quantiles, fitted_inv_gamma_stats_fun)
    , KS_Test = stats::ks.test(data, actuar::pinvgamma, shape = Inverse_Gamma_Fit_MLE_Estimated_shape, scale = Inverse_Gamma_Fit_MLE_Estimated_scale)
    , AD_Test = ADGofTest::ad.test(data, actuar::pinvgamma, shape = Inverse_Gamma_Fit_MLE_Estimated_shape, scale = Inverse_Gamma_Fit_MLE_Estimated_scale)
    , Inverse_Gamma_Fit_Preliminary_Summary = fitdistrplus::fitdist(data, "invgamma")
    , Inverse_Gamma_Fit_Summary = base::summary(Inverse_Gamma_Fit_Preliminary_Summary)
    , alpha_for_CR = value_of_alpha_for_CR
    , Level_of_Confidence.1_T = 1-(alpha_for_CR)
    , Z_Value_at_LoC.1_T = stats::qnorm(Level_of_Confidence.1_T)
    , Number_of_Parameters_Inverse_Gamma_Fit_MLE_Estimates = base::length(stats::coef(Inverse_Gamma_Fit_Summary))
    , Quantile_of_Chi.Square_Dist_Inverse_Gamma_Fit = stats::qchisq(1 - alpha_for_CR, df = Number_of_Parameters_Inverse_Gamma_Fit_MLE_Estimates)
    , Log.Likelihood_Threshold_Inverse_Gamma_Fit = Inverse_Gamma_Fit_Summary$loglik - 0.5 * Quantile_of_Chi.Square_Dist_Inverse_Gamma_Fit
    
    , if(pertinent_graphs == TRUE) {
      
      base::c(
        Graph_of_eCDF_with_Inverse_Gamma_Fitted_CDF_Overlay = graph_of_eCDF_with_Inverse_Gamma_Fitted_CDF_overlay()
        , Graph_of_eLEV_with_Inverse_Gamma_Fitted_LEV_Overlay = graph_of_eLEV_with_Inverse_Gamma_Fitted_LEV_overlay()
        , Graph_of_eMRL_with_Inverse_Gamma_Fitted_MRL_Overlay = graph_of_eMRL_with_Inverse_Gamma_Fitted_MRL_overlay()
        , Graph_of_Histogram_with_Inverse_Gamma_Fitted_PDF_Overlay = graph_of_histogram_with_Inverse_Gamma_Fitted_PDF_overlay()
        , Graph_of_P.P_Plot_with_X.Y_Line_Overlay = graph_of_P.P_Plot_with_X.Y_line_overlay()
        , Q.Q_Plot = Q.Q_Plot.function()
        , Graph_of_Confidence_Region_for_Inverse_Gamma_Fit_Parameters = graph_of_confidence_region_for_Inverse_Gamma_Fit_parameters()
      )
      
      
    } else {}
    
    
  )
  
  
  Useful_List
  
  
}

#' Inverse Weibull Distribution Data Fit
#'
#' This function fits an inverse Weibull distribution to the given data and provides various statistics and graphs.
#'
#' @param data A numeric vector of data to be fitted.
#' @param pertinent_graphs A logical value indicating whether to generate pertinent graphs.
#' @param desired_quantiles A numeric vector of quantiles to be computed.
#' @param value_of_alpha_for_CI A numeric value for the level of significance for confidence intervals.
#' @param value_of_alpha_for_CR A numeric value for the level of significance for confidence regions.
#'
#' @return A list containing various summary statistics, test results, and graphs (if pertinent_graphs = TRUE).
#'
#' @examples
#' \dontrun{
#' install.packages("actuar")    # See below.
#' library(actuar)               # This is done to ensure the required functions are defined.
#' data <- rinvweibull(100)
#' Inverse_Weibull_Distribution_Data_Fit(data, TRUE, c(0.1, 0.5, 0.9), 0.05, 0.05)
#' }
#'
#' @export
Inverse_Weibull_Distribution_Data_Fit <- function(data, pertinent_graphs, desired_quantiles, value_of_alpha_for_CI, value_of_alpha_for_CR) {
  
  
  
  #       ----------------   Actuarial Coding Portfolio  ----------------
  
  #                         Ioannis (Yanni) Papadopoulos
  
  #       =================================================================
  
  
  #   Before we begin, we will load all the necessary libraries.
  
  
  
  
  #   We will now begin by loading the data.
  #   We will proceed with the following method to load the data for reproducibility reasons:
  
  
  
  
  #       ======    Part I:    Model Selection    ======
  
  
  # *           -----   Descriptive Statistics  -----
  
  
  #   Preliminary summary statistics:
  
  preliminary_summary_statistics = base::summary(data)
  
  
  #   Observing the sample size:
  
  sample_size = base::length(data)
  
  
  #   Observing the minimum, maximum, and range:
  
  minimum = base::min(data)
  
  maximum = base::max(data)
  
  
  #   The minimum and maximum of the data can also be observed using the following function:
  
  min_and_max_of_data = base::range(data)
  
  #   To find the range, we will take the difference of the maximum and the minimum:
  
  range_difference = maximum - minimum
  
  
  
  #   Observing the mean of the data:
  
  mean_of_data = base::mean(data)
  
  
  
  #   Observing the median of the data:
  
  median_of_data = stats::median(data)
  
  
  
  #   Observing the variance and standard deviation of the data:
  
  variance_of_data = actuar::var(data)
  
  
  standard_deviation_of_data = actuar::sd(data)
  
  
  
  #   Observing the coefficient of variation for the data:
  
  coefficient_of_variation = standard_deviation_of_data/mean_of_data
  
  
  
  #   Calculation of Quantiles:
  
  #       The following are the required quantiles:
  
  Required_Quantiles = desired_quantiles
  
  
  #      We now observe the required quantile values:
  
  quantile_values = stats::quantile(data, Required_Quantiles)
  
  
  #   Observation of skewness and kurtosis.
  #     For this, the 'moments' package is used.
  
  Skewness_of_data = EnvStats::skewness(data)
  
  
  Kurtosis_of_data = EnvStats::kurtosis(data, excess = FALSE)
  
  
  
  
  
  # *           -----   Supporting Graphs  -----
  
  
  #   In this section, the required exploratory graphs are produced.
  
  #   We will also sort the data in increasing order for when such an operation would be appropriate.
  
  
  data_increasing_order = base::sort(data)
  
  
  #   Boxplot:
  
  
  graph_of_boxplot = function(data) {
    
    graphics::boxplot(data)
    graphics::title(main = "Boxplot", ylab = "Data")
    
  }
  
  
  
  #   Histogram:
  
  
  graph_of_histogram = function(data) {
    
    graphics::hist(data, main="Histogram for Data", xlab="Data", col="green")
    
  }
  
  
  
  
  #   Empirical CDF:
  
  eCDF_of_data = stats::ecdf(data)
  
  graph_of_eCDF = function(data) {
    
    #    eCDF_of_data = ecdf(data)
    graphics::plot(eCDF_of_data, main="Empirical CDF for Data", xlab="Data", ylab="F_n")
    
  }
  
  
  
  
  #   In preparation of the Empirical MRL (Mean Residual Life) function, a CDF function is defined.
  
  CDF = function(x){
    
    a = x
    
    for(i in 1:base::length(x) ) {
      
      a[i] = base::sum(x<=x[i]) / base::length(x)
      
    }
    return(a)
  }
  
  
  #   The Empirical MRL function is as follows:
  
  eMRL = function(x){
    
    a = x
    
    x.Fn = CDF(x)
    
    for(i in 1:base::length(x)){
      
      if(x[i] < base::max(x)){
        
        a[i] = base::sum(x[x - x[i] > 0] - x[i]) / (base::length(x)) / (1 - x.Fn[i])
        
      } else {
        
        a[i] = 0
        
      }
      
    }
    
    return(a)
    
  }
  
  
  
  
  #   The graph of the Empirical MRL is as follows:
  
  data_eMRL = eMRL(data_increasing_order)
  
  graph_of_eMRL = function(data) {
    
    #    data_eMRL = eMRL(data_increasing_order)
    graphics::plot(data_increasing_order, data_eMRL, type="l", xlab="Data", ylab="e_n", main="Empirical MRL")
    
  }
  
  
  
  
  
  #   The graph of the Empirical LEV is as follows:
  
  data_eLEV = base::mean(data_increasing_order) - data_eMRL * (1 - eCDF_of_data(data_increasing_order))
  
  graph_of_eLEV = function(data) {
    
    #    data_eLEV = mean(data_increasing_order) - data_eMRL * (1 - eCDF_of_data(data_increasing_order))
    graphics::plot(data_increasing_order, data_eLEV, type="l", xlab="Data", ylab="Limited Expected Value", main="Empirical LEV", las = 1)
    
  }
  
  
  
  
  
  #   To verify the validity of the assumptions, an IID test is performed:
  
  
  graph_of_IID_Test = function(data) {
    
    IID_Test = stats::acf(x = data, lag.max = sample_size, type = "correlation", main="ACF of Data")
    
  }
  
  
  
  
  
  #       ======    Part II:    Estimation    ======
  
  
  
  # *           -----   Method of Moments (MME) Parameter Estimation  -----
  
  
  #   We will proceed with fitting the model distribution candidates with respect to the Method of Moments (MME) approach:
  
  
  #   Inverse Weibull Distribution
  
  #   Define the system of equations
  eqns = function(params) {
    k = params[1]
    lambda = params[2]
    mean_eqn = lambda*base::gamma(1 + 1/k) - mean_of_data
    var_eqn = lambda^2 * base::gamma(1 + 2/k) - mean_of_data^2 - variance_of_data
    base::c(mean_eqn, var_eqn)
  }
  
  #   Solve the system of equations
  start_params = base::c(10, 10)  # initial guesses for k and lambda
  sol = nleqslv::nleqslv(start_params, eqns)
  
  #   Print the parameters
  #cat("Inverse Weibull Parameters:\n")
  #cat("Shape (k):", sol$x[1], "\n")
  #cat("Scale (lambda):", sol$x[2], "\n")
  
  Inverse_Weibull_Fit_MME = base::c(sol$x[1], sol$x[2])
  
  
  
  
  
  
  # *           -----   Maximum Likelihood Estimator (MLE) Parameter Estimation  -----
  
  
  #   Now, we will proceed with fitting the model distribution candidates with respect to the Maximum Likelihood Estimator (MLE) approach:
  
  
  Inverse_Weibull_Fit_MLE = fitdistrplus::fitdist(data, "invweibull")
  
  
  
  
  
  
  # *           -----   Estimation of Asymptotic Covariance Matrices  -----
  
  
  #   The Information Matrices with respect to each candidate distribution are observed as follows:
  
  
  #     Inverse Weibull:
  
  Inverse_Weibull_Fit_MLE_Estimates = Inverse_Weibull_Fit_MLE$estimate
  
  Inverse_Weibull_Fit_MLE_Estimated_shape = base::as.numeric(Inverse_Weibull_Fit_MLE_Estimates[1])
  Inverse_Weibull_Fit_MLE_Estimated_scale = base::as.numeric(Inverse_Weibull_Fit_MLE_Estimates[2])
  
  Inverse_Weibull_vcov_Matrix = Inverse_Weibull_Fit_MLE$vcov
  
  Inverse_Weibull_Information_Matrix = matlib::Inverse(Inverse_Weibull_vcov_Matrix)
  
  
  
  
  
  
  
  
  
  
  # *           -----   Confidence Regions  -----
  
  
  #   Note:   The Confidence Regions (CR) will be constructed after the Formal Goodness of Fit test have been conducted, which is at the end of the script.
  #               This section will be denoted as "Confidence Regions, revisited"
  
  #   Construction of a Confidence Interval (CI) with respect to each candidate distribution's estimators is observed as follows:
  
  #   The chosen level of confidence is (value_of_alpha_for_CI)% with respect to a two-tailed test.
  
  alpha_for_CI = value_of_alpha_for_CI
  Level_of_Confidence.2_T = 1-(alpha_for_CI/2)
  
  Z_Value_at_LoC.2_T = stats::qnorm(Level_of_Confidence.2_T)
  
  
  #   We now proceed.
  
  
  #     Inverse Weibull:
  
  Lower_bound_of_CI_for_MLE_Inverse_Weibull_shape = Inverse_Weibull_Fit_MLE_Estimated_shape - Z_Value_at_LoC.2_T * base::sqrt(Inverse_Weibull_vcov_Matrix[1, 1])
  Upper_bound_of_CI_for_MLE_Inverse_Weibull_shape = Inverse_Weibull_Fit_MLE_Estimated_shape + Z_Value_at_LoC.2_T * base::sqrt(Inverse_Weibull_vcov_Matrix[1, 1])
  
  CI_of_MLE_Inverse_Weibull_shape = base::as.numeric(base::c(Lower_bound_of_CI_for_MLE_Inverse_Weibull_shape, Upper_bound_of_CI_for_MLE_Inverse_Weibull_shape))
  
  Lower_bound_of_CI_for_MLE_Inverse_Weibull_scale = Inverse_Weibull_Fit_MLE_Estimated_scale - Z_Value_at_LoC.2_T * base::sqrt(Inverse_Weibull_vcov_Matrix[2, 2])
  Upper_bound_of_CI_for_MLE_Inverse_Weibull_scale = Inverse_Weibull_Fit_MLE_Estimated_scale + Z_Value_at_LoC.2_T * base::sqrt(Inverse_Weibull_vcov_Matrix[2, 2])
  
  CI_of_MLE_Inverse_Weibull_scale = base::as.numeric(base::c(Lower_bound_of_CI_for_MLE_Inverse_Weibull_scale, Upper_bound_of_CI_for_MLE_Inverse_Weibull_scale))
  
  
  
  
  
  
  #       ======    Part III:    Goodness-of-fit Tests    ======
  
  
  
  #     Creating the required functions for each distribution before proceeding.
  
  
  #      Inverse Weibull
  
  #   CDF:
  
  CDF_Inverse_Weibull_Fit = function(x){
    
    actuar::pinvweibull(x, shape = Inverse_Weibull_Fit_MLE_Estimated_shape, scale = Inverse_Weibull_Fit_MLE_Estimated_scale)
    
  }
  
  #   Survival Function:
  
  Survival_Function_Inverse_Weibull_Fit = function(x){
    
    1 - CDF_Inverse_Weibull_Fit(x)
    
  }
  
  
  
  #   Probability Density Function:
  
  PDF_Inverse_Weibull_Fit = function(x){
    
    actuar::dinvweibull(x, shape = Inverse_Weibull_Fit_MLE_Estimated_shape, scale = Inverse_Weibull_Fit_MLE_Estimated_scale)
    
  }
  
  
  
  #   Limited Expected Value:
  
  LEV_Inverse_Weibull_Fit = function(t){
    
    stats::integrate(f = Survival_Function_Inverse_Weibull_Fit, lower = 0, upper = t)$value
    
  }
  
  
  
  
  #   Mean Residual Life:
  
  MRL_Inverse_Weibull_Fit = function(t){
    
    ( LEV_Inverse_Weibull_Fit(Inf) - LEV_Inverse_Weibull_Fit(t) ) / Survival_Function_Inverse_Weibull_Fit(t)
    
  }
  
  
  
  
  
  
  # *           -----   Graphical Comparisons  -----
  
  
  #   Now, we will create a graph comparing the Empirical CDF with the Fitted Distributions' CDFs.
  
  #     Firstly, we will create a sequence to aid in the creation of the graph.
  
  Sequence_for_Graphs = base::seq(minimum, maximum, length = sample_size)
  
  #     Also, we will introduce a multiplier for "par("cex")" as well as a line width parameter for flexibility:
  
  Multiplier_for_par.cex = 3
  
  Line_Width = Multiplier_for_par.cex*graphics::par("cex")
  
  
  Empirical_Distribution.text = "Empirical Distribution"
  
  
  Collection_of_Fitted_Distributions = base::c(
    #     "Fitted Normal Distribution"
    #    ,"Fitted Gamma Distribution"
    #    ,"Fitted Weibull Distribution"
    #    ,"Fitted Inverse Gamma Distribution"
    "Fitted Inverse Weibull Distribution"
  )
  
  
  Collection_of_Colours = base::c(
    #     "blue"
    #    ,"orange"
    #    ,"purple"
    #    ,"red"
    "turquoise"
  )
  
  
  
  #   Graphing the Empirical CDF with the Inverse Weibull Fitted CDF overlay:
  
  
  graph_of_eCDF_with_Inverse_Weibull_Fitted_CDF_overlay = function() {
    
    
    #   The Empirical CDF is graphed as follows:
    
    EnvStats::ecdfPlot(data_increasing_order, main = "Graphical Comparison of the eCDF with the Fitted CDFs", xlab = "Data")
    
    
    #   We will now overlay the fitted distributions' CDFs.
    
    Inverse_Weibull_Fitted_CDF = base::sapply(X = Sequence_for_Graphs, CDF_Inverse_Weibull_Fit)
    graphics::lines(x = Sequence_for_Graphs, y = Inverse_Weibull_Fitted_CDF, lwd = Line_Width, col = Collection_of_Colours[1])
    
    
    #   Allowing the legend to be placed outside the plot area:
    graphics::par(xpd = TRUE)
    
    
    # Call locator() in the console and click on the plot to choose the legend position
    legend_position = graphics::locator(1)
    
    # Use the returned coordinates to place the legend
    graphics::legend(x = legend_position$x, y = legend_position$y, base::c(Empirical_Distribution.text, Collection_of_Fitted_Distributions), fill = base::c("black", Collection_of_Colours))
    
    
    # Revoking the allowance of the legend to be placed outside the plot area:
    graphics::par(xpd = FALSE)
    
  }
  
  
  
  
  
  
  
  
  #   Graphing the Empirical LEV with the Inverse Weibull Fitted LEV overlay:
  
  graph_of_eLEV_with_Inverse_Weibull_Fitted_LEV_overlay = function() {
    
    #   The Empirical LEV is graphed as follows:
    
    graphics::plot(actuar::elev(data), type="l", xlab = "Data", ylab = "Limited Expected Value", main="Graphical Comparison of the eLEV with the Fitted LEVs", lwd = Line_Width, las = 1)
    
    
    #   We will now overlay the fitted distributions' LEVs.
    
    Inverse_Weibull_Fitted_LEV = base::sapply(X = Sequence_for_Graphs, LEV_Inverse_Weibull_Fit)
    graphics::lines(x = Sequence_for_Graphs, y = Inverse_Weibull_Fitted_LEV, lwd = Line_Width, col = Collection_of_Colours[1])
    
    
    #   Allowing the legend to be placed outside the plot area:
    graphics::par(xpd = TRUE)
    
    
    # Call locator() in the console and click on the plot to choose the legend position
    legend_position = graphics::locator(1)
    
    # Use the returned coordinates to place the legend
    graphics::legend(x = legend_position$x, y = legend_position$y, base::c(Empirical_Distribution.text, Collection_of_Fitted_Distributions), fill = base::c("black", Collection_of_Colours))
    
    
    # Revoking the allowance of the legend to be placed outside the plot area:
    graphics::par(xpd = FALSE)
    
  }
  
  
  
  
  
  
  
  
  #   Graphing the Empirical MRL with the Inverse Weibull Fitted MRL overlay:
  
  graph_of_eMRL_with_Inverse_Weibull_Fitted_MRL_overlay = function() {
    
    #   The Empirical MRL is graphed as follows:
    
    graphics::plot(data_increasing_order, data_eMRL, type="l", xlab="Data", ylab="e_n", main="Graphical Comparison of the eMRL with the Fitted MRLs", lwd = Line_Width)
    
    
    #   We will now overlay the fitted distributions' MRLs.
    
    Inverse_Weibull_Fitted_MRL = base::sapply(X = Sequence_for_Graphs, MRL_Inverse_Weibull_Fit)
    graphics::lines(x = Sequence_for_Graphs, y = Inverse_Weibull_Fitted_MRL, lwd = Line_Width, col = Collection_of_Colours[1])
    
    
    #   Allowing the legend to be placed outside the plot area:
    graphics::par(xpd = TRUE)
    
    
    # Call locator() in the console and click on the plot to choose the legend position
    legend_position = graphics::locator(1)
    
    # Use the returned coordinates to place the legend
    graphics::legend(x = legend_position$x, y = legend_position$y, base::c(Empirical_Distribution.text, Collection_of_Fitted_Distributions), fill = base::c("black", Collection_of_Colours))
    
    
    # Revoking the allowance of the legend to be placed outside the plot area:
    graphics::par(xpd = FALSE)
    
  }
  
  
  
  
  
  
  
  
  
  #   Graphing the Empirical Histogram with the Inverse Weibull Fitted PDF overlay:
  
  graph_of_histogram_with_Inverse_Weibull_Fitted_PDF_overlay = function() {
    
    #   The histogram of the empirical data is graphed as follows:
    
    graphics::hist(data, prob = T, main="Empirical Histogram compared to the Fitted Empirical Densities", xlab="Data", col="green")
    
    #   We will now overlay the fitted distributions' PDFs altogether.
    
    graphics::lines(data_increasing_order, actuar::dinvweibull(data_increasing_order, shape = Inverse_Weibull_Fit_MLE_Estimated_shape, scale = Inverse_Weibull_Fit_MLE_Estimated_scale), xpd = T, lwd = Line_Width, col = Collection_of_Colours[1])
    
    
    #   Allowing the legend to be placed outside the plot area:
    graphics::par(xpd = TRUE)
    
    
    # Call locator() in the console and click on the plot to choose the legend position
    legend_position = graphics::locator(1)
    
    # Use the returned coordinates to place the legend
    graphics::legend(x = legend_position$x, y = legend_position$y, base::c(Empirical_Distribution.text, Collection_of_Fitted_Distributions), fill = base::c("black", Collection_of_Colours))
    
    
    # Revoking the allowance of the legend to be placed outside the plot area:
    graphics::par(xpd = FALSE)
    
  }
  
  
  
  
  
  
  
  #   Graphing the P-P Plot with the X = Y line overlay:
  
  graph_of_P.P_Plot_with_X.Y_line_overlay = function() {
    
    #     P-P Plots
    
    Legend_Content_for_P.P_Plots = base::c("(Fn, F_Fitted)","X = Y")
    
    Colour_of_X.Y_Line_for_P.P_Plots = "blue"
    Legend_Content_Fill_for_P.P_Plots = base::c("black", Colour_of_X.Y_Line_for_P.P_Plots)
    
    
    
    #   Inverse Weibull P-P Plot:
    
    ppoints_vector = stats::ppoints(sample_size)
    pp_plot = actuar::qinvweibull(ppoints_vector, shape = Inverse_Weibull_Fit_MLE_Estimated_shape, scale = Inverse_Weibull_Fit_MLE_Estimated_scale)
    graphics::plot(data_increasing_order, pp_plot, main = "Inverse Weibull P-P Plot for Data", xlab = "Empirical CDF", ylab = "Fitted Inverse Weibull CDF")
    graphics::abline(0, 1, lwd = Line_Width, col = Colour_of_X.Y_Line_for_P.P_Plots)
    
    
    #   Allowing the legend to be placed outside the plot area:
    graphics::par(xpd = TRUE)
    
    
    # Call locator() in the console and click on the plot to choose the legend position
    legend_position = graphics::locator(1)
    
    # Use the returned coordinates to place the legend
    graphics::legend(x = legend_position$x, y = legend_position$y, Legend_Content_for_P.P_Plots, fill = Legend_Content_Fill_for_P.P_Plots)
    
    
    # Revoking the allowance of the legend to be placed outside the plot area:
    graphics::par(xpd = FALSE)
    
  }
  
  
  
  
  
  
  
  #     Q-Q Plots
  
  
  
  
  #   The following Q-Q Plots compare the graphs with respect to a line of least squares, the 0-1 (y-intercept = 0, slope = 1) line, and the robust line (a line that is fit between the first and third quartiles of the data).
  
  
  Q.Q_Plot.function = function() {
    
    #   The Q-Q Plot for the Fitted Inverse Weibull distribution is as follows:
    inv_weibull_quant <- base::vector()
    for (i in 1:sample_size){
      inv_weibull_quant[i] = actuar::qinvweibull(i/(sample_size+1), shape = Inverse_Weibull_Fit_MLE_Estimated_shape, scale = Inverse_Weibull_Fit_MLE_Estimated_scale)
    }
    graphics::plot(x=data_increasing_order,y=inv_weibull_quant, main="Q-Q Plot", xlab="Empirical Percentile", ylab="Inverse Weibull Percentile" )
    graphics::lines(x = data_increasing_order, y = data_increasing_order, lwd = Line_Width, col = "blue")
    graphics::legend("bottomright", base::c("(Empirical, Inverse Weibull)","X = Y"), cex=1, fill=base::c("black","blue"))
    
  }
  
  
  
  
  
  
  
  
  # Define a function to calculate summary statistics
  Fitted_Distribution_Summary_Statistics = function(dist, par, quantiles, stats_fun) {
    q = base::sapply(quantiles, function(x) base::do.call(dist, base::c(base::list(p = x), par)))
    stats = base::do.call(stats_fun, par)
    base::c(quantiles = q, stats)
  }
  
  # Define parameters
  
  fitted_inv_weibull_par = base::list(shape = Inverse_Weibull_Fit_MLE_Estimated_shape, scale = Inverse_Weibull_Fit_MLE_Estimated_scale)
  
  
  # Define functions to calculate additional statistics
  
  fitted_inv_weibull_stats_fun = function(shape, scale) {
    list(
      mean = scale * base::gamma(1 - 1 / shape),
      sd = base::sqrt(scale^2 * (base::gamma(1 - 2 / shape) - (base::gamma(1 - 1 / shape))^2)),
      cv = base::sqrt(base::gamma(1 - 2 / shape) - (base::gamma(1 - 1 / shape))^2) / base::gamma(1 - 1 / shape),
      skew = ( base::gamma(1 - 3 / shape) - 3 * base::gamma(1 - 2 / shape) * base::gamma(1 - 1 / shape) + 2*(base::gamma(1 - 1 / shape))^3 ) / ( base::sqrt((base::gamma(1 - 2 / shape) - (base::gamma(1 - 1 / shape))^2)^3)),
      kurt = ((base::gamma(1 - 4 / shape) - 4*base::gamma(1 - 3/shape) * base::gamma(1 - 1/shape) + 3*(base::gamma(1 - 2/shape))^2) / (base::gamma(1 - 2/shape) - (base::gamma(1 - 1/shape))^2)^2)-3
    )
  }
  
  
  
  
  
  # Calculate summary statistics
  
  fitted_inv_weibull_stats = Fitted_Distribution_Summary_Statistics(actuar::qinvweibull, fitted_inv_weibull_par, Required_Quantiles, fitted_inv_weibull_stats_fun)
  
  
  
  
  
  
  
  
  # *           -----   Formal Goodness-of-Fit Tests  -----
  
  
  
  
  #   Conducting the Kolmogorov-Smirnov Test on each fitted distribution:
  
  KS_Test = stats::ks.test(data, actuar::pinvweibull, shape = Inverse_Weibull_Fit_MLE_Estimated_shape, scale = Inverse_Weibull_Fit_MLE_Estimated_scale)
  
  
  
  
  #   Conducting the Anderson-Darling Test on each fitted distribution:
  
  AD_Test = ADGofTest::ad.test(data, actuar::pinvweibull, shape = Inverse_Weibull_Fit_MLE_Estimated_shape, scale = Inverse_Weibull_Fit_MLE_Estimated_scale)
  
  
  
  
  
  #   Conducting the AIC / BIC Tests on each fitted distribution:
  
  #   Preparing to conduct the tests:
  
  
  Inverse_Weibull_Fit_Preliminary_Summary = fitdistrplus::fitdist(data, "invweibull")
  
  
  
  
  #   The test results are observable as follows:
  
  Inverse_Weibull_Fit_Summary = base::summary(Inverse_Weibull_Fit_Preliminary_Summary)
  
  
  
  
  
  
  
  
  # *           -----   Confidence Regions, revisited  -----
  
  
  #   Construction of a Confidence Region (CR) with respect to each candidate distribution's estimators is observed as follows:
  
  #   The chosen level of confidence is (value_of_alpha_for_CR)% with respect to the Chi-Square test.
  
  
  alpha_for_CR = value_of_alpha_for_CR
  Level_of_Confidence.1_T = 1-(alpha_for_CR)
  
  Z_Value_at_LoC.1_T = stats::qnorm(Level_of_Confidence.1_T)
  
  
  #   We now proceed.
  
  
  
  graph_of_confidence_region_for_Inverse_Weibull_Fit_parameters = function() {
    
    
    #       Inverse Weibull:
    
    
    #   Generating the confidence ellipse:
    
    Inverse_Weibull_Fit_MLE_Estimates_CR = ellipse::ellipse(x = Inverse_Weibull_vcov_Matrix, center = Inverse_Weibull_Fit_MLE_Estimates, shape = Inverse_Weibull_vcov_Matrix, radius = Z_Value_at_LoC.1_T)
    
    #   Plotting the confidence ellipse:
    
    graphics::plot(Inverse_Weibull_Fit_MLE_Estimates_CR, type = "l", xlab = "Shape", ylab = "Scale")
    graphics::polygon(Inverse_Weibull_Fit_MLE_Estimates_CR, col = grDevices::rgb(0, 0, 1, 0.5))  # Adding the shaded region
    
    graphics::points(Inverse_Weibull_Fit_MLE_Estimated_shape, Inverse_Weibull_Fit_MLE_Estimated_scale, pch = 19)
    
    graphics::title(main = "CR for the ML Estimators of the Fitted Inverse Weibull Distribution")
    
    
  }
  
  
  
  #   Number of parameters in the Inverse Weibull distribution is observed as follows:
  
  Number_of_Parameters_Inverse_Weibull_Fit_MLE_Estimates = base::length(stats::coef(Inverse_Weibull_Fit_Summary))
  
  
  #   Quantile of the Chi-Square distribution:
  
  Quantile_of_Chi.Square_Dist_Inverse_Weibull_Fit = stats::qchisq(1 - alpha_for_CR, df = Number_of_Parameters_Inverse_Weibull_Fit_MLE_Estimates)
  
  # Calculating the threshold on the log-likelihood for the likelihood ratio test at the alpha = 5% significance level:
  
  Log.Likelihood_Threshold_Inverse_Weibull_Fit = Inverse_Weibull_Fit_Summary$loglik - 0.5 * Quantile_of_Chi.Square_Dist_Inverse_Weibull_Fit
  
  
  
  
  
  
  Useful_List = base::list(
    preliminary_summary_statistics = base::summary(data)
    , sample_size = base::length(data)
    , minimum = base::min(data)
    , maximum = base::max(data)
    , min_and_max_of_data = base::range(data)
    , range_difference = maximum - minimum
    , mean_of_data = base::mean(data)
    , median_of_data = stats::median(data)
    , variance_of_data = actuar::var(data)
    , standard_deviation_of_data = actuar::sd(data)
    , coefficient_of_variation = standard_deviation_of_data/mean_of_data
    , Required_Quantiles = desired_quantiles
    , quantile_values = stats::quantile(data, Required_Quantiles)
    , Skewness_of_data = EnvStats::skewness(data)
    , Kurtosis_of_data = EnvStats::kurtosis(data, excess = FALSE)
    , Inverse_Weibull_Fit_MME = base::c(sol$x[1], sol$x[2])
    , Inverse_Weibull_Fit_MLE = fitdistrplus::fitdist(data, "invweibull")
    , Inverse_Weibull_Fit_MLE_Estimates = Inverse_Weibull_Fit_MLE$estimate
    , Inverse_Weibull_Fit_MLE_Estimated_shape = base::as.numeric(Inverse_Weibull_Fit_MLE_Estimates[1])
    , Inverse_Weibull_Fit_MLE_Estimated_scale = base::as.numeric(Inverse_Weibull_Fit_MLE_Estimates[2])
    , Inverse_Weibull_vcov_Matrix = Inverse_Weibull_Fit_MLE$vcov
    , Inverse_Weibull_Information_Matrix = matlib::Inverse(Inverse_Weibull_vcov_Matrix)
    , alpha_for_CI = value_of_alpha_for_CI
    , Level_of_Confidence.2_T = 1-(alpha_for_CI/2)
    , Z_Value_at_LoC.2_T = stats::qnorm(Level_of_Confidence.2_T)
    , Lower_bound_of_CI_for_MLE_Inverse_Weibull_shape = Inverse_Weibull_Fit_MLE_Estimated_shape - Z_Value_at_LoC.2_T * base::sqrt(Inverse_Weibull_vcov_Matrix[1, 1])
    , Upper_bound_of_CI_for_MLE_Inverse_Weibull_shape = Inverse_Weibull_Fit_MLE_Estimated_shape + Z_Value_at_LoC.2_T * base::sqrt(Inverse_Weibull_vcov_Matrix[1, 1])
    , CI_of_MLE_Inverse_Weibull_shape = base::as.numeric(base::c(Lower_bound_of_CI_for_MLE_Inverse_Weibull_shape, Upper_bound_of_CI_for_MLE_Inverse_Weibull_shape))
    , Lower_bound_of_CI_for_MLE_Inverse_Weibull_scale = Inverse_Weibull_Fit_MLE_Estimated_scale - Z_Value_at_LoC.2_T * base::sqrt(Inverse_Weibull_vcov_Matrix[2, 2])
    , Upper_bound_of_CI_for_MLE_Inverse_Weibull_scale = Inverse_Weibull_Fit_MLE_Estimated_scale + Z_Value_at_LoC.2_T * base::sqrt(Inverse_Weibull_vcov_Matrix[2, 2])
    , CI_of_MLE_Inverse_Weibull_scale = base::as.numeric(base::c(Lower_bound_of_CI_for_MLE_Inverse_Weibull_scale, Upper_bound_of_CI_for_MLE_Inverse_Weibull_scale))
    , fitted_inv_weibull_par = base::list(shape = Inverse_Weibull_Fit_MLE_Estimated_shape, scale = Inverse_Weibull_Fit_MLE_Estimated_scale)
    , fitted_inv_weibull_stats = Fitted_Distribution_Summary_Statistics(actuar::qinvweibull, fitted_inv_weibull_par, Required_Quantiles, fitted_inv_weibull_stats_fun)
    , KS_Test = stats::ks.test(data, actuar::pinvweibull, shape = Inverse_Weibull_Fit_MLE_Estimated_shape, scale = Inverse_Weibull_Fit_MLE_Estimated_scale)
    , AD_Test = ADGofTest::ad.test(data, actuar::pinvweibull, shape = Inverse_Weibull_Fit_MLE_Estimated_shape, scale = Inverse_Weibull_Fit_MLE_Estimated_scale)
    , Inverse_Weibull_Fit_Preliminary_Summary = fitdistrplus::fitdist(data, "invweibull")
    , Inverse_Weibull_Fit_Summary = base::summary(Inverse_Weibull_Fit_Preliminary_Summary)
    , alpha_for_CR = value_of_alpha_for_CR
    , Level_of_Confidence.1_T = 1-(alpha_for_CR)
    , Z_Value_at_LoC.1_T = stats::qnorm(Level_of_Confidence.1_T)
    , Number_of_Parameters_Inverse_Weibull_Fit_MLE_Estimates = base::length(stats::coef(Inverse_Weibull_Fit_Summary))
    , Quantile_of_Chi.Square_Dist_Inverse_Weibull_Fit = stats::qchisq(1 - alpha_for_CR, df = Number_of_Parameters_Inverse_Weibull_Fit_MLE_Estimates)
    , Log.Likelihood_Threshold_Inverse_Weibull_Fit = Inverse_Weibull_Fit_Summary$loglik - 0.5 * Quantile_of_Chi.Square_Dist_Inverse_Weibull_Fit
    
    , if(pertinent_graphs == TRUE) {
      
      base::c(
        Graph_of_eCDF_with_Inverse_Weibull_Fitted_CDF_Overlay = graph_of_eCDF_with_Inverse_Weibull_Fitted_CDF_overlay()
        , Graph_of_eLEV_with_Inverse_Weibull_Fitted_LEV_Overlay = graph_of_eLEV_with_Inverse_Weibull_Fitted_LEV_overlay()
        , Graph_of_eMRL_with_Inverse_Weibull_Fitted_MRL_Overlay = graph_of_eMRL_with_Inverse_Weibull_Fitted_MRL_overlay()
        , Graph_of_Histogram_with_Inverse_Weibull_Fitted_PDF_Overlay = graph_of_histogram_with_Inverse_Weibull_Fitted_PDF_overlay()
        , Graph_of_P.P_Plot_with_X.Y_Line_Overlay = graph_of_P.P_Plot_with_X.Y_line_overlay()
        , Q.Q_Plot = Q.Q_Plot.function()
        , Graph_of_Confidence_Region_for_Inverse_Weibull_Fit_Parameters = graph_of_confidence_region_for_Inverse_Weibull_Fit_parameters()
      )
      
      
    } else {}
    
    
  )
  
  
  Useful_List
  
  
}
