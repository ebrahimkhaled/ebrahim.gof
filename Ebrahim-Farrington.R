#' Farrington Goodness-of-Fit Test
#'
#' Computes the Farrington test for goodness-of-fit in binomial logistic regression.
#' This implementation is adapted from the SAS GOFLOGIT macro and is optimized
#' for both grouped and binary data.
#'
#' @param y Response vector. For binary data, this is a vector of 0s and 1s.
#'   For grouped data, this is the vector of success counts.
#' @param predicted_probs Vector of predicted probabilities from the logistic model.
#' @param model The fitted `glm` object.
#' @param m Optional vector of trial counts for each observation (for grouped data).
#'   If NULL (default), the data is assumed to be binary (`m=1` for all observations),
#'   and the test is not applicable, returning a p-value of 1.
#' @param G Optional integer. If provided, binary data will be grouped into G groups
#'   ordered by predicted probability. Each group will have approximately n/G observations.
#'   This converts binary data into grouped data format for the test.
#' @return A data frame containing the test statistic and p-value.
#'
#' @author Converted from SAS GOFLOGIT by Ebrahim
#' @references
#' - Farrington, C. P. (1996). On Assessing Goodness of Fit of Generalized Linear Models to Sparse Data.
#'   Journal of the Royal Statistical Society. Series B (Methodological), 58(2), 349–360.
#' - SAS Institute Inc. "goflogit.sas". SAS/IML Macro.
#' 
farrington_test <- function(y, predicted_probs, model = NULL, m = NULL, G = 10) {

  n <- length(y)
  if (n != length(predicted_probs)) {
    stop("Length mismatch: y and predicted_probs, they must be of the same length")
  }
  
  # If G is provided, group the binary data
  if (!is.null(G)) {
    if (!is.null(model) || !is.null(m)) {
      cat("Note: When using Ebrahim-Farrington grouping (i.e., grouping binary data via G), you do NOT need to provide the number of trials (m) vector or the model specification. These are only required for the original Farrington test on grouped data. Ignoring provided 'model' or 'm' for grouping mode.\n")
    }
    test <- "Ebrahim-Farrington test"
    # Sort data by predicted probabilities
    sorted_indices <- order(predicted_probs)
    sorted_y <- y[sorted_indices]
    sorted_probs <- predicted_probs[sorted_indices]
    
    # Create groups of approximately equal size
    group_size <- floor(n / G)
    remainder <- n %% G
    
    # Create group assignments
    group_assignments <- rep(1:G, each = group_size)
    if (remainder > 0) {
      group_assignments <- c(group_assignments, rep(G, remainder))
    }
    
    # Calculate grouped statistics
    y_grouped <- tapply(sorted_y, group_assignments, sum)
    m_grouped <- tapply(sorted_y, group_assignments, length)
    probs_grouped <- tapply(sorted_probs, group_assignments, mean)
    
    # Update variables for grouped analysis
    y <- as.numeric(y_grouped)
    m <- as.numeric(m_grouped)
    predicted_probs <- as.numeric(probs_grouped)
    n <- length(y)  # Update n to number of groups
    
    # # Create new design matrix for grouped data
    # # For simplicity, use the mean of original design matrix values within each group
    # original_design <- model.matrix(model)
    # design <- matrix(0, nrow = G, ncol = ncol(original_design))
    # for (i in 1:G) {
    #   group_indices <- which(group_assignments == i)
    #   original_group_indices <- sorted_indices[group_indices]
    #   design[i, ] <- colMeans(original_design[original_group_indices, , drop = FALSE])
    # }
    #design <- model.matrix(model)


    # ===== Ebrahim-Farrington Test =====
    #Calculate the expected value of the Ebrahim-Farrington statistic
    # For equal sample size grouped data, the expected value of the Farrington statistic is:
    # NOTE n = number of groups in case of grouped data
    #ewfarr <- (n - 2 * (1 - 1/n))  #more conservative for small n
    ewfarr <- (n - 2) 
    # Variance of the Farrington statistic
    #varfarr <- 2 * (n - 2) * (1 - 1/n)    #more conservative for small n
    varfarr <- 2 * (n - 2) 

    
    # Input validation
    p2_bounded <- pmax(pmin(predicted_probs, 0.999999), 0.000001)  # Numerical stability


  } else {
    test <- "Original-Farrington (1996) test"

    if (!is.null(model)) {
        stop("You didn't enter the model specification. We will use the model you have already fitted. Example: mymodel <- glm(y ~ x1 + x2 + x3, data = dat, family = binomial())\n and then pass mymodel to the function as the model argument. example: ebrahim.gof(y, predicted_probs, model = mymodel)\n")
    }
    # Use original design matrix if no grouping
    design <- model.matrix(model)
    

    #TRIALS (REPEATITIONS) of each covariate pattern (case)
    if (is.null(m)) {
      # If m is not provided, assume binary data (m=1 for all observations)
      cat("WARNING!! . You didn't enter vector m of the trials (repetitions) for each covariate pattern (case). We will set them deliberately to 1: m = c(1,1,1,...). By theory, this will make the Farrington test never reject H0. WARNING!!\n")
      m <- rep(1, n)
    } else {
      if (length(m) != n) stop("Length mismatch: y and m")
    }
        
    #Calculate the expected value of the Farrington statistic
    p <- ncol(design)
    # For the Farrington test, we need the variance of the linear predictor
    # Input validation
    p2_bounded <- pmax(pmin(predicted_probs, 0.999999), 0.000001)  # Numerical stability
    # This requires the inverse of the information matrix (X'WX)^-1
    w_matrix <- diag(as.vector(variance_terms))
    xpxi <- solve(t(design) %*% w_matrix %*% design)
    varlinpr <- diag(design %*% xpxi %*% t(design))
    # Expected value of the Farrington statistic
    # This term uses `pi*(1-pi)`, not `m*pi*(1-pi)`
    ewfarr <- (n - p) + sum(p2_bounded * (1 - p2_bounded) * varlinpr)
    # Variance of the Farrington statistic
    varfarr <- 2 * (1 - (p / n)) * sum((m - 1) / m)
  }


  


  # ===== Pearson Test =====
  # Pre-calculate common components
  variance_terms <- m * p2_bounded * (1 - p2_bounded)
  # Calculate Pearson statistic (needed for Farrington test)
  pearson_components <- ((y - m * p2_bounded)^2) / variance_terms


  # ===== Farrington Test Chi-Square =====
  farring_corrected_components <- ((1 - 2 * p2_bounded) / variance_terms) * (y - m * p2_bounded)
  farring_stat <- sum(pearson_components) - sum(farring_corrected_components)

 



  # For binary data (m=1), the variance is 0, and the test is not applicable.
  # The original SAS code handles this by setting the p-value to 1.
  if (varfarr <= 1e-8 || is.na(varfarr)) {
    zfarr <- 0
    pzfarr <- 1
    if (all(m == 1) && is.null(G)) {
      warning("Farrington test is not applicable to purely binary data (all m=1). Consider using G parameter to group data or provide m vector.")
    }
  } else {
    zfarr <- (farring_stat - ewfarr) / sqrt(varfarr)
    pzfarr <- 1 - pnorm(zfarr)
  }

  # Return results as a data frame
  results <- data.frame(
    Test = test,
    Test_Statistic = zfarr,
    p_value = pzfarr,
    stringsAsFactors = FALSE
  )

  return(results)
}
