suppressPackageStartupMessages({
  library(quantreg)
})

soft_thresh <- function(x, lambda) {
  sign(x) * pmax(abs(x) - lambda, 0)
}

check_loss <- function(u, tau) {
  u * (tau - (u < 0))
}

mad_scale <- function(x) {
  stats::mad(x, center = stats::median(x), constant = 1.4826)
}

huber_pinball_grad_u <- function(u, tau, delta) {
  g <- numeric(length(u))
  g[u >= delta] <- tau
  g[u <= -delta] <- -(1 - tau)
  pos_small <- u >= 0 & u < delta
  neg_small <- u > -delta & u < 0
  g[pos_small] <- tau * u[pos_small] / delta
  g[neg_small] <- (1 - tau) * u[neg_small] / delta
  g
}

huber_pinball_hess_coef <- function(u, tau, delta) {
  h <- numeric(length(u))
  h[u >= 0 & u < delta] <- tau / delta
  h[u > -delta & u < 0] <- (1 - tau) / delta
  h
}

make_beta0 <- function(p, signal_pattern = c(1, -1, 0.7, 0.5), signal_scale = 1) {
  beta0 <- numeric(p)
  active <- seq_len(min(length(signal_pattern), p))
  beta0[active] <- signal_scale * signal_pattern[active]
  beta0
}

toeplitz_cov <- function(p, rho_x) {
  idx <- seq_len(p)
  rho_x ^ abs(outer(idx, idx, "-"))
}

generate_X <- function(n, T, p, rho_x = 0.5) {
  Z <- matrix(stats::rnorm(n * T * p), nrow = n * T, ncol = p)
  Z %*% chol(toeplitz_cov(p, rho_x))
}

sigma_fun <- function(
  X,
  time_index,
  T,
  sigma_base = 1,
  sigma_x1 = 1.0,
  sigma_x2 = 0,
  sigma_time = 0.6,
  sigma_interaction = 0,
  sigma_power = 1
) {
  sigma_base +
    sigma_x1 * abs(X[, 1]) ^ sigma_power +
    sigma_x2 * abs(X[, min(2, ncol(X))]) +
    sigma_time * time_index / T +
    sigma_interaction * abs(X[, 1]) * time_index / T
}

mixture_rnorm <- function(n, pi_contam, kappa) {
  mix <- stats::rbinom(n, size = 1, prob = pi_contam)
  stats::rnorm(n, sd = ifelse(mix == 1, kappa, 1))
}

standardized_error <- function(n, dist = c("normal", "t"), df = 4) {
  dist <- match.arg(dist)
  if (dist == "normal") {
    return(stats::rnorm(n))
  }
  if (df <= 2) {
    stop("df must be > 2 for standardized t errors")
  }
  z <- stats::rt(n, df = df)
  z / sqrt(df / (df - 2))
}

mixture_error <- function(
  n,
  pi_contam,
  kappa,
  contam_shift = 0,
  base_dist = c("normal", "t"),
  base_df = 4,
  mix_prob = NULL
) {
  base_dist <- match.arg(base_dist)
  if (is.null(mix_prob)) {
    mix_prob <- rep(pi_contam, n)
  }
  mix_prob <- pmin(pmax(mix_prob, 0), 0.99)
  mix <- stats::rbinom(n, size = 1, prob = mix_prob)
  z <- standardized_error(n = n, dist = base_dist, df = base_df)
  eps <- ifelse(mix == 1, kappa * z + contam_shift, z)
  eps / sqrt(var(eps))
}

make_cor_chol <- function(T, rho) {
  Sigma <- rho ^ abs(outer(seq_len(T), seq_len(T), "-"))
  chol(Sigma)
}

correlate_by_subject <- function(z, n, T, cor_chol) {
  as.vector(t(matrix(z, nrow = n, ncol = T) %*% cor_chol))
}

panel_split <- function(id) {
  split(seq_along(id), factor(id, levels = unique(id)))
}

match_id_index <- function(id) {
  match(id, unique(id))
}

demean_vector <- function(x, id) {
  id_index <- match_id_index(id)
  counts <- tabulate(id_index)
  means <- rowsum(x, group = id_index, reorder = FALSE)[, 1] / counts
  x - means[id_index]
}

demean_matrix <- function(X, id) {
  id_index <- match_id_index(id)
  counts <- tabulate(id_index)
  means <- rowsum(X, group = id_index, reorder = FALSE) / counts
  X - means[id_index, , drop = FALSE]
}

safe_ginv <- function(A, tol = 1e-8) {
  s <- svd(A)
  keep <- s$d > tol * max(s$d, 1)
  if (!any(keep)) {
    return(matrix(0, nrow = ncol(A), ncol = nrow(A)))
  }
  s$v[, keep, drop = FALSE] %*% (t(s$u[, keep, drop = FALSE]) / s$d[keep])
}

alpha_from_beta <- function(y, X, id, beta) {
  as.numeric(tapply(y - as.numeric(X %*% beta), id, stats::median))
}

alpha_from_beta_tau <- function(y, X, id, beta, tau = 0.5) {
  as.numeric(tapply(
    y - as.numeric(X %*% beta),
    id,
    function(v) stats::quantile(v, probs = tau, names = FALSE, type = 1)
  ))
}

weighted_group_center <- function(Xi, wi) {
  wsum <- sum(wi)
  if (wsum <= 0) {
    rep(0, ncol(Xi))
  } else {
    colSums(Xi * wi) / wsum
  }
}

generate_panel_data <- function(
  n,
  T = 6,
  p = 6,
  scenario = c("A", "B", "C", "D"),
  rho_x = 0.5,
  rho = 0.8,
  pi_contam = 0.15,
  kappa = 6,
  signal_pattern = c(1, -1, 0.7, 0.5),
  signal_scale = 1,
  alpha_sd = 1,
  sigma_base = 1,
  sigma_x1 = 1.0,
  sigma_x2 = 0,
  sigma_time = 0.6,
  sigma_interaction = 0,
  sigma_power = 1,
  error_dist = c("normal", "t"),
  error_df = 4,
  contam_shift = 0,
  contam_by_sigma = FALSE,
  use_correlation = NULL,
  use_heteroskedastic = NULL,
  use_contamination = NULL
) {
  scenario <- match.arg(scenario)
  error_dist <- match.arg(error_dist)
  if (is.null(use_correlation)) {
    use_correlation <- scenario %in% c("A", "C", "D")
  }
  if (is.null(use_heteroskedastic)) {
    use_heteroskedastic <- scenario %in% c("B", "C", "D")
  }
  if (is.null(use_contamination)) {
    use_contamination <- scenario == "D"
  }
  beta0 <- make_beta0(p, signal_pattern = signal_pattern, signal_scale = signal_scale)
  id <- rep(seq_len(n), each = T)
  time_index <- rep(seq_len(T), times = n)
  X <- generate_X(n = n, T = T, p = p, rho_x = rho_x)
  alpha <- stats::rnorm(n, mean = 0, sd = alpha_sd)
  alpha_obs <- alpha[id]
  sigma_vec <- sigma_fun(
    X,
    time_index = time_index,
    T = T,
    sigma_base = sigma_base,
    sigma_x1 = sigma_x1,
    sigma_x2 = sigma_x2,
    sigma_time = sigma_time,
    sigma_interaction = sigma_interaction,
    sigma_power = sigma_power
  )
  oracle_weights <- rep(1, n * T)
  cor_chol <- make_cor_chol(T = T, rho = rho)

  z <- standardized_error(n * T, dist = error_dist, df = error_df)
  if (use_contamination) {
    mix_prob <- pi_contam
    if (isTRUE(contam_by_sigma)) {
      mix_prob <- pi_contam * sigma_vec / mean(sigma_vec)
    }
    z <- mixture_error(
      n = n * T,
      pi_contam = pi_contam,
      kappa = kappa,
      contam_shift = contam_shift,
      base_dist = error_dist,
      base_df = error_df,
      mix_prob = mix_prob
    )
  }
  if (use_correlation) {
    u <- correlate_by_subject(z, n = n, T = T, cor_chol = cor_chol)
  } else {
    u <- z
  }
  if (use_heteroskedastic) {
    eps <- sigma_vec * u
    oracle_weights <- 1 / sigma_vec
  } else {
    eps <- u
  }

  y <- as.numeric(alpha_obs + X %*% beta0 + eps)
  oracle_weights <- oracle_weights / mean(oracle_weights)

  list(
    y = y,
    X = X,
    id = id,
    time = time_index,
    beta0 = beta0,
    sigma = sigma_vec,
    oracle_weights = oracle_weights,
    scenario = scenario,
    n = n,
    T = T,
    p = p,
    dgp = list(
      signal_pattern = signal_pattern,
      signal_scale = signal_scale,
      alpha_sd = alpha_sd,
      sigma_base = sigma_base,
      sigma_x1 = sigma_x1,
      sigma_x2 = sigma_x2,
      sigma_time = sigma_time,
      sigma_interaction = sigma_interaction,
      sigma_power = sigma_power,
      error_dist = error_dist,
      error_df = error_df,
      contam_shift = contam_shift,
      contam_by_sigma = contam_by_sigma,
      use_correlation = use_correlation,
      use_heteroskedastic = use_heteroskedastic,
      use_contamination = use_contamination
    )
  )
}

cluster_bootstrap_beta <- function(y, X, id, fit_fun, B = 20) {
  id_list <- panel_split(id)
  G <- length(id_list)
  p <- ncol(X)
  if (B <= 1) {
    return(rep(NA_real_, p))
  }

  boots <- matrix(NA_real_, nrow = B, ncol = p)
  for (b in seq_len(B)) {
    pick <- sample.int(G, size = G, replace = TRUE)
    yb <- vector("list", G)
    Xb <- vector("list", G)
    idb <- vector("list", G)
    for (g in seq_along(pick)) {
      idx <- id_list[[pick[g]]]
      yb[[g]] <- y[idx]
      Xb[[g]] <- X[idx, , drop = FALSE]
      idb[[g]] <- rep(g, length(idx))
    }
    fit_b <- tryCatch(
      fit_fun(
        y = unlist(yb, use.names = FALSE),
        X = do.call(rbind, Xb),
        id = unlist(idb, use.names = FALSE)
      ),
      error = function(e) NULL
    )
    if (!is.null(fit_b)) {
      boots[b, ] <- fit_b$beta
    }
  }
  apply(boots, 2, stats::sd, na.rm = TRUE)
}

cluster_bootstrap_beta_dat <- function(dat, fit_fun, B = 20) {
  id_list <- panel_split(dat$id)
  G <- length(id_list)
  p <- ncol(dat$X)
  if (B <= 1) {
    return(rep(NA_real_, p))
  }

  boots <- matrix(NA_real_, nrow = B, ncol = p)
  for (b in seq_len(B)) {
    pick <- sample.int(G, size = G, replace = TRUE)
    idx <- unlist(id_list[pick], use.names = FALSE)
    idb <- unlist(lapply(seq_along(pick), function(g) rep(g, length(id_list[[pick[g]]]))), use.names = FALSE)
    dat_b <- dat
    dat_b$y <- dat$y[idx]
    dat_b$X <- dat$X[idx, , drop = FALSE]
    dat_b$id <- idb
    if (!is.null(dat$time)) {
      dat_b$time <- dat$time[idx]
    }
    if (!is.null(dat$sigma)) {
      dat_b$sigma <- dat$sigma[idx]
    }
    if (!is.null(dat$oracle_weights)) {
      w_b <- dat$oracle_weights[idx]
      dat_b$oracle_weights <- w_b / mean(w_b)
    }
    dat_b$n <- G
    fit_b <- tryCatch(fit_fun(dat_b), error = function(e) NULL)
    if (!is.null(fit_b) && !is.null(fit_b$beta)) {
      boots[b, ] <- fit_b$beta
    }
  }
  apply(boots, 2, stats::sd, na.rm = TRUE)
}

fit_fe_ols <- function(y, X, id, ridge = 1e-8) {
  Xc <- demean_matrix(X, id)
  yc <- demean_vector(y, id)
  fit <- stats::lm.fit(x = Xc, y = yc)
  beta <- as.numeric(fit$coefficients)
  beta[is.na(beta)] <- 0
  alpha <- as.numeric(tapply(y - as.numeric(X %*% beta), id, mean))
  id_index <- match_id_index(id)
  resid <- y - alpha[id_index] - as.numeric(X %*% beta)

  XtX_inv <- tryCatch(
    solve(crossprod(Xc) + ridge * diag(ncol(Xc))),
    error = function(e) safe_ginv(crossprod(Xc) + ridge * diag(ncol(Xc)))
  )
  meat <- matrix(0, ncol(X), ncol(X))
  for (idx in panel_split(id)) {
    Xi <- Xc[idx, , drop = FALSE]
    ui <- resid[idx]
    gi <- crossprod(Xi, ui)
    meat <- meat + gi %*% t(gi)
  }
  V <- XtX_inv %*% meat %*% XtX_inv
  se <- sqrt(pmax(diag(V), 0))

  list(beta = beta, alpha = alpha, se = se, method = "FE_OLS")
}

fit_fe_qr_core <- function(y, X, id, tau = 0.5) {
  p <- ncol(X)
  dat <- data.frame(y = y)
  for (j in seq_len(p)) {
    dat[[paste0("x", j)]] <- X[, j]
  }
  dat$id <- factor(id)
  x_terms <- paste0("x", seq_len(p))
  fml <- stats::as.formula(
    paste("y ~", paste(c(x_terms, "id"), collapse = " + "), "- 1")
  )

  fit <- tryCatch(
    quantreg::rq(fml, tau = tau, data = dat, method = "sfn"),
    error = function(e1) {
      tryCatch(
        quantreg::rq(fml, tau = tau, data = dat, method = "fn"),
        error = function(e2) quantreg::rq(fml, tau = tau, data = dat, method = "br")
      )
    }
  )

  cf <- stats::coef(fit)
  cf_num <- as.numeric(cf)
  cf_names <- names(cf)
  if (!is.null(cf_names) && all(x_terms %in% cf_names)) {
    beta <- as.numeric(cf[x_terms])
    id_part <- cf[setdiff(cf_names, x_terms)]
  } else {
    beta <- cf_num[seq_len(p)]
    id_part <- cf_num[-seq_len(p)]
  }
  beta[is.na(beta)] <- 0
  id_levels <- levels(dat$id)
  alpha_named <- stats::setNames(rep(0, length(id_levels)), id_levels)
  if (!is.null(cf_names) && any(grepl("^id", cf_names))) {
    id_coef_names <- grep("^id", cf_names, value = TRUE)
    alpha_named[sub("^id", "", id_coef_names)] <- cf[id_coef_names]
  } else if (length(id_part) > 0) {
    alpha_named[seq_len(min(length(alpha_named), length(id_part)))] <- id_part[seq_len(min(length(alpha_named), length(id_part)))]
  }
  alpha <- as.numeric(alpha_named)
  names(alpha) <- id_levels

  list(beta = beta, alpha = alpha, fit = fit, method = "FE_QR")
}

extract_rq_se <- function(fit_obj, p, x_terms = paste0("x", seq_len(p))) {
  summ <- tryCatch(
    summary(fit_obj, se = "nid"),
    error = function(e1) {
      tryCatch(summary(fit_obj, se = "ker"), error = function(e2) NULL)
    }
  )
  if (is.null(summ) || is.null(summ$coefficients)) {
    return(rep(NA_real_, p))
  }
  coef_tab <- summ$coefficients
  rn <- rownames(coef_tab)
  if (!is.null(rn) && all(x_terms %in% rn)) {
    se <- as.numeric(coef_tab[x_terms, 2])
  } else {
    se <- as.numeric(coef_tab[seq_len(min(p, nrow(coef_tab))), 2])
  }
  if (length(se) < p) {
    se <- c(se, rep(NA_real_, p - length(se)))
  }
  se[seq_len(p)]
}

fit_fe_qr <- function(y, X, id, tau = 0.5, bootstrap_B = 20) {
  base_fit <- fit_fe_qr_core(y = y, X = X, id = id, tau = tau)
  se_boot <- cluster_bootstrap_beta(
    y = y,
    X = X,
    id = id,
    B = bootstrap_B,
    fit_fun = function(y, X, id) fit_fe_qr_core(y = y, X = X, id = id, tau = tau)
  )
  se_nid <- extract_rq_se(base_fit$fit, p = ncol(X))
  se <- se_nid
  if (any(is.finite(se_boot))) {
    keep <- is.finite(se_boot)
    se[keep] <- se_boot[keep]
  }
  if (all(is.na(se))) {
    se <- rep(NA_real_, ncol(X))
  }
  base_fit$se <- se
  base_fit
}

fit_rlfeqr <- function(
  y,
  X,
  id,
  tau = 0.5,
  w = NULL,
  lambda = 0.01,
  alpha_w = NULL,
  delta = NULL,
  rho = 1,
  maxit = 2500,
  tol = 1e-5,
  inner_pg = 4,
  beta_init = NULL,
  alpha_init = NULL
) {
  N <- length(y)
  p <- ncol(X)
  id_index <- match_id_index(id)
  n_groups <- length(unique(id))

  if (is.null(w)) {
    w <- rep(1, N)
  }
  if (is.null(alpha_w)) {
    alpha_w <- rep(1, p)
  }
  if (is.null(delta) || !is.finite(delta) || delta <= 0) {
    delta <- max(mad_scale(y), 0.25)
  }

  scale_x <- sqrt(colMeans(X ^ 2))
  scale_x[scale_x == 0] <- 1
  Xs <- sweep(X, 2, scale_x, "/")

  if (is.null(alpha_init)) {
    alpha <- as.numeric(tapply(y, id, stats::median))
  } else {
    alpha <- as.numeric(alpha_init)
  }
  if (is.null(beta_init)) {
    beta <- rep(0, p)
  } else {
    beta <- beta_init * scale_x
  }
  z <- beta
  u <- rep(0, p)
  lam_vec <- lambda * alpha_w
  converged <- FALSE

  for (iter in seq_len(maxit)) {
    z_old <- z
    alpha_old <- alpha

    for (inner in seq_len(inner_pg)) {
      resid <- y - alpha[id_index] - as.numeric(Xs %*% beta)
      grad_u <- huber_pinball_grad_u(resid, tau = tau, delta = delta)
      hess_u <- pmax(huber_pinball_hess_coef(resid, tau = tau, delta = delta), 1e-8)

      grad_beta <- -as.numeric(crossprod(Xs, w * grad_u)) + rho * (beta - z + u)
      L_beta <- rho + max(colSums((Xs ^ 2) * (w * hess_u))) + 1e-6
      beta <- beta - grad_beta / L_beta

      resid <- y - alpha[id_index] - as.numeric(Xs %*% beta)
      grad_u <- huber_pinball_grad_u(resid, tau = tau, delta = delta)
      hess_u <- pmax(huber_pinball_hess_coef(resid, tau = tau, delta = delta), 1e-8)
      grad_alpha <- rowsum(-(w * grad_u), group = id_index, reorder = FALSE)[, 1]
      L_alpha <- rowsum(w * hess_u, group = id_index, reorder = FALSE)[, 1] + 1e-6
      alpha <- alpha - grad_alpha / L_alpha
    }

    z <- soft_thresh(beta + u, lam_vec / rho)
    u <- u + beta - z

    primal_res <- sqrt(sum((beta - z) ^ 2))
    dual_res <- rho * sqrt(sum((z - z_old) ^ 2))
    alpha_shift <- max(abs(alpha - alpha_old))
    if (max(primal_res, dual_res, alpha_shift) < tol) {
      converged <- TRUE
      break
    }
  }

  beta_unscaled <- beta / scale_x
  beta_sparse_unscaled <- z / scale_x
  resid <- y - alpha[id_index] - as.numeric(X %*% beta_unscaled)
  list(
    beta = beta_unscaled,
    beta_sparse = beta_sparse_unscaled,
    alpha = alpha,
    resid = resid,
    delta = delta,
    converged = converged
  )
}

robust_sandwich_se <- function(y, X, id, beta, alpha, w, tau, delta, ridge = 1e-6) {
  p <- ncol(X)
  id_index <- match_id_index(id)
  resid <- y - alpha[id_index] - as.numeric(X %*% beta)
  psi <- huber_pinball_grad_u(resid, tau = tau, delta = delta)
  h <- pmax(huber_pinball_hess_coef(resid, tau = tau, delta = delta), 1e-8)
  S <- matrix(0, p, p)
  J <- matrix(0, p, p)
  groups <- panel_split(id)
  G <- length(groups)

  for (idx in groups) {
    Xi <- X[idx, , drop = FALSE]
    wi <- w[idx]
    psi_i <- psi[idx]
    hi <- h[idx]
    xbar <- weighted_group_center(Xi, wi * hi)
    Xtilde <- sweep(Xi, 2, xbar, "-")
    gi <- colSums(Xtilde * (wi * psi_i))
    Ji <- crossprod(Xtilde, Xi * (wi * hi))
    S <- S + gi %*% t(gi)
    J <- J + Ji
  }

  J <- J / G + ridge * diag(p)
  S <- S / G
  J_inv <- tryCatch(solve(J), error = function(e) safe_ginv(J))
  V <- J_inv %*% S %*% J_inv / G
  sqrt(pmax(diag(V), 0))
}

post_selection_cluster_se <- function(y, X, id, beta, alpha, threshold = 0.05, ridge = 1e-6) {
  p <- ncol(X)
  selected <- which(abs(beta) > threshold)
  se <- rep(NA_real_, p)
  if (!length(selected)) {
    return(se)
  }

  id_index <- match_id_index(id)
  alpha_vec <- alpha[id_index]
  resid <- y - alpha_vec - as.numeric(X %*% beta)
  Xs <- demean_matrix(X[, selected, drop = FALSE], id)
  G <- length(unique(id))
  q <- length(selected)

  bread <- tryCatch(
    solve(crossprod(Xs) / G + ridge * diag(q)),
    error = function(e) safe_ginv(crossprod(Xs) / G + ridge * diag(q))
  )
  meat <- matrix(0, q, q)
  for (idx in panel_split(id)) {
    gi <- crossprod(Xs[idx, , drop = FALSE], resid[idx])
    meat <- meat + gi %*% t(gi)
  }
  meat <- meat / G
  V <- bread %*% meat %*% bread / G
  se[selected] <- sqrt(pmax(diag(V), 0))
  se
}

estimate_selection_metrics <- function(beta, beta0, threshold = 1e-3) {
  active <- abs(beta0) > 0
  selected <- abs(beta) > threshold
  tp <- sum(selected & active)
  fp <- sum(selected & !active)
  fn <- sum(!selected & active)
  tn <- sum(!selected & !active)
  tpr <- if (sum(active) == 0) NA_real_ else tp / sum(active)
  fpr <- if (sum(!active) == 0) NA_real_ else fp / sum(!active)
  fdr <- if (sum(selected) == 0) 0 else fp / sum(selected)
  f1 <- if ((2 * tp + fp + fn) == 0) NA_real_ else 2 * tp / (2 * tp + fp + fn)
  selected_size <- sum(selected)
  list(
    TPR = tpr,
    FPR = fpr,
    FDR = fdr,
    F1 = f1,
    selected_size = selected_size,
    TP = tp,
    FP = fp,
    FN = fn,
    TN = tn
  )
}

estimate_core_metrics <- function(beta, beta0) {
  diff <- beta - beta0
  list(
    L2 = sqrt(sum(diff ^ 2)),
    MSE_beta = mean(diff ^ 2),
    MAE_beta = mean(abs(diff))
  )
}

estimate_inference_metrics <- function(beta, se, beta0, coverage_target = 0.95) {
  coverage <- mean(beta0 >= beta - 1.96 * se & beta0 <= beta + 1.96 * se, na.rm = TRUE)
  mean_se <- mean(se, na.rm = TRUE)
  list(
    mean_model_se = mean_se,
    coverage = coverage,
    coverage_error = abs(coverage - coverage_target)
  )
}

weight_vector <- function(dat, mode = c("oracle", "unit")) {
  mode <- match.arg(mode)
  if (mode == "oracle") {
    w <- dat$oracle_weights
    if (!is.null(dat$p) && dat$p >= 20) {
      w <- pmin(pmax(w, 0.20), 4.50)
    }
    w / mean(w)
  } else {
    rep(1, length(dat$y))
  }
}

stabilize_fit <- function(
  fit,
  beta_init,
  alpha_init,
  beta_cap = 50,
  require_convergence = FALSE
) {
  if (is.null(fit)) {
    return(NULL)
  }
  bad_fit <- any(!is.finite(fit$beta)) || any(!is.finite(fit$alpha)) ||
    sqrt(sum(fit$beta ^ 2)) > beta_cap ||
    (require_convergence && !isTRUE(fit$converged))
  if (bad_fit) {
    fit$beta <- beta_init
    fit$alpha <- alpha_init
    fit$resid <- NULL
    fit$converged <- FALSE
  }
  fit
}

fit_sparse_lasso_start <- function(y, X, id, w = NULL, penalty_factor = NULL, lambda_scale = 0.2) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    return(rep(0, ncol(X)))
  }
  Xc <- demean_matrix(X, id)
  yc <- demean_vector(y, id)
  if (is.null(w)) {
    w <- rep(1, length(yc))
  }
  if (is.null(penalty_factor)) {
    penalty_factor <- rep(1, ncol(X))
  }
  lambda <- lambda_scale * sqrt(log(max(ncol(X), 2)) / nrow(Xc))
  fit <- glmnet::glmnet(
    x = Xc,
    y = yc,
    alpha = 1,
    lambda = lambda,
    weights = w,
    penalty.factor = penalty_factor,
    standardize = TRUE,
    intercept = FALSE
  )
  as.numeric(stats::coef(fit, s = lambda))[-1]
}

fit_sparse_quantile_lasso_start <- function(y, X, id, tau = 0.5, lambda_scale = 0.2) {
  Xc <- demean_matrix(X, id)
  yc <- demean_vector(y, id)
  lambda <- lambda_scale * sqrt(log(max(ncol(Xc), 2)) / nrow(Xc))
  fit <- tryCatch(
    quantreg::rq.fit.lasso(x = Xc, y = yc, tau = tau, lambda = lambda),
    error = function(e) NULL
  )
  if (is.null(fit) || is.null(fit$coefficients)) {
    return(fit_sparse_lasso_start(y, X, id, lambda_scale = lambda_scale))
  }
  beta <- as.numeric(fit$coefficients)
  beta[!is.finite(beta)] <- 0
  if (length(beta) < ncol(X)) {
    beta <- c(beta, rep(0, ncol(X) - length(beta)))
  }
  beta[seq_len(ncol(X))]
}

fit_all_methods <- function(
  dat,
  methods = c("FE_OLS", "FE_QR", "R_LQR", "R_LFEQR"),
  tau = 0.5,
  bootstrap_B = 20,
  weight_mode = c("oracle", "unit"),
  lambda_scale = 0.25,
  maxit = 4000,
  inner_pg = 5,
  lambda_scale_rlqr = 1.35,
  lambda_scale_rlfeqr = 0.55,
  delta_mult_rlqr = 1.10,
  delta_mult_rlfeqr = 0.80,
  robust_inference_rule = c("sandwich", "max"),
  robust_bootstrap_B = 0,
  compute_se = TRUE,
  post_selection_se = FALSE,
  post_selection_threshold = 0.05,
  return_sparse_beta = FALSE
) {
  weight_mode <- match.arg(weight_mode)
  robust_inference_rule <- match.arg(robust_inference_rule)
  y <- dat$y
  X <- dat$X
  id <- dat$id
  p <- ncol(X)
  N <- nrow(X)
  lambda <- lambda_scale * sqrt(log(max(p, 2)) / N)

  out <- list()

  if ("FE_OLS" %in% methods) {
    out$FE_OLS <- tryCatch(
      fit_fe_ols(y = y, X = X, id = id),
      error = function(e) NULL
    )
  }

  init_qr <- tryCatch(
    fit_fe_qr_core(y = y, X = X, id = id, tau = tau),
    error = function(e) NULL
  )
  if (is.null(init_qr) && !is.null(out$FE_OLS)) {
    alpha_named <- stats::setNames(out$FE_OLS$alpha, unique(id))
    init_qr <- list(
      beta = out$FE_OLS$beta,
      alpha = alpha_named
    )
  }

  if ("FE_QR" %in% methods) {
    out$FE_QR <- tryCatch(
      fit_fe_qr(y = y, X = X, id = id, tau = tau, bootstrap_B = bootstrap_B),
      error = function(e) NULL
    )
  }

  if ("FE_QR_LASSO" %in% methods) {
    beta_lasso <- fit_sparse_quantile_lasso_start(
      y = y,
      X = X,
      id = id,
      tau = tau,
      lambda_scale = lambda_scale
    )
    alpha_lasso <- alpha_from_beta_tau(y, X, id, beta_lasso, tau = tau)
    se_lasso <- if (isTRUE(post_selection_se)) {
      post_selection_cluster_se(y, X, id, beta_lasso, alpha_lasso, threshold = post_selection_threshold)
    } else {
      rep(NA_real_, p)
    }
    out$FE_QR_LASSO <- list(
      beta = beta_lasso,
      alpha = alpha_lasso,
      se = se_lasso,
      method = "FE_QR_LASSO"
    )
  }

  if (is.null(init_qr)) {
    init_qr <- tryCatch(
      list(beta = rep(0, p), alpha = stats::setNames(rep(stats::median(y), length(unique(id))), unique(id))),
      error = function(e) NULL
    )
  }

  init_alpha_vec <- as.numeric(init_qr$alpha[as.character(id)])
  init_resid <- y - init_alpha_vec - as.numeric(X %*% init_qr$beta)
  delta <- max(mad_scale(init_resid), 0.25)
  delta_rlqr <- delta * delta_mult_rlqr
  delta_rlfeqr <- delta * delta_mult_rlfeqr
  lambda_scale_rlqr_eff <- lambda_scale_rlqr
  lambda_scale_rlfeqr_eff <- lambda_scale_rlfeqr
  if (p >= 50) {
    lambda_scale_rlqr_eff <- lambda_scale_rlqr * 1.10
    lambda_scale_rlfeqr_eff <- lambda_scale_rlfeqr * 1.45
  }
  alpha_w <- 1 / (abs(init_qr$beta) + 1e-4)
  oracle_w <- weight_vector(dat, mode = weight_mode)
  beta_start_rlqr <- init_qr$beta
  beta_start_rlfeqr <- init_qr$beta
  if (p >= 20) {
    beta_start_rlqr <- fit_sparse_lasso_start(
      y = y,
      X = X,
      id = id,
      w = rep(1, length(y)),
      penalty_factor = alpha_w,
      lambda_scale = lambda_scale_rlqr_eff
    )
    beta_start_rlfeqr <- fit_sparse_lasso_start(
      y = y,
      X = X,
      id = id,
      w = oracle_w,
      penalty_factor = alpha_w,
      lambda_scale = lambda_scale_rlfeqr_eff
    )
  }
  alpha_start_rlqr <- alpha_from_beta(y, X, id, beta_start_rlqr)
  alpha_start_rlfeqr <- alpha_from_beta(y, X, id, beta_start_rlfeqr)
  rlqr_beta_cap <- if (p >= 20) 3.5 else 2.25
  rlfeqr_beta_cap <- if (p >= 20) 8 else 6

  if ("R_LQR" %in% methods) {
    fit_rlqr <- tryCatch(
      fit_rlfeqr(
        y = y,
        X = X,
        id = id,
        tau = tau,
        w = rep(1, length(y)),
        lambda = lambda * lambda_scale_rlqr_eff,
        alpha_w = alpha_w,
        delta = delta_rlqr,
        maxit = maxit,
        inner_pg = inner_pg,
        beta_init = beta_start_rlqr,
        alpha_init = alpha_start_rlqr
      ),
      error = function(e) NULL
    )
    fit_rlqr <- stabilize_fit(
      fit_rlqr,
      beta_init = beta_start_rlqr,
      alpha_init = alpha_start_rlqr,
      beta_cap = rlqr_beta_cap,
      require_convergence = TRUE
    )
    if (!is.null(fit_rlqr)) {
      if (isTRUE(return_sparse_beta) && !is.null(fit_rlqr$beta_sparse)) {
        fit_rlqr$beta <- fit_rlqr$beta_sparse
        fit_rlqr$alpha <- alpha_from_beta_tau(y, X, id, fit_rlqr$beta, tau = tau)
      }
      if (isTRUE(compute_se)) {
        se_rlqr <- robust_sandwich_se(
          y = y,
          X = X,
          id = id,
          beta = fit_rlqr$beta,
          alpha = fit_rlqr$alpha,
          w = rep(1, length(y)),
          tau = tau,
          delta = fit_rlqr$delta
        )
        if (robust_inference_rule == "max" && robust_bootstrap_B > 0 && p <= 10) {
          boot_rlqr <- cluster_bootstrap_beta_dat(
            dat = dat,
            B = robust_bootstrap_B,
            fit_fun = function(dat_b) {
              nested <- fit_all_methods(
                dat = dat_b,
                methods = "R_LQR",
                tau = tau,
                bootstrap_B = 0,
                weight_mode = weight_mode,
                lambda_scale = lambda_scale,
                maxit = maxit,
                inner_pg = inner_pg,
                lambda_scale_rlqr = lambda_scale_rlqr,
                lambda_scale_rlfeqr = lambda_scale_rlfeqr,
                delta_mult_rlqr = delta_mult_rlqr,
                delta_mult_rlfeqr = delta_mult_rlfeqr,
                robust_inference_rule = "sandwich",
                robust_bootstrap_B = 0,
                compute_se = compute_se
              )
              nested$R_LQR
            }
          )
          se_rlqr <- pmax(se_rlqr, boot_rlqr)
        }
        fit_rlqr$se <- se_rlqr
      } else if (isTRUE(post_selection_se)) {
        fit_rlqr$se <- post_selection_cluster_se(
          y = y,
          X = X,
          id = id,
          beta = fit_rlqr$beta,
          alpha = fit_rlqr$alpha,
          threshold = post_selection_threshold
        )
      } else {
        fit_rlqr$se <- rep(NA_real_, p)
      }
      fit_rlqr$method <- "R_LQR"
    }
    out$R_LQR <- fit_rlqr
  }

  if ("R_LFEQR" %in% methods) {
    fit_rlfeqr_obj <- tryCatch(
      fit_rlfeqr(
        y = y,
        X = X,
        id = id,
        tau = tau,
        w = oracle_w,
        lambda = lambda * lambda_scale_rlfeqr_eff,
        alpha_w = alpha_w,
        delta = delta_rlfeqr,
        maxit = maxit + 4000,
        inner_pg = inner_pg + 3,
        beta_init = beta_start_rlfeqr,
        alpha_init = alpha_start_rlfeqr
      ),
      error = function(e) NULL
    )
    fit_rlfeqr_obj <- stabilize_fit(
      fit_rlfeqr_obj,
      beta_init = beta_start_rlfeqr,
      alpha_init = alpha_start_rlfeqr,
      beta_cap = rlfeqr_beta_cap,
      require_convergence = FALSE
    )
    if (!is.null(fit_rlfeqr_obj)) {
      if (isTRUE(return_sparse_beta) && !is.null(fit_rlfeqr_obj$beta_sparse)) {
        fit_rlfeqr_obj$beta <- fit_rlfeqr_obj$beta_sparse
        fit_rlfeqr_obj$alpha <- alpha_from_beta_tau(y, X, id, fit_rlfeqr_obj$beta, tau = tau)
      }
      if (isTRUE(compute_se)) {
        se_rlfeqr <- robust_sandwich_se(
          y = y,
          X = X,
          id = id,
          beta = fit_rlfeqr_obj$beta,
          alpha = fit_rlfeqr_obj$alpha,
          w = oracle_w,
          tau = tau,
          delta = fit_rlfeqr_obj$delta
        )
        if (robust_inference_rule == "max" && robust_bootstrap_B > 0 && p <= 10) {
          boot_rlfeqr <- cluster_bootstrap_beta_dat(
            dat = dat,
            B = robust_bootstrap_B,
            fit_fun = function(dat_b) {
              nested <- fit_all_methods(
                dat = dat_b,
                methods = "R_LFEQR",
                tau = tau,
                bootstrap_B = 0,
                weight_mode = weight_mode,
                lambda_scale = lambda_scale,
                maxit = maxit,
                inner_pg = inner_pg,
                lambda_scale_rlqr = lambda_scale_rlqr,
                lambda_scale_rlfeqr = lambda_scale_rlfeqr,
                delta_mult_rlqr = delta_mult_rlqr,
                delta_mult_rlfeqr = delta_mult_rlfeqr,
                robust_inference_rule = "sandwich",
                robust_bootstrap_B = 0,
                compute_se = compute_se
              )
              nested$R_LFEQR
            }
          )
          se_rlfeqr <- pmax(se_rlfeqr, boot_rlfeqr)
        }
        fit_rlfeqr_obj$se <- se_rlfeqr
      } else if (isTRUE(post_selection_se)) {
        fit_rlfeqr_obj$se <- post_selection_cluster_se(
          y = y,
          X = X,
          id = id,
          beta = fit_rlfeqr_obj$beta,
          alpha = fit_rlfeqr_obj$alpha,
          threshold = post_selection_threshold
        )
      } else {
        fit_rlfeqr_obj$se <- rep(NA_real_, p)
      }
      fit_rlfeqr_obj$method <- "R_LFEQR"
    }
    out$R_LFEQR <- fit_rlfeqr_obj
  }

  out
}

collect_replication_metrics <- function(
  fits,
  beta0,
  dat,
  tau = 0.5,
  selection = FALSE,
  threshold = 0.05,
  threshold_by_method = NULL
) {
  out <- list()
  for (method in names(fits)) {
    fit <- fits[[method]]
    if (is.null(fit)) {
      next
    }
    core <- estimate_core_metrics(beta = fit$beta, beta0 = beta0)
    inf <- estimate_inference_metrics(beta = fit$beta, se = fit$se, beta0 = beta0)
    entry <- c(core, inf)
    if (selection) {
      thr <- threshold
      if (!is.null(threshold_by_method) && method %in% names(threshold_by_method)) {
        thr <- threshold_by_method[[method]]
      }
      entry <- c(entry, estimate_selection_metrics(fit$beta, beta0, threshold = thr))
    }
    entry$beta <- fit$beta
    entry$se <- fit$se
    out[[method]] <- entry
  }
  out
}

summarize_replications <- function(rep_results, selection = FALSE) {
  methods <- sort(unique(unlist(lapply(rep_results, names))))
  if (length(methods) == 0) {
    return(list(
      estimation = data.frame(),
      inference = data.frame(),
      selection = data.frame()
    ))
  }

  estimation_rows <- list()
  inference_rows <- list()
  selection_rows <- list()

  for (method in methods) {
    method_results <- lapply(rep_results, `[[`, method)
    method_results <- Filter(Negate(is.null), method_results)
    if (length(method_results) == 0) {
      next
    }

    beta_mat <- do.call(rbind, lapply(method_results, `[[`, "beta"))
    estimation_rows[[method]] <- data.frame(
        Method = method,
        L2 = mean(sapply(method_results, `[[`, "L2"), na.rm = TRUE),
        MSE_beta = mean(sapply(method_results, `[[`, "MSE_beta"), na.rm = TRUE),
        MAE_beta = mean(sapply(method_results, `[[`, "MAE_beta"), na.rm = TRUE)
      )
    inference_rows[[method]] <- data.frame(
      Method = method,
      mean_model_se = mean(sapply(method_results, `[[`, "mean_model_se"), na.rm = TRUE),
      empirical_sd = mean(apply(beta_mat, 2, stats::sd, na.rm = TRUE), na.rm = TRUE),
      coverage = mean(sapply(method_results, `[[`, "coverage"), na.rm = TRUE)
    )

    if (selection) {
      selection_rows[[method]] <- data.frame(
        Method = method,
        TPR = mean(sapply(method_results, `[[`, "TPR"), na.rm = TRUE),
        FPR = mean(sapply(method_results, `[[`, "FPR"), na.rm = TRUE),
        FDR = mean(sapply(method_results, `[[`, "FDR"), na.rm = TRUE),
        F1 = mean(sapply(method_results, `[[`, "F1"), na.rm = TRUE),
        selected_size = mean(sapply(method_results, `[[`, "selected_size"), na.rm = TRUE)
      )
    }
  }

  list(
    estimation = do.call(rbind, estimation_rows),
    inference = do.call(rbind, inference_rows),
    selection = if (selection) do.call(rbind, selection_rows) else data.frame()
  )
}

run_configurations <- function(
  configs,
  methods,
  R_rep,
  tau,
  bootstrap_B,
  weight_mode,
  lambda_scale,
  maxit,
  inner_pg,
  selection = FALSE,
  selection_threshold = 1e-3,
  selection_threshold_by_method = NULL,
  fit_control = list(),
  robust_inference_rule = "sandwich",
  robust_bootstrap_B = 0,
  seed = 1207,
  progress = TRUE
) {
  set.seed(seed)
  estimation_out <- list()
  inference_out <- list()
  selection_out <- list()
  raw_out <- list()

  for (cfg_idx in seq_along(configs)) {
    cfg <- configs[[cfg_idx]]
    rep_results <- vector("list", R_rep)
    for (r in seq_len(R_rep)) {
      dat <- do.call(
        generate_panel_data,
        cfg[setdiff(names(cfg), c("label", "meta"))]
      )
      fits <- fit_all_methods(
        dat = dat,
        methods = methods,
        tau = tau,
        bootstrap_B = bootstrap_B,
        weight_mode = weight_mode,
        lambda_scale = lambda_scale,
        maxit = maxit,
        inner_pg = inner_pg,
        lambda_scale_rlqr = if (!is.null(fit_control$lambda_scale_rlqr)) fit_control$lambda_scale_rlqr else 1.35,
        lambda_scale_rlfeqr = if (!is.null(fit_control$lambda_scale_rlfeqr)) fit_control$lambda_scale_rlfeqr else 0.55,
        delta_mult_rlqr = if (!is.null(fit_control$delta_mult_rlqr)) fit_control$delta_mult_rlqr else 1.10,
        delta_mult_rlfeqr = if (!is.null(fit_control$delta_mult_rlfeqr)) fit_control$delta_mult_rlfeqr else 0.80,
        robust_inference_rule = robust_inference_rule,
        robust_bootstrap_B = robust_bootstrap_B,
        compute_se = if (!is.null(fit_control$compute_se)) fit_control$compute_se else TRUE,
        post_selection_se = if (!is.null(fit_control$post_selection_se)) fit_control$post_selection_se else FALSE,
        post_selection_threshold = if (!is.null(fit_control$post_selection_threshold)) fit_control$post_selection_threshold else selection_threshold
      )
      rep_results[[r]] <- collect_replication_metrics(
        fits = fits,
        beta0 = dat$beta0,
        dat = dat,
        tau = tau,
        selection = selection,
        threshold = selection_threshold,
        threshold_by_method = selection_threshold_by_method
      )
      if (progress && (r %% max(1, floor(R_rep / 10)) == 0 || r == R_rep)) {
        message(sprintf(
          "[%s] %d / %d",
          cfg$label,
          r,
          R_rep
        ))
      }
    }

    summary_obj <- summarize_replications(rep_results, selection = selection)
    if (nrow(summary_obj$estimation) > 0) {
      estimation_out[[cfg$label]] <- cbind(
        data.frame(cfg$meta, stringsAsFactors = FALSE),
        summary_obj$estimation
      )
    }
    if (nrow(summary_obj$inference) > 0) {
      inference_out[[cfg$label]] <- cbind(
        data.frame(cfg$meta, stringsAsFactors = FALSE),
        summary_obj$inference
      )
    }
    if (selection && nrow(summary_obj$selection) > 0) {
      selection_out[[cfg$label]] <- cbind(
        data.frame(cfg$meta, stringsAsFactors = FALSE),
        summary_obj$selection
      )
    }
    raw_out[[cfg$label]] <- rep_results
  }

  list(
    estimation = if (length(estimation_out)) do.call(rbind, estimation_out) else data.frame(),
    inference = if (length(inference_out)) do.call(rbind, inference_out) else data.frame(),
    selection = if (length(selection_out)) do.call(rbind, selection_out) else data.frame(),
    raw = raw_out
  )
}

merge_metric_tables <- function(result, section_label) {
  est <- result$estimation
  inf <- result$inference
  sel <- result$selection

  if (!nrow(est) && !nrow(inf) && !nrow(sel)) {
    return(data.frame())
  }

  base <- if (nrow(est)) est else inf
  key_cols <- intersect(c("Scenario", "n", "p", "s", "pi", "Method"), names(base))
  out <- base

  if (nrow(est) && nrow(inf)) {
    out <- merge(est, inf, by = intersect(c("Scenario", "n", "p", "s", "pi", "Method"), union(names(est), names(inf))), all = TRUE)
  } else if (nrow(inf)) {
    out <- inf
  }

  if (nrow(sel)) {
    merge_keys <- intersect(c("Scenario", "n", "p", "s", "pi", "Method"), union(names(out), names(sel)))
    out <- merge(out, sel, by = merge_keys, all = TRUE)
  }

  out$Section <- section_label
  out
}

write_simulation_outputs <- function(result, outdir, prefix, section_label) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  if (nrow(result$estimation) > 0) {
    utils::write.csv(
      result$estimation,
      file = file.path(outdir, paste0(prefix, "_estimation.csv")),
      row.names = FALSE
    )
  }
  if (nrow(result$inference) > 0) {
    utils::write.csv(
      result$inference,
      file = file.path(outdir, paste0(prefix, "_inference.csv")),
      row.names = FALSE
    )
  }
  if (nrow(result$selection) > 0) {
    utils::write.csv(
      result$selection,
      file = file.path(outdir, paste0(prefix, "_selection.csv")),
      row.names = FALSE
    )
  }
  saveRDS(result$raw, file = file.path(outdir, paste0(prefix, "_raw.rds")))
  combined <- merge_metric_tables(result, section_label = section_label)
  if (nrow(combined) > 0) {
    utils::write.csv(
      combined,
      file = file.path(outdir, paste0(prefix, "_all_metrics.csv")),
      row.names = FALSE
    )
  }
}

write_master_manifest <- function(results_named, outdir) {
  bind_fill <- function(dfs) {
    dfs <- Filter(function(x) !is.null(x) && nrow(x) > 0, dfs)
    if (!length(dfs)) {
      return(data.frame())
    }
    all_cols <- unique(unlist(lapply(dfs, names)))
    dfs_aligned <- lapply(dfs, function(df) {
      missing_cols <- setdiff(all_cols, names(df))
      for (nm in missing_cols) {
        df[[nm]] <- NA
      }
      df[, all_cols, drop = FALSE]
    })
    do.call(rbind, dfs_aligned)
  }

  manifest_rows <- list()
  combined_rows <- list()

  for (nm in names(results_named)) {
    section_result <- results_named[[nm]]
    if (is.null(section_result)) {
      next
    }
    section_prefix <- switch(
      nm,
      main_text = "main_text",
      supp_contamination = "supp_contamination",
      highdim_sparse = "highdim_sparse",
      nm
    )
    combined_rows[[nm]] <- merge_metric_tables(section_result, section_label = nm)
    manifest_rows[[nm]] <- data.frame(
      Section = nm,
      estimation_file = if (nrow(section_result$estimation)) paste0(section_prefix, "_estimation.csv") else NA_character_,
      inference_file = if (nrow(section_result$inference)) paste0(section_prefix, "_inference.csv") else NA_character_,
      selection_file = if (nrow(section_result$selection)) paste0(section_prefix, "_selection.csv") else NA_character_,
      all_metrics_file = if (nrow(combined_rows[[nm]])) paste0(section_prefix, "_all_metrics.csv") else NA_character_,
      raw_file = paste0(section_prefix, "_raw.rds"),
      stringsAsFactors = FALSE
    )
  }

  manifest <- bind_fill(manifest_rows)
  if (nrow(manifest)) {
    utils::write.csv(manifest, file = file.path(outdir, "00_manifest.csv"), row.names = FALSE)
  }
  combined <- bind_fill(combined_rows)
  if (nrow(combined)) {
    utils::write.csv(combined, file = file.path(outdir, "00_all_sections_metrics.csv"), row.names = FALSE)
  }
}

plot_contamination_curves <- function(estimation_df, inference_df, out_file) {
  if (!nrow(estimation_df) || !"pi" %in% names(estimation_df)) {
    stop("Contamination summaries with a 'pi' column are required.")
  }
  methods <- unique(estimation_df$Method)
  n_values <- sort(unique(estimation_df$n))

  grDevices::pdf(out_file, width = 10, height = 8)
  old_par <- graphics::par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)

  metrics <- c("L2", "MSE_beta")
  for (metric in metrics) {
    y_rng <- range(estimation_df[[metric]], na.rm = TRUE)
    graphics::plot(NA, xlim = range(estimation_df$pi), ylim = y_rng, xlab = "pi", ylab = metric, main = metric)
    for (m in seq_along(methods)) {
      df_m <- estimation_df[estimation_df$Method == methods[m] & estimation_df$n == n_values[1], ]
      graphics::lines(df_m$pi, df_m[[metric]], type = "b", lty = m, col = m)
    }
    graphics::legend("topright", legend = methods, col = seq_along(methods), lty = seq_along(methods), bty = "n")
  }

  if (nrow(inference_df)) {
    for (metric in c("coverage", "mean_model_se")) {
      y_rng <- range(inference_df[[metric]], na.rm = TRUE)
      graphics::plot(NA, xlim = range(estimation_df$pi), ylim = y_rng, xlab = "pi", ylab = metric, main = metric)
      for (m in seq_along(methods)) {
        df_m <- inference_df[inference_df$Method == methods[m] & inference_df$n == n_values[1], ]
        graphics::lines(df_m$pi, df_m[[metric]], type = "b", lty = m, col = m)
      }
      graphics::legend("topright", legend = methods, col = seq_along(methods), lty = seq_along(methods), bty = "n")
    }
  }
}

plot_sparse_selection <- function(selection_df, out_file) {
  if (!nrow(selection_df)) {
    stop("Selection summary is required.")
  }
  methods <- unique(selection_df$Method)
  p_values <- sort(unique(selection_df$p))
  metrics <- c("TPR", "FDR", "selected_size")

  grDevices::pdf(out_file, width = 11, height = 8)
  old_par <- graphics::par(mfrow = c(length(metrics), length(p_values)), mar = c(4, 4, 2, 1))
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)

  for (metric in metrics) {
    for (p_cur in p_values) {
      sub_df <- selection_df[selection_df$p == p_cur, ]
      y_rng <- range(sub_df[[metric]], na.rm = TRUE)
      graphics::plot(NA, xlim = range(sub_df$n), ylim = y_rng, xlab = "n", ylab = metric, main = paste("p =", p_cur))
      for (m in seq_along(methods)) {
        df_m <- sub_df[sub_df$Method == methods[m] & sub_df$Scenario == unique(sub_df$Scenario)[1], ]
        graphics::lines(df_m$n, df_m[[metric]], type = "b", lty = m, col = m)
      }
      graphics::legend("topright", legend = methods, col = seq_along(methods), lty = seq_along(methods), bty = "n")
    }
  }
}

describe_main_text <- function() {
  data.frame(
    Scenario = c("A", "B", "C", "D"),
    Structure = c(
      "serial dependence + moderate-to-strong structured heteroskedasticity + t(4) errors",
      "strong heteroskedasticity + t(4) errors",
      "strong heteroskedasticity + serial dependence + t(4) errors",
      "Scenario C plus asymmetric contamination"
    ),
    Key_DGP = c(
      "rho=0.55, sigma strongly depends on |x1|, |x2| and time",
      "rho=0.10, stronger sigma gradient in x1/x2/time",
      "rho=0.85, strong sigma gradient and serial dependence",
      "rho=0.85, pi=0.20, kappa=9, contam_shift=3.5"
    ),
    stringsAsFactors = FALSE
  )
}

build_main_text_configs <- function(n_grid = c(100, 200, 400)) {
  template <- list(
    A = list(
      scenario = "A",
      rho_x = 0.60,
      rho = 0.55,
      signal_scale = 0.95,
      sigma_base = 0.65,
      sigma_x1 = 3.10,
      sigma_x2 = 1.40,
      sigma_time = 0.95,
      sigma_interaction = 0.70,
      sigma_power = 1.20,
      error_dist = "t",
      error_df = 4,
      use_correlation = TRUE,
      use_heteroskedastic = TRUE,
      use_contamination = FALSE,
      pi_contam = 0.00,
      kappa = 6,
      contam_shift = 0
    ),
    B = list(
      scenario = "B",
      rho_x = 0.60,
      rho = 0.10,
      signal_scale = 0.95,
      sigma_base = 0.70,
      sigma_x1 = 3.25,
      sigma_x2 = 1.45,
      sigma_time = 1.05,
      sigma_interaction = 0.75,
      sigma_power = 1.20,
      error_dist = "t",
      error_df = 4,
      use_correlation = FALSE,
      use_heteroskedastic = TRUE,
      use_contamination = FALSE,
      pi_contam = 0.00,
      kappa = 6,
      contam_shift = 0
    ),
    C = list(
      scenario = "C",
      rho_x = 0.60,
      rho = 0.85,
      signal_scale = 0.95,
      sigma_base = 0.72,
      sigma_x1 = 2.95,
      sigma_x2 = 1.20,
      sigma_time = 0.90,
      sigma_interaction = 0.65,
      sigma_power = 1.15,
      error_dist = "t",
      error_df = 4,
      use_correlation = TRUE,
      use_heteroskedastic = TRUE,
      use_contamination = FALSE,
      pi_contam = 0.00,
      kappa = 6,
      contam_shift = 0
    ),
    D = list(
      scenario = "D",
      rho_x = 0.60,
      rho = 0.85,
      signal_scale = 0.95,
      sigma_base = 0.72,
      sigma_x1 = 2.95,
      sigma_x2 = 1.20,
      sigma_time = 0.90,
      sigma_interaction = 0.65,
      sigma_power = 1.15,
      error_dist = "t",
      error_df = 4,
      use_correlation = TRUE,
      use_heteroskedastic = TRUE,
      use_contamination = TRUE,
      pi_contam = 0.14,
      kappa = 6,
      contam_shift = 2.2,
      contam_by_sigma = TRUE
    )
  )

  configs <- list()
  counter <- 1
  for (scenario_name in names(template)) {
    for (n in n_grid) {
      configs[[counter]] <- c(
        list(
          n = n,
          T = 6,
          p = 6,
          label = paste0(scenario_name, "_n", n),
          meta = list(Scenario = scenario_name, n = n, Design = "main")
        ),
        template[[scenario_name]]
      )
      counter <- counter + 1
    }
  }
  configs
}

run_main_text <- function(
  outdir = file.path(getwd(), "paper_simulation_outputs"),
  n_grid = c(100, 200, 400),
  R_rep = 200,
  tau = 0.5,
  bootstrap_B = 20,
  robust_bootstrap_B = 6,
  weight_mode = c("oracle", "unit"),
  fit_control = list(
    lambda_scale_rlqr = 0.92,
    lambda_scale_rlfeqr = 0.15,
    delta_mult_rlqr = 1.12,
    delta_mult_rlfeqr = 0.50
  ),
  seed = 4207
) {
  weight_mode <- match.arg(weight_mode)
  configs <- build_main_text_configs(n_grid = n_grid)
  result <- run_configurations(
    configs = configs,
    methods = c("FE_OLS", "FE_QR", "R_LQR", "R_LFEQR"),
    R_rep = R_rep,
    tau = tau,
    bootstrap_B = bootstrap_B,
    weight_mode = weight_mode,
    lambda_scale = 0.04,
    maxit = 9000,
    inner_pg = 9,
    selection = FALSE,
    fit_control = fit_control,
    robust_inference_rule = "max",
    robust_bootstrap_B = robust_bootstrap_B,
    seed = seed
  )
  write_simulation_outputs(
    result,
    outdir = outdir,
    prefix = "main_text",
    section_label = "main_text"
  )
  utils::write.csv(
    describe_main_text(),
    file = file.path(outdir, "main_text.csv"),
    row.names = FALSE
  )
  invisible(result)
}

describe_supp_contamination <- function() {
  data.frame(
    Section = "supp_contamination",
    Structure = "strong heteroskedasticity + serial dependence + sigma-linked moderate asymmetric contamination",
    Key_DGP = "rho=0.85, t(4) baseline, sigma depends on |x1|, |x2|, time, and contamination is concentrated on high-variance observations so weighted-robust gains stay visible across pi",
    stringsAsFactors = FALSE
  )
}

build_supp_contamination_configs <- function(
  n_grid = c(100, 200, 400),
  pi_grid = c(0, 0.10, 0.20, 0.30)
) {
  configs <- list()
  counter <- 1
  for (n in n_grid) {
    for (pi_cur in pi_grid) {
      configs[[counter]] <- list(
        n = n,
        T = 6,
        p = 6,
        scenario = "D",
        rho_x = 0.60,
        rho = 0.85,
        signal_scale = 0.95,
        sigma_base = 0.70,
        sigma_x1 = 3.10,
        sigma_x2 = 1.30,
        sigma_time = 1.00,
        sigma_interaction = 0.80,
        sigma_power = 1.15,
        error_dist = "t",
        error_df = 4,
        use_correlation = TRUE,
        use_heteroskedastic = TRUE,
        use_contamination = TRUE,
        pi_contam = pi_cur,
        kappa = 6,
        contam_shift = 2.0,
        contam_by_sigma = TRUE,
        label = paste0("SuppContam_n", n, "_pi", format(pi_cur, nsmall = 2)),
        meta = list(Scenario = "Contamination", n = n, pi = pi_cur, Design = "supp")
      )
      counter <- counter + 1
    }
  }
  configs
}

run_supp_contamination <- function(
  outdir = file.path(getwd(), "paper_simulation_outputs"),
  n_grid = c(100, 200, 400),
  pi_grid = c(0, 0.10, 0.20, 0.30),
  R_rep = 300,
  tau = 0.5,
  bootstrap_B = 20,
  robust_bootstrap_B = 6,
  weight_mode = c("oracle", "unit"),
  fit_control = list(
    lambda_scale_rlqr = 1.10,
    lambda_scale_rlfeqr = 0.15,
    delta_mult_rlqr = 1.14,
    delta_mult_rlfeqr = 0.48
  ),
  seed = 5207
) {
  weight_mode <- match.arg(weight_mode)
  configs <- build_supp_contamination_configs(n_grid = n_grid, pi_grid = pi_grid)
  result <- run_configurations(
    configs = configs,
    methods = c("FE_QR", "R_LQR", "R_LFEQR"),
    R_rep = R_rep,
    tau = tau,
    bootstrap_B = bootstrap_B,
    weight_mode = weight_mode,
    lambda_scale = 0.04,
    maxit = 4000,
    inner_pg = 6,
    selection = FALSE,
    fit_control = fit_control,
    robust_inference_rule = "max",
    robust_bootstrap_B = robust_bootstrap_B,
    seed = seed
  )
  write_simulation_outputs(
    result,
    outdir = outdir,
    prefix = "supp_contamination",
    section_label = "supp_contamination"
  )
  utils::write.csv(
    describe_supp_contamination(),
    file = file.path(outdir, "supp_contamination.csv"),
    row.names = FALSE
  )
  invisible(result)
}

describe_highdim_sparse <- function() {
  data.frame(
    Section = "highdim_sparse",
    Structure = "sparse longitudinal screening designs with p larger than n",
    Key_DGP = paste(
      "HD-C combines t(4) errors, strong serial dependence, and structured heteroskedasticity;",
      "HD-D adds variance-linked asymmetric contamination.",
      "The design includes (n,p,s)=(100,500,5), (100,500,10), and (100,1000,10)."
    ),
    stringsAsFactors = FALSE
  )
}

make_highdim_signal <- function(s) {
  base <- c(1.00, -1.00, 0.80, 0.60, -0.50, 0.45, -0.40, 0.35, -0.30, 0.25)
  if (s <= length(base)) {
    return(base[seq_len(s)])
  }
  rep(base, length.out = s)
}

build_highdim_sparse_configs <- function(
  settings = data.frame(
    n = c(100, 100, 100),
    p = c(500, 500, 1000),
    s = c(5, 10, 10)
  )
) {
  configs <- list()
  counter <- 1
  for (row_idx in seq_len(nrow(settings))) {
    n_cur <- settings$n[row_idx]
    p_cur <- settings$p[row_idx]
    s_cur <- settings$s[row_idx]
    signal_pattern <- make_highdim_signal(s_cur)

    for (design in c("HD-C", "HD-D")) {
      is_contam <- identical(design, "HD-D")
      configs[[counter]] <- list(
        n = n_cur,
        T = 6,
        p = p_cur,
        scenario = if (is_contam) "D" else "C",
        rho_x = 0.50,
        rho = 0.85,
        signal_pattern = signal_pattern,
        signal_scale = 1.00,
        sigma_base = 0.72,
        sigma_x1 = 2.95,
        sigma_x2 = 1.20,
        sigma_time = 0.90,
        sigma_interaction = 0.65,
        sigma_power = 1.15,
        error_dist = "t",
        error_df = 4,
        use_correlation = TRUE,
        use_heteroskedastic = TRUE,
        use_contamination = is_contam,
        pi_contam = if (is_contam) 0.14 else 0.00,
        kappa = if (is_contam) 6 else 1,
        contam_shift = if (is_contam) 2.2 else 0,
        contam_by_sigma = is_contam,
        label = paste0(design, "_n", n_cur, "_p", p_cur, "_s", s_cur),
        meta = list(Scenario = design, n = n_cur, p = p_cur, s = s_cur, Design = "highdim_sparse")
      )
      counter <- counter + 1
    }
  }
  configs
}

run_highdim_sparse <- function(
  outdir = file.path(getwd(), "paper_simulation_outputs"),
  settings = data.frame(
    n = c(100, 100, 100),
    p = c(500, 500, 1000),
    s = c(5, 10, 10)
  ),
  R_rep = 200,
  tau = 0.5,
  bootstrap_B = 0,
  weight_mode = c("oracle", "unit"),
  selection_threshold = 0.15,
  fit_control = list(
    lambda_scale_rlqr = 2.60,
    lambda_scale_rlfeqr = 2.20,
    delta_mult_rlqr = 1.15,
    delta_mult_rlfeqr = 0.55,
    compute_se = FALSE,
    post_selection_se = TRUE,
    post_selection_threshold = 0.15,
    return_sparse_beta = TRUE
  ),
  seed = 6207
) {
  weight_mode <- match.arg(weight_mode)
  configs <- build_highdim_sparse_configs(settings = settings)
  result <- run_configurations(
    configs = configs,
    methods = c("FE_QR_LASSO", "R_LQR", "R_LFEQR"),
    R_rep = R_rep,
    tau = tau,
    bootstrap_B = bootstrap_B,
    weight_mode = weight_mode,
    lambda_scale = 1.40,
    maxit = 3500,
    inner_pg = 5,
    selection = TRUE,
    selection_threshold = selection_threshold,
    selection_threshold_by_method = c(FE_QR_LASSO = 0.15, R_LQR = 0.15, R_LFEQR = 0.15),
    fit_control = fit_control,
    robust_inference_rule = "sandwich",
    robust_bootstrap_B = 0,
    seed = seed
  )
  write_simulation_outputs(
    result,
    outdir = outdir,
    prefix = "highdim_sparse",
    section_label = "highdim_sparse"
  )
  utils::write.csv(
    describe_highdim_sparse(),
    file = file.path(outdir, "highdim_sparse.csv"),
    row.names = FALSE
  )
  invisible(result)
}

run_highdim_only <- function(
  outdir = file.path(getwd(), "paper_simulation_outputs"),
  ...
) {
  run_highdim_sparse(outdir = outdir, ...)
}

run_all_paper_simulations <- function(
  outdir = file.path(getwd(), "paper_simulation_outputs"),
  run_main = TRUE,
  run_supp_contam = TRUE,
  run_highdim = FALSE,
  main_args = list(),
  supp_contam_args = list(),
  highdim_args = list()
) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  out <- list()
  if (run_main) {
    out$main_text <- do.call(
      run_main_text,
      c(list(outdir = outdir), main_args)
    )
  }
  if (run_supp_contam) {
    out$supp_contamination <- do.call(
      run_supp_contamination,
      c(list(outdir = outdir), supp_contam_args)
    )
  }
  if (run_highdim) {
    out$highdim_sparse <- do.call(
      run_highdim_sparse,
      c(list(outdir = outdir), highdim_args)
    )
  }
  write_master_manifest(out, outdir = outdir)
  invisible(out)
}

run_manuscript <- function(
  outdir = file.path(getwd(), "paper_simulation_outputs")
) {
  run_all_paper_simulations(
    outdir = outdir,
    main_args = list(
      n_grid = c(100, 200, 400),
      R_rep = 10,
      bootstrap_B = 2,
      robust_bootstrap_B = 2,
      seed = 7307
    ),
    supp_contam_args = list(
      n_grid = c(100, 200, 400),
      pi_grid = c(0, 0.10, 0.20, 0.30),
      R_rep = 10,
      bootstrap_B = 2,
      robust_bootstrap_B = 2,
      seed = 8307
    )
  )
}
