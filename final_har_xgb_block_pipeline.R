## ================================================================
## FINAL HAR + XGBoost + Block-Conformal Pipeline (GitHub-ready)
## Asset universe: XOM, EFA, XLB, FDX, BTC-USD, BRK-B
## Period: 2015-01-01 to 2024-12-31
## Output folder: outml_final/
##
## Main empirical outputs
##   T0  Descriptive statistics
##   T1  RMSE comparison
##   T2  Coverage comparison
##   T3  Band-width comparison
##   T4  Tail-region performance
##   T5  XGBoost tuning summary
##   T6  Block-length sensitivity
##   F1  Test-period variance proxy and forecasts
##   F2  Standard vs block-conformal bands
##   F4  Calibration nonconformity score histograms
##
## Simulation outputs
##   S1  Simulation coverage table
##   S2  Simulation width table
##   F5  Simulation trade-off figure
##
## Optional appendix robustness
##   T8  Alternative variance-proxy robustness
## ================================================================

## ================================================================
## 0) Packages
## ================================================================
pkgs <- c(
  "quantmod", "xts", "zoo",
  "dplyr", "tibble", "purrr", "tidyr",
  "yardstick",
  "ggplot2", "readr",
  "xgboost",
  "moments",
  "knitr"
)

to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

set.seed(123)
options(stringsAsFactors = FALSE)
options(knitr.kable.NA = "--")

## ================================================================
## 1) Global settings
## ================================================================
OUT_DIR <- "outml_final"
dir.create(OUT_DIR, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "plots"),  showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "rds"),    showWarnings = FALSE, recursive = TRUE)

ASSETS_MAIN <- c("XOM", "EFA", "XLB", "FDX", "BTC-USD", "BRK-B")

## ================================================================
## 2) Utilities
## ================================================================
detect_nthread <- function() {
  cores <- parallel::detectCores(logical = FALSE)
  if (is.na(cores) || cores <= 1) return(1L)
  max(1L, as.integer(cores - 1L))
}

asset_label <- function(symbol) {
  dplyr::case_when(
    symbol == "XOM"     ~ "XOM",
    symbol == "EFA"     ~ "EFA",
    symbol == "XLB"     ~ "XLB",
    symbol == "FDX"     ~ "FDX",
    symbol == "BTC-USD" ~ "BTC/USD",
    symbol == "BRK-B"   ~ "BRK/B",
    TRUE ~ symbol
  )
}

pretty_model_label <- function(x) {
  dplyr::case_when(
    x == "HAR"         ~ "HAR",
    x == "XGBoost"     ~ "XGBoost",
    x == "HAR+XGBoost" ~ "HAR+XGBoost",
    TRUE ~ x
  )
}

pretty_proxy_label <- function(x) {
  dplyr::case_when(
    x == "squared_return" ~ "Squared-return proxy",
    x == "parkinson"      ~ "Parkinson proxy",
    x == "gk"             ~ "Garman-Klass proxy",
    x == "rs"             ~ "Rogers-Satchell proxy",
    TRUE ~ x
  )
}

base_clean_theme <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title       = ggplot2::element_blank(),
      legend.title     = ggplot2::element_text(face = "bold"),
      legend.position  = "bottom"
    )
}

safe_write_csv <- function(df, path) {
  readr::write_csv(df, path)
}

save_latex_table <- function(df, file, caption, label,
                             longtable = FALSE,
                             align = NULL,
                             escape = TRUE) {
  tex_lines <- capture.output(
    knitr::kable(
      df,
      format    = "latex",
      booktabs  = TRUE,
      longtable = longtable,
      caption   = caption,
      label     = label,
      align     = align,
      escape    = escape
    )
  )
  writeLines(tex_lines, con = file, useBytes = TRUE)
}

conformal_quantile <- function(scores, alpha = 0.05) {
  scores <- sort(scores[is.finite(scores)])
  n <- length(scores)
  if (n == 0L) stop("No finite scores available for conformal quantile.")
  k <- ceiling((n + 1) * (1 - alpha))
  k <- min(max(k, 1L), n)
  scores[k]
}

make_time_split <- function(n, train_frac = 0.60, calib_frac = 0.20) {
  n_train <- floor(train_frac * n)
  n_calib <- floor(calib_frac * n)
  list(
    train = 1:n_train,
    calib = (n_train + 1):(n_train + n_calib),
    test  = (n_train + n_calib + 1):n
  )
}

## ================================================================
## 3) Data retrieval
## ================================================================
get_price_data <- function(symbol,
                           from = as.Date("2015-01-01"),
                           to   = as.Date("2024-12-31")) {
  cat("==> Downloading:", symbol, "\n")
  
  x <- tryCatch(
    quantmod::getSymbols(
      Symbols     = symbol,
      src         = "yahoo",
      from        = from,
      to          = to,
      auto.assign = FALSE
    ),
    error = function(e) {
      stop("Data download failed for ", symbol, ": ", e$message)
    }
  )
  
  adj <- tryCatch(
    quantmod::Ad(x),
    error = function(e) quantmod::Cl(x)
  )
  
  ret <- quantmod::dailyReturn(adj, type = "log")
  idx <- index(ret)
  
  out <- tibble::tibble(
    date     = as.Date(idx),
    open     = as.numeric(quantmod::Op(x)[idx]),
    high     = as.numeric(quantmod::Hi(x)[idx]),
    low      = as.numeric(quantmod::Lo(x)[idx]),
    close    = as.numeric(quantmod::Cl(x)[idx]),
    adjusted = as.numeric(adj[idx]),
    r        = as.numeric(ret)
  ) %>%
    dplyr::arrange(date)
  
  if (nrow(out) == 0L) stop("No rows returned for ", symbol, ".")
  out
}

## ================================================================
## 4) Variance-proxy construction
## ================================================================
add_variance_proxy <- function(df,
                               proxy_type = c("squared_return", "parkinson", "gk", "rs"),
                               eps_proxy = 1e-12,
                               eps_log = 1e-8) {
  proxy_type <- match.arg(proxy_type)
  
  open_  <- pmax(df$open,  eps_proxy)
  high_  <- pmax(df$high,  eps_proxy)
  low_   <- pmax(df$low,   eps_proxy)
  close_ <- pmax(df$close, eps_proxy)
  
  hl <- log(high_ / low_)
  co <- log(close_ / open_)
  
  raw_proxy <- dplyr::case_when(
    proxy_type == "squared_return" ~ df$r^2,
    proxy_type == "parkinson"      ~ (hl^2) / (4 * log(2)),
    proxy_type == "gk"             ~ 0.5 * hl^2 - (2 * log(2) - 1) * co^2,
    proxy_type == "rs"             ~ log(high_ / close_) * log(high_ / open_) +
      log(low_ / close_) * log(low_ / open_),
    TRUE ~ df$r^2
  )
  
  df %>%
    dplyr::mutate(
      proxy_type         = proxy_type,
      variance_proxy     = pmax(raw_proxy, eps_proxy),
      log_variance_proxy = log(variance_proxy + eps_log)
    )
}

compute_descriptives <- function(df_vp, symbol) {
  tibble::tibble(
    symbol  = symbol,
    n_obs   = nrow(df_vp),
    mean_r  = mean(df_vp$r, na.rm = TRUE),
    sd_r    = sd(df_vp$r, na.rm = TRUE),
    skew_r  = moments::skewness(df_vp$r, na.rm = TRUE),
    kurt_r  = moments::kurtosis(df_vp$r, na.rm = TRUE),
    mean_vp = mean(df_vp$variance_proxy, na.rm = TRUE),
    sd_vp   = sd(df_vp$variance_proxy, na.rm = TRUE),
    skew_vp = moments::skewness(df_vp$variance_proxy, na.rm = TRUE),
    kurt_vp = moments::kurtosis(df_vp$variance_proxy, na.rm = TRUE)
  )
}

## ================================================================
## 5) Feature engineering
## ================================================================
make_features <- function(df, vp_col = "variance_proxy", p_lag = 5L) {
  vp <- df[[vp_col]]
  
  vp_d <- vp
  vp_w <- zoo::rollapply(vp, 5, mean, align = "right", fill = NA)
  vp_m <- zoo::rollapply(vp, 22, mean, align = "right", fill = NA)
  
  abs_r    <- abs(df$r)
  r2       <- df$r^2
  r3       <- df$r^3
  roll_sd5 <- zoo::rollapply(df$r, 5, sd, align = "right", fill = NA)
  
  log_vp <- df$log_variance_proxy
  lag_list <- vector("list", p_lag + 1L)
  names(lag_list) <- paste0("log_vp_l", 0:p_lag)
  
  for (k in 0:p_lag) {
    lag_list[[paste0("log_vp_l", k)]] <- dplyr::lag(log_vp, k)
  }
  
  lag_df <- tibble::as_tibble(lag_list)
  
  out <- df %>%
    dplyr::bind_cols(lag_df) %>%
    dplyr::mutate(
      vp_d          = vp_d,
      vp_w          = vp_w,
      vp_m          = vp_m,
      log_vp_d      = log(vp_d + 1e-8),
      log_vp_w      = log(vp_w + 1e-8),
      log_vp_m      = log(vp_m + 1e-8),
      abs_r         = abs_r,
      r2            = r2,
      r3            = r3,
      roll_sd5      = roll_sd5,
      y_log_vp_next = dplyr::lead(log_vp, 1)
    ) %>%
    dplyr::filter(
      !is.na(vp_d),
      !is.na(vp_w),
      !is.na(vp_m),
      !is.na(roll_sd5),
      !is.na(y_log_vp_next),
      dplyr::if_all(dplyr::starts_with("log_vp_l"), ~ !is.na(.x))
    )
  
  if (nrow(out) == 0L) stop("Feature engineering produced zero usable rows.")
  out
}

## ================================================================
## 6) HAR model
## ================================================================
fit_har <- function(df_feat, train_idx) {
  form <- y_log_vp_next ~ log_vp_d + log_vp_w + log_vp_m
  fit  <- stats::lm(form, data = df_feat[train_idx, , drop = FALSE])
  df_feat$har_pred <- as.numeric(stats::predict(fit, newdata = df_feat))
  list(model = fit, df_feat = df_feat)
}

## ================================================================
## 7) XGBoost: time-ordered tuning + final fit
## ================================================================
build_xgb_weights <- function(abs_r_vec, tail_lambda = 0) {
  if (tail_lambda <= 0) return(rep(1, length(abs_r_vec)))
  q90 <- stats::quantile(abs_r_vec, 0.9, na.rm = TRUE)
  if (is.na(q90) || q90 == 0) return(rep(1, length(abs_r_vec)))
  ratio <- pmin(abs_r_vec / q90, 3)
  1 + tail_lambda * (ratio^2)
}

fit_xgb_direct_ts <- function(df_feat,
                              train_idx,
                              feature_cols,
                              param_grid = NULL,
                              inner_valid_frac = 0.20,
                              nrounds_max = 800L,
                              early_stopping_rounds = 25L,
                              tail_lambda = 0) {
  if (is.null(param_grid)) {
    param_grid <- tidyr::crossing(
      eta              = c(0.03, 0.05, 0.10),
      max_depth        = c(2L, 4L, 6L),
      subsample        = c(0.80, 1.00),
      colsample_bytree = c(0.80),
      lambda           = c(1.0, 2.0)
    )
  }
  
  n_train <- length(train_idx)
  n_inner_train <- floor((1 - inner_valid_frac) * n_train)
  
  if (n_inner_train < 50L || n_inner_train >= n_train) {
    stop("Inner validation split is not feasible. Check training size and inner_valid_frac.")
  }
  
  inner_train_idx <- train_idx[1:n_inner_train]
  inner_valid_idx <- train_idx[(n_inner_train + 1):n_train]
  
  X_inner_train <- as.matrix(df_feat[inner_train_idx, feature_cols, drop = FALSE])
  y_inner_train <- df_feat$y_log_vp_next[inner_train_idx]
  X_inner_valid <- as.matrix(df_feat[inner_valid_idx, feature_cols, drop = FALSE])
  y_inner_valid <- df_feat$y_log_vp_next[inner_valid_idx]
  
  if (any(!is.finite(X_inner_train)) || any(!is.finite(X_inner_valid))) {
    stop("Non-finite values detected in XGBoost features.")
  }
  if (any(!is.finite(y_inner_train)) || any(!is.finite(y_inner_valid))) {
    stop("Non-finite values detected in XGBoost targets.")
  }
  
  w_inner_train <- build_xgb_weights(abs(df_feat$r[inner_train_idx]), tail_lambda)
  w_full_train  <- build_xgb_weights(abs(df_feat$r[train_idx]), tail_lambda)
  
  dtrain_inner <- xgboost::xgb.DMatrix(
    data   = X_inner_train,
    label  = y_inner_train,
    weight = w_inner_train
  )
  dvalid_inner <- xgboost::xgb.DMatrix(
    data  = X_inner_valid,
    label = y_inner_valid
  )
  
  tuning_results <- purrr::pmap_dfr(
    param_grid,
    function(eta, max_depth, subsample, colsample_bytree, lambda) {
      params <- xgboost::xgb.params(
        objective        = "reg:squarederror",
        eval_metric      = "rmse",
        learning_rate    = eta,
        max_depth        = as.integer(max_depth),
        subsample        = subsample,
        colsample_bytree = colsample_bytree,
        lambda           = lambda,
        alpha            = 0,
        nthread          = detect_nthread(),
        seed             = 123
      )
      
      bst <- tryCatch(
        xgboost::xgb.train(
          params                = params,
          data                  = dtrain_inner,
          nrounds               = as.integer(nrounds_max),
          evals                 = list(train = dtrain_inner, valid = dvalid_inner),
          early_stopping_rounds = as.integer(early_stopping_rounds),
          verbose               = 0
        ),
        error = function(e) NULL
      )
      
      if (is.null(bst)) {
        return(tibble::tibble(
          eta              = eta,
          max_depth        = as.integer(max_depth),
          subsample        = subsample,
          colsample_bytree = colsample_bytree,
          lambda           = lambda,
          best_nrounds     = NA_integer_,
          validation_rmse  = NA_real_
        ))
      }
      
      best_nrounds <- bst$best_iteration
      if (is.null(best_nrounds) || length(best_nrounds) == 0L ||
          !is.finite(best_nrounds) || is.na(best_nrounds) || best_nrounds < 1) {
        best_nrounds <- as.integer(nrounds_max)
      } else {
        best_nrounds <- as.integer(best_nrounds)
      }
      
      pred_valid <- as.numeric(stats::predict(bst, X_inner_valid))
      val_rmse_manual <- sqrt(mean((y_inner_valid - pred_valid)^2))
      
      best_score <- bst$best_score
      if (is.null(best_score) || length(best_score) == 0L ||
          !is.finite(best_score) || is.na(best_score)) {
        best_score <- val_rmse_manual
      } else {
        best_score <- as.numeric(best_score)
      }
      
      tibble::tibble(
        eta              = eta,
        max_depth        = as.integer(max_depth),
        subsample        = subsample,
        colsample_bytree = colsample_bytree,
        lambda           = lambda,
        best_nrounds     = best_nrounds,
        validation_rmse  = best_score
      )
    }
  )
  
  tuning_results <- tuning_results %>%
    dplyr::filter(is.finite(best_nrounds), best_nrounds >= 1L, is.finite(validation_rmse)) %>%
    dplyr::arrange(validation_rmse, eta, max_depth, subsample, colsample_bytree, lambda)
  
  if (nrow(tuning_results) == 0L) {
    stop("All XGBoost tuning configurations failed.")
  }
  
  best_cfg <- tuning_results[1, ]
  
  X_full_train <- as.matrix(df_feat[train_idx, feature_cols, drop = FALSE])
  y_full_train <- df_feat$y_log_vp_next[train_idx]
  
  if (any(!is.finite(X_full_train)) || any(!is.finite(y_full_train))) {
    stop("Non-finite values detected in full XGBoost training data.")
  }
  
  dtrain_full <- xgboost::xgb.DMatrix(
    data   = X_full_train,
    label  = y_full_train,
    weight = w_full_train
  )
  
  final_params <- xgboost::xgb.params(
    objective        = "reg:squarederror",
    eval_metric      = "rmse",
    learning_rate    = best_cfg$eta,
    max_depth        = as.integer(best_cfg$max_depth),
    subsample        = best_cfg$subsample,
    colsample_bytree = best_cfg$colsample_bytree,
    lambda           = best_cfg$lambda,
    alpha            = 0,
    nthread          = detect_nthread(),
    seed             = 123
  )
  
  bst_final <- xgboost::xgb.train(
    params  = final_params,
    data    = dtrain_full,
    nrounds = as.integer(best_cfg$best_nrounds),
    verbose = 0
  )
  
  X_all <- as.matrix(df_feat[, feature_cols, drop = FALSE])
  df_feat$xgb_pred <- as.numeric(stats::predict(bst_final, X_all))
  
  list(
    model          = bst_final,
    df_feat        = df_feat,
    best_cfg       = best_cfg,
    tuning_results = tuning_results
  )
}

## ================================================================
## 8) Stacking
## ================================================================
stack_weight <- function(y, har_pred, xgb_pred) {
  f_obj <- function(w) {
    mean((y - (w * har_pred + (1 - w) * xgb_pred))^2, na.rm = TRUE)
  }
  stats::optimize(f_obj, interval = c(0, 1))$minimum
}

## ================================================================
## 9) Conformal calibration
## ================================================================
add_conformal_intervals_block <- function(df_feat, calib_idx, test_idx,
                                          alpha = 0.05, block_len = 5L) {
  eps_calib     <- df_feat$y_log_vp_next[calib_idx] - df_feat$hyb_pred[calib_idx]
  abs_eps_calib <- abs(eps_calib)
  
  qhat_iid <- conformal_quantile(abs_eps_calib, alpha = alpha)
  
  n_calib <- length(calib_idx)
  B_calib <- floor(n_calib / block_len)
  if (B_calib < 1L) stop("Calibration set too small for selected block length.")
  
  block_scores <- numeric(B_calib)
  for (b in seq_len(B_calib)) {
    idx_block <- calib_idx[((b - 1) * block_len + 1):(b * block_len)]
    block_eps <- df_feat$y_log_vp_next[idx_block] - df_feat$hyb_pred[idx_block]
    block_scores[b] <- max(abs(block_eps), na.rm = TRUE)
  }
  
  qhat_block <- conformal_quantile(block_scores, alpha = alpha)
  
  df_feat <- df_feat %>%
    dplyr::mutate(
      conf_lo_log_iid = hyb_pred - qhat_iid,
      conf_hi_log_iid = hyb_pred + qhat_iid,
      conf_lo_vp_iid  = exp(conf_lo_log_iid),
      conf_hi_vp_iid  = exp(conf_hi_log_iid),
      conf_lo_log_blk = hyb_pred - qhat_block,
      conf_hi_log_blk = hyb_pred + qhat_block,
      conf_lo_vp_blk  = exp(conf_lo_log_blk),
      conf_hi_vp_blk  = exp(conf_hi_log_blk)
    )
  
  y_test <- df_feat$y_log_vp_next[test_idx]
  
  cov_point_iid <- mean(
    y_test >= df_feat$conf_lo_log_iid[test_idx] &
      y_test <= df_feat$conf_hi_log_iid[test_idx],
    na.rm = TRUE
  )
  
  cov_point_blk <- mean(
    y_test >= df_feat$conf_lo_log_blk[test_idx] &
      y_test <= df_feat$conf_hi_log_blk[test_idx],
    na.rm = TRUE
  )
  
  n_test <- length(test_idx)
  B_test <- floor(n_test / block_len)
  
  if (B_test > 0L) {
    block_cov_iid <- logical(B_test)
    block_cov_blk <- logical(B_test)
    
    for (b in seq_len(B_test)) {
      idx_block <- test_idx[((b - 1) * block_len + 1):(b * block_len)]
      
      block_cov_iid[b] <- all(
        df_feat$y_log_vp_next[idx_block] >= df_feat$conf_lo_log_iid[idx_block] &
          df_feat$y_log_vp_next[idx_block] <= df_feat$conf_hi_log_iid[idx_block]
      )
      
      block_cov_blk[b] <- all(
        df_feat$y_log_vp_next[idx_block] >= df_feat$conf_lo_log_blk[idx_block] &
          df_feat$y_log_vp_next[idx_block] <= df_feat$conf_hi_log_blk[idx_block]
      )
    }
    
    cov_block_iid <- mean(block_cov_iid)
    cov_block_blk <- mean(block_cov_blk)
  } else {
    cov_block_iid <- NA_real_
    cov_block_blk <- NA_real_
  }
  
  list(
    df_feat            = df_feat,
    qhat_iid           = qhat_iid,
    qhat_block         = qhat_block,
    coverage_iid       = cov_point_iid,
    coverage_blk_point = cov_point_blk,
    coverage_block_iid = cov_block_iid,
    coverage_block_blk = cov_block_blk,
    scores_iid         = abs_eps_calib,
    scores_block       = block_scores
  )
}

compute_block_sensitivity <- function(df_feat, calib_idx, test_idx,
                                      alpha = 0.05,
                                      block_grid = c(3L, 5L, 10L)) {
  purrr::map_dfr(block_grid, function(BL) {
    conf_res <- add_conformal_intervals_block(
      df_feat   = df_feat,
      calib_idx = calib_idx,
      test_idx  = test_idx,
      alpha     = alpha,
      block_len = as.integer(BL)
    )
    df_tmp <- conf_res$df_feat
    
    width_iid <- df_tmp$conf_hi_vp_iid[test_idx] - df_tmp$conf_lo_vp_iid[test_idx]
    width_blk <- df_tmp$conf_hi_vp_blk[test_idx] - df_tmp$conf_lo_vp_blk[test_idx]
    
    tibble::tibble(
      block_len        = as.integer(BL),
      cov_point_iid    = conf_res$coverage_iid,
      cov_point_block  = conf_res$coverage_blk_point,
      cov_block_iid    = conf_res$coverage_block_iid,
      cov_block_block  = conf_res$coverage_block_blk,
      mean_width_iid   = mean(width_iid, na.rm = TRUE),
      mean_width_block = mean(width_blk, na.rm = TRUE)
    )
  })
}

## ================================================================
## 10) Evaluation helpers
## ================================================================
compute_rmse_summary <- function(df_feat, test_idx) {
  y_true <- df_feat$y_log_vp_next[test_idx]
  tibble::tibble(
    model = c("HAR", "XGBoost", "HAR+XGBoost"),
    RMSE  = c(
      yardstick::rmse_vec(truth = y_true, estimate = df_feat$har_pred[test_idx]),
      yardstick::rmse_vec(truth = y_true, estimate = df_feat$xgb_pred[test_idx]),
      yardstick::rmse_vec(truth = y_true, estimate = df_feat$hyb_pred[test_idx])
    )
  )
}

compute_bandwidth_summary <- function(df_feat, test_idx) {
  width_iid <- df_feat$conf_hi_vp_iid[test_idx] - df_feat$conf_lo_vp_iid[test_idx]
  width_blk <- df_feat$conf_hi_vp_blk[test_idx] - df_feat$conf_lo_vp_blk[test_idx]
  
  tibble::tibble(
    band_type    = c("Standard split-conformal", "Block-conformal"),
    mean_width   = c(mean(width_iid, na.rm = TRUE), mean(width_blk, na.rm = TRUE)),
    median_width = c(median(width_iid, na.rm = TRUE), median(width_blk, na.rm = TRUE))
  )
}

compute_tail_performance <- function(df_feat, test_idx, tail_prob = 0.90) {
  vp_test  <- exp(df_feat$y_log_vp_next[test_idx])
  thr      <- stats::quantile(vp_test, tail_prob, na.rm = TRUE)
  pos_tail <- which(vp_test >= thr)
  
  if (length(pos_tail) < 5L) {
    return(
      tibble::tibble(
        tail_prob            = tail_prob,
        n_tail               = length(pos_tail),
        RMSE_HAR_tail        = NA_real_,
        RMSE_XGB_tail        = NA_real_,
        RMSE_HYB_tail        = NA_real_,
        cov_iid_tail         = NA_real_,
        cov_block_tail       = NA_real_,
        mean_width_iid_tail  = NA_real_,
        mean_width_blk_tail  = NA_real_
      )
    )
  }
  
  idx_tail <- test_idx[pos_tail]
  y_tail   <- df_feat$y_log_vp_next[idx_tail]
  
  rmse_har_tail <- yardstick::rmse_vec(truth = y_tail, estimate = df_feat$har_pred[idx_tail])
  rmse_xgb_tail <- yardstick::rmse_vec(truth = y_tail, estimate = df_feat$xgb_pred[idx_tail])
  rmse_hyb_tail <- yardstick::rmse_vec(truth = y_tail, estimate = df_feat$hyb_pred[idx_tail])
  
  cov_iid_tail <- mean(
    y_tail >= df_feat$conf_lo_log_iid[idx_tail] &
      y_tail <= df_feat$conf_hi_log_iid[idx_tail],
    na.rm = TRUE
  )
  
  cov_block_tail <- mean(
    y_tail >= df_feat$conf_lo_log_blk[idx_tail] &
      y_tail <= df_feat$conf_hi_log_blk[idx_tail],
    na.rm = TRUE
  )
  
  width_iid_tail <- df_feat$conf_hi_vp_iid[idx_tail] - df_feat$conf_lo_vp_iid[idx_tail]
  width_blk_tail <- df_feat$conf_hi_vp_blk[idx_tail] - df_feat$conf_lo_vp_blk[idx_tail]
  
  tibble::tibble(
    tail_prob            = tail_prob,
    n_tail               = length(idx_tail),
    RMSE_HAR_tail        = rmse_har_tail,
    RMSE_XGB_tail        = rmse_xgb_tail,
    RMSE_HYB_tail        = rmse_hyb_tail,
    cov_iid_tail         = cov_iid_tail,
    cov_block_tail       = cov_block_tail,
    mean_width_iid_tail  = mean(width_iid_tail, na.rm = TRUE),
    mean_width_blk_tail  = mean(width_blk_tail, na.rm = TRUE)
  )
}

## ================================================================
## 11) Plotting
## ================================================================
plot_test_fit <- function(df_feat, test_idx, symbol, out_dir) {
  df_plot <- df_feat[test_idx, ] %>%
    dplyr::transmute(
      date,
      `Observed daily variance proxy` = exp(y_log_vp_next),
      `HAR forecast`                  = exp(har_pred),
      `XGBoost forecast`              = exp(xgb_pred),
      `Hybrid forecast`               = exp(hyb_pred)
    ) %>%
    tidyr::pivot_longer(-date, names_to = "series", values_to = "value")
  
  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = date, y = value, colour = series)) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::labs(
      x      = "",
      y      = "Daily variance proxy",
      colour = "Series"
    ) +
    base_clean_theme()
  
  ggplot2::ggsave(
    filename = file.path(
      out_dir, "plots",
      paste0("F1_test_fit_", gsub("[^A-Za-z0-9]", "_", symbol), ".png")
    ),
    plot = p, width = 8, height = 3.8, dpi = 600
  )
}

plot_conf_band_iid_vs_block <- function(df_feat, test_idx, symbol, out_dir, last_n = 250) {
  sel_idx <- if (length(test_idx) <= last_n) test_idx else tail(test_idx, last_n)
  
  df_plot <- df_feat[sel_idx, ] %>%
    dplyr::mutate(
      observed_vp = exp(y_log_vp_next),
      hybrid_pred = exp(hyb_pred)
    )
  
  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = date)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = conf_lo_vp_iid, ymax = conf_hi_vp_iid,
                   fill = "Standard split-conformal"),
      alpha = 0.18
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = conf_lo_vp_blk, ymax = conf_hi_vp_blk,
                   fill = "Block-conformal"),
      alpha = 0.18
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = hybrid_pred, colour = "Hybrid forecast"),
      linewidth = 0.7
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = observed_vp, colour = "Observed daily variance proxy"),
      linewidth = 0.5
    ) +
    ggplot2::labs(
      x      = "",
      y      = "Daily variance proxy",
      fill   = "Band type",
      colour = "Series"
    ) +
    base_clean_theme()
  
  ggplot2::ggsave(
    filename = file.path(
      out_dir, "plots",
      paste0("F2_conformal_bands_", gsub("[^A-Za-z0-9]", "_", symbol), ".png")
    ),
    plot = p, width = 8, height = 3.8, dpi = 600
  )
}

plot_score_hist <- function(scores_iid, scores_block, symbol, out_dir) {
  df_scores <- tibble::tibble(
    score = c(scores_iid, scores_block),
    type  = c(
      rep("Absolute residuals", length(scores_iid)),
      rep("Block maxima of absolute residuals", length(scores_block))
    )
  )
  
  p <- ggplot2::ggplot(df_scores, ggplot2::aes(x = score, fill = type)) +
    ggplot2::geom_histogram(alpha = 0.50, position = "identity", bins = 30) +
    ggplot2::labs(
      x    = "Nonconformity score",
      y    = "Count",
      fill = "Score type"
    ) +
    base_clean_theme()
  
  ggplot2::ggsave(
    filename = file.path(
      out_dir, "plots",
      paste0("F4_score_hist_", gsub("[^A-Za-z0-9]", "_", symbol), ".png")
    ),
    plot = p, width = 6.8, height = 4.2, dpi = 600
  )
}

## ================================================================
## 12) Per-asset pipeline
## ================================================================
run_for_asset <- function(symbol,
                          from           = as.Date("2015-01-01"),
                          to             = as.Date("2024-12-31"),
                          proxy_type     = "squared_return",
                          alpha          = 0.05,
                          block_len_main = 5L,
                          block_len_grid = c(3L, 5L, 10L),
                          p_lag          = 5L,
                          tail_lambda    = 0,
                          out_dir        = OUT_DIR) {
  df_raw <- get_price_data(symbol, from = from, to = to)
  df_vp  <- add_variance_proxy(df_raw, proxy_type = proxy_type)
  desc_tbl <- compute_descriptives(df_vp, symbol)
  
  df_feat <- make_features(df_vp, p_lag = p_lag)
  n <- nrow(df_feat)
  
  if (n < 300L) {
    warning("Too few observations for ", symbol, " (n = ", n, ").")
    return(NULL)
  }
  
  idx <- make_time_split(n)
  train_idx <- idx$train
  calib_idx <- idx$calib
  test_idx  <- idx$test
  
  har_res <- fit_har(df_feat, train_idx)
  df_feat <- har_res$df_feat
  
  feat_cols_xgb <- c(
    "log_vp_d", "log_vp_w", "log_vp_m",
    "abs_r", "r2", "r3", "roll_sd5",
    paste0("log_vp_l", 0:p_lag)
  )
  
  if (any(!is.finite(as.matrix(df_feat[, feat_cols_xgb, drop = FALSE])))) {
    stop("Non-finite values detected in feat_cols_xgb before XGBoost fitting.")
  }
  if (any(!is.finite(df_feat$y_log_vp_next))) {
    stop("Non-finite values detected in y_log_vp_next before XGBoost fitting.")
  }
  
  xgb_res <- fit_xgb_direct_ts(
    df_feat      = df_feat,
    train_idx    = train_idx,
    feature_cols = feat_cols_xgb,
    tail_lambda  = tail_lambda
  )
  df_feat <- xgb_res$df_feat
  
  valid_calib <- calib_idx[
    !is.na(df_feat$har_pred[calib_idx]) &
      !is.na(df_feat$xgb_pred[calib_idx])
  ]
  
  y_c   <- df_feat$y_log_vp_next[valid_calib]
  har_c <- df_feat$har_pred[valid_calib]
  xgb_c <- df_feat$xgb_pred[valid_calib]
  
  w_hat <- stack_weight(y_c, har_c, xgb_c)
  df_feat$hyb_pred <- w_hat * df_feat$har_pred + (1 - w_hat) * df_feat$xgb_pred
  
  conf_res <- add_conformal_intervals_block(
    df_feat   = df_feat,
    calib_idx = calib_idx,
    test_idx  = test_idx,
    alpha     = alpha,
    block_len = block_len_main
  )
  df_feat <- conf_res$df_feat
  
  rmse_tbl <- compute_rmse_summary(df_feat, test_idx) %>%
    dplyr::mutate(symbol = symbol)
  
  coverage_tbl <- tibble::tibble(
    symbol           = symbol,
    alpha            = alpha,
    block_len        = block_len_main,
    cov_point_iid    = conf_res$coverage_iid,
    cov_point_block  = conf_res$coverage_blk_point,
    cov_block_iid    = conf_res$coverage_block_iid,
    cov_block_block  = conf_res$coverage_block_blk
  )
  
  band_tbl <- compute_bandwidth_summary(df_feat, test_idx) %>%
    dplyr::mutate(symbol = symbol)
  
  tail_tbl <- compute_tail_performance(df_feat, test_idx, tail_prob = 0.90) %>%
    dplyr::mutate(symbol = symbol)
  
  tuning_selected_tbl <- xgb_res$best_cfg %>%
    dplyr::mutate(
      symbol = symbol,
      stacking_weight_har = w_hat
    ) %>%
    dplyr::select(
      symbol, eta, max_depth, subsample, colsample_bytree,
      lambda, best_nrounds, validation_rmse, stacking_weight_har
    )
  
  block_sensitivity_tbl <- compute_block_sensitivity(
    df_feat    = df_feat,
    calib_idx  = calib_idx,
    test_idx   = test_idx,
    alpha      = alpha,
    block_grid = block_len_grid
  ) %>%
    dplyr::mutate(symbol = symbol) %>%
    dplyr::select(
      symbol, block_len,
      cov_point_iid, cov_point_block,
      cov_block_iid, cov_block_block,
      mean_width_iid, mean_width_block
    )
  
  plot_test_fit(df_feat, test_idx, symbol, out_dir)
  plot_conf_band_iid_vs_block(df_feat, test_idx, symbol, out_dir)
  plot_score_hist(conf_res$scores_iid, conf_res$scores_block, symbol, out_dir)
  
  safe_write_csv(
    xgb_res$tuning_results %>%
      dplyr::mutate(symbol = symbol) %>%
      dplyr::select(symbol, dplyr::everything()),
    file.path(
      out_dir, "tables",
      paste0("XGB_tuning_grid_", gsub("[^A-Za-z0-9]", "_", symbol), ".csv")
    )
  )
  
  saveRDS(
    list(
      symbol                = symbol,
      proxy_type            = proxy_type,
      df_feat               = df_feat,
      idx                   = idx,
      desc_tbl              = desc_tbl,
      rmse_tbl              = rmse_tbl,
      coverage_tbl          = coverage_tbl,
      band_tbl              = band_tbl,
      tail_tbl              = tail_tbl,
      tuning_selected_tbl   = tuning_selected_tbl,
      tuning_grid_tbl       = xgb_res$tuning_results,
      block_sensitivity_tbl = block_sensitivity_tbl,
      har_model             = har_res$model,
      xgb_model             = xgb_res$model,
      w_hat                 = w_hat,
      scores_iid            = conf_res$scores_iid,
      scores_block          = conf_res$scores_block
    ),
    file = file.path(
      out_dir, "rds",
      paste0("har_xgb_block_", gsub("[^A-Za-z0-9]", "_", symbol), "_", proxy_type, ".rds")
    )
  )
  
  cat("==== ", symbol, " ====\n", sep = "")
  cat("Proxy type: ", proxy_type, "\n", sep = "")
  cat("Stacking weight (HAR share): ", round(w_hat, 3), "\n", sep = "")
  print(desc_tbl)
  print(rmse_tbl)
  print(coverage_tbl)
  print(band_tbl)
  print(tail_tbl)
  print(tuning_selected_tbl)
  print(block_sensitivity_tbl)
  cat("\n")
  
  list(
    symbol                = symbol,
    proxy_type            = proxy_type,
    desc_tbl              = desc_tbl,
    rmse_tbl              = rmse_tbl,
    coverage_tbl          = coverage_tbl,
    band_tbl              = band_tbl,
    tail_tbl              = tail_tbl,
    tuning_selected_tbl   = tuning_selected_tbl,
    block_sensitivity_tbl = block_sensitivity_tbl
  )
}

## ================================================================
## 13) Run main empirical pipeline
## ================================================================
run_main_pipeline <- function(assets          = ASSETS_MAIN,
                              from            = as.Date("2015-01-01"),
                              to              = as.Date("2024-12-31"),
                              proxy_type      = "squared_return",
                              alpha           = 0.05,
                              block_len_main  = 5L,
                              block_len_grid  = c(3L, 5L, 10L),
                              p_lag           = 5L,
                              tail_lambda     = 0,
                              out_dir         = OUT_DIR) {
  res_list <- purrr::map(
    assets,
    ~ run_for_asset(
      symbol          = .x,
      from            = from,
      to              = to,
      proxy_type      = proxy_type,
      alpha           = alpha,
      block_len_main  = block_len_main,
      block_len_grid  = block_len_grid,
      p_lag           = p_lag,
      tail_lambda     = tail_lambda,
      out_dir         = out_dir
    )
  )
  res_list <- purrr::compact(res_list)
  
  summary_desc_tbl <- dplyr::bind_rows(purrr::map(res_list, "desc_tbl")) %>%
    dplyr::mutate(Asset = asset_label(symbol)) %>%
    dplyr::transmute(
      Asset,
      N                         = as.integer(n_obs),
      `Mean return`             = round(mean_r, 4),
      `SD return`               = round(sd_r, 4),
      `Skewness return`         = round(skew_r, 4),
      `Kurtosis return`         = round(kurt_r, 4),
      `Mean variance proxy`     = round(mean_vp, 4),
      `SD variance proxy`       = round(sd_vp, 4),
      `Skewness variance proxy` = round(skew_vp, 4),
      `Kurtosis variance proxy` = round(kurt_vp, 4)
    ) %>%
    dplyr::arrange(Asset)
  
  safe_write_csv(summary_desc_tbl, file.path(out_dir, "tables", "T0_descriptive_by_asset.csv"))
  save_latex_table(
    summary_desc_tbl,
    file      = file.path(out_dir, "tables", "T0_descriptive_by_asset.tex"),
    caption   = "Descriptive statistics of daily returns and the daily variance proxy.",
    label     = "tab:T0_descriptive",
    longtable = TRUE
  )
  
  summary_rmse_tbl <- dplyr::bind_rows(purrr::map(res_list, "rmse_tbl")) %>%
    dplyr::mutate(
      Asset = asset_label(symbol),
      model = pretty_model_label(model)
    ) %>%
    dplyr::select(Asset, model, RMSE) %>%
    tidyr::pivot_wider(names_from = model, values_from = RMSE) %>%
    dplyr::mutate(
      HAR           = round(HAR, 3),
      XGBoost       = round(XGBoost, 3),
      `HAR+XGBoost` = round(`HAR+XGBoost`, 3)
    ) %>%
    dplyr::arrange(Asset)
  
  safe_write_csv(summary_rmse_tbl, file.path(out_dir, "tables", "T1_rmse_summary_by_asset.csv"))
  save_latex_table(
    summary_rmse_tbl,
    file    = file.path(out_dir, "tables", "T1_rmse_summary_by_asset.tex"),
    caption = "Out-of-sample RMSE for HAR, XGBoost, and hybrid HAR+XGBoost forecasts on the log-variance-proxy scale.",
    label   = "tab:T1_rmse"
  )
  
  summary_cov_tbl <- dplyr::bind_rows(purrr::map(res_list, "coverage_tbl")) %>%
    dplyr::mutate(Asset = asset_label(symbol)) %>%
    dplyr::transmute(
      Asset,
      `Pointwise (standard)` = round(cov_point_iid, 3),
      `Pointwise (block)`    = round(cov_point_block, 3),
      `Blockwise (standard)` = round(cov_block_iid, 3),
      `Blockwise (block)`    = round(cov_block_block, 3)
    ) %>%
    dplyr::arrange(Asset)
  
  safe_write_csv(summary_cov_tbl, file.path(out_dir, "tables", "T2_coverage_by_asset.csv"))
  save_latex_table(
    summary_cov_tbl,
    file    = file.path(out_dir, "tables", "T2_coverage_by_asset.tex"),
    caption = "Empirical pointwise and blockwise coverage of standard and block-conformal prediction bands.",
    label   = "tab:T2_coverage"
  )
  
  summary_band_tbl <- dplyr::bind_rows(purrr::map(res_list, "band_tbl")) %>%
    dplyr::mutate(Asset = asset_label(symbol)) %>%
    dplyr::select(Asset, band_type, mean_width, median_width) %>%
    tidyr::pivot_wider(
      names_from  = band_type,
      values_from = c(mean_width, median_width)
    ) %>%
    dplyr::transmute(
      Asset,
      `Mean width (standard)`   = round(`mean_width_Standard split-conformal`, 4),
      `Median width (standard)` = round(`median_width_Standard split-conformal`, 4),
      `Mean width (block)`      = round(`mean_width_Block-conformal`, 4),
      `Median width (block)`    = round(`median_width_Block-conformal`, 4)
    ) %>%
    dplyr::arrange(Asset)
  
  safe_write_csv(summary_band_tbl, file.path(out_dir, "tables", "T3_bandwidth_by_asset.csv"))
  save_latex_table(
    summary_band_tbl,
    file    = file.path(out_dir, "tables", "T3_bandwidth_by_asset.tex"),
    caption = "Average and median widths of standard and block-conformal bands on the daily variance-proxy scale.",
    label   = "tab:T3_bandwidth"
  )
  
  summary_tail_tbl <- dplyr::bind_rows(purrr::map(res_list, "tail_tbl")) %>%
    dplyr::mutate(Asset = asset_label(symbol)) %>%
    dplyr::transmute(
      Asset,
      `Tail observations`     = as.integer(n_tail),
      `HAR RMSE`              = round(RMSE_HAR_tail, 3),
      `XGBoost RMSE`          = round(RMSE_XGB_tail, 3),
      `HAR+XGBoost RMSE`      = round(RMSE_HYB_tail, 3),
      `Coverage (standard)`   = round(cov_iid_tail, 3),
      `Coverage (block)`      = round(cov_block_tail, 3),
      `Mean width (standard)` = round(mean_width_iid_tail, 4),
      `Mean width (block)`    = round(mean_width_blk_tail, 4)
    ) %>%
    dplyr::arrange(Asset)
  
  safe_write_csv(summary_tail_tbl, file.path(out_dir, "tables", "T4_tail_performance_by_asset.csv"))
  save_latex_table(
    summary_tail_tbl,
    file    = file.path(out_dir, "tables", "T4_tail_performance_by_asset.tex"),
    caption = "Tail-region performance in the upper decile of the daily variance-proxy distribution.",
    label   = "tab:T4_tail"
  )
  
  summary_tuning_tbl <- dplyr::bind_rows(purrr::map(res_list, "tuning_selected_tbl")) %>%
    dplyr::mutate(Asset = asset_label(symbol)) %>%
    dplyr::transmute(
      Asset,
      Eta          = round(eta, 2),
      Depth        = as.integer(max_depth),
      Subsample    = round(subsample, 2),
      `Col. samp.` = round(colsample_bytree, 2),
      Lambda       = round(lambda, 2),
      Rounds       = as.integer(best_nrounds),
      `Val. RMSE`  = round(validation_rmse, 4),
      `HAR share`  = round(stacking_weight_har, 3)
    ) %>%
    dplyr::arrange(Asset)
  
  safe_write_csv(summary_tuning_tbl, file.path(out_dir, "tables", "T5_xgboost_tuning_summary.csv"))
  save_latex_table(
    summary_tuning_tbl,
    file    = file.path(out_dir, "tables", "T5_xgboost_tuning_summary.tex"),
    caption = "Time-ordered inner-validation tuning results for the XGBoost component and the selected stacking weight.",
    label   = "tab:T5_tuning"
  )
  
  summary_block_sens_tbl <- dplyr::bind_rows(purrr::map(res_list, "block_sensitivity_tbl")) %>%
    dplyr::mutate(Asset = asset_label(symbol)) %>%
    dplyr::transmute(
      Asset,
      `Block length`          = as.integer(block_len),
      `Pointwise (standard)`  = round(cov_point_iid, 3),
      `Pointwise (block)`     = round(cov_point_block, 3),
      `Blockwise (standard)`  = round(cov_block_iid, 3),
      `Blockwise (block)`     = round(cov_block_block, 3),
      `Mean width (standard)` = round(mean_width_iid, 4),
      `Mean width (block)`    = round(mean_width_block, 4)
    ) %>%
    dplyr::arrange(Asset, `Block length`)
  
  safe_write_csv(summary_block_sens_tbl, file.path(out_dir, "tables", "T6_block_length_sensitivity.csv"))
  save_latex_table(
    summary_block_sens_tbl,
    file      = file.path(out_dir, "tables", "T6_block_length_sensitivity.tex"),
    caption   = "Sensitivity of coverage and band width to the block length used in block-conformal calibration.",
    label     = "tab:T6_blocklength",
    longtable = TRUE
  )
  
  cat("Main empirical pipeline finished. Outputs saved in: ",
      normalizePath(out_dir), "\n", sep = "")
  
  invisible(
    list(
      res_list             = res_list,
      T0_descriptive       = summary_desc_tbl,
      T1_rmse              = summary_rmse_tbl,
      T2_coverage          = summary_cov_tbl,
      T3_bandwidth         = summary_band_tbl,
      T4_tail              = summary_tail_tbl,
      T5_tuning            = summary_tuning_tbl,
      T6_block_sensitivity = summary_block_sens_tbl
    )
  )
}

## ================================================================
## 14) Simulation study: AR(1)-GARCH(1,1)
## ================================================================
simulate_ar1_garch <- function(n_total = 2500L,
                               burn_in = 500L,
                               phi     = 0.10,
                               omega   = 5e-6,
                               alpha   = 0.05,
                               beta    = 0.90,
                               df_t    = 8L) {
  n <- as.integer(n_total + burn_in)
  
  r      <- numeric(n)
  sigma2 <- numeric(n)
  z      <- stats::rt(n, df = df_t)
  
  sigma2[1] <- omega / (1 - alpha - beta)
  r[1]      <- sqrt(sigma2[1]) * z[1]
  
  for (t in 2:n) {
    sigma2[t] <- omega + alpha * r[t - 1]^2 + beta * sigma2[t - 1]
    r[t]      <- phi * r[t - 1] + sqrt(sigma2[t]) * z[t]
  }
  
  idx_use <- (burn_in + 1):n
  r_use   <- r[idx_use]
  vp_use  <- r_use^2
  
  tibble::tibble(
    t                  = seq_along(r_use),
    open               = NA_real_,
    high               = NA_real_,
    low                = NA_real_,
    close              = NA_real_,
    adjusted           = NA_real_,
    r                  = r_use,
    proxy_type         = "squared_return",
    variance_proxy     = vp_use,
    log_variance_proxy = log(vp_use + 1e-8)
  )
}

sim_one_rep <- function(rep_id     = 1L,
                        n_total    = 2500L,
                        burn_in    = 500L,
                        alpha      = 0.05,
                        block_grid = c(5L, 10L),
                        p_lag      = 5L) {
  df_vp <- simulate_ar1_garch(
    n_total = n_total,
    burn_in = burn_in
  )
  
  df_feat <- make_features(df_vp, p_lag = p_lag)
  n <- nrow(df_feat)
  if (n < 300L) return(NULL)
  
  idx <- make_time_split(n)
  train_idx <- idx$train
  calib_idx <- idx$calib
  test_idx  <- idx$test
  
  har_res <- fit_har(df_feat, train_idx = train_idx)
  df_feat <- har_res$df_feat
  df_feat$hyb_pred <- df_feat$har_pred
  
  purrr::map_dfr(block_grid, function(BL) {
    conf_res <- add_conformal_intervals_block(
      df_feat   = df_feat,
      calib_idx = calib_idx,
      test_idx  = test_idx,
      alpha     = alpha,
      block_len = as.integer(BL)
    )
    
    df_conf <- conf_res$df_feat
    width_iid <- df_conf$conf_hi_vp_iid[test_idx] - df_conf$conf_lo_vp_iid[test_idx]
    width_blk <- df_conf$conf_hi_vp_blk[test_idx] - df_conf$conf_lo_vp_blk[test_idx]
    
    tibble::tibble(
      rep_id = rep_id,
      B      = as.integer(BL),
      scheme_iid_point_cov  = conf_res$coverage_iid,
      scheme_iid_block_cov  = conf_res$coverage_block_iid,
      scheme_iid_mean_width = mean(width_iid, na.rm = TRUE),
      scheme_blk_point_cov  = conf_res$coverage_blk_point,
      scheme_blk_block_cov  = conf_res$coverage_block_blk,
      scheme_blk_mean_width = mean(width_blk, na.rm = TRUE)
    )
  })
}

run_simulation <- function(
    n_rep      = 500L,
    n_total    = 2500L,
    burn_in    = 500L,
    alpha      = 0.05,
    block_grid = c(5L, 10L),
    p_lag      = 5L,
    out_dir    = OUT_DIR
) {
  set.seed(999)
  
  sim_out_dir_tables <- file.path(out_dir, "tables")
  sim_out_dir_plots  <- file.path(out_dir, "plots")
  dir.create(sim_out_dir_tables, showWarnings = FALSE, recursive = TRUE)
  dir.create(sim_out_dir_plots,  showWarnings = FALSE, recursive = TRUE)
  
  all_res <- purrr::map_dfr(seq_len(n_rep), function(b) {
    if (b %% 25 == 0) cat("Simulation replication:", b, "/", n_rep, "\n")
    sim_one_rep(
      rep_id     = b,
      n_total    = n_total,
      burn_in    = burn_in,
      alpha      = alpha,
      block_grid = block_grid,
      p_lag      = p_lag
    )
  })
  
  res_long <- all_res %>%
    tidyr::pivot_longer(
      cols      = -c(rep_id, B),
      names_to  = "metric",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      scheme = dplyr::case_when(
        grepl("^scheme_iid_", metric) ~ "Standard split-conformal",
        grepl("^scheme_blk_", metric) ~ "Block-conformal",
        TRUE ~ NA_character_
      ),
      metric_clean = dplyr::case_when(
        grepl("point_cov$", metric)  ~ "Pointwise coverage",
        grepl("block_cov$", metric)  ~ "Blockwise coverage",
        grepl("mean_width$", metric) ~ "Mean width",
        TRUE                         ~ metric
      )
    ) %>%
    dplyr::filter(!is.na(scheme)) %>%
    dplyr::select(rep_id, B, scheme, metric_clean, value) %>%
    tidyr::pivot_wider(
      names_from  = metric_clean,
      values_from = value
    )
  
  sim_summary <- res_long %>%
    dplyr::group_by(scheme, B) %>%
    dplyr::summarise(
      `Mean pointwise coverage` = mean(`Pointwise coverage`, na.rm = TRUE),
      `SD pointwise coverage`   = sd(`Pointwise coverage`, na.rm = TRUE),
      `Mean blockwise coverage` = mean(`Blockwise coverage`, na.rm = TRUE),
      `SD blockwise coverage`   = sd(`Blockwise coverage`, na.rm = TRUE),
      `Mean width`              = mean(`Mean width`, na.rm = TRUE),
      `SD width`                = sd(`Mean width`, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      `Mean pointwise coverage` = round(`Mean pointwise coverage`, 4),
      `SD pointwise coverage`   = round(`SD pointwise coverage`, 4),
      `Mean blockwise coverage` = round(`Mean blockwise coverage`, 4),
      `SD blockwise coverage`   = round(`SD blockwise coverage`, 4),
      `Mean width`              = round(`Mean width`, 4),
      `SD width`                = round(`SD width`, 4)
    )
  
  sim_cov_tbl <- sim_summary %>%
    dplyr::select(
      Scheme = scheme,
      `Block length` = B,
      `Mean pointwise coverage`,
      `SD pointwise coverage`,
      `Mean blockwise coverage`,
      `SD blockwise coverage`
    )
  
  sim_width_tbl <- sim_summary %>%
    dplyr::select(
      Scheme = scheme,
      `Block length` = B,
      `Mean width`,
      `SD width`
    )
  
  safe_write_csv(sim_cov_tbl, file.path(sim_out_dir_tables, "S1_simulation_coverage.csv"))
  safe_write_csv(sim_width_tbl, file.path(sim_out_dir_tables, "S2_simulation_width.csv"))
  
  save_latex_table(
    sim_cov_tbl,
    file    = file.path(sim_out_dir_tables, "S1_simulation_coverage.tex"),
    caption = "Simulation summary: pointwise and blockwise coverage under the AR(1)-GARCH(1,1) design.",
    label   = "tab:S1_sim_cov"
  )
  
  save_latex_table(
    sim_width_tbl,
    file    = file.path(sim_out_dir_tables, "S2_simulation_width.tex"),
    caption = "Simulation summary: mean band width under the AR(1)-GARCH(1,1) design.",
    label   = "tab:S2_sim_width"
  )
  
  df_plot <- sim_summary %>%
    dplyr::transmute(
      scheme = scheme,
      B      = as.factor(B),
      mean_width = `Mean width`,
      mean_cov_block = `Mean blockwise coverage`
    )
  
  p_tradeoff <- ggplot2::ggplot(
    df_plot,
    ggplot2::aes(
      x      = mean_width,
      y      = mean_cov_block,
      colour = scheme,
      shape  = B
    )
  ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_line(ggplot2::aes(group = scheme)) +
    ggplot2::labs(
      x      = "Mean band width on the variance-proxy scale",
      y      = "Mean blockwise coverage",
      colour = "Scheme",
      shape  = "Block length"
    ) +
    base_clean_theme()
  
  ggplot2::ggsave(
    filename = file.path(sim_out_dir_plots, "F5_simulation_tradeoff.png"),
    plot     = p_tradeoff,
    width    = 7, height = 4.5, dpi = 600
  )
  
  cat("Simulation finished.\n")
  invisible(list(raw_results = res_long, summary_table = sim_summary))
}

## ================================================================
## 15) Optional appendix robustness: alternative variance proxies
## ================================================================
run_proxy_robustness <- function(assets      = ASSETS_MAIN,
                                 proxy_grid  = c("squared_return", "parkinson", "gk", "rs"),
                                 from        = as.Date("2015-01-01"),
                                 to          = as.Date("2024-12-31"),
                                 alpha       = 0.05,
                                 block_len   = 5L,
                                 p_lag       = 5L,
                                 tail_lambda = 0,
                                 out_dir     = OUT_DIR) {
  all_results <- list()
  
  for (px in proxy_grid) {
    cat("\n==============================\n")
    cat("Running proxy robustness for:", px, "\n")
    cat("==============================\n")
    
    proxy_res <- purrr::map(
      assets,
      ~ run_for_asset(
        symbol          = .x,
        from            = from,
        to              = to,
        proxy_type      = px,
        alpha           = alpha,
        block_len_main  = block_len,
        block_len_grid  = c(3L, 5L, 10L),
        p_lag           = p_lag,
        tail_lambda     = tail_lambda,
        out_dir         = out_dir
      )
    )
    proxy_res <- purrr::compact(proxy_res)
    all_results[[px]] <- proxy_res
  }
  
  T8_proxy_tbl <- purrr::imap_dfr(all_results, function(res_list, px) {
    dplyr::bind_rows(purrr::map(res_list, "rmse_tbl")) %>%
      dplyr::filter(model == "HAR+XGBoost") %>%
      dplyr::left_join(
        dplyr::bind_rows(purrr::map(res_list, "coverage_tbl")) %>%
          dplyr::select(symbol, cov_point_block, cov_block_block),
        by = "symbol"
      ) %>%
      dplyr::left_join(
        dplyr::bind_rows(purrr::map(res_list, "band_tbl")) %>%
          dplyr::filter(band_type == "Block-conformal") %>%
          dplyr::select(symbol, mean_width),
        by = "symbol"
      ) %>%
      dplyr::mutate(proxy_type = px)
  }) %>%
    dplyr::mutate(
      Asset = asset_label(symbol),
      Proxy = pretty_proxy_label(proxy_type)
    ) %>%
    dplyr::transmute(
      Asset,
      Proxy,
      `Hybrid RMSE`        = round(RMSE, 3),
      `Pointwise coverage` = round(cov_point_block, 3),
      `Blockwise coverage` = round(cov_block_block, 3),
      `Mean block width`   = round(mean_width, 4)
    ) %>%
    dplyr::arrange(Asset, Proxy)
  
  safe_write_csv(T8_proxy_tbl, file.path(out_dir, "tables", "T8_proxy_robustness.csv"))
  save_latex_table(
    T8_proxy_tbl,
    file      = file.path(out_dir, "tables", "T8_proxy_robustness.tex"),
    caption   = "Alternative variance-proxy robustness for the hybrid model with block-conformal calibration.",
    label     = "tab:T8_proxy",
    longtable = TRUE
  )
  
  cat("Proxy robustness finished.\n")
  invisible(T8_proxy_tbl)
}

## ================================================================
## 16) Example calls
## ================================================================

## --- Main empirical run ---
## main_res <- run_main_pipeline(
##   assets         = ASSETS_MAIN,
##   from           = as.Date("2015-01-01"),
##   to             = as.Date("2024-12-31"),
##   proxy_type     = "squared_return",
##   alpha          = 0.05,
##   block_len_main = 5L,
##   block_len_grid = c(3L, 5L, 10L),
##   p_lag          = 5L,
##   tail_lambda    = 0,
##   out_dir        = OUT_DIR
## )

## --- Simulation run ---
## sim_res <- run_simulation(
##   n_rep      = 500L,
##   n_total    = 2500L,
##   burn_in    = 500L,
##   alpha      = 0.05,
##   block_grid = c(5L, 10L),
##   p_lag      = 5L,
##   out_dir    = OUT_DIR
## )

## --- Optional appendix robustness ---
## proxy_res <- run_proxy_robustness(
##   assets      = ASSETS_MAIN,
##   proxy_grid  = c("squared_return", "parkinson", "gk", "rs"),
##   from        = as.Date("2015-01-01"),
##   to          = as.Date("2024-12-31"),
##   alpha       = 0.05,
##   block_len   = 5L,
##   p_lag       = 5L,
##   tail_lambda = 0,
##   out_dir     = OUT_DIR
## )

cat("All final corrected code blocks are ready.\n")