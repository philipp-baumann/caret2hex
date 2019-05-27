################################################################################
## Functions to tune and evaluate Cubist models based on RMSE and
##   by using a nested resampling approach
################################################################################

## Implement repeated group v-fold cross-validation using additional functions
## to repeat rsample::group_vfold_cv() and create combined resampling `id` for
## `repeat` and `fold` =========================================================

# Add a `id_repeat` column for a list of `rset` objects ------------------------

add_cv_repeat_rset <- function(rset_lst,
                          id_repeat = id_repeat) {
  id_repeat <- rlang::ensym(id_repeat)
  imap(rset_lst,
    ~ dplyr::mutate(.x, !!id_repeat := paste0("Repeat", .y))
  )
}

# Rename "Resample<ii>" to "Fold<ii>" in `id` column ---------------------------

replace_resample_values <- function(rset_lst,
                                    resample_col = id,
                                    fold_col = id_fold) {
  resample_col <- rlang::ensym(resample_col)
  fold_col <- rlang::ensym(fold_col)
  
  map(rset_lst,
    ~ dplyr::mutate(.x,
        !!fold_col := stringr::str_replace(
          string = !!resample_col, pattern = "Resample", replacement = "Fold"))
  )
}

# Create a new combined resample `id` by pasting `repeat_col` and `fold_col` ---
bind_repeat_fold <- function(rset_lst,
                             repeat_col = id_repeat,
                             fold_col = id_fold,
                             repeat_fold_col = id) {
  repeat_col <- rlang::ensym(repeat_col)
  fold_col <- rlang::ensym(fold_col)
  repeat_fold_col <- rlang::ensym(repeat_fold_col)
  
  map(rset_lst,
    ~ dplyr::mutate(.x, !!repeat_fold_col := paste(
      !!repeat_col, !!fold_col, sep = "."))
  )
}

# Row bind repeated `group_vfold_cv` `rset` objects contained in `rset_lst` 
# list of rsplit objects and create a combined resample `id` containing
# a combined "Repeat<i>" and "Fold<ii>" identifier for the outer splits
# (`splits` list-column conating `rsplit` objects) -----------------------------

bind_groupcv_reps <- function(rset_lst, ...) {
                              # resample_col = id,
                              # repeat_col = id_repeat,
                              # fold_col = id_fold,
                              # repeat_fold_col = id) {
  user_symbols <- rlang::ensyms(...)
  vars_chr <- map_chr(user_symbols, rlang::quo_name)
  vars_rm <- unname(vars_chr[c("repeat_col")]) # "fold_col"
  # see partition_cols()
  # https://adv-r.hadley.nz/quasiquotation.html#quasi-case-studies
  vars_rm <- map(vars_rm, ~ rlang::expr(-!!rlang::sym(.x)))

  cv_repeat <- add_cv_repeat_rset(rset_lst = rset_lst,
    id_repeat = !!user_symbols[["repeat_col"]])
  
  cv_fold <- replace_resample_values(
    rset_lst = cv_repeat,
    resample_col = !!user_symbols[["resample_col"]],
    fold_col = !!user_symbols[["fold_col"]])
  
  cv_repeat_fold <- bind_repeat_fold(
    rset_lst = cv_fold,
    repeat_col = !!user_symbols[["repeat_col"]],
    fold_col = !!user_symbols[["fold_col"]],
    repeat_fold_col = !!user_symbols[["repeat_fold_col"]])
  
  cv_repeat_fold <- do.call(rbind, cv_repeat_fold)
  dplyr::select(cv_repeat_fold, !!!vars_rm)
}


# Extract and mutate fold id into a new column ---------------------------------

add_cv_fold <- function(data,
                        repeat_fold_col = id,
                        fold_col = id_fold) {
  repeat_fold_col <- rlang::ensym(repeat_fold_col)
  fold_col <- rlang::ensym(fold_col)
  dplyr::mutate(data,
    !!fold_col := stringr::str_extract(
      string = !!repeat_fold_col, pattern = "Fold[:digit:]{2}"))
}

# Extract and mutate repeat id into a new column -------------------------------

add_cv_repeat <- function(data,
                          repeat_fold_col = id,
                          repeat_col = id_repeat) {
  repeat_fold_col <- rlang::ensym(repeat_fold_col)
  repeat_col <- rlang::ensym(repeat_col)
  dplyr::mutate(data,
    !!repeat_col := stringr::str_extract(
      string = !!repeat_fold_col, pattern = "Repeat[:digit:]"))
}

# Row bind list of tibbles that each contain resample id (`repeat_fold_col`)
# and observed (`obs`) and predicted values (`pred`) and sample metadata
# as columns -------------------------------------------------------------------

rbind_cv_repeat <- function(.data,
                            repeat_fold_col = resample_id,
                            repeat_col = id_repeat) {
  repeat_fold_col <- rlang::ensym(repeat_fold_col)
  repeat_col = rlang::ensym(repeat_col)
  
  dplyr::bind_rows(.data) %>%
    add_cv_repeat(repeat_fold_col = !!repeat_fold_col,
      repeat_col = !!repeat_col)
}


# Extract predicticted and observed values for all list_columns across
# all resamples ----------------------------------------------------------------

extract_predobs_cols <- function(.data,
                                 repeat_fold_col = resample_id,
                                 repeat_col = id_repeat) {
  col_syms <- rlang::syms(names(.data))
  vars_chr <- map_chr(col_syms, rlang::quo_name)
  
  repeat_fold_col <- rlang::ensym(repeat_fold_col)
  repeat_col <- rlang::ensym(repeat_col)
  
  out_cols <- map_at(
    .x = .data, 
    .at = vars_chr,
    .f = ~ rbind_cv_repeat(.x, repeat_fold_col = !!repeat_fold_col,
      repeat_col = !!repeat_col))
  
  tibble::tibble(variable = vars_chr, data = as.list(out_cols)) %>%
    tidyr::unnest(data)
}


# Calculate performance metrics by the resampling repeat -----------------------

metrics_by_repeat <- function(.lcol,
                              repeat_col = id_repeat) {
  
  repeat_col <- rlang::ensym(repeat_col)
  # Bind rows and add repeat id column
  df <- dplyr::bind_rows(.lcol)
  df_repeat <- add_cv_repeat(data = df,
    repeat_fold_col = resample_id,
    repeat_col = id_repeat)
  
  df_repeat %>%
    dplyr::group_by(!!repeat_col) %>%
    tidyr::nest() %>%
    dplyr::mutate(metrics = map(data, ~ simplerspec::evaluate_model(
      data = .x, obs = "obs", pred = "pred"))) %>%
    tidyr::unnest(metrics)
}

# Calculate performance metrics across the resampling repeats ------------------

metrics_across_repeats <- function(.data,
                                   repeat_col = id_repeat) {
  repeat_col <- rlang::ensym(repeat_col)
  
  metrics <- metrics_by_repeat(.lcol = .data)
  metrics %>%
    dplyr::summarize_if(is.numeric,
      funs(mean = mean, sd = sd)) # // pb 20180814: ci = simplerspec:::sem_ci
}


# Calculate metrics across repeats for all tibble list-columns -----------------

metrics_cols <- function(.data,
                         repeat_col = id_repeat) {
  col_syms <- rlang::syms(names(.data))
  vars_chr <- map_chr(col_syms, rlang::quo_name)
  
  repeat_col <- rlang::ensym(repeat_col)
  
  out_cols <- map_at(
    .x = .data, 
    .at = vars_chr,
    .f = ~ metrics_across_repeats(.x, repeat_col = !!repeat_col))
  
  tibble::tibble(variable = vars_chr, data = as.list(out_cols)) %>%
    tidyr::unnest(data)
}


## Helper functions to format and paste columns for model evaluation tables ====

# Determine optimal number of digits after the comma for each element in 
# a numeric atomic vector, for a subsequent call to 
# base::signif() or base::round() and further variants in this family ----------

digits_signif <- function(x, digits = NULL, .digit_lag = 1) {
  if (!is.null(digits)) {
    digits_comma <- rep(digits, length(x))
  } else {
    digits_comma <- ifelse(abs(x) > 100, 1,
      ifelse(abs(x) < 10 & abs(x) > 1, 2,
        ceiling(log10(abs(x))) - .digit_lag)) # // 20180514: round(log10(abs(x)), 0)
  }
  digits_comma
}

# Variant for base::sprintf(), his function will always return 
# positive integer values, which can be later be supplied to sprintf() for nice
# formatting
digits_sprintf <- function(x, digits = NULL, .digit_lag = 1) {
  if (!is.null(digits)) {
    digits_comma <- rep(digits, length(x))
  } else {
    digits_comma <- ifelse(abs(x) > 100, 1,
      ifelse(abs(x) < 10 & abs(x) > 1, 2,
        abs(ceiling(log10(abs(x))) - .digit_lag)))
  }
  digits_comma
}


# sprintf() wrapper that return formated digits as character vector
format_sprintf <- function(x, .digits) {
  # Create character vector of format strings
  fmt_strings <- map_chr(.digits, ~ paste0("%.", .x, "f"))
  # Format numeric vector to character vector, formatting to specified digits
  purrr::imap_chr(x, ~ sprintf(fmt = fmt_strings[.y], .x))
}

# Round and paste two columns into one column that contains string 
# of first column value pasted with second column value, using "$\\pm$" as
# separator --------------------------------------------------------------------

paste_plusminus <- function(.data, ..., digits = NULL, .digit_lag = 1) {
  cols <- rlang::enquos(...)
  
  # Determine name of the new pasted plusminus column
  cols_nm <- map_chr(cols, quo_name)
  nm_notundersc <- map(cols_nm, ~ str_extract_all(.x, "[^_]+"))
  nm_first <- purrr::modify_depth(nm_notundersc, 2, 1) %>% purrr::flatten()

  if (Reduce(all.equal, nm_first) == TRUE) {
    nm_plusminus <- nm_first[[1]]
  } else {
    nm_plusminus <- do.call(paste, c(as.list(cols_nm), sep = "_pm_"))
  }
  
  # Vector of digits (integer); this is derived from first column
  # specified in the `...`. argument
  digits_mean <- digits_sprintf(x = dplyr::pull(.data, !!cols[[1]]),
    digits = digits, .digit_lag = .digit_lag)
  
  new_cols <- map(cols,
    ~ format_sprintf(x = dplyr::pull(.data, !!.x), .digits = digits_mean))
  
  new_plusminus <- pmap(new_cols,
    ~ paste(..1, ..2, sep = "$\\,\\pm{}\\,$")) %>% purrr::flatten_chr()
  
  .data %>%
    dplyr::mutate(!!nm_plusminus := new_plusminus)
}



## Custom wrappers for performing nested resampling with Cubist ================


# Modified from:
# https://topepo.github.io/rsample/articles/Applications/Nested_Resampling.html

rules_per_committee_splits <- function(cubist) {
  splits <- cubist[["splits"]]
  splits %>%
    dplyr::group_by(committee) %>%
    dplyr::summarize(n_rules = dplyr::n_distinct(rule)) %>%
    dplyr::pull(n_rules)
}

rules_per_committee <- function(model) {
  cubist_call <- rlang::set_names(capture.output(model))
  rules_match <- stringr::str_match(cubist_call, "rules") %>% as.vector() %>%
    is.na() %>% magrittr::not()
  rules_chr <- cubist_call[rules_match]
  rules_chr_split <- stringr::str_split(rules_chr, ",") %>% unlist()
  rules <- map_int(rules_chr_split, ~ as.integer(stringr::str_extract(.x, "\\d")))
  rules
}

# Fit a model to an `rsplit` object using a Cubist parameter set
# and predict on hold out using inner or outer resampling hold-outs ------------

# `object` will be an `rsplit` object from our `results` tibble
# `commitees` and `neighbors` are the tuning parameters
cubist_rmse <- function(object, committees = 5, neighbors = 2,
                        x_var, y_var,
                        keep_predictions = FALSE) {
  
  data_analysis <- rsample::analysis(object) %>% na.omit()
  y_values <- data_analysis[[y_var]]
  X_values <- data.table::rbindlist(data_analysis[[x_var]])
  mod <- Cubist::cubist(x = X_values, y = y_values,
    committees = committees, neighbors = neighbors)
  data_assessment <- rsample::assessment(object) %>% na.omit()
  y_assessment <- data_assessment[[y_var]]
  X_assessment <- data.table::rbindlist(data_assessment[[x_var]])
  holdout_pred <- predict(object = mod, newdata = X_assessment,
    neighbors = neighbors)
  rmse <- sqrt(mean((data_assessment[[y_var]] - holdout_pred) ^ 2,
    na.rm = TRUE))
  if (!keep_predictions) {rmse} else {
    tibble::tibble(
      RMSE = rmse,
      rules = list(rules_per_committee(mod)), # rules_per_committee
      predobs = list(tibble::tibble(pred = holdout_pred, obs = y_assessment))
    )
  }
}
  

# In some case, we want to parameterize the function over the tuning parameter
rmse_wrapper <- function(committees, neighbors, object,
                         x_var, y_var) {
  cubist_rmse(object, committees, neighbors, x_var = x_var, y_var = y_var)
}

# For the nested resampling, a model needs to be fit for each tuning parameter
# and each inner resampling split. To do this, a wrapper can be created;
# `object` will be an `rsplit` object for the bootstrap samples
tune_over_cubist_grid <- function(object, x_var, y_var) {
  results <- expand.grid(
    committees = c(5, 10, 20), neighbors = c(2, 5, 7, 9))
  results$RMSE <- map2_dbl(.x = results$committees, .y = results$neighbors,
    ~ rmse_wrapper(committees = .x, neighbors = .y, object = object,
      x_var = x_var, y_var = y_var))
  results
}


# Since this will be called across the set of outer cross-validation splits,
# another wrapper is required:
# `object` is an `rsplit` object in `results$inner_resamples`
summarize_tune_results <- function(object, x_var, y_var) {
  # Return row-bound tibble that has the bootstrap results
  purrr::map_df(object$splits,
    ~ tune_over_cubist_grid(object = .x, x_var = x_var, y_var = y_var)) %>%
  # For each value of the tuning parameter, compute the
  # average RMSE which is the inner bootstrap estimate.
  dplyr::ungroup() %>%
  dplyr::group_by_at(.vars = vars(one_of("committees", "neighbors"))) %>%
  dplyr::mutate( # summarize() does not work in combination with group_by_at() !
    mean_inner_RMSE = mean(RMSE, na.rm = TRUE),
    n = length(RMSE)) %>%
  dplyr::mutate(y_var = y_var)
}


# Parallel wrapper 
# Note that an issue with dplyr::group_by_at() was recently resolved:
# see https://github.com/tidyverse/dplyr/issues/3351
# Update dplyr because issue was fixed on March 12, 2018
# devtools::install_github("krlmlr/bindr")
# devtools::install_github("krlmlr/bindrcpp")
# devtools::install_github("hadley/dplyr")
summarize_tune_results_parallel <- function(object, x_var, y_var) {
  # Return row-bound tibble that has the B bootstrap results
  foreach::foreach(i = seq_along(object$splits), .combine = "rbind",
    .export = c("tune_over_cubist_grid", "rmse_wrapper", "cubist_rmse"),
    .packages = c("purrr", "rsample", "dplyr", "Cubist")) %dopar% { # only first two previously
    tune_over_cubist_grid(
      object = object$splits[[i]], x_var = x_var, y_var = y_var)
    } %>%
  dplyr::ungroup() %>%
  dplyr::group_by_at(.vars = vars(one_of("committees", "neighbors"))) %>%
  dplyr::mutate( # summarize() does not work in combination with group_by_at() !
    mean_inner_RMSE = mean(RMSE, na.rm = TRUE),
    n = length(RMSE)) %>%
  dplyr::mutate(y_var = y_var) %>%
  dplyr::select(-dplyr::one_of("RMSE")) %>%
  dplyr::slice(1L)
}


# Create another wrapper to fit response and predictors from data columns
# within rsample objects -------------------------------------------------------

summarize_tune_results_wrapper <- function(object, x_var, y_var) {
  purrr::map(object, 
    ~ summarize_tune_results(object = .x, x_var = x_var, y_var = y_var))
}

summarize_tune_results_parallel_wrapper <- function(object, x_var, y_var) {
  purrr::map(object,
    ~ summarize_tune_results_parallel(
      object = .x, x_var = x_var, y_var = y_var))
}


# Determine best model parameters for multiple models (one for each response
# variable in <y_var>) by testing the entire tuning grid 
# across the inner resample assessment sets ------------------------------------

extract_inner_resampling <- function(object, x_var, y_var) {
  l <- list(
    "object" = rep(list(object$inner_resamples), length(y_var)),
    "x_vars" = as.list(rep(x_var, length(y_var))),
    "y_vars" = as.list(y_var))
  # names(l["object"]) <- y_var
  names(l[["x_vars"]]) <- rep(x_var, length(y_var))
  names(l[["y_vars"]]) <- y_var
  l
}

tune_inner_multi_parallel <- function(object, x_vars, y_vars) {
  purrr::pmap(
    .l = extract_inner_resampling(object = object, # PB: before spc_prams_nested
      x_var = x_vars, y_var = y_vars),
    .f = ~ summarize_tune_results_parallel_wrapper(
      object = ..1, x_var = ..2, y_var = ..3))
}


## Extract results from Cubist tuning across all combinations of 
## inner resamples and the tuning grid =========================================

best_cubist_tune <- function(dat) dat[which.min(dat$mean_inner_RMSE), ]


# Extract best performing inner Cubist tuning within all training sets of
# outer resamples, and merge results with rsplit object ------------------------

extract_best_inner <- function(object, tuned_inner) {
  pooled_inner <- dplyr::bind_rows(tuned_inner)
  final_params <- tuned_inner %>%
    purrr::map_df(best_cubist_tune) # %>%
    # dplyr::select(
    #   !!! rlang::syms(
    #       c("y_var", "committees", "neighbors", "mean_inner_RMSE", "n"))) %>%
    # dplyr::rowwise()
  final_params
}

merge_best_inner <- function(object, tuned_inner) {
  final_params <- extract_best_inner(object = object, tuned_inner = tuned_inner)
  dplyr::bind_cols(object, final_params)
}

# Merge inner tuning results for multiple y response variables
merge_best_inner_multi <- function(object, tuned_inner) {
  final_params <- purrr::map(tuned_inner, 
    ~ extract_best_inner(object = object, tuned_inner = .x))
  names(final_params) <- lapply(final_params, function(x) unique(x$y_var))
  final_params <- purrr::map(final_params, list)
  final_params <- purrr::map(final_params,
    ~ split(.x[[1]], seq(nrow(.x[[1]]))))
  df_final_params <- tibble::tibble(!!! final_params)
  dplyr::bind_cols(object, df_final_params)
}


## Evaluate performance of multiple y response variable (y_vars) final models on
## assessment sets of outer resamples based on best tuning results across
## inner resamples =============================================================

extract_outer_resampling <- function(object) {
  list(
    "splits" = object$splits,
    "committees" = as.list(object$committees),
    "neighbors" = as.list(object$neighbors))
}


extract_outer_resampling_multi <- function(object, y_var) {
  list(
    "splits" = object$splits,
    "committees" = purrr::map(object[[y_var]], "committees"),
    "neighbors" = purrr::map(object[[y_var]], "neighbors")
  )
}


assess_outer_resampling <- function(object, x_var, y_var) {
  outer_assessments <- purrr::pmap(
    .l = extract_outer_resampling_multi(object = object, y_var = y_var),
    .f = ~ cubist_rmse(
      object = ..1, committees = ..2, neighbors = ..3,
      x_var = x_var, y_var = y_var, keep_predictions = TRUE)
  )
  purrr::map2(.x = object[[y_var]], .y = outer_assessments,
    .f = ~ dplyr::bind_cols(.x, .y))
}

assess_outer_resampling_multi <- function(object, x_vars, y_vars) {
  outer_assessments <- purrr::map2(.x = x_vars, .y = y_vars,
    ~ assess_outer_resampling(object = object, x_var = .x, y_var = .y))
  names(outer_assessments) <- y_vars
  outer_assessments <- tibble::tibble(!!! outer_assessments)
  object %>%
    dplyr::select(-dplyr::one_of(y_vars)) %>%
    dplyr::bind_cols(outer_assessments)
}


## Summarize multiple response models on assessment samples of outer 
## resamples ===================================================================

assess_cubist_nested <- function(object, x_vars, y_vars) {
  # Summarize results from testing Cubist on the entire tuning grid and
  # on all inner resamples assessment sets
  tuned_inner <- tune_inner_multi_parallel(
    object = object,
    x_vars = x_vars, y_vars = y_vars)
  # Merge best inner resampling results with rsplit object
  best_inner <- merge_best_inner_multi(object = object,
    tuned_inner = tuned_inner)
  # Append results from outer assessments using best inner tuning parameters
  outer_assessments <- assess_outer_resampling_multi(object = best_inner,
    x_vars = x_vars, y_vars = y_vars)
  outer_assessments
}


# Helper functions to append <site_id> and <depth_cm> to outer resamples
# assessment predictions of all estimated parameters ---------------------------

extract_outer_assessment_cols <- function(object, id_cols) {
  data_outer_assessment <- purrr::map(object$splits,
    rsample::assessment)
  df_id <- purrr::map(data_outer_assessment, 
    ~ dplyr::select(.x, dplyr::one_of(id_cols)))
  df_id
}


merge_outer_assessment_cols <- function(object, id_cols, y_var) {
  df_id <- extract_outer_assessment_cols(object = object, id_cols = id_cols)
  id_predobs <- purrr::map2(.x = df_id, .y = object[[y_var]],
    ~ dplyr::bind_cols(.x, .y$predobs))
  id_predobs
  y_var_new <- purrr::map(.x = seq_along(object[[y_var]]),
    ~ dplyr::select(object[[y_var]][[.x]], -dplyr::one_of("predobs")) %>%
      dplyr::mutate(predobs = id_predobs[.x])
  )
  y_var_new <- tibble::tibble(!!! list(y_var_new))
  names(y_var_new) <- y_var
  y_var_new
}

# Add IDs for multiple modeled y variables (`y_vars`)
cubist_nested_add_ids <- function(object, x_vars, y_vars, id_cols) {
  # outer_assessments <- assess_cubist_nested(object = object,
  #   x_vars = x_vars, y_vars = y_vars)
  y_vars_new <- purrr::map(y_vars,
    ~ merge_outer_assessment_cols(object = object, id_cols = id_cols,
        y_var = .x))
  names(y_vars_new) <- y_vars
  object %>%
    dplyr::select(-dplyr::one_of(y_vars)) %>%
    dplyr::bind_cols(y_vars_new)
}


## Helper functions to extract outer resampling predictions ====================

# Extract tuned parameters
extract_params_y_var <- function(object, y_var, rep_cv = NULL) {
  committees <- purrr::map(.x = object[[y_var]], "committees")
  neighbors <- purrr::map(.x = object[[y_var]], "neighbors")
  # // pb 2018-05-16: implement tidy evaluation for `id` and `resample_id`
  resample_id <- purrr::map(as.list(object$id),
    ~ tibble::tibble(resample_id = .x))
  # // pb 2018-05-16: implement tidy evaluation for resample_id
  params_id <- purrr::pmap(list(resample_id, committees, neighbors),
    ~ tibble::add_column(..1,
        committees = ..2, neighbors = ..3))
  if (!is.null(rep_cv)) {
    params_id <- purrr::pmap(list(fold_id, committees, neighbors),
      ~ tibble::add_column(..1,
          committees = ..2, neighbors = ..3))
  }
  params_id <- list(params_id)
  names(params_id) <- y_var
  params_id
}

merge_params_y_vars <- function(object, y_vars) {
  params_id_yvars <- purrr::map(.x = y_vars,
    ~ extract_params_y_var(object = object, y_var = .x)) %>%
      purrr::flatten()
  tibble::tibble(!!! params_id_yvars)
}

extract_predobs_y_var <- function(object, y_var, rep_cv = NULL) {
  predobs <- purrr::map(.x = object[[y_var]], "predobs") %>%
    purrr::flatten()
  resample_id <- as.list(object$id)
  # // pb 2018-05-16: implement tidy evaluation for resample_id
  predobs_id <- purrr::map2(.x = predobs, .y = resample_id,
    ~ tibble::add_column(.x,
        resample_id = .y, .before = 1))
  if (!is.null(rep_cv)) {
    predobs_id <- purrr::map2(.x = predobs_id, .y = fold_id,
      ~ tibble::add_column(.x, fold_id = .y, .before = 1))
  }
  predobs_id <- list(predobs_id)
  names(predobs_id) <- y_var
  predobs_id
}

# // pb 2018-05-15: 
merge_predobs_y_vars <- function(object, y_vars) {
  predobs_id_yvars <- purrr::map(.x = y_vars,
    ~ extract_predobs_y_var(object = object, y_var = .x)) %>%
      purrr::flatten()
  tibble::tibble(!!! predobs_id_yvars)
}


# Function to extract and nest resampling prediction results by a column
# variable contained in a data frame
rename_predobs_y_var <- function(object, y_var) {
  # use ensym() instead of enquo() to caputure the y_var argument
  # supplied by the user; enymb check the captured expression is a string or
  # a symbol, and will return a symbol in both cases
  quo_var <- rlang::ensym(y_var)
  new_names <- paste0(rlang::quo_name(quo_var), c("_pred", "_obs"))
  vars <- c("pred", "obs")
  names(vars) <- new_names

  df_renamed <- object %>%
    tidyr::unnest(!! quo_var) %>%
    # note that there will probably be soon a tidyselect feature request
    # to support "bang bang", !!, and
    # the "triple bang", !!!, is needed, see
    # https://github.com/tidyverse/dplyr/issues/3030
    dplyr::rename(!!! vars)
  df_renamed
}


join_predobs_y_vars <- function(..., object) {
  # y_vars <- ensyms(y_vars)
  y_vars <- rlang::ensyms(...)
  # Unquoting `.x`: see 20.6.1 Map-reduce to generate code
  # in https://adv-r.hadley.nz/quasiquotation.html
  dfs_unnested <- purrr::map(.x = y_vars,
    ~ rename_predobs_y_var(object = object, y_var = !!.x))
  # // pb: 20180711: Remove `by = "resample_id"` argument because 
  # there may be other common character or factor variables to join by
  purrr::reduce(dfs_unnested, function(x, y) dplyr::inner_join(x = x, y = y))
}

# Summarize Cubist parameter tuning and rules (majority cases) across resamples 
summarize_outer_tuning <- function(data) {
  
  data_bound <- dplyr::bind_rows(data)
  
  summary_committees <- data_bound %>%
    dplyr::select(y_var, committees) %>%
    dplyr::group_by(y_var) %>%
    dplyr::count(committees) %>%
    dplyr::rename(n_committees = n) %>%
    dplyr::mutate(committees_count =
      paste0(committees, "(", n_committees, ")")) %>% 
    dplyr::arrange(desc(n_committees)) %>%
    select(y_var, committees_count)
  
  summary_neighbors <- data_bound %>%
    dplyr::select(y_var, neighbors) %>%
    dplyr::group_by(y_var) %>%
    dplyr::count(neighbors) %>%
    dplyr::rename(n_neighbors = n) %>%
    dplyr::mutate(neighbors_count = 
      paste0(neighbors, "(", n_neighbors, ")")) %>% 
    dplyr::arrange(desc(n_neighbors)) %>%
    dplyr::select(y_var, neighbors_count)
  
  committees_count_chr <- paste(summary_committees$committees_count,
    sep = "; ", collapse = "; ")
  neighbors_count_chr <- paste(summary_neighbors$neighbors_count,
    sep = "; ", collapse = "; ")
  
  rules_all <- unlist(data_bound$rules)
  rules_ranked <- sort(table(rules_all), decreasing = TRUE)
  rules_majority_chr <- names(rules_ranked[1])
  
  tibble::tibble(
    variable = unique(data_bound$y_var),
    committees = committees_count_chr,
    neighbors = neighbors_count_chr,
    rules = rules_majority_chr
  )
}

# Custom rsample nested resampling assessment indices conversion to
# caret `trainControl` `index` interface
rsample2caret_nested <- function(object, data = "analysis") {
  data <- match.arg(data)
  in_ind <- purrr::map(object$splits, as.integer, data = "analysis")
  names(in_ind) <- paste0("Fold", sprintf("%02d", 1L:nrow(object)))
  in_ind
}

add_repeat <- function(lst,
                       repeat_colname = Repeat,
                       repeat_string = "Repeat") {
  repeat_colnm <- rlang::enquo(repeat_colname)
  purrr::imap(lst, 
    ~ dplyr::mutate(.x, !!repeat_colnm := paste0(repeat_string, .y))
  )
}

# Adapt fold `in` indices because there is a shift of - 1 starting from
shift_fold_idx <- function(x, idx_from) {
  purrr::map(.x = x,
    ~ {.x[.x > idx_from] <- .x[.x > idx_from] - 1L; .x})
}

summarize_outer_tuning_vars <- function(data) {
  col_syms <- rlang::syms(names(data))
  vars_chr <- purrr::map_chr(col_syms, rlang::quo_name)
  
  purrr::map_at(
    .x = data, 
    .at = vars_chr,
    .f = ~ summarize_outer_tuning(data = .x)) %>%
    dplyr::bind_rows()
}


# Extract # of committees, neighbors and rules ---------------------------------

cubist_tuning_params <- function(model, method = "A") {
  rules_all <- rules_per_committee(model$finalModel)
  awc_rules_ranked <- sort(table(rules_all), decreasing = TRUE)
  rule <- names(awc_rules_ranked)[1]
  
  best_tune <- model$bestTune
  
  tibble::tibble(
    method = method,
    committees = best_tune$committees,
    neighbors = best_tune$neighbors,
    rules = rule)
}

freq_highest <- function(x) {
  ranked <- sort(table(x), decreasing = TRUE)
  names(ranked)[1]
}

cubist_tuning_params_frequent <- function(x) {
  tibble::tibble(
    method = x$method[1],
    committees = freq_highest(x$committees),
    neighbors = freq_highest(x$neighbors),
    rules = freq_highest(x$rules))
}
