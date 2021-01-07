# Austin Schumacher
# Aug 20, 2020
# Multivariate modeling of aggregated India MDS data

rm(list = ls())
library(SUMMER);
library(milliondeaths); library(tidyverse);
library(INLA);library(mapmisc);library(spdep);
library(ggpubr); library(broom); library(geofacet);
library(knitr); library(kableExtra); library(cowplot);

# set wd
setwd("/Users/austin/Desktop/india_mds_csmr/sandbox")

#### functions used in this code ####

# fit the INLA model and get results
fitModel <- function(data, reference_cause, formula, alphavar, quantiles) {
  
  # testing
  # data <- dat.df.agg.long
  # reference_cause <- reference_cause
  # formula <- f.ptrick
  # alphavar <- alpha_var
  # quantiles <- res.quantiles
  
  # set causenames
  causenames_v <- unique(data$cause)
  ncause <- length(causenames_v)
  
  # format data
  if (!is.null(data$year)) {
    data$year_num <- data$year - min(data$year) + 1
    data$year_num_iid <- data$year_num
    data$year_num_rw <- data$year_num
    data$year_cent <- data$year - mean(data$year)
  }
  
  # make variables for each variable and cause in order to specify REs
  for (varname in c("year", "year_num", "year_num_iid", "year_num_rw", "state_code", "state.year", "state.int", "year.int")) {
    if (!(varname %in% names(data))) next
    for (causeind in 1:ncause) {
      vnamecause <- paste(varname, causenames_v[causeind], sep = "_")
      data[[vnamecause]] <- data[[varname]]
      data[data$cause != causenames_v[causeind], vnamecause] <- NA
    }
  }
  
  # make cause a numeric vector, with cause 1 equal to the other causes 
  data$cause_factor <- as.factor(data$cause)
  
  fe.prec.pois <- list(prec.intercept = 100000000,
                       prec = list(default = 0, 
                                   year_cent = 100000000))
  cat(paste("\nRunning INLA model...")); flush.console()
  mod.ptrick <- inla(as.formula(formula),
                     control.fixed = fe.prec.pois,
                     family = "poisson",
                     data = data,
                     quantiles = quantiles,
                     control.predictor=list(compute=TRUE),
                     control.compute=list(config = TRUE))
  cat(paste(" done")); flush.console()
  
  # posterior sample to get CSMFs
  cat(paste("\nSampling from posterior...")); flush.console()
  resultsCSMF <- getPostCSMF(data = data, 
                             reference_cause = reference_cause,
                             model = mod.ptrick, 
                             alpha.effect = alphavar, 
                             nsamples = 1000, 
                             quantiles = quantiles)
  cat(paste(" done")); flush.console()
  return(list(model = mod.ptrick,
              posteriorCSMFsummary = resultsCSMF))
}

# get posterior distribution of CSMFs
getPostCSMF <- function(data, reference_cause, model, alpha.effect, nsamples, quantiles) {
  
  # testing
  # nsamples <- 100
  # alpha.effect <- alpha_var
  # data <- dat.df.agg.long
  # model <- results$model
  # quantiles <- res.quantiles
  # reference_cause <- reference_cause
  
  x.samples <- inla.posterior.sample(nsamples, model, use.improved.mean = TRUE, seed = 80085)
  x.contents <- model$misc$configs$contents
  effect <- "Predictor"
  id.effect <- which(x.contents$tag==effect)
  ind.effects <- x.contents$start[id.effect]-1 + (1:x.contents$length[id.effect])
  id.alphas <- which(x.contents$tag==alpha.effect)
  ind.alphas <- x.contents$start[id.alphas]-1 + (1:x.contents$length[id.alphas])
  expbetaresults <- matrix(NA, nrow = nrow(data), ncol = nsamples)
  for (i in 1:nsamples) {
    for (j in 1:ncause) {
      coi <- causenames_v[j]
      expbetaresults[which(data$cause == coi), i] <-
        exp(x.samples[[i]]$latent[ind.effects[which(data$cause == coi)]] - x.samples[[i]]$latent[ind.alphas])
    }
  }
  
  results <- matrix(NA, nrow = nrow(data), ncol = nsamples)
  idlevels <- unique(unlist(data[, alpha.effect]))
  for (i in 1:length(idlevels)) {
    idnum <- idlevels[i]
    denom <- apply(expbetaresults[which(data[, alpha.effect] == idnum),], 2, sum)
    for (j in 1:nsamples) {
      results[which(data[, alpha.effect] == idnum),j] <- 
        expbetaresults[which(data[, alpha.effect] == idnum),j]/denom[j]
    }
  }
  
  res.csmf <- as.data.frame(t(apply(results, 1, quantile, quantiles, na.rm = TRUE)))
  names(res.csmf) <- paste0("csmf", quantiles, "quant")
  res.csmf <- cbind(data[, c("cause", alpha.effect)], res.csmf)
  return(res.csmf)
}

# make table of FEs
makeTableFE <- function(model, format = 'markdown', digits = 5) {
  kable(model$summary.fixed,
        caption = "Post dist of FEs",
        booktabs = TRUE, format = format, digits = digits)
}

# make plot of CSMFs
makePlotCSMFs <- function(data, results, aggregationDimension, modelName, facet.grid = NULL) {
  
  # testing
  # data <- dat.df.agg.long
  # results <- results
  # aggregationDimension <- "year-state"
  # modelName <- model_name
  # facet.grid <- india_grid_MDS

  # merge on cause names for legend
  causecw <- unique(data[, c("cause", "cause_name")])
  res <- merge(results$posteriorCSMFsummary, causecw, by = "cause")
  
  # merge on state codes for plotting via geofacet if applicable
  if (aggregationDimension == "state") {
    crosswalk <- unique(data[, c("state", "state_code")])
    res <- merge(res, crosswalk, by = "state")
  } else if (aggregationDimension == "year-state") {
    crosswalk <- unique(data[, c("state.year", "state_code", "state", "year")])
    res <- merge(res, crosswalk, by = "state.year")
  }
  
  # make plots
  if (aggregationDimension == "year") {
    # plot csmfs
    p <- ggplot(data = res,
           aes(x = year_num, y = `csmf0.5quant`, color = cause_name, fill = cause_name)) +
      geom_line() +
      geom_ribbon(aes(ymin = `csmf0.025quant`, ymax = `csmf0.975quant`), alpha = 0.2, color = NA) +
      labs(x="year",
           y="Cause fraction",
           title = paste0("Estimated CSMF (", ci.level * 100, "% posterior intervals)"),
           subtitle = paste0(modelName, " model, data aggregated to year level")) +
      theme_light() + 
      theme(axis.text.x = element_blank(),
            legend.position = c(0.75, 0.95))
  } else if (aggregationDimension == "state") {
    # plot csmfs
    p <- ggplot(data = res, 
           aes(x = cause_name, y = `csmf0.5quant`, color = cause_name)) +
      geom_errorbar(aes(ymin = `csmf0.025quant`, ymax = `csmf0.975quant`), width = .1) +
      geom_point() +
      facet_geo(~ state, grid = facet.grid, scales = "fixed") +
      labs(x="cause", 
           y="Cause fraction",
           title = paste0("Estimated CSMF (", ci.level * 100, "% posterior intervals)"),
           subtitle = paste0(modelName, " model, data aggregated to state level")) +
      theme_light() + 
      theme(axis.text.x = element_blank(),
            legend.position = c(0.73, 0.95))
    return(p)
    
  } else if (aggregationDimension == "year-state") {
    # plot csmfs
    p <- ggplot(data = res, 
                aes(x = year, y = `csmf0.5quant`, color = cause_name, fill = cause_name)) +
      geom_area(position = "stack") +
      facet_geo(~ state, grid = facet.grid, scales = "fixed") +
      labs(x="year", 
           y="Cause fraction",
           title = paste0("Estimated CSMFs over time"),
           subtitle = paste0(modelName, " model, data aggregated to state-year level")) +
      theme_light() + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            legend.position = c(0.75, 0.95))
  }
  
  return(p)
}

# make plot of posterior median REs
makePlotPostREs <- function(data, results, reference.cause, aggregationDimension, modelName, causenames, facet.grid = NULL, value.y = FALSE, tertile = FALSE) {
  
  # testing
  # data <- dat.df.agg.long
  # results <- results
  # reference.cause <- reference_cause
  # aggregationDimension <- "year-state"
  # modelName <- model_name
  # causenames <- causenames_v
  # facet.grid <- india_grid_MDS
  # tertile <- FALSE
  # value.y = FALSE
  
  if (aggregationDimension == "year") {

  } else if (aggregationDimension == "state") {
    
    # formatting
    causenames_noref <- causenames[causenames != reference.cause]
    states <- unique(results$posteriorCSMFsummary$state_code)
    
    # format RE data
    re.post <- results$posteriorCSMFsummary[results$posteriorCSMFsummary$cause != reference.cause, c("cause", "state_code")]
    re.post$postMedRE <- NA
    for (i in 1:length(causenames_noref)) {
      cc <- causenames_noref[i]
      for (j in 1:length(states)) {
        ss <- states[j]
        re.post$postMedRE[re.post$cause == cc & re.post$state_code == ss] <- results$model$summary.random[[paste0("state_code_", cc)]]$`0.5quant`[j]
      }
    }
    
    # merge on state names
    re.post <- merge(re.post, unique(data[, c("state", "state_code")]))
    
    # make plots
    if (value.y == FALSE) {
      if (tertile == FALSE) {
        # plot REs gradient
        p1 <- ggplot(data = re.post, 
                     aes(x = cause, y = 0, color = `postMedRE`)) +
          scale_color_viridis_c(option = "C") +
          geom_point(cex = 4) +
          facet_geo(~ state, grid = facet.grid, scales = "free_y") +
          labs(x="cause", 
               y="Posterior median RE",
               title = "Posterior median RE",
               subtitle = paste0(modelName, " model, data aggregated to state level")) +
          theme_light() + 
          theme(axis.text.x = element_blank(),
                legend.position = c(0.75, 0.95))
      } else if (tertile == TRUE) {
        # get tertiles
        re.post$tertile <- cut(re.post$postMedRE, c(-Inf, quantile(re.post$postMedRE, c(1/3, 2/3)), Inf), )
        
        # plot REs gradient
        p1 <- ggplot(data = re.post, 
                     aes(x = cause, y = 0, color = tertile)) +
          scale_color_viridis_d(option = "C") +
          geom_point(cex = 4) +
          facet_geo(~ state, grid = facet.grid, scales = "free_y") +
          labs(x="cause", 
               y="Posterior median RE",
               title = "Posterior median RE",
               subtitle = paste0(modelName, " model, data aggregated to state level")) +
          theme_light() + 
          theme(axis.text.x = element_blank(),
                legend.position = c(0.75, 0.95))
      }
      
    } else if (value.y == TRUE) {
      if (tertile == FALSE) {
        # plot REs gradient
        p1 <- ggplot(data = re.post, 
                     aes(x = cause, y = `postMedRE`, color = `postMedRE`)) +
          scale_color_viridis_c(option = "C") +
          geom_point(cex = 4) +
          facet_geo(~ state, grid = facet.grid, scales = "free_y") +
          labs(x="cause", 
               y="Posterior median RE",
               title = "Posterior median RE",
               subtitle = paste0(modelName, " model, data aggregated to state level")) +
          theme_light() + 
          theme(axis.text.x = element_blank(),
                legend.position = c(0.75, 0.95))
      } else if (tertile == TRUE) {
        # get tertiles
        re.post$tertile <- cut(re.post$postMedRE, c(-Inf, quantile(re.post$postMedRE, c(1/3, 2/3)), Inf), )
        
        # plot REs gradient
        p1 <- ggplot(data = re.post, 
                     aes(x = cause, y = `postMedRE`, color = tertile)) +
          scale_color_viridis_d(option = "C") +
          geom_point(cex = 4) +
          facet_geo(~ state, grid = facet.grid, scales = "free_y") +
          labs(x="cause", 
               y="Posterior median RE",
               title = "Posterior median RE",
               subtitle = paste0(modelName, " model, data aggregated to state level")) +
          theme_light() + 
          theme(axis.text.x = element_blank(),
                legend.position = c(0.75, 0.95))
      }
      
    }
    
    return(p1)
    
  } else if (aggregationDimension == "year-state") {
    
    # formatting
    causenames_noref <- causenames[causenames != reference.cause]
    crosswalk <- unique(data[, c("state.year", "state_code", "state", "year")])
    res <- results$posteriorCSMFsummary %>% left_join(crosswalk, by = "state.year")
    
    # format RE data
    re.post <- res[res$cause != reference.cause, c("cause", "state.year", "state_code", "state", "year")]
    re.post$postMedYearIID <- re.post$postMedYearRW2 <- re.post$postMedBYM2i <- re.post$postMedBYM2s <- NA
    states <- unique(re.post$state)
    years <- unique(re.post$year)
    stateyears <- unique(re.post$state.year)
    for (i in 1:length(causenames_noref)) {
      cc <- causenames_noref[i]
      for (j in 1:length(stateyears)) {
        st <- stateyears[j]
        ss <- crosswalk$state[crosswalk$state.year == st]
        yy <- crosswalk$year[crosswalk$state.year == st]
        jy <- which(years == yy)
        js <- which(states == ss)
        
        # year IID re
        iidRE <- results$model$summary.random[[paste0("year_num_iid_", cc)]]$`0.5quant`[jy]
        if (is.null(iidRE)) {
          re.post$postMedYearIID[re.post$year == yy & re.post$cause == cc] <- 0
        } else {
          re.post$postMedYearIID[re.post$year == yy & re.post$cause == cc] <- iidRE
        }
        
        # year RW2 re
        re.post$postMedYearRW2[re.post$year == yy & re.post$cause == cc] <- results$model$summary.random[[paste0("year_num_rw_", cc)]]$`0.5quant`[jy]
        
        # state bym2 spatial
        u <- results$model$summary.random[[paste0("state_code_", cc)]]$`0.5quant`[js + length(states)]
        re.post$postMedBYM2s[re.post$state == ss & re.post$cause == cc] <- u
        
        # state bym2 iid
        tau <- results$model$summary.hyperpar[paste0("Precision for state_code_", cc), "0.5quant"]
        phi <- results$model$summary.hyperpar[paste0("Phi for state_code_", cc), "0.5quant"]
        re_val <- results$model$summary.random[[paste0("state_code_", cc)]]$`0.5quant`[js]
        re.post$postMedBYM2i[re.post$state == ss & re.post$cause == cc] <- (re_val * sqrt(tau) - sqrt(phi) * u)/(sqrt(1 - phi))
      }
      
      if (grepl(modelName, "type1")) {
        stop("no code for spatiotemporal interaction RE plots yet")
        
        # state-year interaction RE
        re.post$postMedSpacetime[re.post$cause == cc] <- results$model$summary.random[[paste0("state.year_", cc)]]$`0.5quant`
        
      } else if (grepl(modelName, "type4")) {
        stop("no code for spatiotemporal interaction RE plots yet")
      }
    }
    
    # make plots
    plot.list <- vector(mode = "list", length = 5)
    names(plot.list) <- c("yearIID", "yearRW2", "bym2spatial", "bym2IID", "spacetime")
    
    plot.list[["spacetime"]] <- NA
    
    plot.list[["yearIID"]] <- ggplot(data = re.post %>% filter(state == state[1]), 
                                     aes(x = year, y = `postMedYearIID`, color = cause)) +
      geom_line() +
      labs(x="year", 
           y="Posterior median RE",
           title = "Year IID RE posterior median",
           subtitle = paste0(modelName, " model, data aggregated to state-year level")) +
      theme_light() + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 4))
    
    plot.list[["yearRW2"]] <- ggplot(data = re.post %>% filter(state == state[1]), 
                                     aes(x = year, y = `postMedYearRW2`, color = cause)) +
      geom_line() +
      labs(x="year", 
           y="Posterior median RE",
           title = "Year RW2 RE posterior median",
           subtitle = paste0(modelName, " model, data aggregated to state-year level")) +
      theme_light() + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 4))
    
    if (value.y == FALSE) {
      if (tertile == FALSE) {
        plot.list[["bym2spatial"]] <- ggplot(data = re.post %>% filter(year == min(year)), 
                                             aes(x = cause, y = 0, color = `postMedBYM2s`)) +
          scale_color_viridis_c(option = "C") +
          geom_point(pch = 15, cex = 4.5) +
          facet_geo(~ state, grid = facet.grid, scales = "free_y") +
          labs(x="cause", 
               y="Posterior median RE",
               title = "BYM2 spatial RE posterior median",
               subtitle = paste0(modelName, " model, data aggregated to state-year level")) +
          theme_light() + 
          theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 4),
                legend.position = c(0.75, 0.95),
                axis.text.y = element_blank())
        plot.list[["bym2IID"]] <- ggplot(data = re.post %>% filter(year == min(year)), 
                                         aes(x = cause, y = 0, color = `postMedBYM2i`)) +
          scale_color_viridis_c(option = "C") +
          geom_point(pch = 15, cex = 4.5) +
          facet_geo(~ state, grid = facet.grid, scales = "free_y") +
          labs(x="cause", 
               y="Posterior median RE",
               title = "BYM2 IID RE posterior median",
               subtitle = paste0(modelName, " model, data aggregated to state-year level")) +
          theme_light() + 
          theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 4),
                legend.position = c(0.75, 0.95),
                axis.text.y = element_blank())
        
      } else if (tertile == TRUE) {
        # get tertiles
        re.post$tertileSpatial <- cut(re.post$postMedBYM2s, c(-Inf, quantile(re.post$postMedBYM2s, c(1/3, 2/3)), Inf), )
        re.post$tertileIID <- cut(re.post$postMedBYM2i, c(-Inf, quantile(re.post$postMedBYM2i, c(1/3, 2/3)), Inf), )
        
        plot.list[["bym2spatial"]] <- ggplot(data = re.post %>% filter(year == min(year)), 
                                             aes(x = cause, y = 0, color = tertileSpatial)) +
          scale_color_viridis_d(option = "C") +
          geom_point(pch = 15, cex = 4.5) +
          facet_geo(~ state, grid = facet.grid, scales = "free_y") +
          labs(x="cause", 
               y="Posterior median RE",
               title = "BYM2 spatial RE posterior median",
               subtitle = paste0(modelName, " model, data aggregated to state-year level")) +
          theme_light() + 
          theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 4),
                legend.position = c(0.75, 0.95),
                axis.text.y = element_blank())
        plot.list[["bym2IID"]] <- ggplot(data = re.post %>% filter(year == min(year)), 
                                         aes(x = cause, y = 0, color = tertileIID)) +
          scale_color_viridis_d(option = "C") +
          geom_point(pch = 15, cex = 4.5) +
          facet_geo(~ state, grid = facet.grid, scales = "free_y") +
          labs(x="cause", 
               y="Posterior median RE",
               title = "BYM2 IID RE posterior median",
               subtitle = paste0(modelName, " model, data aggregated to state-year level")) +
          theme_light() + 
          theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 4),
                legend.position = c(0.75, 0.95),
                axis.text.y = element_blank())
      }
      
    } else if (value.y == TRUE) {
      if (tertile == FALSE) {
        plot.list[["bym2spatial"]] <- ggplot(data = re.post %>% filter(year == min(year)), 
                                             aes(x = cause, y = `postMedBYM2s`, color = `postMedBYM2s`)) +
          scale_color_viridis_c(option = "C") +
          geom_point(pch = 15, cex = 4.5) +
          facet_geo(~ state, grid = facet.grid, scales = "free_y") +
          labs(x="cause", 
               y="Posterior median RE",
               title = "BYM2 spatial RE posterior median",
               subtitle = paste0(modelName, " model, data aggregated to state-year level")) +
          theme_light() + 
          theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 4),
                legend.position = c(0.75, 0.95),
                axis.text.y = element_blank())
        plot.list[["bym2IID"]] <- ggplot(data = re.post %>% filter(year == min(year)), 
                                         aes(x = cause, y = `postMedBYM2i`, color = `postMedBYM2i`)) +
          scale_color_viridis_c(option = "C") +
          geom_point(pch = 15, cex = 4.5) +
          facet_geo(~ state, grid = facet.grid, scales = "free_y") +
          labs(x="cause", 
               y="Posterior median RE",
               title = "BYM2 IID RE posterior median",
               subtitle = paste0(modelName, " model, data aggregated to state-year level")) +
          theme_light() + 
          theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 4),
                legend.position = c(0.75, 0.95),
                axis.text.y = element_blank())
      } else if (tertile == TRUE) {
        # get tertiles
        re.post$tertileSpatial <- cut(re.post$postMedBYM2s, c(-Inf, quantile(re.post$postMedBYM2s, c(1/3, 2/3)), Inf), )
        re.post$tertileIID <- cut(re.post$postMedBYM2i, c(-Inf, quantile(re.post$postMedBYM2i, c(1/3, 2/3)), Inf), )
        
        plot.list[["bym2spatial"]] <- ggplot(data = re.post %>% filter(year == min(year)), 
                                             aes(x = cause, y = `postMedBYM2s`, color = tertileSpatial)) +
          scale_color_viridis_d(option = "C") +
          geom_point(pch = 15, cex = 4.5) +
          facet_geo(~ state, grid = facet.grid, scales = "free_y") +
          labs(x="cause", 
               y="Posterior median RE",
               title = "BYM2 spatial RE posterior median",
               subtitle = paste0(modelName, " model, data aggregated to state-year level")) +
          theme_light() + 
          theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 4),
                legend.position = c(0.75, 0.95),
                axis.text.y = element_blank())
        plot.list[["bym2IID"]] <- ggplot(data = re.post %>% filter(year == min(year)), 
                                         aes(x = cause, y = `postMedBYM2i`, color = tertileIID)) +
          scale_color_viridis_d(option = "C") +
          geom_point(pch = 15, cex = 4.5) +
          facet_geo(~ state, grid = facet.grid, scales = "free_y") +
          labs(x="cause", 
               y="Posterior median RE",
               title = "BYM2 IID RE posterior median",
               subtitle = paste0(modelName, " model, data aggregated to state-year level")) +
          theme_light() + 
          theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 4),
                legend.position = c(0.75, 0.95),
                axis.text.y = element_blank())
      }
      
    }
    
    return(plot.list)
  }
}

#### set up MDS info and directories ####

# MDS username/password
username <- "schumachera"
password <- "schumachera2402"

# dropbox directory to save data
datadir <- "../../../Dropbox/csmr_india/data/mds"
savedir <- "../../../Dropbox/csmr_india/sandbox"

#### parameters ####

# set years and ages for analysis
years <- 2004:2013
ages <- c(0,1,5)

# dimensions for data aggregation
aggregationDims <- c("none", "state", "year", "year-state")
# aggregationDims <- c("year")

# set CI level and quantiles for analysis
ci.level <- 0.95
res.quantiles <- c((1-ci.level) - (1-ci.level)/2, 0.5, ci.level + (1-ci.level)/2)

# set reference cause
reference_cause <- "cause_1A"

##### load aggregated data (created from 'aggregate_format_data_for_exploration.R') #####
# this loads:
# - dat.df.agg.list: named list where each element is an aggregated data frame of MDS deaths by cause
# - allChildCodes: data frame with all child COD codes
# - allChildCodesBroad: data frame with broad child COD codes
# - india: shapefile for Indian states/regions
# - india.adj: neighborhood file for india states/regions for use with INLA spatial models
# - india.toplot: tidy() version of India map shapefile for plotting with ggplot
# - india_grid_MDS: grid file for use with geo_facet plots in ggplot with the geofacet package

load(paste0(datadir, "/aggregated_formatted_data_lists_and_codes_and_mapping_data.RData"))

# set cause name references
cause_ref <- data.frame(cause = c("cause_1A",
                                  "cause_1B",
                                  "cause_1C01",
                                  "cause_1C02",
                                  "cause_1D",
                                  "cause_2",
                                  "cause_3",
                                  "cause_4"),
                        cause_name = c("Pneumonia, Acute bacterial sepsis, severe infections, Meningitis/encephalitis",
                                       "Diarrhoeal, TB, Tetanus, Poliomyelitis, Measles, HIV/AIDS, Malaria, Other infectious, unknown fever",
                                       "Prematurity, low birthweight",
                                       "Birth asphyxia, birth trauma",
                                       "Nutritional, congenital anomalies",
                                       "Other non-communicable",
                                       "Injuries",
                                       "Ill-defined or cause unknown"))

##### single cause binomial models (with and without Poisson Trick) #####

# vector of cause names
causenames_v <- sort(allChildCodesBroad$cause)
ncause <- length(causenames_v)

# create list of models to run for each aggregation dimension
aggDimModels <- vector(mode = "list", length = length(aggregationDims))
names(aggDimModels) <- aggregationDims
for (i in 1:length(aggDimModels)) {
    ad <- aggregationDims[i]
    if (ad == "none") aggDimModels[[i]] <- c("int-only")
    if (ad == "state") aggDimModels[[i]] <- c("iid", "bym2")
    if (ad == "year") aggDimModels[[i]] <- c("int-only", "rw2")
    if (ad == "year-state") aggDimModels[[i]] <- c("bym2_rw2_no-int", "bym2_rw2_type1", "bym2_rw2_type4")
}

############
## Multivariate modeling
############

##### Intercept only (yearly data) #####

# set aggregation dimension
# aggdim <- "year"
# aggdim <- "state"
aggdim <- "year-state"

# set ID variables
if (aggdim == "none") {
  idvariables <- "id"
} else if (aggdim == "year") {
  idvariables <- "year"
} else if (aggdim == "state") {
  idvariables <- c("state", "state_code")
} else if (aggdim == "year-state") {
  idvariables <- c("year", "state", "state_code", "state.year", "state.int", "year.int")
}

# set data to use
mod_data <- dat.df.agg.list[[aggdim]]

# set model
model_name <- aggDimModels[[aggdim]]

# set model to run
# model_name <- "rw2"
# model_name <- "iid"
# model_name <- "int-only"
# model_name <- "bym2"
model_name <- "bym2_rw2_no-int"

# format data
mod_data$id <- 1:nrow(mod_data)

# make long data for Poisson model
dat.df.agg.long <- mod_data %>%
  pivot_longer(cols = c(all_of(causenames_v)),
               names_to = "cause",
               values_to = "deaths") %>%
  filter(cause %in% c(causenames_v)) %>%
  dplyr::select(all_of(idvariables), cause, deaths)

# merge on cause names
dat.df.agg.long <- dat.df.agg.long %>% left_join(cause_ref, by = "cause")

# set index variable for alpha parameters
if (aggdim == "none") {
  alpha_var <- "id"
} else if (aggdim == "year") {
  alpha_var <- "year_num"
} else if (aggdim == "state") {
  alpha_var <- "state_code"
} else if (aggdim == "year-state") {
  alpha_var <- "state.year"
}

# make formula for model
f.ptrick <- "deaths ~ 1 + cause_factor"
if (grepl("iid", model_name) & grepl("state", aggdim)) {
  for (i in 2:ncause) {
    f.ptrick <- paste0(f.ptrick, 
                       "+ f(state_code_", causenames_v[i], ", model = 'iid',
                            hyper = list(prec = list(prior = 'pc.prec', 
                                                     param = c(0.3, 0.01))))")
  }
}
if (grepl("bym2", model_name) & grepl("state", aggdim)) {
  for (i in 2:ncause) {
    f.ptrick <- paste0(f.ptrick, 
                       "+ f(state_code_", causenames_v[i], ", model = 'bym2',
                                     graph=india.adj, 
                                     scale.model=T, 
                                     constr=T,
                                     hyper=list(phi=list(prior='pc', 
                                                         param=c(0.5, 0.5), 
                                                         initial=1),
                                                prec=list(prior='pc.prec', 
                                                          param=c(0.3,0.01), 
                                                          initial=5)))")
  }
}
if (grepl("rw2", model_name) & grepl("year", aggdim)) {
  for (i in 2:ncause) {
    f.ptrick <- paste0(f.ptrick, 
                       "+ f(year_num_iid_", causenames_v[i], ", model = 'iid',
                            hyper = list(prec = list(prior = 'pc.prec', 
                                                     param = c(0.3, 0.01))))")
    f.ptrick <- paste0(f.ptrick, 
                       "+ f(year_num_rw_", causenames_v[i], ", model = 'rw2',
                            constr = TRUE,
                            hyper = list(prec = list(prior = 'pc.prec', 
                                                     param = c(0.3, 0.01))))")
  }
}
if (grepl("type1", model_name) & grepl("state", aggdim) & grepl("year", aggdim)) {
  for (i in 2:ncause) {
    f.ptrick <- paste0(f.ptrick, 
                       "+ f(state.year_", causenames_v[i], ", model = 'iid',
                             hyper = list(prec = list(prior = 'pc.prec',
                                                      param = c(0.5, 0.01))))")
  }
} else if (grepl("type4", model_name) & grepl("state", aggdim) & grepl("year", aggdim)) {
  for (i in 2:ncause) {
    f.ptrick <- paste0(f.ptrick, 
                       "+ f(state.int_", causenames_v[i], ", model = 'besag', graph = india.adj, constr = TRUE,
                             group = year.int_", causenames_v[i], ", control.group = list(model = 'rw2'),
                             hyper = list(prec = list(prior = 'pc.prec',
                                                      param = c(0.5, 0.01))))")
  }
}
f.ptrick <- paste0(f.ptrick, 
                  " + f(", alpha_var, ", model = 'iid',
                         hyper = list(prec = list(initial = log(0.000001),
                                                  fixed = TRUE)))")

# run model
results <- fitModel(data = dat.df.agg.long, 
                    reference_cause = reference_cause, 
                    formula = f.ptrick, 
                    alphavar = alpha_var, 
                    quantiles = res.quantiles)

# make tables
makeTableFE(results$model)

# plot CSMFs
p1 <- makePlotCSMFs(data = dat.df.agg.long,
                    results = results, 
                    aggregationDimension = aggdim, 
                    modelName = model_name,
                    facet.grid = india_grid_MDS)

# save plots
if (aggdim == "state") {
  ggexport(p1, filename = paste0(savedir, 
                                 "/graphs/spatial_only/multiple_causes/csmf_posterior_", 
                                 model_name, 
                                 ".pdf"),
           width = 12, height = 12)
} else if (aggdim == "year") {
  
} else if (aggdim == "year-state") {
  ggexport(p1, filename = paste0(savedir, 
                                 "/graphs/spatiotemporal/multiple_causes/csmf_posterior_", 
                                 model_name, 
                                 ".pdf"),
           width = 12, height = 12)
} 

# plot posterior median REs
p2 <- makePlotPostREs(data = dat.df.agg.long,
                      results = results,
                      reference.cause = reference_cause,
                      aggregationDimension = aggdim,
                      modelName = model_name,
                      causenames = causenames_v,
                      facet.grid = india_grid_MDS,
                      tertile = FALSE,
                      value.y = FALSE)
p3 <- makePlotPostREs(data = dat.df.agg.long,
                      results = results,
                      reference.cause = reference_cause,
                      aggregationDimension = aggdim,
                      modelName = model_name,
                      causenames = causenames_v,
                      facet.grid = india_grid_MDS,
                      tertile = TRUE,
                      value.y = FALSE)

# save plots
if (aggdim == "state") {
  ggexport(p2, filename = paste0(savedir, 
                                 "/graphs/spatial_only/multiple_causes/re_posteriorMedian_", 
                                 model_name, 
                                 ".pdf"),
           width = 12, height = 12)
  ggexport(p3, filename = paste0(savedir, 
                                 "/graphs/spatial_only/multiple_causes/re_posteriorMedian_", 
                                 model_name, 
                                 "_tertile.pdf"),
           width = 12, height = 12)
} else if (aggdim == "year") {
  
} else if (aggdim == "year-state") {
  pdf(paste0(savedir, 
             "/graphs/spatiotemporal/multiple_causes/re_posteriorMedian_", 
             model_name, 
             ".pdf"),
      onefile = TRUE,
      width = 12, height = 12)
  
  for (i in 1:4) {
    print(p2[[i]])
  }
  
  dev.off()
  
  pdf(paste0(savedir, 
             "/graphs/spatiotemporal/multiple_causes/re_posteriorMedian_", 
             model_name, 
             "_tertile.pdf"),
      onefile = TRUE,
      width = 12, height = 12)
  
  for (i in 1:4) {
    print(p3[[i]])
  }
  
  dev.off()
  
  
} 

