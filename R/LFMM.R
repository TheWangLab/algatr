
#' Run LFMM
#'
#' @param gen genotype matrix
#' @param env dataframe with environmental data
#' @param coords dataframe with coordinates (only needed if K selection is performed with TESS)
#' @param K number of latent factors (if left blank, K value selection will be conducted)
#' @param lfmm_method lfmm method (either \code{"ridge"} or \code{"lasso"})
#' @param k_selection method for performing k selection (can either by "tracy.widom", "quick.elbow", or "tess")
#' @param Kvals values of K to test if using "tess" method of K selection
#' @param sig alpha level for determining candidate loci
#'
#' @return
#' @export
#'
#' @examples
lfmm_run <- function(gen, env, coords = NULL, K = NULL, lfmm_method = "ridge", k_selection = "tracy.widom", Kvals = 1:20, sig = 0.05){
  
  # PCA to determine number of latent factors
  # if K is not specified it is calculated based on given K selection method
  if(is.null(K)){
    K <- get_K(gen, coords, k_selection = k_selection, Kvals = Kvals)
  }
  
  # gen matrix
  genmat <- as.matrix(gen)
  # env matrix
  envmat <- as.matrix(env)
  
  # remove NAs
  # FIX THIS
  # NOTE: IF YOU HAVE NAs YOU WILL GET WEIRD ERRORS
  genmat <- genmat[complete.cases(envmat),]
  envmat <- envmat[complete.cases(envmat),]
  
  
  # All Env
  # run model
  if(lfmm_method == "ridge"){lfmm_mod <- lfmm_ridge(genmat, envmat, K = K)}
  if(lfmm_method == "lasso"){lfmm_mod <- lfmm_lasso(genmat, envmat, K = K)}
  
  # performs association testing using the fitted model:
  pv <- lfmm_test(Y = genmat, 
                  X = envmat, 
                  lfmm = lfmm_mod, 
                  calibrate = "gif")
  # adjust pvalues
  pval_df <- as_tibble(pv$calibrated.pvalue)
  pvalues <- map_df(pval_df, p.adjust, method = "fdr")
  
  # check qqplots
  par(mfrow = c(1, ncol(pvalues)))
  print(lfmm_qqplot(pvalues))

  # get list of candidate loci and make manhattan plot
  # FIGURE OUT HOW TO ADD PLOT TITLES
  cand_loci <- map(pvalues, get_loci, sig = sig)
  # make plots
  print(lfmm_manhattanplot(pvalues, sig))
  
  # make list of results 
  results <- list(K = K,
                  loci = cand_loci,
                  pvalues = pvalues)
  
  return(results)
}


#' Function to Select K value
#'
#' @param gen genotype matrix
#' @param k_selection method for performing k selection (can either by "tracy.widom", "quick.elbow", or "tess")
#' @param coords coordinates for TESS based K selection
#'
#' @return
#' @export
#'
#' @examples
get_K <- function(gen, coords = NULL, k_selection = "tracy.widom", Kvals = Kvals, ...){
  
  if(k_selection == "tracy.widom"){K <- get_K_tw(gen)}
  
  if(k_selection == "quick.elbow"){K <- get_K_elbow(gen)}
  
  if(k_selection == "tess"){K <- get_K_tess(gen, coords, Kvals)}
  
  return(K)
}

#' @describeIn get_K Determine K using Tracy-Widom Test
#' @param gen genotype matrix
#'
#' @return
#' @export
#'
#' @examples
get_K_tw <- function(gen){
  # run pca
  pc <- prcomp(gen)
  
  # get eig
  eig <- pc$sdev^2
  
  # run tracy widom test
  # NOTE: 	
  # the critical point is a numeric value corresponding to the significance level. 
  # If the significance level is 0.05, 0.01, 0.005, or 0.001, 
  # the criticalpoint should be set to be 0.9793, 2.0234, 2.4224, or 3.2724, accordingly. 
  # The default is 2.0234.
  tw_result <- AssocTests::tw(eig, eigenL = length(eig), criticalpoint = 0.9793)
  
  # get K based on number of significant eigenvalues
  K <- tw_result$SigntEigenL
  
  return(K)
}

#' @describeIn get_K Determine K using PCA and \code{quick.elbow}
#' 
#' @param gen genotype matrix
#'
#' @return
#' @export
#'
#' @examples
get_K_elbow <- function(gen){
  # run pca
  pc <- prcomp(gen)
  
  # get eig
  eig <- pc$sdev^2
  # estimate number of latent factors using quick.elbow (see general functions for description of how this function works)
  # this is a crude way to determine the number of latent factors that is based on an arbitrary "low" value 
  K <- quick.elbow(eig, low = 0.08, max.pc = 0.9)
  
  par(pty = "s",mfrow = c(1,1))
  plot(eig, xlab = 'PC', ylab = "Variance explained")
  abline(v = K, col = "red", lty = "dashed")
  
  return(K)
}


#' @describeIn get_K Determine K using TESS and \code{bestK}
#' 
#' @param gen genotype matrix
#' @param coords coordinates for "tess"
#' @param Kvals values of K to test for "tess"
#' @param tess_method method to use for "tess"
#' @param ploidy ploidy for "tess"
#'
#' @return
#' @export
#'
#' @examples
get_K_tess <- function(gen, coords, Kvals = 1:10, tess_method = "projected.ls", ploidy = 2){
  # run tess for all K values
  tess3_obj <- tess3(X = gen, coord = coords, K = Kvals, method = tess_method, ploidy = ploidy) 
  
  # plot CV results and mark the K-value automatically selected  
  plot(tess3_obj, pch = 19, col = "blue",
       xlab = "Number of ancestral populations",
       ylab = "Cross-validation score")
  
  # get best K value
  K <-  bestK(tess3_obj, Kvals)
  
  return(K)
}

#' LFMM QQplot
#'
#' @param pvalues pvalue dataframe 
#'
#' @return
#' @export
#'
#' @examples
lfmm_qqplot <- function(pvalues){
  
  pvalues <- -log10(pvalues)
  pvalues$loci <- 1:nrow(pvalues)
  pvalues_tidy <- pvalues %>% gather(env, p, -loci)
  
  plt <- ggplot(pvalues_tidy, aes(sample = p)) +
    stat_qq() +
    geom_abline(intercept = 0, slope = 1) +
    facet_wrap( ~ env, nrow = 1) +
    labs(x = NULL, 
         y = "-log10(p)") + 
    theme_bw() +
    theme(line = element_blank())
  
  return(plt)
}

#' LFMM Manhattan Plot
#'
#' @param pvalues pvalue dataframe
#' @param sig alpha level
#'
#' @return
#' @export
#'
#' @examples
lfmm_manhattanplot <- function(pvalues, sig){
  pvalues10 <- -log10(pvalues)
  pvalues10$loci <- 1:nrow(pvalues10)
  pvalues_tidy <- pvalues10 %>% gather(env, p, colnames(pvalues))
  
  plt <- ggplot(pvalues_tidy, aes(x = loci, y = p)) +
    geom_hline(yintercept = -log10(sig), color = "red", linetype = "dashed") + 
    geom_point(alpha = 0.75) +
    facet_wrap( ~ env, nrow = 1) +
    labs(x = NULL, 
         y = "-log10(p)") + 
    theme_bw() +
    theme( 
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank())
  
  return(plt)
}

#' Get candidate loci row numbers
#'
#' @param pvec vector of pvalues
#' @param sig alpha level
#'
#' @return
#' @export
#'
#' @examples
get_loci <- function(pvec, sig){
  loci <- which(pvec < sig) 
  return(loci)
}


# quickly choose an elbow for a PC. 
# at variance below 5% per component, choose the largest % drop
# designed for variance percentages, but will also work given a full set of Evalues
#' Quickly estimate the 'elbow' of a scree plot (PCA)
#' 
#' This function uses a rough algorithm to estimate a sensible 'elbow' to
#' choose for a PCA scree plot of eigenvalues. The function looks at an initial arbitrarily 'low'
#' level of variance and looks for the first eigenvalue lower than this. If the very first eigenvalue
#' is actually lower than this (i.e, when the PCs are not very explanatory) then this 'low' value is
#' iteratively halved until this is no longer the case. After starting below this arbitrary threshold
#' the drop in variance explained by each pair of consecutive PCs is standardized by dividing over the 
#' larger of the pair. The largest percentage drop in the series below 'low' % is selected as the 'elbow'.
#' @param varpc numeric, vector of eigenvalues, or 'percentage of variance' explained datapoints for
#'  each principle component. If only using a partial set of components, should first pass to 
#'  estimate.eig.vpcs() to estimate any missing eigenvalues.
#' @param low numeric, between zero and one, the threshold to define that a principle component
#'  does not explain much 'of the variance'.
#' @param max.pc maximum percentage of the variance to capture before the elbow (cumulative sum to PC 'n')
#' @return The number of last principle component to keep, prior to the determined elbow cutoff
#' @export
#' @seealso \code{\link{estimate.eig.vpcs}}
#' @author Nicholas Cooper 
#' @examples
#' # correlated data
#' mat <- sim.cor(100,50)
#' result <- princomp(mat)
#' eig <- result$sdev^2
#' elb.a <- quick.elbow(eig)
#' pca.scree.plot(eig,elbow=elb.a,M=mat) 
#' elb.b <- quick.elbow(eig,low=.05) # decrease 'low' to select more components
#' pca.scree.plot(eig,elbow=elb.b,M=mat) 
#' # random (largely independent) data, usually higher elbow #
#' mat2 <- generate.test.matrix(5,3)
#' result2 <- princomp(mat2)
#' eig2 <- result2$sdev^2
#' elb2 <- quick.elbow(result2$sdev^2)
#' pca.scree.plot(eig2,elbow=elb2,M=mat2)
quick.elbow <- function(varpc,low=.08,max.pc=.9) {
  ee <- varpc/sum(varpc) # ensure sums to 1
  #print(round(log(ee),3))
  while(low>=max(ee)) { low <- low/2 } # when no big components, then adjust 'low'
  lowie <- (ee<low) ; highie <- ee>low/8
  low.ones <- which(lowie & highie)
  others <- length(which(!lowie))
  if(length(low.ones)>0) {
    if(length(low.ones)==1) {
      elbow <- low.ones 
    } else {
      set <- ee[low.ones]
      pc.drops <- abs(diff(set))/(set[1:(length(set)-1)])
      infz <- is.infinite(pc.drops)
      #print(pc.drops)
      elbow <- which(pc.drops==max(pc.drops[!infz],na.rm=T))[1]+others
    }
  } else { 
    # if somehow there are no small eigenvalues, just choose the elbow as the second last
    cat("no eigenvalues weresignificantly smaller than the previous\n")
    elbow <- length(ee) 
  }
  if(tail(cumsum(ee[1:elbow]),1)>max.pc) {
    elbow <- which(cumsum(ee)>max.pc)[1]-1
  }
  if(elbow<1) {
    warning("elbow calculation failed, return zero")
    return(0)
  }
  names(elbow) <- NULL
  return(elbow)
}

#' Best K Selection based on cross entropy
#'
#' @param tess3_obj list produced by \code{\link{tess3}}
#' @param Kvals vector of K values for testing
#'
#' @note  (source: https://chazhyseni.github.io/NALgen/post/determining_bestk/)
#' @return
#' @export
#'
#' @examples
bestK <- function(tess3_obj, Kvals){
  ce <- list()
  for(k in Kvals) ce[[k]] <- tess3_obj[[k]]$crossentropy
  ce.K <- c()
  for(k in Kvals) ce.K[k] <- min(ce[[k]])
  diff <- ce.K[-1] - ce.K[-max(Kvals)]
  slope <- exp(-diff) - 1
  #K is selected based on the smallest slope value in the upper quartile
  K <- min(which(slope <= quantile(slope)[4]))
  return(K)
}

