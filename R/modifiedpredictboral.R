## Intervals from marginal predictions will typically be wider than those from conditional predictions, because the former takes into account the additional uncertainty due to the lv being unknown. 
## predict.type = "conditional": Predictions conditional on the latent variables observed and everything else e.g., random row effects. 
## predict.type = "marginal": marginalizes from the linear predictor down i.e., LVs and random row effects if included. Does not marginalize over species random effects as it is not clear that you need to take into account their uncertainty i.e., you are not predicting to new species! 

## Modification for LHR:
## - added noroweff argument to allow one to "turn off" row effects" by not adding them into the linear predictor, or equivalently set them all to zero. Defaults to FALSE.
## - added return.alllinpred argument to return all MCMC samples from prediction. Defaults to FALSE. This will actually be in the next version of boral.
## - Allow new trait matrix, but the problem is that you don't know the loadings for these new species. After some discussion with colleagues, the best solution we have is ad-hoc: either take the loadings as zero or some median or mean across the observed species. 
###################


function() 
     {
object <- R.nb.phyloonX.boral
newX = age.new
#newtraits <- matrix(colMeans(traits), nrow = 1) ## For a hypothetical taxa
newtraits = NULL
lvcoefs.type <- "zero"
noroweff = TRUE
newrow.ids = NULL
distmat = NULL
predict.type = "marginal"
est = "median"
prob = 0.95
lv.mc = 1000
return.alllinpred = FALSE
     }
     
     
	
predict.boral <- function(object, newX = NULL, newtraits = NULL, lvcoefs.type = "zero", noroweff = FALSE, newrow.ids = NULL, 
	distmat = NULL, predict.type = "conditional", est = "median", prob = 0.95, lv.mc = 1000, 
	return.alllinpred = FALSE, ...) 
     {  
     
     num.lv <- object$num.lv
     predict.type <- match.arg(predict.type, choices = c("conditional","marginal"))
     lvcoefs.type <- match.arg(lvcoefs.type, choices = c("zero","mean","median"))
    
     if(predict.type == "marginal") 
          {
          message("Marginal predictions take a long time, because there is a lot of (Monte-Carlo) integration involved. Apologies in advance!")
          if(num.lv == 0) 
               {
               message("Please note if there are no latent variables in the model, then marginal and conditional predictions are equivalent")
               predict.type <- "conditional"
               }
          }
     if(predict.type == "conditional" & !is.null(newrow.ids)) 
          { 
          message("For conditional predictions, newrow.ids is ignored since predictions are made conditional on the set of row effects i.e., on the same set of sites")
          newrow.ids <- NULL
          }
        
     ## Check properties of X
     if(!is.null(newX)) 
          {
          X <- as.matrix(newX) 
          if(is.null(object$X.coefs.mean)) 
               stop("Cannot find coefficients for X in object, even though you supplied newX")
          if(ncol(object$X.coefs.mean) != ncol(newX)) 
               stop("Number of columns in newX does not match number of columns in object$X.coefs.mean")
          if(predict.type == "conditional") 
               {
               if(object$row.eff != "none" & noroweff == FALSE) 
                    { 
                    if(nrow(object$row.ids) != nrow(X)) 
                         stop("For conditional predictions, the number of rows in newX must be equal to number of rows in object$row.ids")
                    }
               if(num.lv > 0) 
                    { 
                    if(nrow(object$lv.mean) != nrow(X)) 
                         stop("For conditional predictions, the number of rows in newX must be equal to number of rows in object$lv.mean") 
                    }
               }
          }
     if(is.null(newX)) 
          X <- object$X 
     n <- nrow(X) 

     ## Check properties of newrow.ids; this should only be activated once the predictions are marginal 
     if(!is.null(newrow.ids) & noroweff == FALSE) 
          { 
          if(object$row.eff == "none") 
               stop("Cannot find row effects parameters in object, even though you supplied in newrow.ids")
    
          newrow.ids <- as.matrix(newrow.ids)
          if(is.null(colnames(newrow.ids))) 
               colnames(newrow.ids) <- paste0("ID", 1:ncol(newrow.ids))
          if(ncol(object$row.ids) != ncol(newrow.ids)) 
               stop("The number of columns in newrow.ids must be equal to number of columns in object$row.ids")
          for(k in 1:ncol(newrow.ids)) 
               {
               if(!all(unique(newrow.ids[,k]) %in% unique(object$row.ids[,k])))
                    stop(paste0("Not all levels of newrow.ids[,",k,"] can be found in object$row.ids[,",k,"]. This is a problem as then thre will be some IDs which are unknown (as based on object)"))
               }
          }     
     if(is.null(newrow.ids) & noroweff == FALSE) 
          newrow.ids <- object$row.ids 
     if(!is.null(X) & !is.null(newrow.ids) & noroweff == FALSE) 
          {
          if(n != nrow(newrow.ids)) 
               stop("Number of rows in X does not match number of rows in newrow.ids")
          }

          
     ## Check properties of newtraits
     if(!is.null(newtraits))
          {
          if(is.null(object$traits))
               stop("The fitted boral model does not involve traits; please set newtraits = NULL instead")
          if(!is.null(object$traits)) 
               {  
               if(!is.matrix(newtraits)) 
                    stop("newtraits should be a matrix with the number of rows equal to the number of \"new\" species")
          
               message("Overriding object$traits with newtraits...")
               object$p <- nrow(newtraits)
               }
          }
    if(is.null(object$jags.model)) 
        stop("MCMC samples not found")

        
     #################
     ## Checks done
     #################
    
     combined_fit_mcmc <- get.mcmcsamples(object) 
     mcmc_names <- colnames(combined_fit_mcmc)
     all_linpred <- array(NA, dim = c(n, object$p, nrow(combined_fit_mcmc)))
     pt_pred <- lower_linpred <- upper_linpred <- matrix(NA, nrow = n, ncol = object$p)

     if(predict.type == "conditional") 
          {
          for(k0 in 1:nrow(combined_fit_mcmc)) 
               {
               if(object$row.eff != "none" & noroweff == FALSE) 
                    cw_row_coefs <- vector("list", ncol(newrow.ids))

               cw_lv_coefs <- matrix(combined_fit_mcmc[k0, grep("lv.coefs", mcmc_names)], nrow = nrow(object$lv.coefs.mean)) ## Not using object$p since this might be changed by newtraits
               if(is.null(newtraits)) 
                    cw_eta <- matrix(cw_lv_coefs[,1,drop=FALSE], nrow = n, ncol = object$p, byrow = TRUE) 
            
               if(!is.null(X)) 
                    {
                    if(is.null(object$traits)) 
                         {
                         cw_X_coefs <- matrix(combined_fit_mcmc[k0, grep("X.coefs", mcmc_names)], nrow = object$p) 
                         cw_eta <- cw_eta + tcrossprod(X, cw_X_coefs)
                         }
                    if(!is.null(object$traits) & is.null(newtraits)) 
                         { ## Since the MCMC samples in X.coefs already account for residual noise given by trait.sigma
                         cw_X_coefs <- matrix(combined_fit_mcmc[k0, grep("X.coefs", mcmc_names)], nrow = object$p) 
                         cw_eta <- cw_eta + tcrossprod(X, cw_X_coefs)
                         }
                    if(!is.null(newtraits)) 
                         {
                         cw_traits_coefs <- cbind(combined_fit_mcmc[k0, grep("traits.int",colnames(combined_fit_mcmc))], matrix(combined_fit_mcmc[k0, grep("traits.coefs",colnames(combined_fit_mcmc))], nrow = ncol(X)+1))
                         rownames(cw_traits_coefs) <- c("beta0", colnames(object$X))
                         trait_X_coefs <- tcrossprod(cbind(1,newtraits), cw_traits_coefs) ## beta = intercept + trait %*% trait.coefs
                         cw_trait_sigma <- as.vector(combined_fit_mcmc[k0, grep("trait.sigma", mcmc_names)])
                         cw_X_coefs <- matrix(0, nrow = nrow(newtraits), ncol = ncol(X)+1)
                         for(k1 in 1:nrow(newtraits)) 
                         cw_X_coefs[k1,] <- rmvnorm(1, mean = trait_X_coefs[k1,], sigma = diag(cw_trait_sigma^2, nrow = ncol(X)+1)) 
                         
                         cw_eta <- tcrossprod(cbind(1,X), cw_X_coefs) ## Reset cw_eta here since the intercept is also new!                                         
                         }
                    }

               if(num.lv > 0) 
                    {
                    cw.lv <- matrix(combined_fit_mcmc[k0, grep("lvs", mcmc_names)], nrow = n)
                    if(is.null(newtraits))
                         cw_eta <- cw_eta + tcrossprod(cw.lv, cw_lv_coefs[,2:(num.lv+1),drop=FALSE])
                    if(!is.null(newtraits)) 
                         {
                         if(lvcoefs.type == "zero")
                              cw_eta <- cw_eta
                         if(lvcoefs.type == "median") 
                              {
                              new_lv_coefs <- matrix(apply(cw_lv_coefs, 2, median), nrow = nrow(newtraits), byrow = TRUE)
                              cw_eta <- cw_eta + tcrossprod(cw.lv, new_lv_coefs[,2:(num.lv+1),drop=FALSE])
                              }
                         if(lvcoefs.type == "mean") 
                              {
                              new_lv_coefs <- matrix(apply(cw_lv_coefs, 2, mean), nrow = nrow(newtraits), ncol = ncol(cw_lv_coefs), byrow = TRUE)
                              cw_eta <- cw_eta + tcrossprod(cw.lv, new_lv_coefs[,2:(num.lv+1),drop=FALSE])
                              }
                         }
                    }
                
               if(!is.null(object$offset)) 
                    cw_eta <- cw_eta + object$offset
               if(object$row.eff != "none" & noroweff == FALSE) 
                    { 
                    for(k in 1:ncol(newrow.ids)) 
                         {
                         cw_row_coefs[[k]] <- combined_fit_mcmc[k0, grep(paste0("row.coefs.ID",k,"\\["), mcmc_names)] 
                         cw_eta <- cw_eta + cw_row_coefs[[k]][newrow.ids[,k]] 
                         }
                    }

               all_linpred[,,k0] <- cw_eta
               }
          }


     if(predict.type == "marginal") 
          {
          mc_lv <- rmvnorm(n*lv.mc, mean = rep(0,num.lv))
          mc_lv <- array(c(mc_lv), dim = c(lv.mc, n, num.lv))
        
          for(k0 in 1:nrow(combined_fit_mcmc)) 
               {
        
               if(k0 %% 100 == 0) 
                    message("Onto MCMC sample ", k0)              

               cw_lv_coefs <- matrix(combined_fit_mcmc[k0, grep("lv.coefs", mcmc_names)], nrow = nrow(object$lv.coefs.mean)) ## Not using object$p since this might be changed by newtraits
               all_linpredmc <- array(NA, dim = c(n, object$p, lv.mc))

               if(!is.null(X)) 
                    {
                    if(is.null(newtraits))
                         cw_X_coefs <- matrix(combined_fit_mcmc[k0, grep("X.coefs", mcmc_names)], nrow = object$p) 
                    if(!is.null(newtraits)) 
                         {
                         cw_traits_coefs <- cbind(combined_fit_mcmc[k0, grep("traits.int",colnames(combined_fit_mcmc))], matrix(combined_fit_mcmc[k0, grep("traits.coefs",colnames(combined_fit_mcmc))], nrow = ncol(X)+1))
                         rownames(cw_traits_coefs) <- c("beta0", colnames(object$X))
                         trait_X_coefs <- tcrossprod(cbind(1,newtraits), cw_traits_coefs) ## beta = intercept + trait %*% trait.coefs
                         cw_trait_sigma <- as.vector(combined_fit_mcmc[k0, grep("trait.sigma", mcmc_names)])
                         cw_X_coefs <- matrix(0, nrow = nrow(newtraits), ncol = ncol(X)+1)
                         for(k1 in 1:nrow(newtraits)) 
                              cw_X_coefs[k1,] <- rmvnorm(1, mean = trait_X_coefs[k1,], sigma = diag(cw_trait_sigma^2, nrow = ncol(X)+1)) 
                         }
                    }
                
          if(object$row.eff == "fixed" & noroweff == FALSE) 
               { 
               cw_row_coefs <- vector("list", ncol(newrow.ids))
               for(k in 1:ncol(newrow.ids)) 
                    cw_row_coefs[[k]] <- combined_fit_mcmc[k0, grep(paste0("row.coefs.ID",k,"\\["), mcmc_names)] 
               } 
          if(object$row.eff == "random" & noroweff == FALSE) 
               {
               cw_row_coefs <- vector("list", ncol(newrow.ids))
               ## Need to generate from length(unique(object$row.ids[,k])) to account for fact that we may have less levels in newrow.ids (which is OK)
               for(k in 1:ncol(newrow.ids)) 
                    cw_row_coefs[[k]] <- matrix(rnorm(length(unique(object$row.ids[,k]))*lv.mc, mean = 0, sd = combined_fit_mcmc[k0, grep(paste0("row.sigma.ID",k,"$"), mcmc_names)]), ncol = lv.mc) 
               }
            
#            cw_lv_covparams <- combined_fit_mcmc[k0, grep("lv.covparams", mcmc_names)]
#             if(object$lv.control$type == "exponential")
#                 covmat_chol <- chol(exp(-distmat/cw_lv_covparams[1]))
#             if(object$lv.control$type == "squared.exponential")
#                 covmat_chol <- (chol(exp(-(distmat/cw_lv_covparams[1])^2)))
#            if(object$lv.control$type == "cauchy")
#                 covmat_chol <- (chol((1 + (distmat/cw_lv_covparams[1])^2)^(-cw_lv_covparams[2])))
#            if(object$lv.control$type == "spherical")
#                 covmat_chol <- (chol((distmat < cw_lv_covparams[1])*(1 - 1.5*distmat/cw_lv_covparams[1] + 0.5*(distmat/cw_lv_covparams[1])^3)))


          for(b in 1:lv.mc) 
               {
#                 if(object$lv.control$type != "independent")
#                     mc_lv[b,,] <- crossprod(covmat_chol, mc_lv[b,,])                    
                
               if(is.null(newtraits))
                    cw_eta <- matrix(cw_lv_coefs[,1,drop=FALSE], nrow = n, ncol = object$p, byrow = TRUE)                 
                
               if(!is.null(X)) 
                    {
                    if(is.null(newtraits))                    
                         cw_eta <- cw_eta + tcrossprod(X,cw_X_coefs)
                    if(!is.null(newtraits))
                         cw_eta <- tcrossprod(cbind(1,X), cw_X_coefs) ## Reset cw_eta here since the intercept is also new!                                                                 
                    }
                
               if(is.null(newtraits)) 
                    cw_eta <- cw_eta + tcrossprod(mc_lv[b,,], cw_lv_coefs[,2:(num.lv+1),drop=FALSE]) ## LVs to marginalize over
               if(!is.null(newtraits)) 
                    {
                    if(lvcoefs.type == "zero")
                         cw_eta <- cw_eta
                    if(lvcoefs.type == "median") 
                         {
                        new_lv_coefs <- matrix(apply(cw_lv_coefs, 2, median), nrow = nrow(newtraits), ncol = ncol(cw_lv_coefs), byrow = TRUE)
                        cw_eta <- cw_eta + tcrossprod(mc_lv[b,,], new_lv_coefs[,2:(num.lv+1),drop=FALSE])
                        }
                    if(lvcoefs.type == "mean") 
                         {
                         new_lv_coefs <- matrix(apply(cw_lv_coefs, 2, mean), nrow = nrow(newtraits), ncol = ncol(cw_lv_coefs), byrow = TRUE)
                         cw_eta <- cw_eta + tcrossprod(mc_lv[b,,], new_lv_coefs[,2:(num.lv+1),drop=FALSE])
                         }
                    }

               if(!is.null(object$offset)) 
                    cw_eta <- cw_eta + object$offset
               if(object$row.eff == "fixed" & noroweff == FALSE) 
                    { 
                    for(k in 1:ncol(newrow.ids)) 
                        cw_eta <- cw_eta + cw_row_coefs[[k]][newrow.ids[,k]] 
                    }
               if(object$row.eff == "random" & noroweff == FALSE) 
                    { 
                    for(k in 1:ncol(newrow.ids)) 
                        cw_eta <- cw_eta + cw_row_coefs[[k]][newrow.ids[,k],b] 
                    }
               all_linpredmc[,,b] <- cw_eta
               rm(cw_eta)
               }
                
          all_linpred[,,k0] <- apply(all_linpredmc, c(1,2), mean)
          rm(all_linpredmc)
          }
     }

        
     for(i in 1:n) { for(j in 1:object$p) 
          {
          if(est == "mean") 
               pt_pred[i,j] <- mean(all_linpred[i,j,]) ## Posterior mean
          if(est == "median") 
               pt_pred[i,j] <- median(all_linpred[i,j,]) ## Posterior median
          lower_linpred[i,j] <- quantile(all_linpred[i,j,], probs = (1-prob)/2)
          upper_linpred[i,j] <- quantile(all_linpred[i,j,], probs = 1-(1-prob)/2)
          } }
    
     out <- list(linpred = pt_pred, lower = lower_linpred, upper = upper_linpred)
     if(return.alllinpred) 
          out$all.linpred <- all_linpred

     return(out)
     }
    

## Simple extraction of MCMC samples from fitted boral object
get.mcmcsamples <- function(object) {
    fit.mcmc <- object$jags.model$BUGSoutput
    if(is.null(fit.mcmc)) 
    stop("MCMC samples not found. Please use save.model = TRUE to save MCMC samples when using boral")
    fit.mcmc <- mcmc(fit.mcmc$sims.matrix, start = 1, thin = object$mcmc.control$n.thin) ## Thanks to Guilliaume Blanchet for the original formatting!

    return(fit.mcmc)
    }
