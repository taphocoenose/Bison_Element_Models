# Plots for output generated from the first two scripts

# REMOVED GLOBAL ALPHA AND BETA

# Load libraries and data
library(ggplot2)
library(patchwork)
library(parallel)
library(MASS)

# Load data
load("01_Out.RData")

# Get output file names for LOO outputs and bind
# all LOO outputs into one table
dirfiles <- list.files()
dirfiles <- dirfiles[which(grepl("Out_LOO", dirfiles))]
LOO_table <- lapply(dirfiles, function(x){
  load(x)
  return(LOO_comp)
})
LOO_table <- do.call("rbind", LOO_table)

# Print stat summaries across predictors for each mode
for(i in unique(LOO_table$model)){
  cat(paste0(i, ": med min_neff ", 
             round(median(LOO_table$min_neff[which(LOO_table$model==i)]), 3),
             ", med max_rhat ",
             round(median(LOO_table$max_rhat[which(LOO_table$model==i)]), 3),
             "\n"))
}


# Define inverse logit function
inv_logit <- function(x){
  exp(x)/(1+exp(x))
}

# Get output file names and find which ones contain 
# model output data
dirfiles <- list.files()
dirfiles <- dirfiles[which(grepl("Out_Models", dirfiles))]
dirfiles <- dirfiles[which(grepl("m1", dirfiles))]


# Start cluster and export packages and variables
corecount <- detectCores()-1
if(corecount > 2){
  corecount <- 2
}
cl <- makeCluster(corecount)
clusterExport(cl, varlist=c("densities", "dm", "uc",
                            "inv_logit"))
clusterEvalQ(cl, library("ggplot2"))
clusterEvalQ(cl, library("patchwork"))
clusterEvalQ(cl, library("MASS"))


# Sequentially load files for plotting
for(q in 1:length(dirfiles)){
  
  # Display progress text
  cat(paste("\nPlotting", dirfiles[q], Sys.time()))
  starttime <- Sys.time()
  
  # Load models saved in file q
  load(dirfiles[q])
  Model_List <- output_models
  rm(output_models)
  
  #######################################################################
  #################### PLOTS OF COVARIATE VALUES #######################
  #####################################################################
  
  # Find which model names start with 'm1', indicating that only
  # one covariate is in the model
  mnames <- sapply(Model_List, function(x) x[[1]])
  mID <- which(grepl("m1", mnames))
  
  # Subset list of model for relevant variables
  m1_List <- lapply(Model_List[mID], function(x){
    list(model_name=x$model_name,
         post_samples=x$post_samples,
         loo=x$loo)
  })
  rm(Model_List)
  
  # CREATE A SERIES OF COVARIATE VALUE PLOTS FOR EVERY FITTED MODEL
  Successful_Value_Plots <- parSapply(cl, m1_List, function(z){
    
    # Get model name
    mname <- z$model_name
    predictor <- substr(mname, nchar(mname)-2, nchar(mname))
    predictor <- paste0(substr(predictor, 1, 1), 
                        tolower(substr(predictor, 2, 3)))
    
    # Set x-axis title for plotting
    x_title <- paste0(toupper(predictor), " value")
    if(grepl("_", x_title)){
      x_title <- gsub("_", "", x_title)
    }
    
    # Create directory name
    dirname <- paste0("Plots/", mname, "_plots")
    
    # If subdirectory for plots does not exist, create it
    if(!dir.exists("Plots")) dir.create("Plots")
    if(!dir.exists(dirname)) dir.create(dirname)
    
    # Extract posterior samples from the first model
    ps <- z$post_samples
    
    # Get unique indices for element utilities or densities and
    # select out the elements that are combinations of scan site
    # or density values
    if(predictor=="Den"){
      uss <- sort(unique(dm$I_De))
      uss <- uss[which(uss %in% densities$I)]
      
      # Data frame of known predictors
      pDF <- densities
      pDF$o <- pDF$VD_adj
      
      dm$temp_I <- dm$I_De
      
    }else{
      uss <- sort(unique(dm$I_Ut))
      uss <- uss[which(uss %in% uc$I)]
      
      # Data frame of known predictors
      pDF <- uc
      pDF$o <- pDF[,which(colnames(pDF)==toupper(predictor))]
      
      dm$temp_I <- dm$I_Ut
    }
    
    # Remove imputed data rows from covariates data frame
    d2 <- pDF[which(pDF$o<99),]
    
    # Assemble plots for the posterior value of each element
    # with accompanying data and site-wise value estimates.
    p_plots <- lapply(uss, function(x){
      
      # Get values of observed values for element x
      obs_d <- d2$o[which(d2$I==x)]
      # Get indices for x in archaeological dataset
      est_d_ind <- which(dm$temp_I==x)
      # Get posterior distributions for element x
      # values for the archaeological dataset
      est_d_dist <- ps$PRED1[, est_d_ind]
      if(is.matrix(est_d_dist)){
        est_d <- apply(est_d_dist, 2, median)
      }else{
        est_d <- median(est_d_dist)
      }
      
      # Get 100 samples of beta distribution
      # parameters for element x
      shp <- as.vector(ps$shape[1:100])
      sgma <- ps$sigma[1:100]
      loc <- ps$loc[1:100, x]
      if(grepl("exp", mname)){
        shp <- exp(shp) + 2
        sgma <- exp((sgma-4)*2)
      }else{
        shp <- shp + 2
      }
      if(grepl("nc", mname)){
        loc <- ps$mu[1:100] + loc*sgma
      }
      loc <- inv_logit(loc)
      
      # Get x-axis scaling parameters for scaled beta distribution
      raCol <- which(grepl(paste0(toupper(predictor), "ran"),
                           colnames(pDF)))[1]
      miCol <- which(grepl(paste0(toupper(predictor), "min"),
                           colnames(pDF)))[1]
      ra <- pDF[which(pDF$I==x), raCol][1]
      mi <- pDF[which(pDF$I==x), miCol][1]
      
      # Sequence of x values over which to get the density
      # of the scaled beta distribution
      betaseq <- seq(0, 1, 0.01)
      
      # Create data frame to plot beta distributions
      bDF <- lapply(1:length(shp), function(u){
        grp <- rep(as.character(u), length(betaseq))
        X <- betaseq*ra + mi
        Y <- dbeta(betaseq, loc[u]*shp[u], (1-loc[u])*shp[u])
        return(data.frame(grp=grp, X=X, Y=Y,
                          stringsAsFactors=FALSE))
      })
      bDF <- do.call("rbind", bDF)
      bDF <- bDF[!is.infinite(bDF$Y), ]
      
      # Get y axis heights
      ymax <- max(bDF$Y)*1.05
      ymaxT <- ymax*1.1
      ymax0 <- -0.03*ymax
      ymax1 <- -0.10*ymax
      ymax2 <- -0.17*ymax
      ymax3 <- -0.20*ymax
      
      # Get x element name and scan site name
      if(predictor=="Den"){
        s_name <- d2$ScanSite[which(d2$I==x)][1]
        s_element <- dm$Element[which(dm$temp_I==x)][1]
        s_element <- gsub("E_", "", s_element)
        s_element <- gsub("_", " ", s_element)
        p_title <- paste0(s_element, " (", s_name, ")")
      }else{
        s_element <- dm$Element[which(dm$temp_I==x)][1]
        s_element <- gsub("E_", "", s_element)
        p_title <- gsub("_", " ", s_element) 
      }
      
      # Abbreviate longer common words in plot_title
      if(grepl("Undifferentiated", p_title)){
        p_title <- gsub("Undifferentiated", "Undiff.", p_title)
      }
      if(grepl("Unspecified", p_title)){
        p_title <- gsub("Unspecified", "Unspec.", p_title)
      }
      if(grepl(".3.7.", p_title)){
        p_title <- gsub(".3.7.", "3-7", p_title)
      }
      if(grepl("phalanges", p_title)){
        p_title <- gsub("phalanges", "phal.", p_title)
      }
      
      # X positions for text plotting
      xp <- quantile(bDF$X, probs=c(0.01, 0.80, 0.99))
      
      # Initialize plot
      p <- ggplot(bDF, aes(x=X, y=Y, group=grp))+
        annotate("text", label=p_title, x=xp[1], y=ymax, vjust=0, hjust=0)+
        annotate("text", label=paste0("italic(n)==", length(est_d)),
                 x=xp[2], y=ymax, vjust=0, hjust=1, parse=TRUE, color="red")+
        annotate("text", label=paste0("italic(n)==", length(obs_d)),
                 x=xp[3], y=ymax, vjust=0, hjust=1, parse=TRUE, color="purple")+
        annotate("segment", x=est_d, xend=est_d, y=rep(ymax0, length(est_d)),
                 yend=rep(ymax1, length(est_d)), color="red", alpha=0.6)+
        annotate("segment", x=obs_d, xend=obs_d, y=rep(ymax1, length(obs_d)),
                 yend=rep(ymax2, length(obs_d)), color="purple", alpha=0.6)+
        geom_line(color="blue", alpha=0.2)
      # If predictor is utility values, rescale x-axis text to 1-100.
      if(predictor=="Den"){
        p <- p +
          scale_x_continuous(limits=c(min(bDF$X), max(bDF$X)), 
                             breaks=c(min(bDF$X), max(bDF$X)),
                             labels=round(c(min(bDF$X), max(bDF$X)), 2), 
                             expand=c(0, 0))
      }else{
        p <- p +
        scale_x_continuous(limits=c(min(bDF$X), max(bDF$X)), 
                           breaks=c(min(bDF$X), max(bDF$X)),
                           labels=round(c(min(bDF$X)*100, 
                                          max(bDF$X)*100), 2), 
                           expand=c(0, 0))
      }
      # If in the last five plots, add x-axis title
      if(x %in% uss[(length(uss)-4):length(uss)]){
        p <- p + labs(x=x_title) +
          scale_y_continuous(expand=c(0.02,0.02), limits=c(ymax3, ymaxT))+
          theme(panel.background=element_blank(), panel.grid=element_blank(),
                axis.ticks.y=element_blank(), axis.text.y=element_blank(),
                axis.title.y=element_blank(), axis.line.x=element_line(color="black"),
                axis.text.x=element_text(hjust=c(0, 1)))
        
        # Otherwise, exclude x-axis title
      }else{
        p <- p +
          scale_y_continuous(expand=c(0.02,0.02), limits=c(ymax3, ymaxT))+
          theme(panel.background=element_blank(), panel.grid=element_blank(),
                axis.ticks.y=element_blank(), axis.text.y=element_blank(),
                axis.title=element_blank(), axis.line.x=element_line(color="black"),
                axis.text.x=element_text(hjust=c(0, 1)))
      }
      
      # Return plot
      return(p)
      
    })
    
    # Combine plots
    p_plotsA <- p_plots[[1]]
    for(i in 2:length(p_plots)){
      p_plotsA <- p_plotsA + p_plots[[i]]
    }
    p_plotsA <- p_plotsA + plot_layout(ncol=5)
    
    # Save plots
    height_scale <- 7*(ceiling(length(p_plots)/5)/4)
    ggsave(paste0(dirname, "/Vals_", mname, ".pdf"), 
           device="pdf", plot=p_plotsA, width=14, height=height_scale,
           units="in")
    
    
    #######################################################################
    #################### PLOTS OF LINEAR MODELS ##########################
    #####################################################################
    
    # Calculate median and 95% percentile intervals for
    # the posterior values of each element.
    dm$dv_med <- sapply(1:nrow(dm), function(x){
      median(ps$PRED1[,x])
    })
    dm$l95 <- sapply(1:nrow(dm), function(x){
      quantile(ps$PRED1[,x], p=0.025)
    })
    dm$u95 <- sapply(1:nrow(dm), function(x){
      quantile(ps$PRED1[,x], p=0.975)
    })
    
    # Calculate MAU
    dm$MAU <- dm$MNE/dm$AF
    # Assign an ID to each component that assigns which
    # plot it will be assigned to.
    dm$PlotID <- sapply(dm$I_Co, function(x) ceiling(x/12))
    
    # Max and min x value for plotting
    maxx2 <- max(dm$u95)
    minx2 <- 0
    
    # Indices for drawn posterior samples
    si <- 1:5000
    # Vector of covariate values over which
    # to simulate posterior predictions
    cvval <- seq(minx2, maxx2, length.out=200)
    # Data frame of percentile intervals
    # for posterior predictions
    spi <- data.frame(perc=seq(0.98, 0.02, -0.02))
    spi$IntervalID <- 1:nrow(spi)
    spi$Lower <- (1-spi$perc)/2
    spi$Upper <- spi$perc + spi$Lower
    
    Model_Plots <- lapply(1:max(dm$PlotID), function(x){
      
      dm2 <- dm[which(dm$PlotID==x), ]
      
      inds <- min(dm2$I_Co):max(dm2$I_Co)
      
      Plotx <- lapply(inds, function(y){
        
        dm3 <- dm2[which(dm2$I_Co==y), ]
        
        p_title <- paste0(gsub("S_", "", dm3$Comp[1]),
                          " (MNE=", sum(dm3$MNE), ")")
        
        # Shorten title for plotting
        if(grepl("ocality", p_title)){
          p_title <- gsub("ocality", "oc.", p_title)
        }
        if(grepl("evel", p_title)){
          p_title <- gsub("evel", "vl", p_title)
        }
        if(grepl("ample", p_title)){
          p_title <- gsub("ample", "mpl", p_title)
        }
        if(grepl("onebed", p_title)){
          p_title <- gsub("onebed", "nbd", p_title)
        }
        if(grepl("omponent", p_title)){
          p_title <- gsub("omponent", "mpnt", p_title)
        }
        
        # Calculate predictive distributions
        # for model y
        if(grepl("v", mname)){
          alpha <- ps$t_alpha[si, dm3$I_Ty[1]] + 
            ps$c_alpha[si, dm3$I_Ty[1], dm3$I_Ty_Co[1]]
          
          b <- ps$t_b[si, dm3$I_Ty[1]] + 
            ps$c_b[si, dm3$I_Ty[1], dm3$I_Ty_Co[1]]
          
        }else{
          alpha <- ps$t_alpha[si, dm3$I_Ty[1]] + 
            ps$c_alpha[si, y]*ps$scale_c[si, 1, dm3$I_Ty[1]]
          
          b <- ps$t_b[si, dm3$I_Ty[1]] + 
            ps$c_b[si, y]*ps$scale_c[si, 2, dm3$I_Ty[1]]
        }
        
        # Simulate posterior lambda values
        # for each value in cvval
        s_lambda <- t(sapply(cvval, function(w){
          pd <- exp(alpha + w*b)
          return(c(median(pd), 
                   quantile(pd, prob=0.025),
                   quantile(pd, prob=0.975)))
        }))
        s_lambdaM <- data.frame(X=cvval, med=s_lambda[,1])
        s_lambdaL <- data.frame(X=cvval, l95=s_lambda[,2])
        s_lambdaU <- data.frame(X=cvval, u95=s_lambda[,3])
        
        
        # Simulate posterior prediction offsets for
        # each value in cvval. This will vary depending on
        # if offsets are drawn from a normal or t distribution
        # and whether the dispersion of this distribution
        # varies by component.
        if(grepl("v2", mname)){
          I_s <- ps$I_scale[si, dm3$I_Co[1]]
        }else{
          I_s <- ps$I_scale[si]
        }

        post_scalings <- sapply(I_s, function(q) rnorm(1, 0, q))
        
        
        s_pred <- lapply(cvval, function(w){
          
          pd <- exp(alpha + post_scalings + w*b)
          
          preds <- suppressWarnings(sapply(pd, function(u) rpois(1, u)))
          predsmax <- max(preds[!is.na(preds)])
          preds[is.na(preds)] <- predsmax
          
          percentiles <- data.frame(X=rep(w, nrow(spi)),
                                    PercentileID=spi$IntervalID,
                                    Lower=sapply(spi$Lower, 
                                                 function(b) quantile(preds, prob=b)),
                                    Upper=sapply(spi$Upper, 
                                                 function(b) quantile(preds, prob=b)),
                                    stringsAsFactors=FALSE)
          
          return(percentiles)
        })
        
        s_pred <- do.call("rbind", s_pred)
        if(dm3$I_Ty[1]==1){
          c_fill <- "red"
        }else if(dm3$I_Ty[1]==2){
          c_fill <- "darkorchid1"
        }else{
          c_fill <- "turquoise3"
        }
        
        # Yaxis scale
        miny <- min(dm3$MAU)
        minlab <- floor(miny)
        maxy <- max(dm3$MAU)
        maxlab <- ceiling(maxy)
        
        mody <- (maxy - miny)*0.1
        maxy <- maxy + mody
        miny <- miny - mody
        
        if(maxlab > maxy) maxy <- maxlab
        if(minlab < miny) miny <- minlab
        
        
        # Reformat percentile interval geometry for 
        # y-axis limits
        s_pred$Lower[which(s_pred$Lower < miny)] <- miny
        s_pred$Upper[which(s_pred$Upper > maxy)] <- maxy
        
        # If simulated lambda values are outside of
        # y-axis range, trim values.
        s_lambdaM <- s_lambdaM[which((s_lambdaM$med >= miny) &
                                       (s_lambdaM$med <= maxy)), ]
        s_lambdaL <- s_lambdaL[which((s_lambdaL$l95 >= miny) &
                                       (s_lambdaL$l95 <= maxy)), ]
        s_lambdaU <- s_lambdaU[which((s_lambdaU$u95 >= miny) &
                                       (s_lambdaU$u95 <= maxy)), ]
        
        p <- ggplot(data=s_pred, aes(x=X))+
          geom_ribbon(aes(group=factor(PercentileID),
                          ymin=Lower, ymax=Upper),
                      fill=c_fill, alpha=1/nrow(spi))+
          annotate("line", x=s_lambdaM$X, y=s_lambdaM$med,
                   color=c_fill, size=1.5)+
          annotate("line", x=s_lambdaL$X, y=s_lambdaL$l95,
                   color=c_fill, size=1)+
          annotate("line", x=s_lambdaU$X, y=s_lambdaU$u95,
                   color=c_fill, size=1)+
          annotate("line", x=s_lambdaM$X, y=s_lambdaM$med,
                   color="white", size=1)+
          annotate("line", x=s_lambdaL$X, y=s_lambdaL$l95,
                   color="white", size=0.5)+
          annotate("line", x=s_lambdaU$X, y=s_lambdaU$u95,
                   color="white", size=0.5)+
          geom_segment(data=dm3, aes(x=l95, xend=u95, 
                                     y=MAU, yend=MAU),size=0.5)+
          geom_point(data=dm3, aes(x=dv_med, y=MAU), size=1)+
          labs(x=x_title, y="MAU")+
          ggtitle(p_title)+
          scale_y_continuous(breaks=c(minlab, maxlab),
                             labels=c(minlab, maxlab),
                             limits=c(miny, maxy),
                             expand=c(0.02,0.02))
          
          if(predictor=="Den"){
            p <- p +
              scale_x_continuous(limits=c(0, maxx2), 
                                 breaks=c(0, maxx2),
                                 labels=c(0, round(maxx2, 2)),
                                 expand=c(0.02,0.02))

          }else{
            p <- p +
            scale_x_continuous(limits=c(0, maxx2), 
                               breaks=c(0, maxx2),
                               labels=c(0, round(maxx2*100, 2)),
                               expand=c(0.02,0.02))
          }
        
        if((y %in% seq(1, 5e3, 4)) & 
           (y %in% inds[(length(inds)-3):(length(inds))])){
          p <- p +
            theme(panel.background=element_blank(), 
                  panel.grid=element_blank(),
                  axis.line=element_line(color="grey"),
                  axis.ticks=element_line(color="grey"),
                  plot.title=element_text(size=11, hjust=0.5))   
        }else if(y %in% inds[(length(inds)-3):(length(inds))]){
          p <- p +
            theme(panel.background=element_blank(), 
                  panel.grid=element_blank(),
                  axis.title.y=element_blank(),
                  axis.line=element_line(color="grey"),
                  axis.ticks=element_line(color="grey"),
                  plot.title=element_text(size=11, hjust=0.5))   
        }else if(y %in% seq(1, 5e3, 4)){
          p <- p +
            theme(panel.background=element_blank(), 
                  panel.grid=element_blank(),
                  axis.line=element_line(color="grey"),
                  axis.title.x=element_blank(),
                  axis.ticks=element_line(color="grey"),
                  plot.title=element_text(size=11, hjust=0.5))   
        }else{
          p <- p +
            theme(panel.background=element_blank(), 
                  panel.grid=element_blank(),
                  axis.title=element_blank(),
                  axis.line=element_line(color="grey"),
                  axis.ticks=element_line(color="grey"),
                  plot.title=element_text(size=11, hjust=0.5))   
        }
        
        
        return(p)
        
      })
      
      # Combine plots
      p_A <- Plotx[[1]]
      for(i in 2:length(Plotx)){
        p_A <- p_A + Plotx[[i]]
      }
      p_A <- p_A + plot_layout(ncol=4)
      
      # Save plots
      height_scale <- 7*(ceiling(length(Plotx)/4)/3)
      ggsave(paste0(dirname, "/Model_", x, "_", mname, ".pdf"), 
             device="pdf", plot=p_A, width=14, height=height_scale,
             units="in")
      
      return(p_A)
      
    })
    
    #######################################################################
    #################### PLOTS OF SIMULATED SITES ########################
    #####################################################################
    
    SimPlots <- lapply(1:max(dm$I_Ty), function(y){
      
      # Get type name and component indices
      tname <- dm$Type[which(dm$I_Ty==y)][1]
      CoInds <- unique(dm$I_Co[which(dm$I_Ty==y)])
      TyCoInds <- unique(dm$I_Ty_Co[which(dm$I_Ty==y)])
      
      # Calculate alpha and b parameters for every 
      # component for site type b, with calculations
      # varying by model type.
      if(grepl("v", mname)){
        alpha <- sapply(TyCoInds, function(a){
          ps$t_alpha[ , y] + ps$c_alpha[ , y, a]
        })

        
        b <- sapply(TyCoInds, function(a){
          ps$t_b[ , y] + ps$c_b[ , y, a]
        })
        
      }else{
        alpha <- sapply(CoInds, function(a){
          ps$t_alpha[ , y] + ps$c_alpha[ , a]*ps$scale_c[ , 1, y]
        })
        
        b <- sapply(CoInds, function(a){
          ps$t_b[ , y] + ps$c_b[ , a]*ps$scale_c[ , 2, y]
        })
      }
      
      # Calculate mean alphas and bs for each component
      alpha <- colMeans(alpha)
      b <- colMeans(b)
      
      # Get values over which to simulate
      raCol <- which(grepl(paste0(toupper(predictor), "ran"),
                           colnames(pDF)))[1]
      miCol <- which(grepl(paste0(toupper(predictor), "min"),
                           colnames(pDF)))[1]
      ra <- pDF[ , raCol] + pDF[ , miCol]
      
      # Vector of values over which to simulate
      xsim <- seq(0, max(ra), length.out=200)
      
      # Simulate relationships for observed sites
      obs_sites <- lapply(1:length(b), function(v){
        
        site <- paste0(rep(v, length(xsim)), "o")
        dtype <- rep("o", length(xsim))
        ysim <- exp(alpha[v] + xsim*b[v])
        
        simdata <- data.frame(site=site, dtype=dtype,
                              x=xsim, y=ysim,
                              stringsAsFactors=FALSE)
        return(simdata)
      })
      
      obs_sites <- do.call("rbind", obs_sites)
      
      
      # Simulate covariance matrices from the
      # posterior distributions, depending on
      # the model type.
      s <- 500 # number of values to simulate
      if(grepl("v", mname)){
        vcvm <- lapply(1:s, function(g){
          
          s1 <- ps$scale_c[g, y, 1]
          s2 <- ps$scale_c[g, y, 2]
          rho <- ps$Rho_c[g, y, 1, 2]
          cov_s <- s1*s2*rho
          
          return(matrix(c(s1^2, cov_s, 
                        cov_s, s2^2), ncol=2))
        })
      }else{
        vcvm <- lapply(1:s, function(g){
          
          s1 <- ps$scale_c[g, 1, y]
          s2 <- ps$scale_c[g, 2, y]
          rho <- ps$Rho_c[g, 1, 2]
          cov_s <- s1*s2*rho
          
          return(matrix(c(s1^2, cov_s, 
                          cov_s, s2^2), ncol=2))
        })
      }
      
      # Simulate alphas and bs
      simpars <- t(sapply(1:s, function(a){
        mvrnorm(1, c(0, 0), vcvm[[a]])
      }))
      
      # Calculate general parameters for
      # simulated sites
      alphaSim <- ps$t_alpha[1:s, y] + simpars[,1]
      bSim <- ps$t_b[1:s, y] + simpars[,2]
      
      # Simulate relationships for simulated sites
      sim_sites <- lapply(1:length(bSim), function(v){
        
        site <- paste0(rep(v, length(xsim)), "s")
        dtype <- rep("s", length(xsim))
        ysim <- exp(alphaSim[v] + xsim*bSim[v])
        
        simdata <- data.frame(site=site, dtype=dtype,
                              x=xsim, y=ysim,
                              stringsAsFactors=FALSE)
        return(simdata)
      })
      
      sim_sites <- do.call("rbind", sim_sites)
      
      # Bind data frame of simulated and observed sites
      all_sites <- rbind(sim_sites, obs_sites)
      
      # Trim data frame points with y values over 150 to 
      # avoid plot distortion for extreme values
      ymax <- 100
      all_sites <- all_sites[which(all_sites$y < ymax), ]
      
      # Set plotting colors specific to site type. These are
      # the same colors as the site specific plots.
      if(y==1){
        c_line <- "red"
      }else if(y==2){
        c_line <- "darkorchid1"
      }else{
        c_line <- "turquoise3"
      }
      
      SitePlot <- ggplot(all_sites, aes(x=x, y=y))+
        geom_line(aes(group=site, color=dtype,
                      alpha=dtype), size=0.35)+
        scale_alpha_manual(values=c("s"=0.1, 
                                    "o"=1))+
        scale_color_manual(values=c("s"="grey60", 
                                    "o"=c_line))+
        labs(x=x_title, y="MAU") + ggtitle(tname)+
        scale_y_continuous(limits=c(0, ymax))
      
      if(predictor=="Den"){
        SitePlot <- SitePlot +
          scale_x_continuous(limits=c(0, max(ra)),
                             breaks=seq(0, 3, 0.2),
                             expand=c(0,0))
        
      }else{
        SitePlot <- SitePlot +
          scale_x_continuous(limits=c(0, max(ra)),
                             breaks=seq(0, 3, 0.2),
                             labels=seq(0, 300, 20),
                             expand=c(0,0))
      }

      if(y!=max(dm$I_Ty)){
        SitePlot <- SitePlot +
          theme(panel.background=element_blank(), 
                panel.grid=element_blank(),
                axis.title.x=element_blank(),
                axis.line=element_line(color="grey"),
                axis.ticks=element_line(color="grey"),
                plot.title=element_text(size=11),
                legend.position="none")
      }else{
        SitePlot <- SitePlot +
          theme(panel.background=element_blank(), 
                panel.grid=element_blank(),
                axis.line=element_line(color="grey"),
                axis.ticks=element_line(color="grey"),
                plot.title=element_text(size=11),
                legend.position="none")
      }
      
      return(SitePlot)
      
    })
    
    # Combine plots
    p_S <- SimPlots[[1]]
    for(i in 2:length(SimPlots)){
      p_S <- p_S + SimPlots[[i]]
    }
    p_S <- p_S + plot_layout(ncol=1)
    
    # Save plots
    height_scale <- 6*(length(SimPlots)/4)
    ggsave(paste0(dirname, "/Sites_", mname, ".pdf"), 
           device="pdf", plot=p_S, width=3, height=height_scale,
           units="in")
    
    
    return(mname)
    
  })
  
  gc()
  
  # Display progress statement
  cat(paste0("\n", dirfiles[q], " plotting complete."))
  cat(paste("\nTime:", round(Sys.time() - starttime, 2), "min.\n"))
}

# End cluster
stopCluster(cl)


