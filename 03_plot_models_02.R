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
  mID <- which(grepl("m2", mnames))
  
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
    predictor <- c(substr(mname, nchar(mname)-2, nchar(mname)),
                   substr(mname, nchar(mname)-5, nchar(mname)-3))
    predictor <- paste0(substr(predictor, 1, 1), 
                        tolower(substr(predictor, 2, 3)))
    
    # Set x-axis title for plotting
    x_title <- paste0(toupper(predictor), " value")
    for(i in 1:length(x_title)){
      if(grepl("_", x_title[i])){
        x_title[i] <- gsub("_", "", x_title[i])
      }
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
    uss <- rep(list(NA), length(x_title))
    pDF <- rep(list(NA), length(x_title))
    temp_I <- rep(list(NA), length(x_title))
    
    for(i in 1:length(x_title)){
      if(predictor[i]=="Den"){
        uss[[i]] <- sort(unique(dm$I_De))
        uss[[i]] <- uss[[i]][which(uss[[i]] %in% densities$I)]
        
        # Data frame of known predictors
        pDF[[i]] <- densities
        pDF[[i]]$o <- pDF[[i]]$VD_adj
        
        temp_I[[i]] <- dm$I_De
        
      }else{
        uss[[i]] <- sort(unique(dm$I_Ut))
        uss[[i]] <- uss[[i]][which(uss[[i]] %in% uc$I)]
        
        # Data frame of known predictors
        pDF[[i]] <- uc
        pDF[[i]]$o <- pDF[[i]][,which(colnames(pDF[[i]])==
                                        toupper(predictor[i]))]
        
        temp_I[[i]] <- dm$I_Ut
      }
    }

    
    # Remove imputed data rows from covariates data frame
    d2 <- lapply(1:length(x_title), function(i){
      pDF[[i]][which(pDF[[i]]$o<99),]
    })
    
    # Apply over each predictor, making a plot of values
    # for each one.
    All_p_plots <- lapply(1:length(x_title), function(i){
      
      # Assign plotting variables from list of variables
      # created int he parent environment
      uss <- uss[[i]]
      d2 <- d2[[i]]
      temp_I <- temp_I[[i]]
      pDF <- pDF[[i]]
      predictor <- predictor[i]
      x_title <- x_title[i]
      
      # Assemble plots for the posterior value of each element
      # with accompanying data and site-wise value estimates.
      p_plots <- lapply(uss, function(x){
        
        # Get values of observed values for element x
        obs_d <- d2$o[which(d2$I==x)]
        # Get indices for x in archaeological dataset
        est_d_ind <- which(temp_I==x)
        # Get posterior distributions for element x
        # values for the archaeological dataset
        est_d_dist <- ps$PRED[, i, est_d_ind]
        if(is.matrix(est_d_dist)){
          est_d <- apply(est_d_dist, 2, median)
        }else{
          est_d <- mean(est_d_dist)
        }
        
        # Get 100 samples of beta distribution parameters
        # for element x
        shp <- as.vector(ps$shape[1:100, i])
        sgma <- ps$sigma[1:100, i]
        if(i==1){
          loc <- ps$loc1[1:100, x]
        }else if(i==2){
          loc <- ps$loc2[1:100, x]
        }
        
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
          s_element <- dm$Element[which(temp_I==x)][1]
          s_element <- gsub("E_", "", s_element)
          s_element <- gsub("_", " ", s_element)
          p_title <- paste0(s_element, " (", s_name, ")")
        }else{
          s_element <- dm$Element[which(temp_I==x)][1]
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
      ggsave(paste0(dirname, "/Vals_", predictor, "_", mname, ".pdf"), 
             device="pdf", plot=p_plotsA, width=14, height=height_scale,
             units="in")
      
      # Return plot panels for each predictors
      return(p_plotsA)
      
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
      # component for site type b.
      alpha <- sapply(TyCoInds, function(a){
        ps$t_alpha[ , y] + ps$c_alpha[ , y, a]
      })
      b1 <- sapply(TyCoInds, function(a){
        ps$t_b[ , 1, y] + ps$c_b1[ , y, a]
      })
      b2 <- sapply(TyCoInds, function(a){
        ps$t_b[ , 2, y] + ps$c_b2[ , y, a]
      })
      b12 <- sapply(TyCoInds, function(a){
        ps$t_b[ , 3, y] + ps$c_b12[ , y, a]
      })
      
      # Calculate mean alphas and bs for each component
      alpha <- colMeans(alpha)
      b1 <- colMeans(b1)
      b2 <- colMeans(b2)
      b12 <- colMeans(b12)
      
      # Get values over which to simulate
      raCol <- sapply(1:length(x_title), function(i){
        which(grepl(paste0(toupper(predictor[i]), "ran"),
                           colnames(pDF[[i]])))[1]
      })
      miCol <- sapply(1:length(x_title), function(i){
        which(grepl(paste0(toupper(predictor[i]), "min"),
                           colnames(pDF[[i]])))[1]
      })
      ra <- lapply(1:length(x_title), function(i){
        pDF[[i]][ , raCol[i]] + pDF[[i]][ , miCol[i]]
      })
      
      # Vector of values over which to simulate for
      # each predictor
      xsim <- sapply(1:length(x_title), function(i){
        seq(0, max(ra[[i]]), length.out=200)
      })
      
      # Simulate relationships for observed sites for
      # each predictor
      obs_sites <- lapply(1:length(x_title), function(i){
        
        xsimt <- xsim
        
        if(i==1){
          xsimt[,2] <- rep(0.5, nrow(xsim))
        }else{
          xsimt[,1] <- rep(0.5, nrow(xsim))
        }
        
        obs_sites <- lapply(1:length(b1), function(v){
          
          site <- paste0(rep(v, nrow(xsim)), "o")
          dtype <- rep("o", nrow(xsim))
          ysim <- exp(alpha[v] + xsimt[,1]*b1[v] +
                      xsimt[,2]*b2[v] + xsimt[,1]*xsimt[,2]*b12[v])
          
          simdata <- data.frame(site=site, dtype=dtype, 
                                x1=xsimt[,1], x2=xsimt[,2], y=ysim,
                                stringsAsFactors=FALSE)
          return(simdata)
        })
        
        obs_sites <- do.call("rbind", obs_sites)
        
        return(obs_sites)
      })
      
      # Simulate covariance matrices from the
      # posterior distributions.
      s <- 500 # number of values to simulate
      vcvm <- lapply(1:s, function(g){
        
        s1 <- ps$scale_c[g, y, 1]
        s2 <- ps$scale_c[g, y, 2]
        s3 <- ps$scale_c[g, y, 3]
        s4 <- ps$scale_c[g, y, 4]
        rho12 <- ps$Rho_c[g, y, 1, 2]
        rho13 <- ps$Rho_c[g, y, 1, 3]
        rho14 <- ps$Rho_c[g, y, 1, 4]
        rho23 <- ps$Rho_c[g, y, 2, 3]
        rho24 <- ps$Rho_c[g, y, 2, 4]
        rho34 <- ps$Rho_c[g, y, 3, 4]
        cov_s12 <- s1*s2*rho12
        cov_s13 <- s1*s3*rho13
        cov_s14 <- s1*s4*rho14
        cov_s23 <- s2*s3*rho23
        cov_s24 <- s2*s4*rho24
        cov_s34 <- s3*s4*rho34
        
        return(matrix(c(s1^2, cov_s12, cov_s13, cov_s14,
                        cov_s12, s2^2, cov_s23, cov_s24,
                        cov_s13, cov_s23, s3^2, cov_s34,
                        cov_s14, cov_s24, cov_s34, s4^2),
                      ncol=4))
      })
      
      # Simulate alphas and bs
      simpars <- t(sapply(1:s, function(a){
        mvrnorm(1, rep(0, 4), vcvm[[a]])
      }))
      
      # Calculate general parameters for
      # simulated sites
      alphaSim <- ps$t_alpha[1:s, y] + simpars[,1]
      b1Sim <- ps$t_b[1:s, 1, y] + simpars[,2]
      b2Sim <- ps$t_b[1:s, 2, y] + simpars[,3]
      b12Sim <- ps$t_b[1:s, 3, y] + simpars[,4]
      
      # Simulate relationships for simulated sites
      sim_sites <- lapply(1:length(x_title), function(i){
        
        xsimt <- xsim
        
        if(i==1){
          xsimt[,2] <- rep(0.5, nrow(xsim))
        }else{
          xsimt[,1] <- rep(0.5, nrow(xsim))
        }
        
        sim_sites <- lapply(1:length(b1Sim), function(v){
          
          site <- paste0(rep(v, nrow(xsim)), "s")
          dtype <- rep("s", nrow(xsim))
          ysim <- exp(alphaSim[v] + xsimt[,1]*b1Sim[v] +
                        xsimt[,2]*b2Sim[v] + xsimt[,1]*xsimt[,2]*b12Sim[v])
          
          simdata <- data.frame(site=site, dtype=dtype, 
                                x1=xsimt[,1], x2=xsimt[,2], y=ysim,
                                stringsAsFactors=FALSE)
          return(simdata)
        })
        
        sim_sites <- do.call("rbind", sim_sites)
        
        return(sim_sites)
      })
      
      # Bind data frame of simulated and observed sites
      all_sites <- lapply(1:length(x_title), function(i){
        rbind(sim_sites[[i]], obs_sites[[i]])
      })
      
      # Trim data frame points with y values over 100 to 
      # avoid plot distortion for extreme values
      ymax <- 100
      all_sites <- lapply(all_sites, function(i) i[which(i$y < ymax), ])
        
      
      # Set plotting colors specific to site type. These are
      # the same colors as the site specific plots.
      if(y==1){
        c_line <- "red"
      }else if(y==2){
        c_line <- "darkorchid1"
      }else{
        c_line <- "turquoise3"
      }
      
      SitePlot <- rep(list(NA), 2)
      
      for(i in 1:length(x_title)){

        if(i==1){
          SitePlot[[i]] <- ggplot(all_sites[[i]], aes(y=y, x=x1))
        }else if(i==2){
          SitePlot[[i]] <- ggplot(all_sites[[i]], aes(y=y, x=x2))
        }
        SitePlot[[i]] <- SitePlot[[i]] +
          geom_line(aes(group=site, color=dtype,
                        alpha=dtype), size=0.35)+
          scale_alpha_manual(values=c("s"=0.1, 
                                      "o"=1))+
          scale_color_manual(values=c("s"="grey60", 
                                      "o"=c_line))+
          labs(x=x_title[i], y="MAU") + ggtitle(tname)+
          scale_y_continuous(limits=c(0, 100))
        
        if(predictor[i]=="Den"){
          SitePlot[[i]] <- SitePlot[[i]] +
            scale_x_continuous(limits=c(0, max(ra[[i]])),
                               breaks=seq(0, 3, 0.25),
                               expand=c(0,0))
          
        }else{
          SitePlot[[i]] <- SitePlot[[i]] +
            scale_x_continuous(limits=c(0, max(ra[[i]])),
                               breaks=seq(0, 3, 0.25),
                               labels=seq(0, 300, 25),
                               expand=c(0,0))
        }
        
        if(y!=max(dm$I_Ty)){
          SitePlot[[i]] <- SitePlot[[i]] +
            theme(panel.background=element_blank(), 
                  panel.grid=element_blank(),
                  axis.title.x=element_blank(),
                  axis.line=element_line(color="grey"),
                  axis.ticks=element_line(color="grey"),
                  plot.title=element_text(size=11),
                  legend.position="none")
        }else{
          SitePlot[[i]] <- SitePlot[[i]] +
            theme(panel.background=element_blank(), 
                  panel.grid=element_blank(),
                  axis.line=element_line(color="grey"),
                  axis.ticks=element_line(color="grey"),
                  plot.title=element_text(size=11),
                  legend.position="none")
        }      
        
      }

      SitePlot <- SitePlot[[1]] + SitePlot[[2]]
      
      return(SitePlot)
      
    })
    
    # Combine plots
    p_S <- SimPlots[[1]]
    for(i in 2:length(SimPlots)){
      p_S <- p_S / SimPlots[[i]]
    }
    
    # Save plots
    height_scale <- 6*(length(SimPlots)/4)
    ggsave(paste0(dirname, "/Sites_", mname, ".pdf"), 
           device="pdf", plot=p_S, width=6, height=height_scale,
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


