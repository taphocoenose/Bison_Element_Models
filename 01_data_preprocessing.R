# IMPORT ARCHAEOLOGICAL BISON ELEMENT DATASETS AND FIT STAN
# MODELS TO ESTIMATE THE EFFECTS OF VARIOUS DIETARY UTILITY,
# TAPHONOMIC, AND SITE TYPE COVARIATES ON ELEMENT FREQUENCIES.

# Load libraries
library(reshape2)

############################################################
########## Import and format site-level data ##############
##########################################################

d <- read.csv("Data/MNE_Data.csv", header=TRUE, stringsAsFactors=FALSE)
AF <- d[which(d$Component_or_Index=="AF"), ]
AF <- AF[1,!is.na(AF[1,])][,-which(colnames(AF)=="Component_or_Index")]
d <- d[which(d[,1]!="AF"),]

# Aggregate phalanges into E_Unspecified_phalanges
d$E_Unspecified_phalanges <- sapply(1:nrow(d), function(x){
  
  if(is.na(d$E_Unspecified_phalanges[x])){
    
    # If E_unspecified phalanges is missing an MNE value, create a
    # vector of MNE values for the individual phalanges
    p <- c(d$E_Phalanx_1[x], d$E_Phalanx_2[x], d$E_Phalanx_3[x])
    
    if(length(p[!is.na(p)])==0){
      
      # If thre are no MNE values for specific phalanges, return NA
      # for E_Unspecified_pahalnges
      return(NA)
    }else{
      
      # If there are valid MNE values for specific phalanges,
      # return maximum MNE value across these values
      return(sum(p[!is.na(p)]))
    }
  }else{
    
    # If E_unspecified_phalanges already contains an MNE value,
    # return this value.
    return(d$E_Unspecified_phalanges[x])
  }
  
})
# Reset specific phalanx MNE values to NA
d$E_Phalanx_1 <- d$E_Phalanx_2 <- d$E_Phalanx_3 <- NA

dm <- melt(d, id.vars=colnames(d)[which(substr(colnames(d), 1, 2)!="E_")])
colnames(dm)[which(colnames(dm)=="variable")] <- "Element"
dm$Element <- as.vector(dm$Element)
colnames(dm)[which(colnames(dm)=="Component_or_Index")] <- "Comp"
colnames(dm)[which(colnames(dm)=="value")] <- "MNE"
dm <- dm[,c(1:4, ncol(dm)-1, ncol(dm), seq(5, ncol(dm)-2, 1))]
dm <- dm[!is.na(dm$MNE),]

# Remove mandibles and crania from dataset and remove extra columns
dm <- dm[!(dm$Element %in% c("E_Cranium", "E_Mandible", "E_Caudal_vert")),1:12]

# Store anatomical frequency values in dm
dm$AF <- sapply(dm$Element, function(x) AF[1,which(colnames(AF)==x)])
dm$logAF <- log(dm$AF)
rm(AF)

############################################################
## These blocks of code remove any sites from the dataset
## with specific issues, such as anomalous period or regions
## with only one observation, or random data quality issues.
## Here, I have selected only early/late Paleoindian sites
## from the Great Plains with reliable datasets.Any site types
## other than "kill", "camp", and "unknown" are excluded.
dm <- dm[which(dm$Type %in% c("kill", "camp", "unknown")),]
#dm <- dm[which(dm$Region %in% c("southern plains", "northern plains")),]
dm <- dm[which(dm$Period %in% c("early paleoindian", "late paleoindian")),]
#dm <- dm[which(dm$Data_issues == "no"),]

# Create index values
dm$I_Co <- sapply(dm$Comp, function(x) which(unique(dm$Comp)==x))
dm$I_Ty <- sapply(dm$Type, function(x) which(unique(dm$Type)==x))
dm$I_Pe  <- sapply(dm$Period, function(x) which(unique(dm$Period)==x))
dm$I_Re  <- sapply(dm$Region, function(x) which(unique(dm$Region)==x))
dm$I_El <- sapply(dm$Element, function(x) which(unique(dm$Element)==x))
dm$I <- 1:nrow(dm)



############################################################
########## Import and format utility data #################
##########################################################

u <- read.csv("Data/Utility_Data.csv", header=TRUE, stringsAsFactors=FALSE)
# Remove reference column for data frame processing
u <- u[, which(colnames(u)!="Reference")]
# Remove crania and caudal verts from dataset
u <- u[!is.na(u$FYM),]
u <- u[which(u$Element!="E_Cranium"),]

# Find element-index pairings with 0s. These are below a detection
# threshold.
zeroes <- sapply(1:nrow(u), 
                 function(x) ifelse(sum(u[x,3:6])==0, x, NA))
zeroes <- zeroes[!is.na(zeroes)]
# Replace with random values below detection threshold
u2 <- u
for(i in zeroes){
  
  j <- u$Index[i]
  minv <- as.matrix(u2[which(u2$Index==j), 3:6])
  minv <- min(minv[which(minv>0)])
  
  u[i,3:6] <- minv*rbeta(4, 2, 2)
  
}
rm(i, j, minv, zeroes, u2)

# Melt data frame
um <- melt(u, id.vars=c("Element", "Index"))
colnames(um)[which(colnames(um)=="variable")] <- "BisonID"
um$BisonID <- as.vector(um$BisonID)
# Scale values within indices to [0, max(value)]
um$value_st <- sapply(1:nrow(um), function(x){
  m <- max(um$value[which(um$Index==um$Index[x])])
  return(um$value[x]/m)
})
      
# Cast a version of the data frame in wide format
uc <- dcast(um, Element + BisonID ~ Index, value.var="value_st")
rm(um)
# Create index values
uc$I <- sapply(uc$Element, function(x) which(unique(uc$Element)==x))

# Elements in in dm that are not in uc
missing_elements <- dm$Element[!(dm$Element %in% uc$Element)]
missing_elements <- unique(missing_elements)

# Append index values to main data frame (dm)
dm$I_Ut <- sapply(dm$Element, function(x){
  if(x %in% uc$Element){
    return(uc$I[which(uc$Element==x)][1])
  }else{
    return(which(missing_elements==x)[1]+max(uc$I))  
  }
})

# For elements which are combinations of existing elements,
# retrieve vectors of indices for the elements that comprise
# those combinations.
index_matrix_utility <- t(sapply(dm$Element, function(x){
  if(!(x %in% missing_elements)){
    return(rep(0, 5))
  }else if(x=="E_Cervical_.1.7._vert"){
    b1 <- uc$I[which(uc$Element=="E_Atlas")][1]
    b2 <- uc$I[which(uc$Element=="E_Axis")][1]
    b3 <- uc$I[which(uc$Element=="E_Cervical_.3.7._vert")][1]
    b45 <- c(0, 0)
    return(c(b1, b2, b3, b45))
  }else if(x=="E_Unspecified_vert"){
    b1 <- uc$I[which(uc$Element=="E_Atlas")][1]
    b2 <- uc$I[which(uc$Element=="E_Axis")][1]
    b3 <- uc$I[which(uc$Element=="E_Cervical_.3.7._vert")][1]
    b4 <- uc$I[which(uc$Element=="E_Thoracic_vert")][1]
    b5 <- uc$I[which(uc$Element=="E_Lumbar_vert")][1]
    return(c(b1, b2, b3, b4, b5))
  }else if(x=="E_Sacrum"){
    b1 <- uc$I[which(uc$Element=="E_Pelvis")][1]
    b25 <- rep(0, 4)
    return(c(b1, b25))
  }else if(x=="E_Innominate"){
    b1 <- uc$I[which(uc$Element=="E_Pelvis")][1]
    b25 <- rep(0, 4)
    return(c(b1, b25))    
  }else if(x=="E_Undifferentiated_carpals"){
    b1 <- uc$I[which(uc$Element=="E_Carpal_units")][1]
    b25 <- rep(0, 4)
    return(c(b1, b25))
  }else if(x=="E_Undifferentiated_tarsals"){
    b1 <- uc$I[which(uc$Element=="E_Tarsal_units")][1]
    b25 <- rep(0, 4)
    return(c(b1, b25))
  }else if(x=="E_Unspecified_metapodial"){
    b1 <- uc$I[which(uc$Element=="E_Metacarpal")][1]
    b2 <- uc$I[which(uc$Element=="E_Metatarsal")][1]
    b35 <- rep(0, 3)
    return(c(b1, b2, b35))
  }else if(x=="E_OlsChub_CarpTars"){
    b1 <- uc$I[which(uc$Element=="E_Tarsal_units")][1]
    b2 <- uc$I[which(uc$Element=="E_Carpal_units")][1]
    b35 <- rep(0, 3)
    return(c(b1, b2, b35))
  }
  
}))

# Rescale values for each index within each index
uc_temp <- lapply(1:nrow(uc), function(x){
  I <- uc[which(uc$I==uc$I[x]),]
  
  rngs <- t(sapply(3:(ncol(I)-1), function(y){
    mi <- min(I[,y])
    ma <- max(I[,y])
    
    if(mi < ma){
      r <- ma-mi
      mi <- mi-0.1*r
      ma <- ma+0.1*r
      if(mi < 0){mi <- 0}
      r <- ma-mi
    }else{
      mi <- 0.98*mi
      ma <- 1.02*ma
      r <- ma-mi
    }

    return(c(mi, r))
  }))
  
  df <- data.frame(GREsc=(uc$GRE[x]-rngs[1,1])/rngs[1,2],
                   MARsc=(uc$MAR[x]-rngs[2,1])/rngs[2,2],
                   PROsc=(uc$PRO[x]-rngs[3,1])/rngs[3,2],
                   SF_sc=(uc$SF_[x]-rngs[4,1])/rngs[4,2],
                   SKFsc=(uc$SKF[x]-rngs[5,1])/rngs[5,2],
                   TF_sc=(uc$TF_[x]-rngs[6,1])/rngs[6,2],
                   TP_sc=(uc$TP_[x]-rngs[7,1])/rngs[7,2],
                   GREmin=rngs[1,1], MARmin=rngs[2,1], 
                   PROmin=rngs[3,1], SF_min=rngs[4,1],
                   SKFmin=rngs[5,1], TF_min=rngs[6,1],
                   TP_min=rngs[7,1], 
                   GREran=rngs[1,2], MARran=rngs[2,2], 
                   PROran=rngs[3,2], SF_ran=rngs[4,2],
                   SKFran=rngs[5,2], TF_ran=rngs[6,2],
                   TP_ran=rngs[7,2])
  return(df)
  
})
uc_temp <- do.call("rbind", uc_temp)
uc <- cbind(uc, uc_temp)
rm(uc_temp)

# Make a data from of unique element utilities for later 
# formatting code
u_unique <- uc[,c(10, 18:31)]
u_unique <- u_unique[!duplicated(u_unique), ]


############################################################
########## Organize density data ##########################
##########################################################

# Read table of bison bone densities and calculate volume
# density (VD).
densities <- read.csv("Data/Density_Data.csv", header=TRUE,
                      stringsAsFactors=FALSE)
densities$VD <- densities$MEAN_BMC/densities$Volume

# Remove referene column for data frame processing
densities <- densities[,which(colnames(densities)!="Reference")]

# Read table of density values for wildebeest. These scale
# density estimates that do not account for internal cavities
# (BMD1) to density estimates that do account for internal
# cavities (BMD2).
densities_adj <- read.csv("Data/AdjDensity_Data.csv", header=TRUE,
                  stringsAsFactors=FALSE)

# Remove referene column for data frame processing
densities_adj <- densities_adj[,which(colnames(densities_adj)!="Reference")]

# Adjust bison bone densities to account for internal cavities.
densities$VD_adj <- sapply(1:nrow(densities), function(x){
  s <- which(densities_adj$ScanSite==densities$ScanSite[x])
  adj <- densities_adj$adjs[s][1]
  return(adj*densities$VD[x])
})


# Assign index values to scan sites
densities$I <- sapply(densities$ScanSite, function(x){
  which(unique(densities$ScanSite)==x)
})

# Rescale density values within scan sites to be fitted
# to beta distributions
adj_df <- lapply(1:nrow(densities), function(x){
  I <- densities[which(densities$I==densities$I[x]), ]
  DENran <- max(I$VD_adj) - min(I$VD_adj)
  DENmin <- min(I$VD_adj) - 0.1*DENran
  VDmax <- max(I$VD_adj) + 0.1*DENran
  if(DENmin < 0) DENmin <- 0
  DENran <- VDmax - DENmin
  if(DENran==0){
    DENmin <- 0.98*DENmin
    VDmax <- 1.02*VDmax
    DENran <- VDmax - DENmin
  }
  DENsc <- (densities$VD_adj[x]-DENmin)/DENran
  return(data.frame(DENsc=DENsc, DENran=DENran, DENmin=DENmin))
})
adj_df <- do.call("rbind", adj_df)
densities <- cbind(densities, adj_df)
rm(adj_df)

# Create a list of scan sites for each element or element group.
Elmnt_ScnSts <- list(E_Cranium=NA,
                    E_Mandible=c("DN1", "DN2", "DN3", "DN4", 
                                 "DN5", "DN6", "DN7", "DN8"),
                    E_Atlas=c("AT1", "AT2", "AT3"),
                    E_Axis=c("AX1", "AX2", "AX3"),
                    E_Cervical_.3.7._vert=c("CE1", "CE2"),
                    E_Thoracic_vert=c("TH1", "TH2"),
                    E_Lumbar_vert=c("LU1", "LU2", "LU3"),
                    E_Caudal_vert=NA,
                    E_Sacrum=c("SC1", "SC2"),
                    E_Pelvis=c("IL1", "IL2", "AC1", "PU1", "PU2",
                               "IS1", "IS2", "SC1", "SC2"),
                    E_Rib=c("RI1", "RI2", "RI3", "RI4", "RI5"),
                    E_Scapula=c("SP1", "SP2", "SP3", "SP4", "SP5"),
                    E_Humerus=c("HU1", "HU2", "HU3", "HU4", "HU5"),
                    E_Radius_ulna=c("UL1", "UL2", "RA1", "RA2", "RA3",
                                    "RA4", "RA5"),
                    E_Carpal_units=c("2&3", "LUN", "TM", "UNC", "SCAPH", "FIF"),
                    E_Undifferentiated_carpals=c("2&3", "LUN", "TM", 
                                                 "UNC", "SCAPH", "FIF"),
                    E_Metacarpal=c("MC1", "MC2", "MC3", "MC4", "MC5", "MC6"),
                    E_Innominate=c("IL1", "IL2", "AC1", "PU1", "PU2",
                                   "IS1", "IS2"),
                    E_Femur=c("FE1", "FE2", "FE3", "FE4", "FE5", "FE6"),
                    E_Tibia=c("TI1", "TI2", "TI3", "TI4", "TI5", "TI6"),
                    E_Tarsal_units=c("CA1", "CA2", "CA3", "CA4", "AS1", "AS2",
                                     "AS3", "NC1", "NC2", "NC3", "CUN"),
                    E_Undifferentiated_tarsals=c("CA1", "CA2", "CA3", "CA4", 
                                                 "AS1", "AS2",  "AS3", "NC1", 
                                                 "NC2", "NC3", "CUN"),
                    E_Metatarsal=c("MR1", "MR2", "MR3", "MR4", "MR5", "MR6"),
                    E_Unspecified_phalanges=c("P11", "P12", "P13", "P21", 
                                              "P23", "P31"))

# Elements with scan sites
elements_wss <- sapply(1:length(Elmnt_ScnSts), function(z) ls(Elmnt_ScnSts[z]))
# Elements which are comprised of scan site combinations
elements_nss <- c("E_OlsChub_CarpTars", "E_Unspecified_metapodial", 
                  "E_Unspecified_vert", "E_Cervical_.1.7._vert")

# Find densest scan site index for each element or element group
dm$I_De <- sapply(dm$Element, function(x){
  
  if(x %in% elements_wss){
    
    # Get scan sites for element x
    pos <- which(elements_wss==x)
    scansites <- Elmnt_ScnSts[[pos]]
    
    # Find median value for each scan site for element x
    medians <- sapply(scansites, function(y){
      median(densities$VD_adj[which(densities$ScanSite==y)])
    })
    
    # Get scan site with highest median density for most elements,
    # unless undifferentiated carpals or tarsals, in which case, 
    # get the site with the median density.
    if(!(x %in% c("E_Undifferentiated_carpals", 
                  "E_Undifferentiated_tarsals",
                  "E_Unspecified_phalanges"))){
      sc <- scansites[which.max(medians)]
    }else{
      ssmed <- which.min(abs(medians - median(medians)))
      sc <- scansites[ssmed]
    }
    
    
    # Return index value for that scan site
    return(densities$I[which(densities$ScanSite==sc)][1])
    
  }else{
    
    # Return index for element that is represented by a combination
    # of scan sites.
    return(max(densities$I)+which(elements_nss==x)[1])
  }
  
})

# Make a data from of unique element utilities for later 
# formatting code
d_unique <- densities[,c(7, 9:10)]
d_unique <- d_unique[!duplicated(d_unique), ]

#################################################################
############## FINAL DATA FORMATTING STEPS #####################
###############################################################

# For elements which are combinations of existing elements,
# retrieve vectors of density indices for the elements that
# comprise those combinations.
index_matrix_density <- t(sapply(dm$Element, function(x){
  
  indices <- rep(0, 5)
  
  if(x %in% elements_nss){

    if(x=="E_Cervical_.1.7._vert"){
      el <- c("E_Atlas", "E_Axis", "E_Cervical_.3.7._vert")
    }else if(x=="E_Unspecified_vert"){
      el <- c("E_Atlas", "E_Axis", "E_Cervical_.3.7._vert",
              "E_Thoracic_vert", "E_Lumbar_vert")
    }else if(x=="E_Unspecified_metapodial"){
      el <- c("E_Metacarpal", "E_Metatarsal")
    }else if(x=="E_OlsChub_CarpTars"){
      el <- c("E_Undifferentiated_carpals", 
              "E_Undifferentiated_tarsals") 
    }
    
    for(i in 1:length(el)){
      
      # Get scan sites for element x
      pos <- which(elements_wss==el[i])
      scansites <- Elmnt_ScnSts[[pos]]
      
      # Find median value for each scan site for element x
      medians <- sapply(scansites, function(y){
        median(densities$VD_adj[which(densities$ScanSite==y)])
      })
      
      sc <- scansites[which.max(medians)]
      indices[i] <- densities$I[which(densities$ScanSite==sc)][1]
    }
  }

    return(indices)
    
}))


# Number of element combination utility and density values
# for the model
u_comb <- sapply(1:ncol(index_matrix_utility), function(x){
  length(which(index_matrix_utility[,x]>0))
})
d_comb <- sapply(1:ncol(index_matrix_density), function(x){
  length(which(index_matrix_density[,x]>0))
})

# Remove row names from matrix rows
index_matrix_utility <- unname(index_matrix_utility)
index_matrix_density <- unname(index_matrix_density)


# Append unknown density values to density data frame
densities$I_dm <- rep(0, nrow(densities))
densities$I_combo <- rep(0, nrow(densities))
ddf1 <- data.frame(Specimen=dm$Element,
                   ScanSite=rep("Unkn", nrow(dm)),
                   Volume=rep(99, nrow(dm)),
                   MEAN_BMC=rep(99, nrow(dm)),
                   VD=rep(99, nrow(dm)),
                   VD_adj=rep(99, nrow(dm)),
                   I=dm$I_De, 
                   DENsc=rep(99, nrow(dm)),
                   DENran=rep(99, nrow(dm)),
                   DENmin=rep(99, nrow(dm)),
                   I_dm=1:nrow(dm),
                   I_combo=rep(0, nrow(dm)),
                   stringsAsFactors=FALSE)

# Append unknown density values for combinations
combospecs <- lapply(1:ncol(index_matrix_density), function(i){
  
  inds <- which(index_matrix_density[, i]>0)
  tdf <- index_matrix_density[which(index_matrix_density[,i]>0), i]
  tdf2 <- data.frame(Specimen=dm$Element[inds],
                     ScanSite=rep(paste0("Unkn", i), length(tdf)),
                     Volume=rep(99, length(tdf)),
                     MEAN_BMC=rep(99, length(tdf)),
                     VD=rep(99, length(tdf)),
                     VD_adj=rep(99, length(tdf)),
                     I=tdf, 
                     DENsc=rep(99, length(tdf)),
                     DENran=rep(99, length(tdf)),
                     DENmin=rep(99, length(tdf)),
                     I_dm=inds,
                     I_combo=1:length(tdf),
                     stringsAsFactors=FALSE)
  
  return(tdf2)
})
combospecs <- do.call("rbind", combospecs)
densities <- do.call("rbind", list(ddf1, combospecs, densities))
# Remove combo rows
densities <- densities[-which(index_matrix_density[,1]>0),]

rm(combospecs, ddf1, densities_adj, index_matrix_density)

# Append unknown utility values to utility data frame
uc$I_dm <- rep(0, nrow(uc))
uc$I_combo <- rep(0, nrow(uc))
udf1 <- data.frame(Element=dm$Element,
                   BisonID=rep("Unkn", nrow(dm)),
                   GRE=rep(99, nrow(dm)),
                   MAR=rep(99, nrow(dm)),
                   PRO=rep(99, nrow(dm)),
                   SF_=rep(99, nrow(dm)),
                   SKF=rep(99, nrow(dm)),
                   TF_=rep(99, nrow(dm)),
                   TP_=rep(99, nrow(dm)),
                   I=dm$I_Ut, 
                   GREsc=rep(99, nrow(dm)),
                   MARsc=rep(99, nrow(dm)),
                   PROsc=rep(99, nrow(dm)),
                   SF_sc=rep(99, nrow(dm)),
                   SKFsc=rep(99, nrow(dm)),
                   TF_sc=rep(99, nrow(dm)),
                   TP_sc=rep(99, nrow(dm)),
                   GREmin=rep(99, nrow(dm)),
                   MARmin=rep(99, nrow(dm)),
                   PROmin=rep(99, nrow(dm)),
                   SF_min=rep(99, nrow(dm)),
                   SKFmin=rep(99, nrow(dm)),
                   TF_min=rep(99, nrow(dm)),
                   TP_min=rep(99, nrow(dm)),
                   GREran=rep(99, nrow(dm)),
                   MARran=rep(99, nrow(dm)),
                   PROran=rep(99, nrow(dm)),
                   SF_ran=rep(99, nrow(dm)),
                   SKFran=rep(99, nrow(dm)),
                   TF_ran=rep(99, nrow(dm)),
                   TP_ran=rep(99, nrow(dm)),
                   I_dm=1:nrow(dm),
                   I_combo=rep(0, nrow(dm)),
                   stringsAsFactors=FALSE)

# Append unknwon utility values for combinations
combospecs <- lapply(1:ncol(index_matrix_utility), function(i){
  
  inds <- which(index_matrix_utility[,i]>0)
  tdf <- index_matrix_utility[which(index_matrix_utility[,i]>0), i]
  tdf2 <- data.frame(Element=dm$Element[inds],
                     BisonID=rep(paste0("Unkn", i), length(tdf)),
                     GRE=rep(99, length(tdf)),
                     MAR=rep(99, length(tdf)),
                     PRO=rep(99, length(tdf)),
                     SF_=rep(99, length(tdf)),
                     SKF=rep(99, length(tdf)),
                     TF_=rep(99, length(tdf)),
                     TP_=rep(99, length(tdf)),
                     I=tdf, 
                     GREsc=rep(99, length(tdf)),
                     MARsc=rep(99, length(tdf)),
                     PROsc=rep(99, length(tdf)),
                     SF_sc=rep(99, length(tdf)),
                     SKFsc=rep(99, length(tdf)),
                     TF_sc=rep(99, length(tdf)),
                     TP_sc=rep(99, length(tdf)),
                     GREmin=rep(99, length(tdf)),
                     MARmin=rep(99, length(tdf)),
                     PROmin=rep(99, length(tdf)),
                     SF_min=rep(99, length(tdf)),
                     SKFmin=rep(99, length(tdf)),
                     TF_min=rep(99, length(tdf)),
                     TP_min=rep(99, length(tdf)),
                     GREran=rep(99, length(tdf)),
                     MARran=rep(99, length(tdf)),
                     PROran=rep(99, length(tdf)),
                     SF_ran=rep(99, length(tdf)),
                     SKFran=rep(99, length(tdf)),
                     TF_ran=rep(99, length(tdf)),
                     TP_ran=rep(99, length(tdf)),
                     I_dm=inds,
                     I_combo=1:length(tdf),
                     stringsAsFactors=FALSE)
  
  return(tdf2)
})
combospecs <- do.call("rbind", combospecs)
uc <- do.call("rbind", list(udf1, combospecs, uc))
# Remove combo rows
uc <- uc[-which(index_matrix_utility[,1]>0),]

rm(combospecs, udf1, index_matrix_utility)

# Find the indices in densities and uc that correspond
# to dm and append these to dm.
De_cols <- max(data.frame(table(densities$I_dm[which(densities$I_dm>0)]))[,2])
I_De_Rev <- sapply(1:De_cols,function(y){
    r <- sapply(1:nrow(dm), function(x){
       if(length(which(densities$I_dm==x)) > (y-1)){
         return(which(densities$I_dm==x)[y])
       }else{
         return(0)
       }
    })
    return(r)
 })

Ut_cols <- max(data.frame(table(uc$I_dm[which(uc$I_dm>0)]))[,2])
I_Ut_Rev <- sapply(1:Ut_cols,function(y){
  r <- sapply(1:nrow(dm), function(x){
    if(length(which(uc$I_dm==x)) > (y-1)){
      return(which(uc$I_dm==x)[y])
    }else{
      return(0)
    }
  })
  return(r)
})


### Replace minimum and range 99 values in densities and uc with
### with known values
densities$DENran <- sapply(1:nrow(densities), function(x){
  if(densities$DENran[x] != 99){
    return(densities$DENran[x])
  }else{
    I <- densities$I[x]
    return(d_unique$DENran[which(d_unique$I==I)[1]])
  }
})

densities$DENmin <- sapply(1:nrow(densities), function(x){
  if(densities$DENmin[x] != 99){
    return(densities$DENmin[x])
  }else{
    I <- densities$I[x]
    return(d_unique$DENmin[which(d_unique$I==I)[1]])
  }
})

uc$GREmin <- sapply(1:nrow(uc), function(x){
  if(uc$GREmin[x] != 99){
    return(uc$GREmin[x])
  }else{
    I <- uc$I[x]
    return(u_unique$GREmin[which(u_unique$I==I)[1]])
  }
})
uc$GREran <- sapply(1:nrow(uc), function(x){
  if(uc$GREran[x] != 99){
    return(uc$GREran[x])
  }else{
    I <- uc$I[x]
    return(u_unique$GREran[which(u_unique$I==I)[1]])
  }
})

uc$MARmin <- sapply(1:nrow(uc), function(x){
  if(uc$MARmin[x] != 99){
    return(uc$MARmin[x])
  }else{
    I <- uc$I[x]
    return(u_unique$MARmin[which(u_unique$I==I)[1]])
  }
})
uc$MARran <- sapply(1:nrow(uc), function(x){
  if(uc$MARran[x] != 99){
    return(uc$MARran[x])
  }else{
    I <- uc$I[x]
    return(u_unique$MARran[which(u_unique$I==I)[1]])
  }
})

uc$PROmin <- sapply(1:nrow(uc), function(x){
  if(uc$PROmin[x] != 99){
    return(uc$PROmin[x])
  }else{
    I <- uc$I[x]
    return(u_unique$PROmin[which(u_unique$I==I)[1]])
  }
})
uc$PROran <- sapply(1:nrow(uc), function(x){
  if(uc$PROran[x] != 99){
    return(uc$PROran[x])
  }else{
    I <- uc$I[x]
    return(u_unique$PROran[which(u_unique$I==I)[1]])
  }
})

uc$SF_min <- sapply(1:nrow(uc), function(x){
  if(uc$SF_min[x] != 99){
    return(uc$SF_min[x])
  }else{
    I <- uc$I[x]
    return(u_unique$SF_min[which(u_unique$I==I)[1]])
  }
})
uc$SF_ran <- sapply(1:nrow(uc), function(x){
  if(uc$SF_ran[x] != 99){
    return(uc$SF_ran[x])
  }else{
    I <- uc$I[x]
    return(u_unique$SF_ran[which(u_unique$I==I)[1]])
  }
})

uc$SKFmin <- sapply(1:nrow(uc), function(x){
  if(uc$SKFmin[x] != 99){
    return(uc$SKFmin[x])
  }else{
    I <- uc$I[x]
    return(u_unique$SKFmin[which(u_unique$I==I)[1]])
  }
})
uc$SKFran <- sapply(1:nrow(uc), function(x){
  if(uc$SKFran[x] != 99){
    return(uc$SKFran[x])
  }else{
    I <- uc$I[x]
    return(u_unique$SKFran[which(u_unique$I==I)[1]])
  }
})

uc$TF_min <- sapply(1:nrow(uc), function(x){
  if(uc$TF_min[x] != 99){
    return(uc$TF_min[x])
  }else{
    I <- uc$I[x]
    return(u_unique$TF_min[which(u_unique$I==I)[1]])
  }
})
uc$TF_ran <- sapply(1:nrow(uc), function(x){
  if(uc$TF_ran[x] != 99){
    return(uc$TF_ran[x])
  }else{
    I <- uc$I[x]
    return(u_unique$TF_ran[which(u_unique$I==I)[1]])
  }
})

uc$TP_min <- sapply(1:nrow(uc), function(x){
  if(uc$TP_min[x] != 99){
    return(uc$TP_min[x])
  }else{
    I <- uc$I[x]
    return(u_unique$TP_min[which(u_unique$I==I)[1]])
  }
})
uc$TP_ran <- sapply(1:nrow(uc), function(x){
  if(uc$TP_ran[x] != 99){
    return(uc$TP_ran[x])
  }else{
    I <- uc$I[x]
    return(u_unique$TP_ran[which(u_unique$I==I)[1]])
  }
})

# Get number of components for each site type
t_N <- sapply(1:max(dm$I_Ty), function(x){
  df <- dm[which(dm$I_Ty==x),]
  return(length(unique(df$I_Co)))
})
dm$I_Ty_Co <- sapply(1:nrow(dm), function(x){
  t <- dm$I_Ty[x]
  c <- dm$I_Co[x]
  cs <- unique(dm$I_Co[which(dm$I_Ty==t)])
  return(which(cs==c))
})

# Get unique preditor variables
predictors <- c(colnames(uc)[which(nchar(colnames(uc))==5)],
                colnames(densities)[which(nchar(colnames(densities))==5)])

ModelData <- lapply(predictors, function(x){
  
  if(x=="DENsc"){
    df <- densities
  }else{
    df <- uc
  }
  
  ModelData <- list(model_name=substr(x, 1, 3),
                    N=nrow(dm), 
                    I=dm$I, MNE=dm$MNE, 
                    logAF=dm$logAF,
                    n_c=max(dm$I_Co),
                    I_Component=dm$I_Co,
                    I_Type=dm$I_Ty,
                    n_t=max(dm$I_Ty),
                    t_N=t_N,
                    I_Ty_Co=dm$I_Ty_Co,
                    
                    n_=max(df$I), N_=nrow(df),
                    
                    pred1=df[,which(colnames(df)==x)],
                    pred1range=df[,which(colnames(df)==
                                         paste0(substr(x, 1, 3), 
                                                "ran"))],
                    pred1min=df[,which(colnames(df)==
                                         paste0(substr(x, 1, 3), 
                                                "min"))],
                    I_=df$I, 
                    N_impute=length(which(df[,4]==99)),
                    i_impute=1:length(which(df[,4]==99)))
  
  if(x=="DENsc"){
    ModelData[[length(ModelData)+1]] <- dm$I_De
    ModelData[[length(ModelData)+1]] <- I_De_Rev
    ModelData[[length(ModelData)+1]] <- De_cols
  }else{
    ModelData[[length(ModelData)+1]] <- dm$I_Ut
    ModelData[[length(ModelData)+1]] <- I_Ut_Rev
    ModelData[[length(ModelData)+1]] <- Ut_cols
  } 
  
  names(ModelData)[(length(ModelData)-2):length(ModelData)] <-
    c("I_Element", "I_Rev", "I_cols")
  
  return(ModelData)
  
})

# Create another set of models with two predictors, the first
# of which is density.
predictors2 <- predictors[which(predictors!="DENsc")]
ModelData2 <- lapply(predictors2, function(x){
  
  mn <- paste0("DEN", substr(x, 1, 3))
  
  ModelData2 <- list(model_name=mn,
                    N=nrow(dm), 
                    I=dm$I, MNE=dm$MNE, 
                    logAF=dm$logAF,
                    n_c=max(dm$I_Co),
                    I_Component=dm$I_Co,
                    I_Type=dm$I_Ty,
                    n_t=max(dm$I_Ty),
                    t_N=t_N,
                    I_Ty_Co=dm$I_Ty_Co,
                    
                    n_1=max(uc$I), N_1=nrow(uc),
                    
                    pred1=uc[,which(colnames(uc)==x)],
                    pred1range=uc[,which(colnames(uc)==
                                           paste0(substr(x, 1, 3), 
                                                  "ran"))],
                    pred1min=uc[,which(colnames(uc)==
                                         paste0(substr(x, 1, 3), 
                                                "min"))],
                    I_1=uc$I, 
                    N_impute1=length(which(uc[,4]==99)),
                    i_impute1=1:length(which(uc[,4]==99)),
                    
                    n_2=max(densities$I), N_2=nrow(densities),
                    
                    pred2=densities$DENsc,
                    pred2range=densities$DENran,
                    pred2min=densities$DENmin,
                    I_2=densities$I, 
                    N_impute2=length(which(densities[,4]==99)),
                    i_impute2=1:length(which(densities[,4]==99)),
                    
                    I_Element1=dm$I_Ut,
                    I_Rev1=I_Ut_Rev,
                    I_cols1=Ut_cols,
                    
                    I_Element2=dm$I_De,
                    I_Rev2=I_De_Rev,
                    I_cols2=De_cols)
  
  return(ModelData2)
  
})


# Save output
save(ModelData, ModelData2, dm, densities, uc, file="01_Out.RData")
