// This code fits a model for bison element frequencies
// using two predictor variables.

data{
 // Total observations
 int<lower=1> N;
 // Number of components
 int<lower=1> n_c;
 // Index of observations
 int I[N];
 // Count outcome
 int MNE[N];
 // Exposure
 real logAF[N];
 // Index of element types across observations. Element
 // types can vary depending on the type of predictor
 // variable (e.g., bone density vs dietary utility),
 // so these index variables are assumed here to differ.
 int I_Element1[N];
 int I_Element2[N];
 
 // Site type index
 int I_Type[N];
  // Component indices
 int I_Component[N];
 // Component indices nested within site types; In other
 // words, within each site type, components recieve an
 // index value independent of their I_Component value,
 // which indexes all components in a sequence, regardless
 // of the site type to which they belong.
 int I_Ty_Co[N];
 // Number of site types
 int n_t;
 // Number of components per site type
 int t_N[n_t];
 
 // Number of total predictor 1 values (estimated + 
 // observed) across all observations for predictor 1. 
 // This can exceed N, because some some predictor 
 // values are combinations of values (e.g.,
 // metapodial = (metatarsal + metacarpal)/2).
 int<lower=1> N_1;
 // Number unique element types for predictor 1
 int<lower=1> n_1;
 // Vector of predictor 1 values
 real<lower=0> pred1[N_1];  
 // Index of element types in predictor 1 data
 int I_1[N_1]; 
 // Range of predictor 1 values per element type
 real<lower=0> pred1range[N_1];
 // Minimum predictor 1 value per element type
 real<lower=0> pred1min[N_1];
 
 ///// Predictor 1 imputation variables ////
 ///// (i.e., N_1 - observed) ////
 // Number of imputed predictor 1 values
 int<lower=1> N_impute1;
 // Vector of positions for values to impute within
 // pred1
 int<lower=1> i_impute1[N_impute1];
 // Maximum number of pred1 values in a combination
 // for a single predictor 1 value
 int I_cols1;
 // Matrix that specifies indices in pred1 for combi-
 // nations of predictor 1 values across N observations
 int I_Rev1[N, I_cols1];
 
 // Number of total predictor 2 values (estimated + 
 // observed) across all observations for predictor 2. 
 // This can exceed N, because some some predictor 2
 // values are combinations of values (e.g.,
 // metapodial = (metatarsal + metacarpal)/2).
 int<lower=1> N_2;
 // Number unique element types for predictor 2
 int<lower=1> n_2;
 // Vector of predictor 2 values
 real<lower=0> pred2[N_2];  
// Index of element types in predictor 2 data
 int I_2[N_2]; 
 // Range of predictor values 2 per element type
 real<lower=0> pred2range[N_2];
 // Minimum predictor 2 value per element type
 real<lower=0> pred2min[N_2];
 
 ///// Predictor 2 imputation variables ////
 ///// (i.e., N_2 - observed) ////
 // Number of imputed predictor 2 values
 int<lower=1> N_impute2;
 // Vector of positions for values to impute
 // within pred2
 int<lower=1> i_impute2[N_impute2];
 // Maximum number of pred2 values in a combination
 // for a single predictor 2 value
 int I_cols2;
 // Matrix that specifies indices in pred2 for combi-
 // nations of predictor 2 values across N observations
 int I_Rev2[N, I_cols2];

}

parameters{
  
 // Global intercept parameters, one for
 // each site type.
 vector[n_t] t_alpha;
 // Beta parameters for predictor 1, 
 // predictor 2, and an interaction, as
 // a vector of site types for each beta.
 vector[n_t] t_b[3];
 
 // Observation level overdispersion parameters,
 // one value for each observation.
 vector[N] I_offset;   
 
 // Z-scores for component-level varying effects,
 // c_alpha, c_b1, c_b2, c_b12. There is one
 // z-score matrix for each site type.
 matrix[4, n_c] z_c[n_t];
 cholesky_factor_corr[4] LRho_c[n_t];
 // Scaling parameters for z_c.
 vector<lower=0>[4] scale_c[n_t];
 
 // Mean location parameters for element type
 // values in predictor 1 and predictor 2.
 vector[2] mu;
 // SD parameters for element type values in
 // predictor 1 and predictor 2.
 vector<lower=0>[2] sigma;
 // Shape parameters for distributions of predictor 1
 // and predictor 2 values; beta distributed.
 vector<lower=0>[2] shape;
 // Mu offsets for each mean element type value within
 // predictor 1 and predictor 2. These are in z-score
 // units, as they are scaled by sigma[1] and sigma[2].
 vector[n_1] loc1;
 vector[n_2] loc2;
 
 // Predictor 1 and predictor 2 are modelled as beta
 // distributed: beta(LOC*shape, (1-LOC)*shape)
 // LOC1/2 is on the (0, 1) scale, but the parameter 
 // loc1/2 specified here is unbounded, being inv_logit
 // tranformed in the linear model to obtain LOC1/2.
 
 /// The next two parameters relate to the overdispersion
 /// of observations within each component. This over-
 /// dispersion is controlled by an SD parameter that
 /// varies between components and is distributed:
 ///  SD ~ gamma(gloc/gshape, 1/gshape)
 real<lower=0> gloc;
 real<lower=0> gshape;
 // Vector of scaling parameters for overdispersed
 // observations in components
 vector<lower=0>[n_c] I_scale;
 
 // Imputed predictor 1 and predictor 2 values
 vector<lower=0, upper=1>[N_impute1] val1;
 vector<lower=0, upper=1>[N_impute2] val2;
 
}
transformed parameters{
 
 // Initialize component level parameters
 // for alpha and b parameters.
 matrix[n_c, 4] v_c[n_t];
 vector[n_c] c_alpha[n_t];
 vector[n_c] c_b1[n_t];
 vector[n_c] c_b2[n_t];
 vector[n_c] c_b12[n_t];
 matrix[4,4] Rho_c[n_t];

 // Cholesky decomposition for a matrix of component
 // level parameters for each site type (the loop
 // traverses matrixes for each site type).
 for (j in 1:n_t) {
   
  v_c[j][1:t_N[j], ] = (diag_pre_multiply(scale_c[j], 
                        LRho_c[j]) * z_c[j][ , 1:t_N[j]])';
  
  // Since the number of components is uneven between site
  // types, there are unfilled bins in each matrix. Fill these
  // bins with zeros.
  v_c[j][(t_N[j] + 1):n_c, ] = rep_matrix(0, n_c - t_N[j], 4);     
  
  c_alpha[j] = col(v_c[j], 1);
  c_b1[j] = col(v_c[j], 2);
  c_b2[j] = col(v_c[j], 3);
  c_b12[j] = col(v_c[j], 4);
  
  Rho_c[j] = LRho_c[j] * LRho_c[j]';
 
 }
 
}
model{
 
 // Vector of merged imputed + observed values for
 // predictor 1 and predictor 2.
 vector[N_1] merged1;
 vector[N_2] merged2;
 // Vector of merged values put back on predictor
 // scale from the (0, 1) scale that was used to model
 // values on the beta scale.
 vector[N_1] pred1sc;
 vector[N_2] pred2sc;
 
 // Lambdas for Poisson outcome
 vector[N] lambda; 
 // Predictor 1 and predictor 2 beta distribution
 // location values, (0, 1) scale.
 vector[N_1] LOC1;
 vector[N_2] LOC2;
 // Linear combinations of alpha parameters
 vector[N] ALPHA;    
 // Linear combinations of each beta parameter
 vector[N] BETA[3];    
 // Across N observations: Estimated predictor 1
 // and predictor 2 values. Some of these are 
 // combinations of estimated values.
 vector[N] PRED[2];

 // Prior distributions for global alpha and beta
 // parameters for each site type.
 t_alpha ~ normal(0, 1);
 for (j in 1:3) {
   t_b[j] ~ normal(0, 1);
 }
 
 
 // Priors for hyper parameters for components 
 // nested within site types.
 for (j in 1:n_t) {
  LRho_c[j] ~ lkj_corr_cholesky(4);
  scale_c[j] ~ exponential(2);
  to_vector(z_c[j]) ~ normal(0, 1);
 }
 
 // Priors and hyper priors for observation
 // level over dispersion parameters.
 I_offset ~ normal(0, 1);
 I_scale ~ gamma(gloc/gshape, 1/gshape);
 gloc ~ exponential(1);
 gshape ~ exponential(2);
 
 // Priors and hyper priors for predictor 1
 // and predictor 2 values.
 shape ~ exponential(2);
 mu ~ normal(0, 0.5);
 sigma ~ exponential(2);
 loc1 ~ normal(0, 1);
 loc2 ~ normal(0, 1);
 
 // Merge imputed and observed predictor 1 values.
 for (k in 1:N_1) {
  if (k > N_impute1) {
   merged1[k] = pred1[k];
  }else{
   merged1[k] = val1[k];
  }
  // Put merged (0, 1) values back on predictor 1 scale
  pred1sc[k] = pred1range[k] * merged1[k] + pred1min[k];
 }
 
 // Model for predictor 1 values, where beta
 // location parameter varies by element type.
 for (k in 1:N_1) {
  LOC1[k] = inv_logit(mu[1] + loc1[I_1[k]] * sigma[1]);
 }
 
 // Likelihood for merged predictor 1 values
 merged1 ~ beta(LOC1 * (shape[1] + 2), 
              (1 - LOC1) * (shape[1] + 2));
              
 // Merge imputed and observed predictor 2 values.
 for (k in 1:N_2) {
  if (k > N_impute2) {
   merged2[k] = pred2[k];
  }else{
   merged2[k] = val2[k];
  }
  // Put merged (0, 1) values back on predictor 2 scale
  pred2sc[k] = pred2range[k] * merged2[k] + pred2min[k];
 }
 
 // Model for predictor 2 values, where beta
 // location parameter varies by element type.
 for (k in 1:N_2) {
  LOC2[k] = inv_logit(mu[2] + loc2[I_2[k]] * sigma[2]);
 }
 
 // Likelihood for merged predictor 2 values
 merged2 ~ beta(LOC2 * (shape[2] + 2), 
               (1 - LOC2) * (shape[2] + 2));
 
 // LOOP FOR LINEAR MODEL
 for (j in 1:N) {
 
   // OBSERVATION j PREDICTOR 1 VALUE:
   // If the element index indicates that observation j
   // corresponds to an element associated with a predictor 1
   // value that was directly imputed, obtain the imputed
   // value that has been rescaled from (0, 1) back to the
   // predictor 1 scale. Otherwise, observation j corresponds
   // to an element that is some combination of predictor 1
   // values. In this case, combine the necessary rescaled
   // predictor 1 values.
   if (I_Element1[j] <= n_1) {
    PRED[1][j] = pred1sc[I_Rev1[j, 1]];
    
   } else if (I_Rev1[j, 2] == 0) {
     // Element type corresponds to another already specified
     // element type, e.g., "pelvis" = "sacrum"
    PRED[1][j] = pred1sc[I_Rev1[j, 1]];
    
   } else if (I_Rev1[j, 3] == 0) {
     // Various element types that are aggregates ot two element
     // types, e.g., "metapodial" = ("metacarpal" + "metatarsal")/2.
    PRED[1][j] = pred1sc[I_Rev1[j, 1]] * 0.5 + pred1sc[I_Rev1[j, 2]] * 0.5;

   } else if (I_Rev1[j, 4] == 0) {
     // "cervical 1-7 vert"
    PRED[1][j] = pred1sc[I_Rev1[j, 1]]/7 + pred1sc[I_Rev1[j, 2]]/7
             + pred1sc[I_Rev1[j, 3]] * 5/7;

   } else if (I_Rev1[j, 5] > 0) {
     // "unspecified vert"
    PRED[1][j] = pred1sc[I_Rev1[j, 1]]/26 + pred1sc[I_Rev1[j, 2]]/26 
             + pred1sc[I_Rev1[j, 3]] * 5/26 + pred1sc[I_Rev1[j, 4]] * 7/13 
             + pred1sc[I_Rev1[j, 5]] * 5/26;

   }
   
   // OBSERVATION j PREDICTOR 2 VALUE:
   // If the element index indicates that observation j
   // corresponds to an element associated with a predictor 2
   // value that was directly imputed, obtain the imputed
   // value that has been rescaled from (0, 1) back to the
   // predictor 2 scale. Otherwise, observation j corresponds
   // to an element that is some combination of predictor 2
   // values. In this case, combine the necessary rescaled
   // predictor 2 values.
   if (I_Element2[j] <= n_2) {
    PRED[2][j] = pred2sc[I_Rev2[j, 1]];

   } else if (I_Rev2[j, 2] == 0) {
     // Element type corresponds to another already specified
     // element type, e.g., "pelvis" = "sacrum"
    PRED[2][j] = pred2sc[I_Rev2[j, 1]];

   } else if (I_Rev2[j, 3] == 0) {
     // Various element types that are aggregates ot two element
     // types, e.g., "metapodial" = ("metacarpal" + "metatarsal")/2.
    PRED[2][j] = pred2sc[I_Rev2[j, 1]] * 0.5 + pred2sc[I_Rev2[j, 2]] * 0.5;

   } else if (I_Rev2[j, 4] == 0) {
     // "cervical 1-7 vert"
    PRED[2][j] = pred2sc[I_Rev2[j, 1]]/7 + pred2sc[I_Rev2[j, 2]]/7
             + pred2sc[I_Rev2[j, 3]] * 5/7;

   } else if (I_Rev2[j, 5] > 0) {
     // "unspecified vert"
    PRED[2][j] = pred2sc[I_Rev2[j, 1]]/26 + pred2sc[I_Rev2[j, 2]]/26 
             + pred2sc[I_Rev2[j, 3]] * 5/26 + pred2sc[I_Rev2[j, 4]] * 7/13 
             + pred2sc[I_Rev2[j, 5]] * 5/26;

   }
 
   // Intercept linear model
   ALPHA[j] = t_alpha[I_Type[j]] + c_alpha[I_Type[j]][I_Ty_Co[j]];
   // Predictor effects linear beta models
   BETA[1][j] = t_b[1][I_Type[j]] + c_b1[I_Type[j]][I_Ty_Co[j]];
   BETA[2][j] = t_b[2][I_Type[j]] + c_b2[I_Type[j]][I_Ty_Co[j]];
   BETA[3][j] = t_b[3][I_Type[j]] + c_b12[I_Type[j]][I_Ty_Co[j]];
   // Lambda linear model
   lambda[j] = logAF[j] + I_offset[I[j]] * I_scale[I_Component[j]] + ALPHA[j] 
               + BETA[1][j] * PRED[1][j] + BETA[2][j] * PRED[2][j]
               + BETA[3][j] * PRED[1][j] * PRED[2][j];
 }
 
 // Likelihood for MNE counts
 MNE ~ poisson_log(lambda);
 
}
generated quantities{
 
 // Log likelihoods
 vector[N] log_lik;
 // Across N observations: Estimated predictor 1 and
 // predictor 2 values. Some of these are combinations 
 // of estimated values.
 vector[N] PRED[2];
 // Vector of imputed values put back on predictor 1
 // and predictor 2 scales from the (0, 1) scale that 
 // was used to model values as beta distributed.
 vector[N_impute1] pred1sc;
 vector[N_impute2] pred2sc;
 
 // Put imputed (0, 1) values back on the predictor 1 
 // and predictor 2 scales.
 for (j in 1:N_impute1) {
  pred1sc[j] = pred1range[j] * val1[j] + pred1min[j];
 }
 for (j in 1:N_impute2) {
  pred2sc[j] = pred2range[j] * val2[j] + pred2min[j];
 }
 
 // LOOP FOR LINEAR MODEL
 for (j in 1:N) {
 
   // See comments in analogous 'if else' code chunk
   // in the model block for details about these lines.
   if (I_Element1[j] <= n_1) {
    PRED[1][j] = pred1sc[I_Rev1[j, 1]];
   } else if (I_Rev1[j, 2] == 0) {
    PRED[1][j] = pred1sc[I_Rev1[j, 1]];
   } else if (I_Rev1[j, 3] == 0) {
    PRED[1][j] = pred1sc[I_Rev1[j, 1]] * 0.5 + pred1sc[I_Rev1[j, 2]] * 0.5;
   } else if (I_Rev1[j, 4] == 0) {
    PRED[1][j] = pred1sc[I_Rev1[j, 1]] / 7 + pred1sc[I_Rev1[j, 2]] / 7
             + pred1sc[I_Rev1[j, 3]] * 5 / 7;
   } else if (I_Rev1[j, 5] > 0) {
    PRED[1][j] = pred1sc[I_Rev1[j, 1]] / 26 + pred1sc[I_Rev1[j, 2]] / 26
             + pred1sc[I_Rev1[j, 3]] * 5 / 26 + pred1sc[I_Rev1[j, 4]] * 7 / 13 
             + pred1sc[I_Rev1[j, 5]] * 5 / 26;
   }
   
   // See comments in analogous 'if else' code chunk
   // in the model block for details about these lines.
   if (I_Element2[j] <= n_2) {
    PRED[2][j] = pred2sc[I_Rev2[j, 1]];
   } else if (I_Rev2[j, 2] == 0) {
    PRED[2][j] = pred2sc[I_Rev2[j, 1]];
   } else if (I_Rev2[j, 3] == 0) {
    PRED[2][j] = pred2sc[I_Rev2[j, 1]] * 0.5 + pred2sc[I_Rev2[j, 2]] * 0.5;
   } else if (I_Rev2[j, 4] == 0) {
    PRED[2][j] = pred2sc[I_Rev2[j, 1]] / 7 + pred2sc[I_Rev2[j, 2]] / 7
             + pred2sc[I_Rev2[j, 3]] * 5 / 7;
   } else if (I_Rev2[j, 5] > 0) {
    PRED[2][j] = pred2sc[I_Rev2[j, 1]] / 26 + pred2sc[I_Rev2[j, 2]] / 26
             + pred2sc[I_Rev2[j, 3]] * 5 / 26 + pred2sc[I_Rev2[j, 4]] * 7 / 13 
             + pred2sc[I_Rev2[j, 5]] * 5 / 26;
   }
 
   // Generated log likelihood
   log_lik[j] = poisson_log_lpmf(MNE[j] | logAF[j]
               + I_offset[I[j]] * I_scale[I_Component[j]]
               + t_alpha[I_Type[j]] + c_alpha[I_Type[j]][I_Ty_Co[j]]
               + (t_b[1][I_Type[j]] + c_b1[I_Type[j]][I_Ty_Co[j]])
               * PRED[1][j]
               + (t_b[2][I_Type[j]] + c_b2[I_Type[j]][I_Ty_Co[j]])
               * PRED[2][j]
               + (t_b[3][I_Type[j]] + c_b12[I_Type[j]][I_Ty_Co[j]])
               * PRED[1][j] * PRED[2][j]);
 }

}
