// This code fits a model for bison element frequencies
// using one predictor variable.

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
 // Index of element types across observations
 int I_Element[N];
 
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
 
 // Number of total predictor values (estimated + 
 // observed) across all observations for predictor. 
 // This can exceed N, because some some predictor 
 // values are combinations of values (e.g.,
 // metapodial = (metatarsal + metacarpal)/2).
 int<lower=1> N_;
 // Number unique element types for predictor
 int<lower=1> n_;
 // Vector of predictor values
 real<lower=0> pred1[N_];  
 // Index of element types in predictor data
 int I_[N_]; 
 // Range of predictor values per element type
 real<lower=0> pred1range[N_];
 // Minimum predictor value per element type
 real<lower=0> pred1min[N_];
 
 ///// Predictor imputation variables ////
 ///// (i.e., N_ - observed) ////
 // Number of imputed predictor values
 int<lower=1> N_impute;
 // Vector of positions for values to impute
 // within pred1
 int<lower=1> i_impute[N_impute];
 // Maximum number of pred1 values in a combination
 // for a single predictor value
 int I_cols;
 // Matrix that specifies indices in pred1 for combi-
 // nations of predictor values across N observations
 int I_Rev[N, I_cols];

}

parameters{
  
 // Global intercept parameters, one for
 // each site type.
 vector[n_t] t_alpha;
 // Beta parameter for predictor, one for
 // each site type
 vector[n_t] t_b;
 
 // Observation level overdispersion parameters,
 // one value for each observation.
 vector[N] I_offset;   
 
 // Z-scores for component-level varying effects,
 // c_alpha and c_b. There is one z-score matrix
 // for each site type.
 matrix[2, n_c] z_c[n_t];
 cholesky_factor_corr[2] LRho_c[n_t];
 // Scaling parameters for z_c.
 vector<lower=0>[2] scale_c[n_t];
 
 // Mean location parameter for element type
 // values in predictor.
 real mu;
 // SD parameter for element type values in
 // predictor.
 real<lower=0> sigma;
 // Shape parameter for distribution of predictor
 // values; beta distributed.
 real<lower=0> shape;
 // Mu offsets for each mean element type value within
 // predictor. These are in z-score units, as they
 // are scaled by sigma.
 vector[n_] loc;
 
 // Predictor is modelled as beta distribution:
 // beta(LOC*shape, (1-LOC)*shape)
 // LOC is on the (0, 1) scale, but the parameter loc
 // specified here is unbounded, being inv_logit
 // tranformed in the linear model to obtain LOC.
 
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
 
 // Imputed predictor values
 vector<lower=0, upper=1>[N_impute] val;
 
}
transformed parameters{
 
 // Initialize component level parameters
 // for alpha and b parameters. 
 matrix[n_c, 2] v_c[n_t];
 vector[n_c] c_alpha[n_t];
 vector[n_c] c_b[n_t];
 matrix[2,2] Rho_c[n_t];
 
 // Cholesky decomposition for a matrix of component
 // level parameters for each site type (the loop
 // traverses matrixes for each site type).
 for (j in 1:n_t) {
   
  v_c[j][1:t_N[j], ] = (diag_pre_multiply(scale_c[j], 
                        LRho_c[j]) * z_c[j][ , 1:t_N[j]])';

  // Since the number of components is uneven between site
  // types, there are unfilled bins in each matrix. Fill these
  // bins with zeros. 
  v_c[j][(t_N[j] + 1):n_c, ] = rep_matrix(0, n_c - t_N[j], 2);     
  
  c_alpha[j] = col(v_c[j], 1);
  c_b[j] = col(v_c[j], 2);
  
  Rho_c[j] = LRho_c[j] * LRho_c[j]';
 
 }
 
}
model{
 
 // Vector of merged imputed + observed values for
 // predictor.
 vector[N_] merged;
 // Vector of merged values put back on predictor
 // scale from the (0, 1) scale that was used to model
 // values on the beta scale.
 vector[N_] pred1sc;
 
 // Lambdas for Poisson outcome
 vector[N] lambda; 
 // Predictor beta distribution location values, 
 // (0, 1) scale.
 vector[N_] LOC;
 // Linear combinations of alpha parameters
 vector[N] ALPHA;    
 // Linear combinations of each beta parameter
 vector[N] BETA;    
 // Across N observations: Estimated predictor
 // values. Some of these are combinations of
 // estimated values.
 vector[N] PRED1;
 
 // Prior distributions for global alpha and beta
 // parameters for each site type.
 t_alpha ~ normal(0, 1);
 t_b ~ normal(0, 1);
 
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

 // Priors and hyper priors for predictor values. 
 shape ~ exponential(2);
 mu ~ normal(0, 0.5);
 sigma ~ exponential(2);
 loc ~ normal(0, 1);
 
 // Merge imputed and observed predictor values.
 for (k in 1:N_) {
  if (k > N_impute) {
   merged[k] = pred1[k];
  }else{
   merged[k] = val[k];
  }
  // Put merged (0, 1) values back on predictor scale
  pred1sc[k] = pred1range[k] * merged[k] + pred1min[k];
 }
 
 // Model for predictor values, where beta
 // location parameter varies by element type.
 for (k in 1:N_) {
  LOC[k] = inv_logit(mu + loc[I_[k]] * sigma);
 }
 
 // Likelihood for merged predictor values
 merged ~ beta(LOC * (shape + 2), 
              (1 - LOC) * (shape + 2));
 
 // LOOP FOR LINEAR MODEL
 for (j in 1:N) {

   // OBSERVATION j PREDICTOR VALUE:
   // If the element index indicates that observation j
   // corresponds to an element associated with a predictor
   // value that was directly imputed, obtain the imputed
   // value that has been rescaled from (0, 1) back to the
   // predictor scale. Otherwise, observation j corresponds
   // to an element that is some combination of predictor 
   // values. In this case, combine the necessary rescaled
   // predictor values. 
   if (I_Element[j] <= n_) {
    PRED1[j] = pred1sc[I_Rev[j, 1]];

   } else if (I_Rev[j, 2] == 0) {
     // Element type corresponds to another already specified
     // element type, e.g., "pelvis" = "sacrum"
    PRED1[j] = pred1sc[I_Rev[j, 1]];

   } else if (I_Rev[j, 3] == 0) {
     // Various element types that are aggregates ot two element
     // types, e.g., "metapodial"" = ("metacarpal" + "metatarsal")/2.
    PRED1[j] = pred1sc[I_Rev[j, 1]] * 0.5 + pred1sc[I_Rev[j, 2]] * 0.5;

   } else if (I_Rev[j, 4] == 0) {
     // "cervical 1-7 vert"
    PRED1[j] = pred1sc[I_Rev[j, 1]]/7 + pred1sc[I_Rev[j, 2]]/7
             + pred1sc[I_Rev[j, 3]] * 5/7;

   } else if (I_Rev[j, 5] > 0) {
    // "unspecified vert"
    PRED1[j] = pred1sc[I_Rev[j, 1]]/26 + pred1sc[I_Rev[j, 2]]/26 
             + pred1sc[I_Rev[j, 3]] * 5/26 + pred1sc[I_Rev[j, 4]] * 7/13 
             + pred1sc[I_Rev[j, 5]] * 5/26;

   }
   
   // Intercept linear model
   ALPHA[j] = t_alpha[I_Type[j]] + c_alpha[I_Type[j]][I_Ty_Co[j]];
   // Predictor effects linear beta model
   BETA[j] = t_b[I_Type[j]] + c_b[I_Type[j]][I_Ty_Co[j]];
   // Lambda linear model
   lambda[j] = logAF[j] + I_offset[I[j]] * I_scale[I_Component[j]] 
               + ALPHA[j] + BETA[j] * PRED1[j];
 }
 
 // Likelihood for MNE counts
 MNE ~ poisson_log(lambda);
 
}
generated quantities{
 
 // Log likelihoods
 vector[N] log_lik;
 // Across N observations: Estimated predictor
 // values. Some of these are combinations of 
 // estimated values.
 vector[N] PRED1;
 // Vector of imputed values put back on predictor
 // scale from the (0, 1) scale that was used to model
 // values as beta distributed.
 vector[N_impute] pred1sc;
 
 // Put imputed (0, 1) values back on the predictor scale.
 for (j in 1:N_impute) {
  pred1sc[j] = pred1range[j] * val[j] + pred1min[j];
 }
 
 // LOOP FOR LINEAR MODEL
 for (j in 1:N) {
 
   // See comments in analogous 'if else' code chunk
   // in the model block for details about these lines.
   if (I_Element[j] <= n_) {
    PRED1[j] = pred1sc[I_Rev[j, 1]];
   } else if (I_Rev[j, 2] == 0) {
    PRED1[j] = pred1sc[I_Rev[j, 1]];
   } else if (I_Rev[j, 3] == 0) {
    PRED1[j] = pred1sc[I_Rev[j, 1]] * 0.5 + pred1sc[I_Rev[j, 2]] * 0.5;
   } else if (I_Rev[j, 4] == 0) {
    PRED1[j] = pred1sc[I_Rev[j, 1]] / 7 + pred1sc[I_Rev[j, 2]] / 7
             + pred1sc[I_Rev[j, 3]] * 5 / 7;
   } else if (I_Rev[j, 5] > 0) {
    PRED1[j] = pred1sc[I_Rev[j, 1]] / 26 + pred1sc[I_Rev[j, 2]] / 26
             + pred1sc[I_Rev[j, 3]] * 5 / 26 + pred1sc[I_Rev[j, 4]] * 7 / 13 
             + pred1sc[I_Rev[j, 5]] * 5 / 26;
   }
   
   // Generated log likelihood
   log_lik[j] = poisson_log_lpmf(MNE[j] | logAF[j] + t_alpha[I_Type[j]]
                                 + c_alpha[I_Type[j]][I_Ty_Co[j]]
                                 + I_offset[I[j]] * I_scale[I_Component[j]] 
                                 + (t_b[I_Type[j]] + c_b[I_Type[j]][I_Ty_Co[j]])
                                 * PRED1[j]);
 }

}
