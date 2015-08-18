#' @include unmarkedExtra-internal.R
NULL

#' Fit the MacKenzie et al. (2002) Occupancy Model
#'
#' This function fits the single season occupancy model of MacKenzie et al (2002). 
#'
#' @param formula Double right-hand side formula describing covariates of detection and occupancy in that order.
#' @param data An \code{\link[unmarked]{UnmarkedFrameOccu}} object
#' @param kwnownOcc Vector of sites that are known to be occupied. These should be supplied as row numbers of the y matrix, eg, c(3,8) if sites 3 and 8 were known to be occupied a priori.
#' @param method Method used to solve the model. Either "stan", or a method used by \code{\link[nloptr]{nloptr}}.
#' @param control List of options to specify model fitting procedure. See Details.
#' @export
#' @return \code{\link[unmarked]{unmarkedFitOccu}} object.
#' @author Jeffrey O. Hanson
occu<-function(formula, data, test.data=NULL, knownOcc=numeric(0), method='BFGS', control=list()) {
	# check data for valid inputs
	if (!is(data, "unmarkedFrameOccu")) 
		stop("data argument is not an unmarkedFrameOccu object.")
	if (!is(test.data, "unmarkedFrameOccu") && !is.null(test.data)) 
		stop("test argument is not an unmarkedFrameOccu object.")
		
	# solve model and return results
	if (method=='stan') {
		return(
			occu.stan(
				formula, data, test.data, knownOcc, method, control
			)
		)
	} else {
		return(
			occu.nloptr(
				formula, data, knownOcc, method, control
			)
		)
	}
}

occu.stan=function(formula, data, test.data, knownOcc, method, control) {
	## set default controls
	if (is.null(control$gp))
		control$gp=FALSE
	if (is.null(control$horseshoe))
		control$horseshoe=FALSE
	control$data=list(
		train.data=convertUMF(formula, data, knownOcc, control),
		test.data=convertUMF(formula, test.data, knownOcc, control)
	)
	
	## main processing	
	if (control$gp) {
		return(
			occu.stan.gp(control)
		)
	} else {
		if (is.null(test.data)) {
			control$data=control$data$train.data
			return(occu.stan.lin(control))
		} else {
			if (control$horseshoe) {
				return(occu.stan.test.horseshoe.lin(control))
			} else {
				return(occu.stan.test.ols.lin(control))
			}
		}
	}
}

occu.stan.test.horseshoe.lin=function(control) {
	## set priors
	if (is.null(control$priors)) {
		control$priors=list(
			dpars=Normal(0,10)
		)
	}
	
	# set pars
	names(control$data$train.data)=paste0(names(control$data$train.data), '_train')
	names(control$data$test.data)=paste0(names(control$data$test.data), '_test')
	control$data=append(control$data, control$data$train.data)
	control$data=append(control$data, control$data$test.data)
	control$data$nopars=control$data$nopars_train
	control$data$ndpars=control$data$ndpars_train
	control$data=control$data[-which(names(control$data) %in% c('train.data', 'test.data', 'nopars_train', 'nopars_test','ndpars_train', 'ndpars_test'))]
	
	## set parameters to return
	if (is.null(control$pars))
		control$pars=c(
			'log_lik',
			'dpars',
			'opars',
			'psi_test',
			'p_test',
			'sites_occupied_test',
			'number_sites_occupied_test',
			'fraction_sites_occupied_test'
		)
	
	## parse code
	control$model_code=paste0('
	data {
		int<lower=0> nobs_train; // number of total training observations
		int<lower=0> nsites_train; // number of training sites

		int<lower=0> nobs_test; // number of total test observations
		int<lower=0> nsites_test; // number of test sites

		int<lower=0> nopars; // number of parameters for observation model
		int<lower=0> ndpars; // number of parameters for detection model
		
		int<lower=0,upper=1> y_train[nobs_train]; // observations
		int<lower=1,upper=nobs_train> site_starts_train[nsites_train]; // index of first observation for i\'th site
		int<lower=1,upper=nobs_train> site_visits_train[nsites_train]; // number of observations for i\'th site
		
		matrix[nsites_train,nopars] X_train; // design matrix for observation model
		matrix[nobs_train,ndpars] V_train; // design matrix for detection model

		matrix[nsites_test,nopars] X_test; // design matrix for observation model
		matrix[nobs_test,ndpars] V_test; // design matrix for detection model
	}
	
	transformed data {
		// declare variables
		int<lower=0> site_detections_train[nsites_train]; // number of detections per training site		

//		vector[nopars] X_train_means;
//		vector[nopars] X_train_sds;
//		vector[ndpars] V_train_means;
//		vector[ndpars] V_train_sds;
//		matrix[nsites_train,nopars] X_train_std;
//		matrix[nobs_train,ndpars] V_train_std;
//		matrix[nsites_test,nopars] X_test_std;
//		matrix[nobs_test,ndpars] V_test_std;
		
		// calculate number of detections per training site
		for (i in 1:nsites_train)
			site_detections_train[i] <- sum(
				segment(
					y_train,
					site_starts_train[i],
					site_visits_train[i]
				)
			);
		
		/// standardise site-level covariates
		// first column is assumed to contain the intercept
//		for (i in 1:nsites_train) X_train_std[i,1] <- X_train[i,1];
//		for (i in 1:nsites_test) X_test_std[i,1] <- X_test[i,1];
//		X_train_means[1] <- 1;
//		X_train_sds[1] <- 1;
//		z-score remaining columns
//		if (nopars > 1) {
// 			for (i in 2:nopars) {
// 				// calculate means + sds
// 				X_train_means[i] <- mean(col(X_train, i));
// 				X_train_sds[i] <- sd(col(X_train, i));
// 				
//				// calculate z-scored values
// 				for (j in 1:nsites_train) X_train_std[j,i] <- (X_train[j,i] - X_train_means[i]) / X_train_sds[i];
// 				for (j in 1:nsites_test) X_test_std[j,i] <- (X_test[j,i] - X_train_means[i]) / X_train_sds[i];
// 			}
// 		}
//		
//		/// standardise observation-level covariates
//		// first column is assumed to contain the intercept
// 		for (i in 1:nobs_train) V_train_std[i,1] <- V_train[i,1];
// 		for (i in 1:nobs_test) V_test_std[i,1] <- V_test[i,1];
// 		V_train_means[1] <- 1;
// 		V_train_sds[1] <- 1;
//
//		// z-score remaining columns
// 		if (ndpars > 1) {
// 			for (i in 2:ndpars) {
// 				// calculate means + sds
// 				V_train_means[i] <- mean(col(V_train, i));
// 				V_train_sds[i] <- sd(col(V_train, i));
//
// 				// calculate z-scored values
// 				for (j in 1:nobs_train) V_train_std[j,i] <- (V_train[j,i] - V_train_means[i]) / V_train_sds[i];
// 				for (j in 1:nobs_test) V_test_std[j,i] <- (V_test[j,i] - V_train_means[i]) / V_train_sds[i];
// 			}
// 		}

 	}
	
	parameters {
		vector[ndpars] dpars;
		vector[nopars] ornorm;
		real<lower=0> otau;
		vector<lower=0>[nopars] olambda;
	}
	
	transformed parameters {
		vector[nopars] opars;
		vector[nsites_train] log1m_psi_train;
		vector[nsites_train] log_psi_train;
		vector[nobs_train] logit_p_train;
		
		{
			// declare internal variables
			vector[nsites_train] logit_psi_train;
			vector[nsites_train] psi_train;

			// calculate opars using matt trick
			opars <- ornorm // .* olambda * otau;
			
			logit_psi_train <- X_train * opars;
			logit_p_train <- V_train * dpars;
			
			for (i in 1:nsites_train) {
				psi_train[i] <- inv_logit(logit_psi_train[i]);
				log1m_psi_train[i] <- log1m(psi_train[i]);
			}
			log_psi_train <- log(psi_train);
		}
	}
	
	model {
		// priors
		dpars ~ ',repr(control$priors$dpars),';
		ornorm ~ normal(0, 1);
		olambda ~ cauchy(0, 1);
		otau ~ cauchy(0, 1);
		
		// likelihood
		for (i in 1:nsites_train) {
			if (site_detections_train[i] > 0)
				increment_log_prob(
					log_psi_train[i] + bernoulli_logit_log(
						segment(
							y_train, 
							site_starts_train[i],
							site_visits_train[i]
						),
						segment(
							logit_p_train,
							site_starts_train[i],
							site_visits_train[i]
						)
					)
				); 
			else 
				increment_log_prob(
					log_sum_exp(
						log_psi_train[i] + 
						bernoulli_logit_log(
							segment(
								y_train, 
								site_starts_train[i],
								site_visits_train[i]
							),
							segment(
								logit_p_train,
								site_starts_train[i],
								site_visits_train[i]
							)
						), 
						log1m_psi_train[i]
					)
				);
		}
	}
	
	generated quantities {
		// global variables
		real sites_occupied_test[nsites_test];
		real number_sites_occupied_test;
		real fraction_sites_occupied_test;

		vector[nsites_test] psi_test;
		vector[nobs_test] p_test;

		vector[nsites_train] log_lik;

		{
			// local variables
			vector[nsites_test] logit_psi_test;
			vector[nobs_test] logit_p_test;	
			
			// calculate psi_test and p_test
			logit_psi_test <- X_test * opars;
			logit_p_test <- V_test * dpars;
			
			for (i in 1:nsites_test) psi_test[i] <- inv_logit(logit_psi_test[i]);
			for (i in 1:nobs_test) p_test[i] <- inv_logit(logit_p_test[i]);
		
			// site-level summaries		
			for (i in 1:nsites_test) sites_occupied_test[i] <- bernoulli_rng(psi_test[i]);
			number_sites_occupied_test <- sum(sites_occupied_test);
			fraction_sites_occupied_test <- number_sites_occupied_test / nsites_test;
				
			// log-likelihood
			for (i in 1:nsites_train) {
				if (site_detections_train[i] > 0) 
					log_lik[i] <- log_psi_train[i] + bernoulli_logit_log(
							segment(
								y_train, 
								site_starts_train[i],
								site_visits_train[i]
							),
							segment(
								logit_p_train,
								site_starts_train[i],
								site_visits_train[i]
							)
						); 
				else 
					log_lik[i] <- log_sum_exp(
							log_psi_train[i] + 
							bernoulli_logit_log(
								segment(
									y_train, 
									site_starts_train[i],
									site_visits_train[i]
								),
								segment(
									logit_p_train,
									site_starts_train[i],
									site_visits_train[i]
								)
							), 
							log1m_psi_train[i]
						);
			}
		}
	}
	
	'	
	)

	if (is.null(options()$occu.stan.test.horseshoe.lin)) {
		options(occu.stan.test.horseshoe.lin=stan_model(model_code=control$model_code))
	}
	
	
	## run model
	return(
		do.call(
			sampling,
			append(
				list(object=options()$occu.stan.test.horseshoe.lin),
				control[!names(control) %in% c('gp','priors','model_code','horseshoe')]
			)
		)
	)
}



occu.stan.test.ols.lin=function(control) {
	## set priors
	if (is.null(control$priors)) {
		control$priors=list(
			opars=Normal(0,10),
			dpars=Normal(0,10)
		)
	}
	
	# set pars
	names(control$data$train.data)=paste0(names(control$data$train.data), '_train')
	names(control$data$test.data)=paste0(names(control$data$test.data), '_test')
	control$data=append(control$data, control$data$train.data)
	control$data=append(control$data, control$data$test.data)
	control$data$nopars=control$data$nopars_train
	control$data$ndpars=control$data$ndpars_train
	control$data=control$data[-which(names(control$data) %in% c('train.data', 'test.data', 'nopars_train', 'nopars_test','ndpars_train', 'ndpars_test'))]
	
	## set parameters to return
	if (is.null(control$pars))
		control$pars=c(
			'log_lik',
			'dpars',
			'opars',
			'psi_test',
			'p_test',
			'sites_occupied_test',
			'number_sites_occupied_test',
			'fraction_sites_occupied_test'
		)
	
	## parse code
	control$model_code=paste0('
	data {
		int<lower=0> nobs_train; // number of total training observations
		int<lower=0> nsites_train; // number of training sites

		int<lower=0> nobs_test; // number of total test observations
		int<lower=0> nsites_test; // number of test sites

		int<lower=0> nopars; // number of parameters for observation model
		int<lower=0> ndpars; // number of parameters for detection model
		
		int<lower=0,upper=1> y_train[nobs_train]; // observations
		int<lower=1,upper=nobs_train> site_starts_train[nsites_train]; // index of first observation for i\'th site
		int<lower=1,upper=nobs_train> site_visits_train[nsites_train]; // number of observations for i\'th site
		
		matrix[nsites_train,nopars] X_train; // design matrix for observation model
		matrix[nobs_train,ndpars] V_train; // design matrix for detection model

		matrix[nsites_test,nopars] X_test; // design matrix for observation model
		matrix[nobs_test,ndpars] V_test; // design matrix for detection model
	}
	
	transformed data {
		// declare variables
		int<lower=0> site_detections_train[nsites_train]; // number of detections per training site		
		vector[nopars] X_train_means;
		vector[nopars] X_train_sds;
		vector[ndpars] V_train_means;
		vector[ndpars] V_train_sds;
		matrix[nsites_train,nopars] X_train_std;
		matrix[nobs_train,ndpars] V_train_std;
		matrix[nsites_test,nopars] X_test_std;
		matrix[nobs_test,ndpars] V_test_std;
		
		// calculate number of detections per training site
		for (i in 1:nsites_train)
			site_detections_train[i] <- sum(
				segment(
					y_train,
					site_starts_train[i],
					site_visits_train[i]
				)
			);
		
		/// standardise site-level covariates
		// first column is assumed to contain the intercept
		for (i in 1:nsites_train) X_train_std[i,1] <- X_train[i,1];
		for (i in 1:nsites_test) X_test_std[i,1] <- X_test[i,1];
		X_train_means[1] <- 1;
		X_train_sds[1] <- 1;
		// z-score remaining columns
		if (nopars > 1) {
			for (i in 2:nopars) {
				// calculate means + sds
				X_train_means[i] <- mean(col(X_train, i));
				X_train_sds[i] <- sd(col(X_train, i));
				
				// calculate z-scored values
				for (j in 1:nsites_train) X_train_std[j,i] <- (X_train[j,i] - X_train_means[i]) / X_train_sds[i];
				for (j in 1:nsites_test) X_test_std[j,i] <- (X_test[j,i] - X_train_means[i]) / X_train_sds[i];
			}
		}
		
		/// standardise observation-level covariates
		// first column is assumed to contain the intercept
		for (i in 1:nobs_train) V_train_std[i,1] <- V_train[i,1];
		for (i in 1:nobs_test) V_test_std[i,1] <- V_test[i,1];
		V_train_means[1] <- 1;
		V_train_sds[1] <- 1;

		// z-score remaining columns
		if (ndpars > 1) {
			for (i in 2:ndpars) {
				// calculate means + sds
				V_train_means[i] <- mean(col(V_train, i));
				V_train_sds[i] <- sd(col(V_train, i));

				// calculate z-scored values
				for (j in 1:nobs_train) V_train_std[j,i] <- (V_train[j,i] - V_train_means[i]) / V_train_sds[i];
				for (j in 1:nobs_test) V_test_std[j,i] <- (V_test[j,i] - V_train_means[i]) / V_train_sds[i];
			}
		}
	}
	
	parameters {
		vector[nopars] opars;
		vector[ndpars] dpars;
	}
	
	transformed parameters {
		vector[nsites_train] log1m_psi_train;
		vector[nsites_train] log_psi_train;
		vector[nobs_train] logit_p_train;
		
		{
			// declare internal variables
			vector[nsites_train] logit_psi_train;
			vector[nsites_train] psi_train;
			
			logit_psi_train <- X_train_std * opars;
			logit_p_train <- V_train_std * dpars;
			
			for (i in 1:nsites_train) {
				psi_train[i] <- inv_logit(logit_psi_train[i]);
				log1m_psi_train[i] <- log1m(psi_train[i]);
			}
			log_psi_train <- log(psi_train);
		}
	}
	
	model {
		// priors
		dpars ~ ',repr(control$priors$dpars),';
		opars ~ ',repr(control$priors$opars),';
		
		// likelihood
		for (i in 1:nsites_train) {
			if (site_detections_train[i] > 0)
				increment_log_prob(
					log_psi_train[i] + bernoulli_logit_log(
						segment(
							y_train, 
							site_starts_train[i],
							site_visits_train[i]
						),
						segment(
							logit_p_train,
							site_starts_train[i],
							site_visits_train[i]
						)
					)
				); 
			else 
				increment_log_prob(
					log_sum_exp(
						log_psi_train[i] + 
						bernoulli_logit_log(
							segment(
								y_train, 
								site_starts_train[i],
								site_visits_train[i]
							),
							segment(
								logit_p_train,
								site_starts_train[i],
								site_visits_train[i]
							)
						), 
						log1m_psi_train[i]
					)
				);
		}
	}
	
	generated quantities {
		// global variables
		real sites_occupied_test[nsites_test];
		real number_sites_occupied_test;
		real fraction_sites_occupied_test;

		vector[nsites_test] psi_test;
		vector[nobs_test] p_test;

		vector[nsites_train] log_lik;

		{
			// local variables
			vector[nsites_test] logit_psi_test;
			vector[nobs_test] logit_p_test;	
			
			// calculate psi_test and p_test
			logit_psi_test <- X_test_std * opars;
			logit_p_test <- V_test_std * dpars;
			
			for (i in 1:nsites_test) psi_test[i] <- inv_logit(logit_psi_test[i]);
			for (i in 1:nobs_test) p_test[i] <- inv_logit(logit_p_test[i]);
		
			// site-level summaries		
			for (i in 1:nsites_test) sites_occupied_test[i] <- bernoulli_rng(psi_test[i]);
			number_sites_occupied_test <- sum(sites_occupied_test);
			fraction_sites_occupied_test <- number_sites_occupied_test / nsites_test;
				
			// log-likelihood
			for (i in 1:nsites_train) {
				if (site_detections_train[i] > 0) 
					log_lik[i] <- log_psi_train[i] + bernoulli_logit_log(
							segment(
								y_train, 
								site_starts_train[i],
								site_visits_train[i]
							),
							segment(
								logit_p_train,
								site_starts_train[i],
								site_visits_train[i]
							)
						); 
				else 
					log_lik[i] <- log_sum_exp(
							log_psi_train[i] + 
							bernoulli_logit_log(
								segment(
									y_train, 
									site_starts_train[i],
									site_visits_train[i]
								),
								segment(
									logit_p_train,
									site_starts_train[i],
									site_visits_train[i]
								)
							), 
							log1m_psi_train[i]
						);
			}
		}
	}
	
	'	
	)

	if (is.null(options()$occu.stan.test.ols.lin)) {
		options(occu.stan.test.ols.lin=stan_model(model_code=control$model_code))
	}
	
	
	## run model
	return(
		do.call(
			sampling,
			append(
				list(object=options()$occu.stan.test.ols.lin),
				control[!names(control) %in% c('gp','priors','model_code','horseshoe')]
			)
		)
	)
}



occu.stan.lin=function(control) {
	## set priors
	if (is.null(control$priors)) {
		control$priors=list(
			opars=Normal(0,10),
			dpars=Normal(0,10),
			psi_mean=Uniform(0,1),
			p_mean=Uniform(0,1)
		)
	}
	
	## set parameters to return
	control$pars=c('dpars', 'opars', 'psi', 'p', 'log_lik','sites_occupied', 'number_sites_occupied','fraction_sites_occupied')
	
	## parse code
	control$model_code=paste0('
	data {
		int<lower=0> nobs; // number of total observations
		int<lower=0> nsites; // number of sites
		
		int<lower=0> nopars; // number of parameters for observation model
		int<lower=0> ndpars; // number of parameters for detection model
		
		int<lower=0,upper=1> y[nobs]; // observations
		int<lower=1,upper=nobs> site_starts[nsites]; // index of first observation for i\'th site
		int<lower=1,upper=nobs> site_visits[nsites]; // number of observations for i\'th site
		
		matrix[nsites,nopars] X; // design matrix for observation model
		matrix[nobs,ndpars] V; // design matrix for detection model
	}
	
	transformed data {
		int<lower=0> site_detections[nsites]; // number of detections per site
		for (i in 1:nsites)
			site_detections[i] <- sum(
				segment(
					y,
					site_starts[i],
					site_visits[i]
				)
			);		
	}
	
	parameters {
		vector[ndpars] dpars;
		vector[nopars] opars;
	}
	
	transformed parameters {
		vector[nsites] logit_psi;
		vector[nobs] logit_p;
		vector[nsites] psi;
				
		vector[nsites] log1m_psi;
		vector[nsites] log_psi;
	
		logit_psi <- X * opars;
		logit_p <- V * dpars;
		
		for (i in 1:nsites) {
			psi[i] <- inv_logit(logit_psi[i]);
			log1m_psi[i] <- log1m(psi[i]);
		}
		log_psi <- log(psi);		
	}
	
	model {
		// priors
		dpars ~ ',repr(control$priors$dpars),';
		opars ~ ',repr(control$priors$opars),';

		// likelihood
		for (i in 1:nsites) {
			if (site_detections[i] > 0)
				increment_log_prob(
					log_psi[i] + bernoulli_logit_log(
						segment(
							y, 
							site_starts[i],
							site_visits[i]
						),
						segment(
							logit_p,
							site_starts[i],
							site_visits[i]
						)
					)
				); 
			else 
				increment_log_prob(
					log_sum_exp(
						log_psi[i] + 
						bernoulli_logit_log(
							segment(
								y, 
								site_starts[i],
								site_visits[i]
							),
							segment(
								logit_p,
								site_starts[i],
								site_visits[i]
							)
						), 
						log1m_psi[i]
					)
				);
		}
	}
	
	generated quantities {
		real sites_occupied[nsites];
		real number_sites_occupied;
		real fraction_sites_occupied;
		vector[nobs] p;
		vector[nsites] log_lik;

		// calculate p
		for (i in 1:nobs) 
			p[i] <- inv_logit(logit_p[i]);
				
		// site-level summaries
		for (i in 1:nsites) 
			sites_occupied[i] <- bernoulli_rng(psi[i]);		
		number_sites_occupied <- sum(sites_occupied);
		fraction_sites_occupied <- number_sites_occupied / nsites;
		
		// log-likelihood
		for (i in 1:nsites) {
			if (site_detections[i] > 0) 
				log_lik[i] <- log_psi[i] + bernoulli_logit_log(
						segment(
							y, 
							site_starts[i],
							site_visits[i]
						),
						segment(
							logit_p,
							site_starts[i],
							site_visits[i]
						)
					); 
			else 
				log_lik[i] <- log_sum_exp(
						log_psi[i] + 
						bernoulli_logit_log(
							segment(
								y, 
								site_starts[i],
								site_visits[i]
							),
							segment(
								logit_p,
								site_starts[i],
								site_visits[i]
							)
						), 
						log1m_psi[i]
					);
		}
	}
	
	'	
	)
	
	## run model
	return(
		do.call(
			stan,
			control[!names(control) %in% c('gp','priors')]
		)
	)
}

occu.stan.gp=function(control) {

	## remove intercept columns from design matrices
	control$data$X=control$data$X[,-1,drop=FALSE]
	control$data$V=control$data$V[,-1,drop=FALSE]
	control$data$nocovs=control$data$nopars-1
	control$data$ndcovs=control$data$ndpars-1
	

	## set priors
	if (is.null(control$priors)) {
		# init list
		control$priors=list()
		control$pars=c('psi','p','number_sites_occupied','fraction_sites_occupied')
		# set detection model priors
		if (ncol(control$data$V)==0) {
			control$priors$d_intercept=Normal(0,10)
			control$pars=c(control$pars, 'd_intercept')
		} else {
			control$priors$d_eta_sq=Cauchy(0,5)
			control$priors$d_sigma_sq=Cauchy(0,5)
			control$priors$d_inv_rho_sq=Cauchy(0,5)	
			control$pars=c(control$pars, 'd_eta_sq', 'd_sigma_sq', 'd_inv_rho_sq')
		}		
		# set occupancy model priors
		if (ncol(control$data$X)==0) {
			control$priors$o_intercept=Normal(0,10)
			control$pars=c(control$pars, 'o_intercept')
		} else {
			control$priors$o_eta_sq=Cauchy(0,5)
			control$priors$o_sigma_sq=Cauchy(0,5)
			control$priors$o_inv_rho_sq=Cauchy(0,5)
			control$pars=c(control$pars, 'o_eta_sq', 'o_sigma_sq', 'o_inv_rho_sq')
		}		
	}

	## parse code
	control$model_code=paste0('
	data {
		int<lower=0> nobs; // number of total observations
		int<lower=0> nocovs; // number of covariates for observation model
		int<lower=0> ndcovs; // number of covariates for detection model
		int<lower=0> nsites; // number of sites
				
		int<lower=0,upper=1> y[nobs]; // observations
		int<lower=1,upper=nobs> site_starts[nsites]; // index of first observation for i\'th site
		int<lower=1,upper=nobs> site_visits[nsites]; // index of last observation for i\'th site

		vector[nocovs] X[nsites]; // covariates for observation model
		vector[ndcovs] V[nobs]; // covariates for detection model
	}
	
	transformed data {
		int<lower=0> site_detections[nsites]; // number of detections per site
		
',
	if (ncol(control$data$V)>0) {
'		vector[nobs] d_mu;\n'
	},
	if (ncol(control$data$X)>0) {
'		vector[nsites] o_mu;\n'
	},
	if (ncol(control$data$V)>0) {
'		d_mu <- rep_vector(0,nobs);\n'
	},
	if (ncol(control$data$X)>0) {
'		o_mu <- rep_vector(0, nsites);\n'
	},	
'
		for (i in 1:nsites)
			site_detections[i] <- sum(segment(y, site_starts[i], site_visits[i]));
',
	'
	}
	
	parameters {
',
	if (ncol(control$data$V)==0) {
'		real d_intercept;'
	} else {
'
		real<lower=0> d_eta_sq;
		real<lower=0> d_inv_rho_sq[ndcovs];
		real<lower=0> d_sigma_sq;
		vector[nsites] d_rnorm;	
'	
	},
	if (ncol(control$data$X)==0) {
'		real o_intercept;'
	} else {
'
		real<lower=0> o_eta_sq;
		real<lower=0> o_inv_rho_sq[nocovs];
		real<lower=0> o_sigma_sq;
		vector[nsites] o_rnorm;
'
	},	
	'
	}
	
	transformed parameters{
		vector[nsites] psi;
		vector[nobs] p;
		
		
',
	if (ncol(control$data$V)>0) {
'
		vector<lower=0>[ndcovs] d_rho_sq;
		vector[ndcovs] d_exp_term;
		vector[nobs] logit_p;
		cov_matrix[nobs] d_Sigma;
		matrix[nsites,nsites] d_L;
'
	},
	if (ncol(control$data$X)>0) {
'
		vector<lower=0>[nocovs] o_rho_sq;
		vector[nocovs] o_exp_term;
		vector[nsites] logit_psi;
		cov_matrix[nsites] o_Sigma;
		matrix[nsites,nsites] o_L;
'
	},
	if (ncol(control$data$V)==0) {
'
		for (i in 1:nobs)
			p[i] <- inv_logit(d_intercept);
'	
	} else {
		'
		for (i in 1:ndcovs) {
			d_rho_sq[i] <- inv(d_inv_rho_sq[i]);
		}
		
		// off-diagonal elements
		for (i in 1:(nobs-1)) {
			for (j in (i+1):nobs) {
				d_exp_term <- d_rho_sq .* (V[i] - V[j]);
				d_Sigma[i,j] <- d_eta_sq * exp(-dot_self(d_exp_term));
				d_Sigma[j,i] <- d_Sigma[i,j];
			}
		}
		// diagonal elements
		for (i in 1:nobs)
			d_Sigma[i,i] <- d_eta_sq + d_sigma_sq;
		
		// implies:	logit_p ~ multi_normal(d_mu, d_Sigma);
		d_L <- cholesky_decompose(d_Sigma);
		logit_p <- d_mu + d_L * d_rnorm;
		
		for (i in 1:nobs)
			p[i] <- inv_logit(d_intercept);
'
	},
	if (ncol(control$data$X)==0) {
'			
		for (i in 1:nsites)
			psi[i] <- inv_logit(o_intercept);
'
	} else {
'
		for (i in 1:nocovs) {
			o_rho_sq[i] <- inv(o_inv_rho_sq[i]);
		}
		
		// off-diagonal elements
		for (i in 1:(nsites-1)) {
			for (j in (i+1):nsites) {			
				o_exp_term <- o_rho_sq .* (X[i] - X[j]);
				o_Sigma[i,j] <- o_eta_sq * exp(-dot_self(o_exp_term));
				o_Sigma[j,i] <- o_Sigma[i,j];
			}
		}
		// diagonal elements
		for (i in 1:nsites)
			o_Sigma[i,i] <- o_eta_sq + o_sigma_sq;
		
		// implies:	logit_psi ~ multi_normal(o_mu, o_Sigma);
		o_L <- cholesky_decompose(o_Sigma);
		logit_psi <- o_mu + o_L * o_rnorm;
		
		for (i in 1:nsites)
			psi[i] <- inv_logit(logit_psi[i]);
'
	},
'		
	}
	
	model {
		// local variables
		vector[nsites] log_psi;
		vector[nsites] log1m_psi;	


		//  cache log psi calculations
		log_psi <- log(psi);
		for (i in 1:nsites)
			log1m_psi[i]<-log1m(psi[i]);

		// priors
',
	if (ncol(control$data$V)==0) {
paste0('		d_intercept ~ ',repr(control$priors$d_intercept), ';\n')
	} else {
		paste0(
'		d_eta_sq ~ ', repr(control$priors$d_eta_sq), ';\n',
'		d_inv_rho_sq ~ ', repr(control$priors$d_inv_rho_sq), ';\n',
'		d_sigma_sq ~ ', repr(control$priors$d_sigma_sq), ';\n',
'		d_rnorm ~ normal(0,1);\n'

		)
	},
	if (ncol(control$data$X)==0) {
paste0('		o_intercept ~ ',repr(control$priors$o_intercept), ';\n')
	} else {
		paste0(
'		o_eta_sq ~ ', repr(control$priors$o_eta_sq), ';\n',
'		o_inv_rho_sq ~ ', repr(control$priors$o_inv_rho_sq), ';\n',
'		o_sigma_sq ~ ', repr(control$priors$o_sigma_sq), ';\n',
'		o_rnorm ~ normal(0,1);\n'
		)
	},	
	'
	
		// likelihood
		for (i in 1:nsites) {
			if (site_detections[i] > 0) 
				increment_log_prob(
					log_psi[i] + bernoulli_log(
						segment(
							y, 
							site_starts[i],
							site_visits[i]
						),
						segment(
							p,
							site_starts[i],
							site_visits[i]
						)
					)
				); 
			else 
				increment_log_prob(
					log_sum_exp(
						log_psi[i] + 
						bernoulli_log(
							segment(
								y, 
								site_starts[i],
								site_visits[i]
							),
							segment(
								p,
								site_starts[i],
								site_visits[i]
							)
						), 
						log1m_psi[i]
					)
				);
		}
	}
	
	generated quantities {
		real sites_occupied[nsites];
		real nsites2;
		real number_sites_occupied;
		real fraction_sites_occupied;
		
		nsites2 <- nsites;
		for (i in 1:nsites)
			sites_occupied[i] <- bernoulli_rng(psi[i]);
		number_sites_occupied <- sum(sites_occupied);
		fraction_sites_occupied <- number_sites_occupied / nsites;
	}


	'	
	)
		
	# debugging
	assign('control', control, envir=globalenv())
	
	
	## run model
	return(
		do.call(
			stan,
			control
		)
	)

}

occu.nloptr=function(formula, data, knownOcc, method, control) {
	# set default controls
	if (is.null(control$engine))
		control$engine<-'C'
	if (!is.null(control$se))
		control$hessian<-control$se
	if (is.null(control$se) && is.null(control$hessian))
		control$hessian<-TRUE
	
	# validate inputs
	control$engine<-match.arg(control$engine, c("C", "R"))
	
	# preliminary processing
	designMats<-unmarked:::getDesign(data, formula)
	X<-designMats$X # design matrix for observation model
	V<-designMats$V # design matrix for detection model
	y<-unmarked:::truncateToBinary(designMats$y)
	removed<-designMats$removed.sites
	X.offset<-designMats$X.offset
	V.offset<-designMats$V.offset
	if (is.null(X.offset))
		X.offset<-rep(0, nrow(X))
	if (is.null(V.offset))
		V.offset<-rep(0, nrow(V))
	J <- ncol(y) # maximum number of repeat visits
	M <- nrow(y) # number of sites
	knownOccLog <- rep(FALSE, numSites(data))
	knownOccLog[knownOcc] <- TRUE
	if (length(removed) > 0) 
		knownOccLog <- knownOccLog[-removed]
	occParms <- colnames(X) # names of parameters for occupancy model
	detParms <- colnames(V) # names of parameters for detection model
	nOP <- ncol(X) # number of parameters for occupancy model (inc. intercept)
	nDP <- ncol(V) # number of parameters for detection model (inc. intercept)
	nP <- nDP + nOP # total number of parameters
	if (!is.null(control$starts) && length(control$starts) != nP) 
		stop(paste("The number of starting values should be", nP))
	# convert y matrix to vector and replace missing values with -1
	yvec <- as.numeric(t(y))
	navec <- is.na(yvec) # logical vector which values are NA
	nd <- ifelse(rowSums(y, na.rm = TRUE) == 0, 1, 0) # number of sites where spp not detected
	# main processing
	if (identical(control$engine, "C")) {
		control$eval_f <- function(params) {
			beta.psi <- params[1:nOP]
			beta.p <- params[(nOP + 1):nP]
			.Call("nll_occu", yvec, X, V, beta.psi, beta.p, nd, knownOccLog, navec, X.offset, V.offset, PACKAGE = "unmarked")
		}
	} else {
		control$eval_f <- function(params) {
			psi <- plogis(X %*% params[1:nOP] + X.offset)
			psi[knownOccLog] <- 1
			pvec <- plogis(V %*% params[(nOP + 1):nP] + V.offset)
			cp <- (pvec^yvec) * ((1 - pvec)^(1 - yvec))
			cp[navec] <- 1
			cpmat <- matrix(cp, M, J, byrow = TRUE)
			loglik <- log(rowProds(cpmat) * psi + nd * (1 - psi))
			-sum(loglik)
		}
	}
	if (is.null(control$starts))
		control$starts <- rep(0, nP)
	control$opts=list(algorithm=method, xtol_rel=1.0e-10, maxeval=100000)
	control$x0=control$starts
	fm <- do.call(
		nloptr,
		control[!names(control) %in% c('engine', 'hessian', 'starts')]
	)
		
	if (control$hessian) {
	
		assign('fm', fm, envir=globalenv())
		assign('control', control, envir=globalenv())
		
		
		tryCatch(covMat <- solve(optimHess(fm$solution, control$eval_f)), error = function(x) stop(simpleError("Hessian is singular.  Try providing starting values or using fewer covariates.")))
	} else {
		covMat <- matrix(NA, nP, nP)
	}
	ests <- fm$solution
	fmAIC <- 2 * fm$objective + 2 * nP
	names(ests) <- c(occParms, detParms)
	state <- unmarked:::unmarkedEstimate(
		name = "Occupancy",
		short.name = "psi", 
		estimates = ests[1:nOP],
		covMat = as.matrix(covMat[1:nOP, 1:nOP]), 
		invlink = "logistic", 
		invlinkGrad = "logistic.grad"
	)
	det <- unmarked:::unmarkedEstimate(
		name = "Detection", 
		short.name = "p", 
		estimates = ests[(nOP + 1):nP],
		covMat = as.matrix(covMat[(nOP + 1):nP, (nOP + 1):nP]),
		invlink = "logistic",
		invlinkGrad = "logistic.grad"
	)
	estimateList <- unmarked:::unmarkedEstimateList(list(state = state, det = det))
	umfit <- new(
		"unmarkedFitOccu",
		fitType = "occu",
		call = match.call(), 
		formula = formula,
		data = data,
		sitesRemoved = designMats$removed.sites, 
		estimates = estimateList,
		AIC = fmAIC,
		opt = list(nloptr=fm, convergence=fm$status),
		negLogLike = fm$objective, 
		nllFun = control$eval_f,
		knownOcc = knownOccLog
	)
	return(umfit)
}





