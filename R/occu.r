#' @include unmarkedExtra-internal.R
NULL

#' Fit the MacKenzie et al. (2002) Occupancy Model
#'
#' This function fits the single season occupancy model of MacKenzie et al (2002). 
#'
#' @param formula Double right-hand side formula describing covariates of detection and occupancy in that order.
#' @param data An \code{\link[unmarked]{UnmarkedFrameOccu}} object
#' @param kwnownOcc Vector of sites that are known to be occupied. These should be supplied as row numbers of the y matrix, eg, c(3,8) if sites 3 and 8 were known to be occupied a priori.
#' @param method Method used to solve the model. Either "stan", or a method used by \code{\link[optimx]{optimx}}.
#' @param control List of options to specify model fitting procedure. See Details.
#' @export
#' @return \code{\link[unmarked]{unmarkedFitOccu}} object.
#' @author Jeffrey O. Hanson
occu<-function(formula, data, knownOcc=numeric(0), method='BFGS', control=list()) {
	# check data for valid inputs
	if (!is(data, "unmarkedFrameOccu")) 
		stop("Data is not an unmarkedFrameOccu object.")
	# solve model and return results
	if (method=='stan') {
		return(
			occu.stan(
				formula, data, knownOcc, method, control
			)
		)
	} else {
		return(
			occu.optimx(
				formula, data, knownOcc, method, control
			)
		)
	}
}

occu.stan=function(formula, data, knownOcc, method, control) {
	## set default controls
	if (is.null(control$gp))
		control$gp=FALSE
	control$model_name='occupancy-detection model'
	control$data=list()
	## preliminary processing
	# prepare data
	control$data$y<-c(t(unmarked:::truncateToBinary(data@y)))
	if (identical(as.formula(paste("~", formula[3], sep="")), ~ 1)) {
		control$data$X <- model.matrix(as.formula(paste("~", formula[3], sep="")), data.frame(x=seq_len(nrow(data@y))))
	} else {
		control$data$X<-model.matrix(as.formula(paste("~", formula[3], sep="")), data@siteCovs)
	}
	if (identical(as.formula(formula[[2]]), ~ 1)) {
		control$data$V <- model.matrix(as.formula(formula[[2]]), data.frame(x=seq_len(length(data@y))))
	} else {
		control$data$V<-model.matrix(as.formula(formula[[2]]), data@obsCovs)
	}
	# handle NA values
	if (any(is.na(control$data$y))) {
		# remove sites which have never been visited
		is.dets.missing<-rowSums(!is.na(data@y))
		if (any(is.dets.missing==0)) {
			control$data$X<-control$X[which(is.dets.missing>0),]
		}
		# remove missing observations from data
		control$data$V<-control$data$V[which(!is.na(control$data$y)),,drop=FALSE]
		control$data$y<-control$data$y[which(!is.na(control$data$y))]
	}
	# calculate numbers for stan
	control$data$nobs=length(control$data$y) # number observations total
	control$data$nsites=nrow(control$data$X) # number of sites
	obssites=rep(seq_len(nrow(data@y)), each=obsNum(data))[which(!is.na(c(t(data@y))))] # number observations per site
	control$data$site_starts=as.integer(unname(tapply(X=seq_along(obssites), INDEX=obssites, FUN=min)))
	control$data$site_visits=as.vector(table(obssites))
	control$data$nopars=ncol(control$data$X) # number observation parameters
	control$data$ndpars=ncol(control$data$V) # number detection parameters
	## main processing	
	if (control$gp) {
		return(
			occu.stan.gp(control)
		)
	} else {
		return(
			occu.stan.lin(control)
		)
	}
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
		int<lower=0> nopars; // number of parameters for observation model
		int<lower=0> ndpars; // number of parameters for detection model
		int<lower=0> nsites; // number of sites
		
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
	
	assign('control', control, envir=globalenv())
		
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

occu.optimx=function(formula, data, knownOcc, method, control) {
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
		control$fn <- function(params, ...) {
			beta.psi <- params[1:nOP]
			beta.p <- params[(nOP + 1):nP]
			.Call("nll_occu", yvec, X, V, beta.psi, beta.p, nd, knownOccLog, navec, X.offset, V.offset, PACKAGE = "unmarked")
		}
	} else {
		control$fn <- function(params, ...) {
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
	control$method=method
	control$par=control$starts
	fm <- do.call(
		optimx,
		control
	)
	if (control$hessian) {
		tryCatch(covMat <- solve(fm$hessian), error = function(x) stop(simpleError("Hessian is singular.  Try providing starting values or using fewer covariates.")))
	} else {
		covMat <- matrix(NA, nP, nP)
	}
	ests <- fm$par
	fmAIC <- 2 * fm$value + 2 * nP
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
		opt = fm,
		negLogLike = fm$value, 
		nllFun = control$fn,
		knownOcc = knownOccLog
	)
	return(umfit)
}





