#' @include unmarkedExtra-internal.R
NULL

#' Fit the MacKenzie et al. (2002) Occupancy Model
#'
#' This function fits the single season occupancy model of MacKenzie et al (2002). 
#'
#' @param formula Double right-hand side formula describing covariates of detection and occupancy in that order.
#' @param data An \link[unmarked][UnmarkedFrameOccu] object
#' @param kwnownOcc Vector of sites that are known to be occupied. These should be supplied as row numbers of the y matrix, eg, c(3,8) if sites 3 and 8 were known to be occupied a priori.
#' @param method Method used to solve the model. Either "stan", or a method used by \link[optimx][optimx].
#' @param control List of options to specify model fitting procedure. See Details.
#' @export
#' @return \link[unmarkedFitOccu] object.
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
	control$data$list=list()
	
	## preliminary processing
	# prepare data
	control$data$y<-c(t(unmarked:::truncateToBinary(data@y)))
	control$data$X<-model.matrix(as.formula(paste("~", formula[3], sep="")), data@siteCovs)
	control$data$V<-model.matrix(as.formula(formula[[2]]), data@obsCovs)

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
	control$data$site_starts=unname(tapply(X=seq_along(obssites), INDEX=obssites, FUN=min))
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
	control$pars=c('dpars', 'opars', 'psi', 'p')
	
	## parse code
	control$model_code=paste0('
	data {
		int<lower=0> nobs; // number of total observations
		int<lower=0> nopars; // number of parameters for observation model
		int<lower=0> ndpars; // number of parameters for detection model
		int<lower=0> nsites; // number of sites
		
		int<lower=0,upper=1> y[nobs]; // observations
		int<lower=1,upper=nobs> site_starts[nsites]; // index of first observation for i\'th site
		int<lower=1,upper=nobs> site_visits[nsites]; // index of last observation for i\'th site
		
		matrix[nsites,nopars] X; // design matrix for observation model
		matrix[nobs,ndpars] V; // design matrix for detection model
	}
	
	transformed data {
		int<lower=0> site_detections[nsites]; // number of detections per site
		for (i in 1:nsites)
			site_detections[i] <- sum(segment(y, site_starts[i], site_visits[i]));		

	}
	
	parameters {
		vector[ndpars] dpars;
		vector[nopars] opars;
	}
	
	transformed parameters {
		vector[nsites] psi;
		vector[nobs] p;
		
		psi <- X * opars;
		for (i in 1:nsites)
			psi[i] <- inv_logit(psi[i]);
				
		p <- V * dpars;
		for (i in 1:nobs)
			p[i] <- inv_logit(p[i]);
			
	}
	
	model {
		// local variables
		vector[nsites] log_psi;
		vector[nsites] log1m_psi;
		
		log_psi <- log(psi);
		for (i in 1:nsites)
			log1m_psi[i]<-log1m(psi[i]);
		
		// priors
		dpars ~ ',repr(control$priors$dpars),';
		opars ~ ',repr(control$priors$opars),';

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
			
	'	
	)
		
	## run model
	return(
		do.call(
			stan,
			control
		)
	)
}

occu.stan.gp=function(control) {
	## set priors
	if (is.null(control$priors)) {
		# init list
		control$priors=list()
		control$pars=c()
		# set detection model priors
		if (ncol(control$data$V)==1) {
			control$priors$d_intercept=Normal(0,10)
			control$pars=c(control$pars, 'd_intercept')
		} else {
			control$priors$d_eta_sq=Cauchy(0,5)
			control$priors$d_sigma_sq=Cauchy(0,5)
			control$priors$d_inv_rho_sq=Cauchy(0,5)
			control$pars=c(control$pars, 'd_eta_sq', 'd_sigma_sq', 'd_inv_rho_sq')
		}		
		# set occupancy model priors
		if (ncol(control$data$X)==1) {
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
		int<lower=0> nopars; // number of parameters for observation model
		int<lower=0> ndpars; // number of parameters for detection model
		int<lower=0> nsites; // number of sites
				
		int<lower=0,upper=1> y[nobs]; // observations
		int<lower=1,upper=nobs> site_starts[nsites]; // index of first observation for i\'th site
		int<lower=1,upper=nobs> site_visits[nsites]; // index of last observation for i\'th site
		
		matrix[nsites,nopars] X; // design matrix for observation model
		matrix[nobs,ndpars] V; // design matrix for detection model
	}
	
	transformed data {
		int<lower=0> site_detections[nsites]; // number of detections per site
		int<lower=0> ndcovs;
		int<lower=0> nocovs;
		
',
	if (ncol(control$data$V)>1) {
'		vector[nobs] d_mu;'
	},
	if (ncol(control$data$X)>1) {
'		vector[nobs] o_mu;'
	},'
		ndcovs <- ndpars-1;
		nocovs <- nopars-1;
		',
	if (ncol(control$data$V)>1) {
'		for (i in 1:nobs) d_mu <- 0;\n'
	},
		'
		for (i in 1:nsites)  {
			site_detections[i] <- sum(segment(y, site_starts[i], site_visits[i]));
		',
	if (ncol(control$data$V)>1) {
'			o_mu <- 0;\n'
	},
	'
		}
	}
	
	parameters {
	',
	if (ncol(control$data$V)==1) {
	'	real d_intercept;'
	} else {
		'
		real<lower=0> d_eta_sq;
		real<lower=0> d_inv_rho_sq[ndcovs];
		real<lower=0> d_sigma_sq;
		vector[nobs] logit_p;
		'
	},
	if (ncol(control$data$X)==1) {
'		real o_intercept;'
	} else {
		'
		real<lower=0> o_eta_sq;
		real<lower=0> o_inv_rho_sq[nocovs];
		real<lower=0> o_sigma_sq;
		vector[nsites] logit_psi;
		'
	},	
	'
	}
	
	transformed parameters{
',
	if (ncol(control$data$V)>1) {
'		real<lower=0> d_rho_sq[ndpars];\n'
	},
	if (ncol(control$data$X)>1) {
'		real<lower=0> o_rho_sq[nopars];\n'
	},
	if (ncol(control$data$V)>1) {
		'
		for (i in 1:ndcovs) {
			d_rho_sq[i] <- inv(d_inv_rho_sq[i]);
		}
		'
	},
	if (ncol(control$data$X)>1) {
		'
//		for (i in 1:nocovs) {
//			o_rho_sq[i] <- inv(o_inv_rho_sq[i]);
//		}		
		o_rho_sq[1] <- inv(o_inv_rho_sq[1]);
		print(o_rho_sq[1]);
'
	},
	'	
	}
	
	model {
		// local variables
		vector[nsites] psi;
		vector[nobs] p;		
		vector[nsites] log_psi;
		vector[nsites] log1m_psi;	
	',
	if (ncol(control$data$V)>1) {
'		matrix[nobs,nobs] d_Sigma;\n'
	},
	if (ncol(control$data$X)>1) {
'	matrix[nsites,nsites] o_Sigma;\n'
	},
	'		print (1234);		// priors\n',
	if (ncol(control$data$V)==1) {
paste0('		d_intercept ~ ',repr(control$priors$d_intercept), ';\n')
	} else {
		paste0(
'		d_eta_sq ~ ', repr(control$priors$d_eta_sq), ';\n',
'		d_inv_rho_sq ~ ', repr(control$priors$d_inv_rho_sq), ';\n',
'		d_sigma_sq ~ ', repr(control$priors$d_sigma_sq), ';\n'
		)
	},

	if (ncol(control$data$X)==1) {
paste0('		o_intercept ~ ',repr(control$priors$o_intercept), ';\n')
	} else {
		paste0(
'		o_eta_sq ~ ', repr(control$priors$o_eta_sq), ';\n',
'		o_inv_rho_sq ~ ', repr(control$priors$o_inv_rho_sq), ';\n',
'		o_sigma_sq ~ ', repr(control$priors$o_sigma_sq), ';\n'
		)
	},
	'		/// calculate p\n',
	if (ncol(control$data$V)==1) {
'		for (i in 1:nobs)
			p[i] <- inv_logit(d_intercept);
'
	} else {
'
		// off-diagonal elements
		for (i in 1:(nobs-1)) {
			for (j in (i+1):nobs) {
				d_Sigma[i,j] <- 0;
				for (k in 2:ndpars) {
					d_Sigma[i,j] <- d_Sigma[i,j] + d_rho_sq[k-1] * pow(V[i,k] - V[V,k], 2);
				}
				d_Sigma[i,j] <- d_eta_sq * exp(d_Sigma[i,j]);
				d_Sigma[j,i] <- d_Sigma[i,j];
			}
		}
		// diagonal elements
		for (i in 1:nobs)
			d_Sigma[i,i] <- d_eta_sq + d_sigma_sq;
		
		// sample parameters
		logit_p ~ multi_normal(d_mu, d_Sigma);
		for (i in 1:nobs)
			p[i] <- inv_logit(logit_p[i]);
		'
	},	
	'		/// calculate psi',
	if (ncol(control$data$X)==1) {
		'
		for (i in 1:nsites)
			psi[i] <- inv_logit(o_intercept);
		'
	} else {
		'
		// off-diagonal elements
		for (i in 1:(nsites-1)) {
			for (j in (i+1):nsites) {
				o_Sigma[i,j] <- 0;
				for (k in 2:nopars) {
					// o_Sigma[i,j] <- o_Sigma[i,j] + o_rho_sq[k-1] * pow(X[i,k] - X[j,k], 2);
					o_Sigma[i,j] <- o_Sigma[i,j] + o_rho_sq[1] * pow(X[i,k] - X[j,k], 2);
				}
				o_Sigma[i,j] <- o_eta_sq * exp(o_Sigma[i,j]);
				o_Sigma[j,i] <- o_Sigma[i,j];
			}
		}
		// diagonal elements
		for (i in 1:nsites)
			o_Sigma[i,i] <- o_eta_sq + o_sigma_sq;
			
		// sample parameters
		logit_psi ~ multi_normal(o_mu, o_Sigma);
		for (i in 1:nsites)
			psi[i] <- inv_logit(logit_psi[i]);
		'
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
			
	'	
	)
	
	assign('control', control, envir=globalenv())
	# stop()
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
		optim,
		control
	)
	assign('fm', fm, envir=globalenv())
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





