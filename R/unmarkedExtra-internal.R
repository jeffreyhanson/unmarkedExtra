convertUMF=function(formula, data, knownOcc, method, control) {
	if (!is.null(data)) {
		control=list()
		control$y<-c(t(unmarked:::truncateToBinary(data@y)))
		if (identical(as.formula(paste("~", formula[3], sep="")), ~ 1)) {
			control$X <- model.matrix(as.formula(paste("~", formula[3], sep="")), data.frame(x=seq_len(nrow(data@y))))
		} else {
			control$X<-model.matrix(as.formula(paste("~", formula[3], sep="")), data@siteCovs)
		}
		if (identical(as.formula(formula[[2]]), ~ 1)) {
			control$V <- model.matrix(as.formula(formula[[2]]), data.frame(x=seq_len(length(data@y))))
		} else {
			control$V<-model.matrix(as.formula(formula[[2]]), data@obsCovs)
		}
		# handle NA values
		if (any(is.na(control$y))) {
			# remove sites which have never been visited
			is.dets.missing<-rowSums(!is.na(data@y))
			if (any(is.dets.missing==0)) {
				control$X<-control$X[which(is.dets.missing>0),]
			}
			# remove missing observations from data
			control$V<-control$V[which(!is.na(control$y)),,drop=FALSE]
			control$y<-control$y[which(!is.na(control$y))]
		}
		# calculate numbers for stan
		control$nobs=length(control$y) # number observations total
		control$nsites=nrow(control$X) # number of sites
		obssites=rep(seq_len(nrow(data@y)), each=obsNum(data))[which(!is.na(c(t(data@y))))] # number observations per site
		control$site_starts=as.integer(unname(tapply(X=seq_along(obssites), INDEX=obssites, FUN=min)))
		control$site_visits=as.vector(table(obssites))
		control$nopars=ncol(control$X) # number observation parameters
		control$ndpars=ncol(control$V) # number detection parameters
	} else {
		control=NULL
	}
	return(control)
}