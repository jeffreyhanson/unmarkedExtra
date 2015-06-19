
distribution=function(name, p1=NULL, p2=NULL) {
	return(
		structure(
			list(
				name=name,
				p1=p1,
				p2=p2
			),
			.Names=c('name','p1','p2'),
			class='distribution'
		)
	)
}

repr<-function(x, ...) {UseMethod('repr')}

repr.distribution<-function(x, ...) {
	return(
		paste0(x$name, '(', x$p1, ', ', x$p2, ')')
	)
}

#' @export
Normal<-function(mean, sd) {
	return(distribution(name='normal', p1=mean, p2=sd))
}

#' @export
Uniform<-function(lower, upper) {
	return(distribution(name='uniform', p1=lower, p2=upper))
}

#' @export
Cauchy<-function(mean, sd) {
	return(distribution(name='cauchy', p1=mean, p2=sd))
}

