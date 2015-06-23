library(boot)
library(rstan)

# test_that("distribution functions don't work", {
	x1=Uniform(20,50)
	expect_identical(repr(x1), 'uniform(20, 50)')
	x2=Normal(0,10)
	expect_identical(repr(x2), 'normal(0, 10)')
	x3=Cauchy(0, 5)
	expect_identical(repr(x3), 'cauchy(0, 5)')
# })


# test_that("optimx solver doesn't work", {
	# init data
	data(frogs)
	pferUMF <- unmarkedFrameOccu(pfer.bin)
	siteCovs(pferUMF) <- data.frame(sitevar1 = rnorm(numSites(pferUMF)))
	obsCovs(pferUMF) <- data.frame(obsvar1 = rnorm(numSites(pferUMF) * obsNum(pferUMF)))
		
	# run using unmarked
	set.seed(500)
	m1=unmarked::occu(~ obsvar1 ~ 1, pferUMF, method='BFGS')
	
	# run using unmarkedExtra
	set.seed(500)
	m2=unmarkedExtra::occu(~ obsvar1 ~ 1, pferUMF, method='BFGS')
	
	# test if answers are the same
	expect_identical(m1@estimates,m2@estimates)
	expect_identical(m1@opt,m2@opt)
# })


# test_that("stan linear solver doesn't work", {
	## simulate data
	# set params
	n.sites=500
	n.visits=20
	o.int=2
	o.slope=2
	d.int=1
	
	# generate design matrices
	o.covs=data.frame(sitevar=rnorm(n.sites))
	o.form=~sitevar
	d.covs=data.frame(obsvar=rnorm(n.sites*n.visits))
	d.form=~1
	
	# generate prob. predictions
	psi=inv.logit(model.matrix(o.form, o.covs) %*% c(o.int, o.slope))
	site.obs=rbinom(n.sites,1,psi)
	p=inv.logit(model.matrix(d.form, d.covs) %*% c(d.int))
	visit.obs=t(sapply(seq_len(n.sites), function(i) {
		sapply(seq_len(n.visits), function(j) {
			rbinom(1,1,site.obs[i] * p[j])
		})
	}))
	
	# generate binary predictions
	testUMF=unmarkedFrameOccu(
			y=visit.obs,
			siteCovs=o.covs,
			obsCovs=d.covs
	)

	# run using optim
	m1=occu(~1 ~sitevar, data=testUMF, method='BFGS')
	m1.coef=coef(m1)
	
	# run using stan
	m2=occu(~1 ~sitevar, data=testUMF, method='stan')
	m2.samples=extract(m2, c('dpars','opars'))
	m2.coef=structure(
		c(colMeans(m2.samples[[2]]), mean(m2.samples[[1]][[1]])),
		.Names=names(m1.coef)
	)
	
	# test that simulated parameters have been derived using optim and stan
	expect_equal(round(m1.coef), round(m2.coef), c(o.int, o.slope, d.int))
# })

test_that("stan gp solver doesn't work", {
	## simulate data
	# set params
	n.sites=20
	n.visits=10
	o.int=2
	o.slope=2
	d.int=1
	
	# generate design matrices
	o.covs=data.frame(sitevar=rnorm(n.sites))
	o.form=~sitevar
	d.covs=data.frame(obsvar=rnorm(n.sites*n.visits))
	d.form=~1
	
	# generate prob. predictions
	psi=inv.logit(model.matrix(o.form, o.covs) %*% c(o.int, o.slope))
	site.obs=rbinom(n.sites,1,psi)
	p=inv.logit(model.matrix(d.form, d.covs) %*% c(d.int))
	visit.obs=t(sapply(seq_len(n.sites), function(i) {
		sapply(seq_len(n.visits), function(j) {
			rbinom(1,1,site.obs[i] * p[j])
		})
	}))
	
	# generate binary predictions
	testUMF=unmarkedFrameOccu(
			y=visit.obs,
			siteCovs=o.covs,
			obsCovs=d.covs
	)

	# run using stan
	m1=occu(~1 ~sitevar, data=testUMF, method='stan', control=list(gp=TRUE, iter = 5000))
		
	# compare predictions
	# m1.samples=rstan::extract(m1, c('psi','p'))
	# expect_equal(round(colMeans(m1.samples$psi),1),round(psi,1))
	# expect_equal(round(colMeans(m1.samples$p),1), round(p,1))
	
})



