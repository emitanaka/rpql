#################
## Auxilary functions
##################


#' Constructs a start fit for use in the \code{rpql} function
#' 
#' Takes a GLMM fitted using the \code{lme4} package i.e., using either the
#' \code{lmer} or \code{glmer} functions, and construct a list containing
#' starting values for use in the \code{start} argument in main fitting
#' function \code{rpql}. It also constructs adaptive lasso weights, which can
#' subsequently be used in the \code{pen.weights} arguments in the \code{rpql}
#' function, if the adaptive lasso penalty is used for variable selection.
#' 
#' This function is mainly used when: 1) you want to produce good starting
#' values for the main fitting function \code{rpql}, and so you fit a saturated
#' (full) GLMM using \code{lme4} and use the estimates from there as starting
#' values, and/or 2) you want to obtain adaptive lasso weights of the form
#' \eqn{weight_k = |\tilde{parameter}_k|^{-\gamma}}, where \eqn{\gamma > 0} is
#' the power parameter and \eqn{\tilde{parameter}_k} is the parameter estimate
#' from the saturated model fit. For regularized PQL specifically, this
#' function will construct adaptive lasso weights from the \code{lme4} fit as
#' follows: Let \eqn{w^F} and \eqn{w^R} denote fixed and random effect adaptive
#' weights respectively. Then we have,
#' 
#' \deqn{w^F_k = |\tilde{\beta}_k|^{-\gamma_1}} \deqn{w^R_l =
#' |\tilde{\Sigma}_{ll}|^{-\gamma_2},}
#' 
#' where \eqn{\tilde{\beta}_k} is the estimated coefficient for the
#' \eqn{k^{th}} fixed effect, \eqn{\tilde{\Sigma}_{ll}} is the \eqn{l^{th}}
#' diagonal element from the estimated random effects covariance matrix, and
#' \eqn{\gamma} is a vector of two power parameters; see Zou (2006) for the
#' adaptive lasso, and Hui et al. (2016) for regularized PQL selection in GLMMs
#' using on adaptive lasso type penalties.
#' 
#' If \code{cov.groups} is supplied, this means that some of the fixed effects
#' coefficients should be treated and penalized collectively as a group. The
#' most common cases where this is used is when you have factor or categorical
#' variables with more than two levels, or when you have polynomial terms that
#' should be dealt with together. For instance, suppose you have a model matrix
#' consisting of six columns, where first three columns correspond to separate
#' covariates (including the intercept) and the last three columns all
#' correspond to dummy variables created for a factor variable with four levels
#' , e.g. soil moisture with levels dry, moderately moist, very moist, wet. The
#' coefficients from the last three columns should then be penalized together,
#' and so we can set \code{cov.groups = c(1,2,3,4,4,4)}.
#' 
#' In doing so, the adaptive lasso weights for the grouped coefficients are
#' then constructed differently. Following on from the example above, we have
#' the fixed effect weight for soil moisture defined as
#' 
#' \deqn{w^F = \|\tilde{\beta}\|^{-\gamma_1},}
#' 
#' where \eqn{\| \cdot \|} corresponds to the L2-norm and \eqn{\tilde{\beta}}
#' are the fixed effect coefficients belonging in the group (three in this
#' case). When entered into the \code{rpql} function, an adaptive group lasso
#' (Wang and Leng, 2008) is applied to these set of coefficients, such that
#' they are all encouraged to be shrunk to zero at the same time.
#' 
#' Of course, after construction the adaptive lasso weights can be manually
#' altered before entering into the main \code{rpql} function e.g., if one
#' wants certain fixed and/or random effects to not be penalized.
#' 
#' @param lme4.fit An object of class "lmerMod" or "glmerMod", obtained when
#' fitting a (G)LMM using the \code{lmer} and \code{glmer} functions in the
#' \code{lme4} package.
#' @param id A optional list with each element being a vector of IDs that
#' reference the model matrix in the corresponding element in the list
#' \code{Z}. Each vector of IDs \emph{must} be integers (but not factors). Note
#' this is optional argument as it is only use for non-compulsory formatting
#' purposes in the function.
#' @param gamma A vector of power parameters, \eqn{\gamma}, for use in
#' constructing adaptive lasso weights. Can be a vector of one or two elements.
#' If two elements, then the first and second elements are the power parameter
#' for the fixed and random effect weights respectively. If one element, the
#' same power parameter is used for both fixed and random effect weights.
#' Defaults to 0, in which case the weights are all equal to 1 i.e., it reverts
#' to the unweighted lasso penalty.
#' @param cov.groups A vector specifying if fixed effect coefficients
#' (including the intercept) should be regarded and therefore penalized in
#' groups. For example, if one or more of the fixed effect covariates are
#' factors, then \code{lme4} will automatically create dummy variables in the
#' model matrix and estimate coefficients for each level, using one level as
#' the reference. \code{cov.groups} is then used to identify all the
#' coefficients that corresponds to that factor, such that all of these
#' coefficients are penalized collectively as a group. Defaults to NULL, in
#' which case it is assumed all coefficients should be treated independently.
#' Please see the details and examples for more details.
#' @return A list containing the following elements \item{fixef}{Fixed effect
#' coefficient estimates from \code{lme4.fit}.} \item{ranef}{A list of random
#' effect predicted coefficients from \code{lme4.fit}.} \item{ran.cov}{A list
#' of random effects covariance matrices from \code{lme4.fit}.}
#' \item{cov.groups}{The argument \code{cov.groups}. Defaults to \code{NULL}.}
#' \item{pen.weights}{A list of adaptive lasso weights constructed from
#' \code{lme4.fit}. Contains elements \code{pen.weight$fixed} and
#' \code{pen.weights$random}, which are the weights for the fixed and random
#' effects respectively. Please see details above as to their construction.}
#' @section Warnings: \itemize{ \item In order to construct sensible starting
#' values and weights, this function should really only be used when
#' \code{lme4.fit} is a fit of the saturated GLMM, i.e. all fixed and random
#' effects included.  }
#' @author
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_author(\"#1\")}",
#' "rpql")\Sexpr{tools:::Rd_package_author("rpql")}
#' 
#' Maintainer:
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_maintainer(\"#1\")}",
#' "rpql")\Sexpr{tools:::Rd_package_maintainer("rpql")}
#' @seealso \code{\link{rpql}} for fitting and performing model selection in
#' GLMMs using regularized PQL, which may use the values obtained from
#' \code{build.start.fit} for starting values and adaptive lasso weights.
#' @references \itemize{ \item Hui, F.K.C., Mueller, S., and Welsh, A.H.
#' (2016). Joint Selection in Mixed Models using Regularized PQL. Journal of
#' the American Statistical Association: accepted for publication.  \item Wang,
#' H., and Leng, C. (2008). A note on adaptive group lasso. Computational
#' Statistics & Data Analysis, 52, 5277-5286.  \item Zou, H. (2006). The
#' adaptive lasso and its oracle properties. Journal of the American
#' statistical association, 101, 1418-1429.  }
#' @examples
#' 
#' 
#' ##################
#' ## Example 1: Bernoulli GLMM with grouped covariates. 
#' ## Independent cluster model with 50 clusters and equal cluster sizes of 10
#' ## Nine covariates where the last covariate (soil type) is a factor with four levels
#' n <- 50; p <- 8; m <- 10
#' set.seed(123)
#' X <- data.frame(matrix(rnorm(n*m*p),n*m,p), soil=sample(1:4,size=m*n,replace=TRUE))
#' X$soil <- factor(X$soil)
#' X <- model.matrix(~ ., data = X)
#' colnames(X) <- paste("X",1:ncol(X),sep="")
#' 
#' Z <- X[,1:5] ## Random effects model matrix taken as first five columns
#' true_betas <- c(-0.1,1,-1,1,-1,1,-1,0,0,0,0,0) 
#' true_D <- matrix(0,ncol(Z),ncol(Z))
#' true_D[1:3,1:3] <- matrix(c(9,4.8,0.6,4.8,4,1,0.6,1,1),
#' 	3,3,byrow=TRUE) ## 3 important random effects 
#' 
#' simy <- gendat.glmm(id = list(cluster = rep(1:n,each=m)), X = X, beta = true_betas, 
#' 	Z = list(cluster = Z), D = list(cluster = true_D), family = binomial())
#' 
#'   
#' \dontrun{
#' library(lme4)
#' dat <- data.frame(y = simy$y, simy$X, simy$Z$cluster, simy$id)
#' fit_satlme4 <- glmer(y ~ X - 1 + (Z - 1 | cluster), data = dat, 
#' 	family = "binomial")
#' fit_sat <- build.start.fit(fit_satlme4, id = simy$id, gamma = 2, 
#' 	cov.groups = c(1:9,10,10,10)) 
#' 
#' new.fit <- rpql(y = simy$y, X = simy$X, Z = simy$Z, id = simy$id, lambda = 0.01, 
#' 	pen.type = "adl", pen.weights = fit_sat$pen.weights,
#' 	cov.groups = fit_sat$cov.groups, start = fit_sat, family = binomial())  
#' 	}
#' 
#' 
build.start.fit <- function(lme4.fit, id = NULL, gamma = 0, cov.groups = NULL) 
     {
	if(!(class(lme4.fit)[1] %in% c("glmerMod","lmerMod"))) 
          stop("lme4.fit must be fit from the lme4 package. Thanks")

	out <- list(fixef = fixef(lme4.fit)+1e-4)
	get_ranefs <- as.list(ranef(lme4.fit))
	for(k in 1:length(get_ranefs)) 
          get_ranefs[[k]] <- as.matrix(get_ranefs[[k]])
	if(!is.null(id)) 
          names(get_ranefs) <- names(id)
	out$ranef <- get_ranefs
	out$cov.groups <- cov.groups
	
	get_varcors <- vector("list", length(get_ranefs))
     names(get_varcors) <- names(get_ranefs)
	for(k in 1:length(get_ranefs)) 
		get_varcors[[k]] <- VarCorr(lme4.fit)[[k]]+1e-4
	out$ran.cov <- get_varcors

	
	if(length(gamma) == 1) 
		gamma <- rep(gamma,2)
	#cat("Building adaptive lasso weights...\n")
	if(is.null(cov.groups)) 
          out$pen.weights <- list(fixed = 1/abs(out$fixef)^gamma[1])
	if(!is.null(cov.groups)) 
          {
		make_adlweights <- (1/unlist(lapply(split(out$fixef,cov.groups), l2.norm)))[cov.groups] ## Split by cov.groups, calc L2 norm for each, convert to vector, then expand back...awesome!
		names(make_adlweights) <- names(out$fixef)
		out$pen.weights <- list(fixed = make_adlweights^gamma[1]) 
		}

	out$pen.weights$random <- vector("list", length(get_ranefs))
	names(out$pen.weights$random) <- names(get_ranefs)
	for(k in 1:length(get_ranefs)) 
          out$pen.weights$random[[k]] <- 1/diag(out$ran.cov[[k]]^gamma[2])
	
	return(out)
	}


## The marginal logL is calculated using Monte Carlo integration, exploiting the fact that conditional on all the random effects, all obs are independent. So generate a whole bunch of b's, then for b calculate the likelihood, then average over the b's


#' Calculate the marginal log-likelihood for a GLMM fitted using \code{rpql}
#' 
#' After fitting and performing joint (fixed and random effects) using
#' regularized PQL, one may then (for one reason or another) want to calculate
#' the marginal likelihood for the (sub)model, possibly on a test dataset for
#' prediction. This is the main purpose of \code{calc.marglogL}.
#' 
#' Regularized PQL performs penalized joint (fixed and random effects)
#' selection for GLMMs, where the penalized quasi-likelihood (PQL, Breslow and
#' Clayton, 1993) is used the loss function. After fitting, one may then wish
#' to calculate the marginal log-likelihood for the (sub)model, defined as
#' 
#' \deqn{\ell = \log\left(\int f(\bm{y}; \bm{\beta}, \bm{b}, \phi) f(\bm{b};
#' \bm{\Sigma}) d\bm{b}\right),}
#' 
#' where \eqn{f(\bm{y}; \bm{\beta}, \bm{b}, \phi)} denotes the conditional
#' likelihood of the responses \eqn{\bm{y}} given the fixed effects
#' \eqn{\bm{\beta}}, random effects \eqn{\bm{b}}, and nuisance parameters
#' \eqn{\phi} if appropriate, and \eqn{f(\bm{b}; \bm{\Sigma})} is the
#' multivariate normal distribution for the random effects, with covariance
#' matrix \eqn{\bm{\Sigma}}. \code{calc.marglogL} calculates the above marginal
#' likelihood using Monte-Carlo integration.
#' 
#' Admittedly, this function is not really useful for fitting the GLMM
#' \emph{per-se}: it is never called by the main function \code{rpql}, and the
#' marginal likelihood is (approximately) calculated anyway if \code{hybrid.est
#' = TRUE} and the final submodel is refitted using \code{lme4}. Where the
#' function comes in handy is if you have a validation or test dataset, and you
#' want to calculated the predicted (log) likelihood of the test data given the
#' regularized PQL fit.
#' 
#' @param new.data A list containing the elements \code{new.data$y},
#' \code{new.data$X}, and \code{new.data$Z}. These correspond respectively to
#' the responses, fixed effects model matrix, and random effects model matrices
#' that the marginal log-likelihood is be calculated on. No check is made
#' against the elements in \code{fit} to ensure that these are of the correct
#' dimensions compared, and furthermore it is assumed that \code{new.data$Z} is
#' a list in the same order as the \code{Z} used when fitting the original
#' model via \code{rpql}.
#' @param fit An object of class \code{pqrl}. In the least, \code{fit} should
#' be a list containing the elements \code{fit$family} for the family, e.g.
#' gaussian(), poisson(), \code{fit$fixef} for the estimated vector of fixed
#' effects, \code{fit$ran.cov} which is a list of estimated random effects
#' covariance matrices. If appropriate, \code{fit} may also contain the
#' elements \code{fit$phi} for the estimated variance parameter in normal,
#' lognormal, and negative binomial GLMMs, \code{fit$shape} for the estimated
#' shape parameter used in Gamma GLMMs, \code{fit$trial.size} for the trial
#' size(s) for binomial GLMMs, and \code{fit$zeroprob} for the estimated
#' probability of a structural zero in ZIP GLMMs.
#' @param B A positive integer for the number of random effects examples to
#' generate, when performing Monte-Carlo integration. Defaults to 1000.
#' @return The marginal log-likelihood of \code{new.data} given the GLMM in
#' \code{fit}.
#' @section Warnings: \itemize{ \item No check is made to see if the dimensions
#' of the elements \code{new.data} and \code{fit} match, e.g. the number of
#' columns in \code{new.data$X} is equal to the number of elements in
#' \code{fit$fixef}. Please ensure they are!  \item Monte-Carlo integration is
#' computationally intensive especially if \eqn{\bm{y}} is long!  }
#' @author
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_author(\"#1\")}",
#' "rpql")\Sexpr{tools:::Rd_package_author("rpql")}
#' 
#' Maintainer:
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_maintainer(\"#1\")}",
#' "rpql")\Sexpr{tools:::Rd_package_maintainer("rpql")}
#' @seealso \code{\link{rpql}} for fitting and performing model selection in
#' GLMMs using regularized PQL. \code{lme4} also approximately calculates the
#' marginal log-likelihood when fitting a GLMM.
#' @references \itemize{ \item Breslow, N. E., & Clayton, D. G. (1993).
#' Approximate inference in generalized linear mixed models. Journal of the
#' American Statistical Association, 88, 9-25.  }
#' @examples
#' 
#' ## Make your own =D
#' 
calc.marglogL <- function(new.data, fit, B = 1000) 
     {
	if(!all(c("y","X","Z") %in% attributes(new.data)$names)) 
		stop("new.data must be a list containing y, X, Z")
	new.data$X <- as.matrix(new.data$X)
	n <- nrow(new.data$X)
	family <- fit$family

	if(any(sapply(new.data$Z,nrow) != n)) 
		stop("The number of rows in each element of new.data$Z should be equal to the number of rows in new.data$X. Thanks.")

	get_ncolZ <- sapply(new.data$Z,ncol)
	get_ncolranef_covmat <- sapply(fit$ran.cov,ncol)
	if(any(get_ncolZ != get_ncolranef_covmat)) 
		stop("The number of columns in each in the list new.data$Z should equal to dimension of the corresponding element in the list fit$ran.cov. Thanks.")
	rm(get_ncolZ, get_ncolranef_covmat)
	for(k in 1:length(fit$ran.cov)) 
          new.data$Z[[k]] <- as.matrix(new.data$Z[[k]])

	
	loglik <- numeric(B); 
	## Generate random effects
	get_nonzeroranef_covmat <- genb <- vector("list",length(fit$ran.cov)); 
	for(k in 1:length(fit$ran.cov)) 
          { 
		get_nonzeroranef_covmat[[k]] <- which(diag(fit$ran.cov[[k]]) > 0)
		genb[[k]] <- matrix(0, B, ncol(new.data$Z[[k]]))
		genb[[k]][,get_nonzeroranef_covmat[[k]]] <- rmvnorm(B, mean = rep(0,length(get_nonzeroranef_covmat[[k]])), as.matrix(fit$ran.cov[[k]][get_nonzeroranef_covmat[[k]],get_nonzeroranef_covmat[[k]]]))
		}


	eta <- new.data$X %*% fit$fixef + fit$offset
	for(k0 in 1:B) 
          {
		eta2 <- eta
		for(k in 1:length(fit$ran.cov)) 
               eta2 <- eta2 + new.data$Z[[k]] %*% genb[[k]][k0,]
		
		if(family$family[1] == "gaussian") 
               tmp_loglik <- dnorm(new.data$y, mean = eta2, sd=sqrt(fit$phi), log = TRUE)
		if(family$family[1] == "Gamma") 
               tmp_loglik <- (dgamma(new.data$y, shape = fit$shape, scale = family$linkinv(eta2)/fit$shape, log = TRUE))
		if(family$family[1] == "poisson") 
               tmp_loglik <- (dpois(new.data$y, lambda=family$linkinv(eta2), log = TRUE))
		if(family$family[1] == "binomial") 
               tmp_loglik <- (dbinom(new.data$y, size=fit$trial.size, prob=family$linkinv(eta2), log = TRUE))
		if(family$family[1] == "negative.binomial") 
               tmp_loglik <- (dnbinom(new.data$y, mu=family$linkinv(eta2), size=1/fit$phi, log = TRUE))
		if(family$family[1] == "LOGNO") 
               tmp_loglik <- (dlnorm(new.data$y, meanlog=eta, sdlog=sqrt(fit$phi), log=TRUE))
		if(family$family[1] == "ZIP") 
               tmp_loglik <- (dZIP(new.data$y, mu=family$mu.linkinv(eta2), sigma=fit$zeroprob, log=TRUE))

		tmp_loglik[!is.finite(tmp_loglik)] <- NA
		loglik[k0] <- exp(sum(tmp_loglik,na.rm=TRUE))
		} 
	
	return(mean(loglik[is.finite(loglik)]))
	}


	
## Dataset generation for GLMM
## Intercept must be manually included in X and Z if desired
#id = list(cluster = rep(1:n,each=m), cluster2 = id2); beta = truebeta; D = list(cluster = true.D, cluster2 = true.D2); trial.size = 1; family = "binomial"; phi = NULL; upper.count = Inf


#' Simulates datasets based on a Generalized Linear Mixed Model (GLMM).
#' 
#' Datasets are simulated from a GLMM given a series of inputs including: model
#' matrices \code{X} and \code{Z} for the fixed and random effects
#' respectively, a set of true fixed effect coefficients \code{beta}, a list of
#' true random effect covariance matrices \code{D}, the family of response, and
#' some other nusiance parameters if appropriate.
#' 
#' The relationship between the mean of the responses and covariates in a GLMM
#' is given as follows: For \eqn{i = 1,\ldots,n}, where \eqn{n} is,
#' equivalently, the number of rows in \code{X}, the length of each element in
#' \code{id}, and the number of rows in each element of \code{Z}, we have
#' 
#' \deqn{g(\mu_{i}) = \bm{x}^T_i \bm{\beta} + \bm{z}^T_{i1} \bm{b}_{i1} +
#' \bm{z}^T_{i2} \bm{b}_{i2} + \ldots,}
#' 
#' where \eqn{g(\cdot)} is the link function, \eqn{\mu_i} is the mean of the
#' distribution for observation \eqn{i}, \eqn{\bm{x}_i} is row \eqn{i} of the
#' fixed effects model matrix \code{X}, and \eqn{\bm{\beta}} is the fixed
#' effects coefficients. For the random effects, \eqn{\bm{z}_{i1}} is row
#' \eqn{i} of the random effects model matrix in the first element of \code{Z},
#' while \eqn{\bm{b}_{i1}} is the vector of random effects generated for
#' observation \eqn{i} based on the first element of \code{D}. The remaining
#' parameters \eqn{\bm{z}_{i2}}, \eqn{\bm{b}_{i2}} and so on, are defined
#' similarly.
#' 
#' Having lists for \code{id, Z}, and \code{D} allows for multiple sets of
#' random effects to be included in the true GLMM. This is analogous to the
#' \code{lme4} package, where multiple random effects are permitted in the
#' formula, e.g., \code{(1|creek) + (1|creek:sample)}. If the true GLMM
#' contains only one set of random effects, e.g., in longitudinal data, then
#' the three lists will all contain only one element. Cases with multiple sets
#' of random effects include nested and crossed designs, in which case
#' \code{id, Z}, and \code{D} will have two or more elements.
#' 
#' It is recommended that the user think through and design these lists
#' carefully to ensure that they are actually constructing a true GLMM that
#' they want to simulated data from. Yes it takes some getting use too, and we
#' apologize for this =( Please see examples below for some ideas.
#' 
#' Finally, note that some of the elements of \code{beta} can be zero, i.e.
#' truly unimportant fixed effects. Likewise, each element of \code{D} can be a
#' random effects covariance matrix containing zero rows and columns, i.e.
#' truly unimportant random effects.
#' 
#' @param id A list with each element being a vector of IDs that reference the
#' model matrix in the corresponding element in the list \code{Z}. Each vector
#' of IDs \emph{must} be integers (but not factors).
#' @param X A model matrix of corresponding to the fixed effects. A column of
#' ones should be included if a fixed intercept is to be included in the model.
#' @param beta A vector of true fixed effect parameters, with the same length
#' as the number of columns in \code{X}.
#' @param Z A list with each element being a model matrix for a set of random
#' effects. Each element of \code{Z} is referenced by a vector of IDs given by
#' the corresponding element in the list \code{id}. Each model matrix (element
#' of \code{Z}) should have the same number of rows as the length of \code{y}.
#' @param D A list with each element being a symmetric random effects
#' covariance matrix which is used to generate random effects. These random
#' effects are then applied to the corresponding element in the list \code{Z},
#' and are referenced by the corresponding element in the list \code{id}.
#' @param trial.size The trial size if \code{family = binomial()}. Either takes
#' a single non-zero value or a vector of non-zero values with length the same
#' as the number of rows in \code{X}. The latter allows for differing trial
#' sizes across responses. Defaults to 1.
#' @param family The distribution for the responses in GLMM. The argument must
#' be applied as a object of class "family". Currently supported arguments
#' include: \code{gaussian()}, \code{poisson()}, \code{binomial()},
#' \code{Gamma()}, \code{nb2()} for negative binomial, \code{LOGNO()} for
#' log-normal, and \code{ZIP()} for zero-inflated Poisson.
#' @param phi A non-zero value for the true variance parameter \eqn{\sigma^2}
#' if \code{family = gaussian()}, the true variance parameter \eqn{\sigma^2} on
#' the log scale if \code{family = LOGNO()}, or the overdispersion parameter if
#' \code{family = nb2()}, where the negative binomial variance is parameterized
#' as \eqn{V = \mu + \phi\mu^2}. Defaults to NULL.
#' @param shape A non-zero value for the shape parameter \eqn{a} if
#' \code{family = Gamma()}, where the variance is parameterized as \eqn{V =
#' \mu^2/a}. Defaults to NULL.
#' @param zeroprob A value between 0 and 1 for the probability of a structural
#' zero if \code{family = ZIP()} for zero-inflated Poisson. Defaults to NULL.
#' @param upper.count A non-zero integer which allows the user to control the
#' maximum value of the counts generates for datasets when \code{family =
#' poisson()} or \code{nb2()}. When the responses are simulated, a \code{while}
#' loop is run to ensure that all responses generated are less than or equal to
#' \code{upper.count}. Default to \code{Inf}.
#' @return A list containing the following elements \item{y}{The vector
#' simulated responses.} \item{b}{A list with each element being a matrix of
#' random effects simulated from a multivariate normal distribution with mean
#' zero and covariance matrix equal to the corresponding element in the list
#' \code{D}. For each element in \code{b}, the number of columns of the matrix
#' equals the dimension of corresponding covariance matrix element in \code{D},
#' while the number of rows equals to the number of unique IDs in the
#' corresponding element of the list \code{id}.} \item{id, X, Z, beta, D, phi,
#' shape, zeroprob, trial.size, family}{Some of the arguments entered into
#' \code{gendat.glmm}.} \item{nonzero.beta}{A vector indexing the non-zero
#' values of \code{beta}, i.e. the truly important fixed effects.}
#' \item{nonzero.b}{A list with each element being a vector indexing the
#' non-zero diagonal variances in the corresponding element of the list
#' \code{D}, i.e. the truly important random effects.}
#' @author
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_author(\"#1\")}",
#' "rpql")\Sexpr{tools:::Rd_package_author("rpql")}
#' 
#' Maintainer:
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_maintainer(\"#1\")}",
#' "rpql")\Sexpr{tools:::Rd_package_maintainer("rpql")}
#' @seealso \code{\link{rpql}} for fitting and performing model selection in
#' GLMMs using regularized PQL.
#' @references \itemize{ \item Schielzeth, H., & Nakagawa, S. (2013). Nested by
#' design: model fitting and interpretation in a mixed model era. Methods in
#' Ecology and Evolution, 4, 14-24.  }
#' @examples
#' 
#' 
#' ##################
#' ## Example 1: Linear Mixed Models 
#' ## Independent cluster model with 50 clusters
#' ## Nine covariates including a fixed and random intercept
#' library(mvtnorm)
#' library(lme4)
#' 
#' n <- 50; m <- 10; p <- 8; 
#' ## Generate rows of a model matrix from a multivariate normal distribution with 
#' ## AR1 covariance structure. 
#' 
#' H <- abs(outer(1:p, 1:p, "-")) 
#' X <- cbind(1,rmvnorm(n*m,rep(0,p),sigma=0.5^H)); 
#' 
#' Z <- X 
#' true_betas <- c(1,3,2,1.5,-1,0,0,0,0) ## 5 important fixed effects 
#' true_D <- matrix(0,p+1,p+1) ## 3 important random effects
#' true_D[1:3,1:3] <- matrix(c(9,4.8,0.6,4.8,4,1,0.6,1,1),3,3,byrow=TRUE)
#' 
#' simy <- gendat.glmm(id = list(cluster = rep(1:n,each=m)), X = X, beta = true_betas, 
#' 	Z = list(cluster = Z), D = list(cluster = true_D), phi = 1, family = gaussian()) 
#' ## Notice how id, Z, and D all are lists with one element, and that 
#' ## the name of the first element (a generic name "cluster") is the 
#' ## same for all three lists. 
#' ## id is where the action takes place. In particular, id$cluster is 
#' ## designed so that the first 12 elements correspond to cluster 1, 
#' ## the second 12 elements correspond to cluster 2, and so forth. 
#' ## In turn, the first 12 rows of X and Z$cluster correspond 
#' ## to cluster 1, and so on. 
#' 
#' \dontrun{
#' dat <- data.frame(y = simy$y, simy$X, simy$Z$cluster, simy$id)
#' fit_satlme4 <- lmer(y ~ X - 1 + (Z - 1 | cluster), data = dat,
#' 	REML = FALSE)
#' fit_sat <- build.start.fit(fit_satlme4, gamma = 2)
#' 
#' 
#' lambda_seq <- lseq(1e-4,1,length=100)
#' fit <- rpqlseq(y = simy$y, X = simy$X, Z = simy$Z, id = simy$id, 
#' 	family = gaussian(), lambda = lambda_seq, pen.type = "adl", 
#' 	pen.weights = fit_sat$pen.weights, start = fit_sat)
#' 
#' summary(fit$best.fit[[3]])  
#' # apply(fit$collect.ics, 2, which.min) ## Look at best fit chosen by different ICs
#' }
#' 
#' 
#' ##################
#' ## Example 2: Bernoulli GLMMs on simulated data
#' ## Nested data with 200 observations in total: split into 10 creeks, 
#' ## 5 samples nested within each creek
#' 
#' mn <- 200; 
#' X <- as.matrix(rep(1,mn)); 
#' ids <- list(samples = rep(1:50,each=4), creek = rep(1:10,each=20)) 
#' ## We have two sets of random intercepts only, one for creek and one 
#' ## 	for samples nested within creek.
#' Zs <- list(samples = X, creek = X) 
#' 
#' true_betas <- -0.1 
#' ## Please ensure each element of true_D is a matrix
#' true_D <- list(samples = as.matrix(0.001), creek = as.matrix(1)) 
#' 
#' simy <- gendat.glmm(id = ids, X = X, beta = true_betas, Z = Zs, D = true_D, 
#' 	trial.size = 1, family = binomial())
#' 
#' \dontrun{
#' 
#' ## Construct a solution path use adaptive LASSO for selection
#' ## Here is another way of constructing the adaptive weights:
#' ## Use the fact that rpql can do a final fit based on maximum likelihood
#' ## to obtain a good saturated fit.
#' fit_sat <- rpql(y = simy$y, X = simy$X, Z = simy$Z, id = simy$id, 
#' 	family = binomial(), lambda = 0, hybrid = TRUE)
#' fit_sat <- build.start.fit(fit_sat$hybrid, gamma = 2)
#' 	
#' lambda_seq <- lseq(1e-6,1,length=100)
#' fit <- rpqlseq(y = simy$y, X = simy$X, Z = simy$Z, id = simy$id, 
#' 	family = binomial(), lambda = lambda_seq, pen.type = "adl", 
#' 	pen.weights = fit_sat$pen.weights, start = fit_sat)
#' 
#' summary(fit$best.fit[[3]])  
#' # apply(fit$collect.ics, 2, which.min) ## Look at best fit chosen by different ICs
#' }
#' 
gendat.glmm <- function(id, X, beta, Z, D, trial.size = 1, family = gaussian(), phi = NULL, shape = NULL, zeroprob = NULL, upper.count = Inf) 
     {
	X <- as.matrix(X)
	n <- nrow(X)
	if(ncol(X) != length(beta)) 
		stop("The number of columns in X should be equal to the length of beta. Thanks")
	
	if(sum(X[,1]) != nrow(X)) 
		warning("Has an intercept column being included in X?")
	if(!is.list(id) || !is.list(D) || !is.list(Z)) 
		stop("Please supply id, Z, D, as lists. The elements in the three lists corresponding to a vector of IDs, and the model matrix and the random effects covariance matrix for those IDs, respectively. Please see the help file for further details. ")
	if(length(unique(length(id),length(D),length(Z))) != 1) 
		stop("The number of elements in the lists id, D, and Z should be the same. Thanks.")

	get_dimD <- sapply(D, ncol)
	get_ncolZ <- sapply(Z, ncol)
	if(any(get_dimD!=get_ncolZ)) 
		stop("The dimension of each element in the list D should be equal to the number of columns in the corresponding element of list Z. Thanks.")

		
	get_lengths_id <- sapply(id,length)
	if(any(get_lengths_id!=nrow(X))) 
		stop("The length of each element in the list id should be equal to the number of rows in X and Z. Thanks.")
	rm(get_lengths_id, get_dimD, get_ncolZ)
 
	
	if(!(family$family[1] %in% c("gaussian","poisson","binomial","negative.binomial","Gamma","LOGNO","ZIP"))) 
		stop("Current version does not permit specified family. Sorry!")
	if(family$family[1] == "Gamma" & is.null(shape)) 
		stop("Please supply shape parameter for the Gamma family (Variance = mu^2/shape).") 
	if(family$family[1] == "negative.binomial" & is.null(phi)) 
		stop("Please supply the overdispersion parameter phi for the negative binomial family; note the variance is parameterized as V = mu + phi*mu^2.") 
	if(family$family[1] == "gaussian" & is.null(phi)) 
		stop("Please supply the variance parameter phi for the Gaussian family; note the variance is parameterized as V = phi.") 
	if(family$family[1] == "LOGNO" & is.null(phi)) 
		stop("Please supply the variance parameter phi for the lognormal family; note the variance is parameterized as V = phi on the log scale).") 
	if(family$family[1] == "ZIP" & is.null(zeroprob)) 
		stop("Please supply zeroprob, the probability for obtaining a structural zero, for the ZIP family.") 
	
	if(is.null(colnames(X))) 
          colnames(X) <- paste("x",1:ncol(X),sep="")
	for(k in 1:length(Z)) 
          colnames(Z[[k]]) <- paste("z",k,1:ncol(Z[[k]]),sep="")

	get.n.id <- numeric(length(id))
	for(k in 1:length(id)) 
          get.n.id[k] <- length(unique(id[[k]]))
	
	for(k in 1:length(D)) 
          { 
          D[[k]] <- as.matrix(D[[k]])
          Z[[k]] <- as.matrix(Z[[k]]) 
          }
	for(k in 1:length(id)) 
          { 
          id[[k]] <- as.integer(id[[k]]) 
          }

	
	## Generate random effects
	get_nonzeroD <- trueb <- vector("list",length(D)); 
	for(k in 1:length(trueb)) 
          { 
		get_nonzeroD[[k]] <- which(diag(D[[k]]) > 0)
		trueb[[k]] <- matrix(0, nrow = get.n.id[k], ncol = ncol(Z[[k]]))
		trueb[[k]][,get_nonzeroD[[k]]] <- rmvnorm(get.n.id[k], rep(0,length(get_nonzeroD[[k]])), as.matrix(D[[k]][get_nonzeroD[[k]],get_nonzeroD[[k]]]))
		}

		
	## Generate response	
	sim_y <- numeric(n)
	eta <- X %*% beta 
	for(k in 1:length(trueb)) 
          eta <- eta + rowSums(Z[[k]]*trueb[[k]][id[[k]],])

	if(family$family[1] == "gaussian") 
          sim_y <- rnorm(n, mean = eta, sd=sqrt(phi))
	if(family$family[1] == "LOGNO") 
          sim_y <- rlnorm(n, meanlog=eta, sdlog=sqrt(phi))
	if(family$family[1] == "Gamma") 
          sim_y <- rgamma(n, shape=shape, scale=family$linkinv(eta)/shape)
	if(family$family[1] == "binomial") 
          sim_y <- rbinom(n, size=trial.size, prob=family$linkinv(eta))
	for(i in 1:n) 
          {
		if(family$family[1] %in% c("poisson","negative.binomial","ZIP")) 
               { 
			if(family$family[1] == "ZIP") 
                    testsim_y <- rZIP(1, mu=family$mu.linkinv(eta[i]), sigma=zeroprob)
			if(family$family[1] == "poisson") 
                    testsim_y <- rpois(1, lambda=family$linkinv(eta[i]))
			if(family$family[1] == "negative.binomial") 
                    testsim_y <- rnbinom(1, mu=family$linkinv(eta[i]), size=1/phi)
			try_counter <- 0
			
			while(testsim_y > upper.count & try_counter < 5000) 
                    {
				if(family$family[1] == "ZIP") 
                         testsim_y <- rZIP(1, mu=family$mu.linkinv(eta[i]), sigma=zeroprob)
				if(family$family[1] == "poisson") 
                         testsim_y <- rpois(1, lambda=family$linkinv(eta[i]))
				if(family$family[1] == "negative.binomial") 
                         testsim_y <- rnbinom(1, mu=family$linkinv(eta[i]), size=1/phi)
				try_counter <- try_counter + 1
				}
			sim_y[i] <- testsim_y
			}		
		}

		
	out <- list(y = sim_y, id = id, X = X, Z = Z, beta = beta, b = trueb, D = D, phi = phi, shape = shape, zeroprob = zeroprob, trial.size = trial.size, family = family, nonzero.beta = which(beta!=0), nonzero.b = get_nonzeroD)
	names(out$b) <- names(out$nonzero.b) <- names(D)
	
	return(out)
	}

	
	
l2.norm <- function(x) 
     sqrt(sum(x^2)) 




#' Generates a sequence of tuning parameters on the log scale
#' 
#' Generates a sequence of tuning parameters \eqn{\lambda} that are equally
#' spaced on the log-scale. It may be used as part of constructing a solution
#' path for the main fitting function \code{rpql}.
#' 
#' For joint selection of fixed and random effects in GLMMs, regularized PQL
#' (Hui et al., 2016) works taking the penalized quasi-likelihood (PQL, Breslow
#' and Clayton, 1993) as a loss function, and then sticking on some penalties
#' in order to model variable. The penalties will depend upon one or more
#' tuning parameters \eqn{\lambda > 0}, and the typical way this is chosen is
#' to construct a sequence of \eqn{\lambda} values, fit the regularized PQL to
#' each one value, and then use a method like information criterion to select
#' the best \eqn{\lambda} and hence the best model. Please see the help file
#' for \code{\link{rpql}} for more details, and \code{glmnet} (Friedman et al.,
#' 2010) and \code{ncvreg} (Breheny, and Huang, 2011) as examples of other
#' packages that do penalized regression and involve tuning parameter
#' selection.
#' 
#' The idea of equally spacing the sequence of \eqn{\lambda}'s on the log (base
#' 10) scale may not necessary be what you want to do, and one is free to use
#' the standard \code{seq()} function for constructing sequences. By equaling
#' spacing them on log-scale, it means that there will be a large concentration
#' of small tuning parameter values, with less large tuning parameter values
#' (analogous to a right skewed distribution). This may be useful if you
#' believe the that most of the penalization/variable selection action takes
#' place on smaller values of \eqn{\lambda}.
#' 
#' It is somewhat of an art form to construct a good sequence of tuning
#' parameter values: the smallest \eqn{\lambda} should produce the saturated
#' model if possible, and the largest \eqn{\lambda} should shrink most if not
#' all covariates to zero i.e., the null model. Good luck!
#' 
#' @param from The minimum tuning parameter to start the sequence from.
#' @param to The maximum tuning parameter to go to.
#' @param length The length of the sequence.
#' @param decreasing Should the sequence be in ascending or descending order?
#' @return A sequence of tuning parameter values of length equal to
#' \code{length}.
#' @author
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_author(\"#1\")}",
#' "rpql")\Sexpr{tools:::Rd_package_author("rpql")}
#' 
#' Maintainer:
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_maintainer(\"#1\")}",
#' "rpql")\Sexpr{tools:::Rd_package_maintainer("rpql")}
#' @seealso \code{\link{rpql}} for fitting and performing model selection in
#' GLMMs using regularized PQL.
#' @references \itemize{ \item Breheny, P. and Huang, J. (2011) Coordinate
#' descent algorithms fof nonconvex penalized regression, with applications to
#' biological feature selection. The Annals of Appliedv Statistics, 5, 232-253.
#' \item Breslow, N. E., and Clayton, D. G. (1993). Approximate inference in
#' generalized linear mixed models. Journal of the American Statistical
#' Association, 88, 9-25. \item Friedman, J., Hastie T., and Tibshirani, R.
#' (2010). Regularization Paths for Generalized Linear Models via Coordinate
#' Descent. Journal of Statistical Software, 33, 1-22. URL:
#' http://www.jstatsoft.org/v33/i01/.  \item Hui, F.K.C., Mueller, S., and
#' Welsh, A.H. (2016). Joint Selection in Mixed Models using Regularized PQL.
#' Journal of the American Statistical Association: accepted for publication.
#' }
#' @examples
#' 
#' ## Please see examples in help file for the rpql function
#' 
lseq <- function (from, to, length, decreasing = FALSE) 
     {
	stopifnot(from > 0)
	out <- 10^(seq(log10(from), log10(to), length.out = length))
	out <- out[order(out, decreasing = decreasing)]; 
	return(out)
	}

	
## First deriative of MC+ penalty
mcp.deriv <- function(x, lambda, gamma = 2) 
     {
	x2 <- abs(x)
	help1 <- sapply(x2, function(theta) { lambda*(1 - theta/(gamma*(lambda+1e-8))) } )
	out <- help1*(x2 < gamma*lambda) + 0*(x2 >= gamma*lambda)
	return(out)
	}
	

## Create a negative binomial family


#' A negative binomial family
#' 
#' Since the negative binomial is not a family in base \code{R}, an
#' \code{nb2()} family has been created which establishes the negative binomial
#' as a family for use in the main \code{rpql} function. Only the log link is
#' available at the moment, with the variance parameterized as \eqn{V = \mu +
#' \phi\mu^2} where \eqn{\phi} is the overdispersion parameter.
#' 
#' Used in the form \code{rpql(y, ..., family = nb2(), ...)}.
#' 
#' @return An object of class "family"
#' @author
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_author(\"#1\")}",
#' "rpql")\Sexpr{tools:::Rd_package_author("rpql")}
#' 
#' Maintainer:
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_maintainer(\"#1\")}",
#' "rpql")\Sexpr{tools:::Rd_package_maintainer("rpql")}
#' @seealso \code{\link{negative.binomial}} in the \code{MASS} package for
#' another example of a negative.binomial family.
#' @examples
#' 
#' \dontrun{
#' ## The function is currently defined as follows
#' nb2 <- function () {
#'     link <- "log"
#'     linkfun <- function(mu) log(mu)
#'     linkinv <- function(eta) pmax(exp(eta), .Machine$double.eps)
#'     mu.eta <- function(eta) pmax(exp(eta), .Machine$double.eps)
#'     variance <- function(mu, phi) mu + phi * mu^2
#'     valideta <- function(eta) TRUE
#'     validmu <- function(mu) all(mu > 0)
#'     structure(list(family = "negative.binomial", link = "log", 
#'         linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, 
#'         variance = variance, valideta = valideta, validmu = validmu, 
#'         name = link), class = "family")
#'   }
#' }  
#' 
nb2 <- function() 
     {
	link <- "log"
	linkfun <- function(mu) log(mu)
	linkinv <- function(eta) pmax(exp(eta), .Machine$double.eps)
	mu.eta <- function(eta) pmax(exp(eta), .Machine$double.eps)
	variance <- function(mu,phi) mu+phi*mu^2
	valideta <- function(eta) TRUE
	validmu <- function(mu) all(mu > 0)
  
	structure(list(family = "negative.binomial", link = "log", linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, variance = variance, valideta = valideta, validmu = validmu, name = link), class = "family")
	}


	    
## First deriative of SCAD penalty
scad.deriv <- function(x, lambda, a = 3.7) {
	x2 <- abs(x)
	help1 <- sapply(x2, function(theta) { max(a*lambda-theta,0)/((a-1)*(lambda+1e-8)) } )
	out <- lambda*((x2 <= lambda) + help1*(x2 > lambda))
	return(out)
	}


	
	
