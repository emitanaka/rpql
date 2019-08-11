############
## Regularized PQL (with a penalty on fixed effects and a group penalty on the random effects
############

## Version 0.8 (HELP FILES TO BE WRITTEN!)
## - Fixed issue with trial.size > 1 seems to run into issues
## - A new control argument lasso.lambda.scale is used, which scales the penalty on the group random effects (by 1/sqrt(# of clusters) so that is on the same scale as the fixed effects penalty; Thanks to Sarbesh Pandeya for pointing this out! Needs to be tested properly though...


## TODO: 1) consider using threshold operators, including LLA and threshold; 
## 2) Very slow for Gaussian responses!
##############
# library(mvtnorm)
# library(Matrix)
# library(MASS)
# library(lme4)
# library(gamlss.dist)
# library(Rcpp)
# library(inline)
# 
# source("rpql08/R/hiddenfunctions.R")
# source("rpql08/R/auxilaryfunctions.R")
# sourceCpp(file = "rpql08/src/RcppExports.cpp")
# sourceCpp(file = "rpql08/src/newcoeffn.cpp")


## pen.weights is a list with up to two elements for adl; first contains weights for fixed, and second contains a list of weights for random. 
# y = simy$y; X = simy$X; Z = simy$Z; id = simy$id; family = binomial(); lambda = 0.1; pen.type = "mcp"; hybrid.est = FALSE; offset = NULL; trial.size = 1; pen.weights = NULL; cov.groups = NULL; trace = TRUE; control$restarts <- 10; control$seed = NULL; intercept = TRUE; save.data = FALSE; tol = 1e-4; control$maxit = 100; start = NULL



#' Joint effects selection in GLMMs using regularized PQL.
#' 
#' \code{rpql} offers fast joint selection of fixed and random effects in
#' Generalized Linear Mixed Model (GLMMs) via regularization. The penalized
#' quasi-likelihood (PQL) is used as a loss function, and penalties are added
#' on to perform fixed and random effects selection. This method of joint
#' selection in GLMMs, referred to regularized PQL, is fast compared to
#' information criterion and hypothesis testing (Hui et al., 2016).
#' 
#' Please note \code{rpql} is the core workshops function that performed
#' regularized PQL on a single set of tuning parameters. \code{rpqlseq} is a
#' wrapper to permit a sequence of tuning parameter values. The latter is often
#' what users want to use.
#' 
#' \bold{Intro}
#' 
#' Generalized Linear Mixed Models (GLMMs) are an extension of Generalized
#' Linear Models (GLM, see the \code{glm} function) to include one or more sets
#' of random effects. For \eqn{i = 1,\ldots,n}, where \eqn{n} is the length of
#' \code{y}, we have
#' 
#' \deqn{g(\mu_{i}) = \bm{x}^T_i \bm{\beta} + \bm{z}^T_{i1} \bm{b}_{i1} +
#' \bm{z}^T_{i2} \bm{b}_{i2} + \ldots,}
#' 
#' where \eqn{g(\cdot)} is the link function, \eqn{\mu_i} is the mean of the
#' distribution for observation \eqn{i}, \eqn{\bm{x}_i} is row \eqn{i} of the
#' fixed effects model matrix \code{X}, and \eqn{\bm{\beta}} is the fixed
#' effects coefficients. For the random effects, \eqn{\bm{z}_{i1}} is row
#' \eqn{i} of the random effects model matrix in the first element of \code{Z},
#' \eqn{\bm{z}_{i2}} is from the second element of \code{Z} and so forth. The
#' random effects \eqn{\bm{b}_{i1}}, \eqn{\bm{b}_{i2}, \ldots} are drawn from a
#' multivariate normal distribution with mean zero and differing covariance
#' matrices \eqn{\bm{D}_1, \bm{D}_2, \ldots}.
#' 
#' Note that having lists for \code{id, Z}, allows for multiple sets of random
#' effects to be included in the GLMM. This is analogous to the \code{lme4}
#' package, where multiple random effects are permitted in the formula e.g.,
#' \code{(1|creek) + (1|creek:sample)}. If the GLMM contains only one set of
#' random effects, e.g., in longitudinal data, then the two lists will all
#' contain only one element. Cases where multiple sets of random effects may be
#' used include nested and crossed designs, in which case \code{id, Z}, will
#' have two or more elements. It is recommended that the user think through and
#' design these lists carefully to ensure that they are actually constructing
#' the appropriate GLMM of interest. Yes it takes some getting use too, and we
#' apologize for this =( Please see examples below for some ideas.
#' 
#' %While they are a flexible class of models for handling non-normal,
#' correlated data (see for instance Bolker et al., 2009, and Verbeke and
#' Molenberghs, 2009), estimating and performing variable selection in GLMMs
#' presents some major challenges. In the former, the marginal log-likelihood
#' requires integration over the unobserved random effects, and this is
#' computationally challenging. The \code{lme4} package, for instance,
#' typically uses the Laplace's approximation to approximate the integral. In
#' the latter, one may be interested in performing joint selection over the
#' fixed and random effects selection. Methods like hypothesis testing using
#' \code{anova} and information criterion may become computationally burdensome
#' in such case, and also may not necessarily perform that well (see Mueller et
#' al., 2013, for a review of variable selection in linear mixed models).
#' 
#' \bold{Regularized PQL}
#' 
#' Regularized PQL is designed as a fast approach to joint selection to GLMMs
#' (Hui et al., 2016). It works by taking the penalized quasi-likelihood (PQL,
#' Breslow and Clayton, 1993) and adding on penalties to perform selection of
#' the fixed and random effects. That is, maximize the regularized PQL function
#' 
#' \deqn{\ell = \sum\limits_{i=1}^n \log(f(y_i | \bm{\beta}, \bm{b}_{i1},
#' \bm{b}_{i2}, \ldots)) - \frac{1}{2} \sum\limits_{i=1}^n
#' \bm{b}^T_{i1}\bm{D}^{-1}_1 \bm{b}_{i1} - \frac{1}{2} \sum\limits_{i=1}^n
#' \bm{b}^T_{i2}\bm{D}^{-1}_2 \bm{b}_{i2} - \ldots - P_{\lambda}}
#' 
#' where \eqn{P_{\lambda}} denotes penalties to shrink the fixed effect
#' \eqn{\bm{\beta}} and random effect \eqn{\bm{b}_{i1}}, \eqn{\bm{b}_{i2},
#' \ldots} coefficients, which depend on one or more tuning parameters
#' \eqn{\lambda}. Like the PQL itself, regularized PQL is a fast approach for
#' estimating GLMMs because it treats the random effects as "fixed"
#' coefficients, and therefore no integration is required. Penalties are then
#' used to shrunk one or more \eqn{\bm{\beta}}'s and \eqn{\bm{b}}'s to zero,
#' the latter done so in a group-based manner, in order to perform joint
#' selection (see Hui et al., 2016, for details). In short, regularized PQL is
#' able to fit many GLMMs in a relatively short period of time, which in turn
#' facilitates the construction of a solution or regularization path ranging
#' from the null (intercept-only) to the full (saturated) model. A tuning
#' parameter selection method such as information criterion can then be used to
#' pick the select the final subset of fixed and random effects. A few penalty
#' types are available in the package, from which we prefer to use the adaptive
#' LASSO (with weights based on the full model, Zou, 2006) mainly because by
#' having weights, we can avoids have to search through a two-dimensional grid
#' of tuning parameter values.
#' 
#' Note that if one only wanted to penalize the fixed effects and leave the
#' random effects unpenalized, this can be achieved by setting the second
#' element/s of lambda equal to to e.g., \code{lambda = c(1,0)}. Note though
#' that in longitudinal studies, for covariates included as both fixed and
#' random effects, if the random effects is not penalized then neither should
#' the fixed effect. This ensures that no covariates end up being selected in
#' the model as a purely random effects (non-hierarchical shrinkage, Hui et
#' al., 2016). This can be accounted for also setting the corresponding
#' elements of \code{pen.weights$fixed} to zero.
#' 
#' \bold{AN IMPORTANT NOTE}
#' 
#' While regularized PQL is relatively fast, it will produce biased estimates
#' of the fixed and random effects parameters for non-normal responses,
#' especially if the amount of data to estimate each random effect is not large
#' e.g., if the number of time points or cluster size is not large. We envision
#' regularized PQL as a method of joint variable selection ONLY, and strongly
#' encourage the user to adopt a hybrid estimation approach (using
#' \code{hybrid.est = TRUE}, for instance). That is, once model selection is
#' performed using regularized PQL, the final submodel should be re-estimated
#' using more exact methods like quadrature or MCMC.
#' 
#' Because regularized PQL treats the random effects as ``fixed" coefficients
#' and therefore penalizes these, then the random effects covariance matrices
#' \eqn{\bm{D}_1, \bm{D}_2, \ldots} are regarded more as nuisance parameters.
#' This is in contrast to traditional maximum likelihood estimation where the
#' random effect coefficients \eqn{\bm{b}_{i1}}, \eqn{\bm{b}_{i2}, \ldots} are
#' integrated over. As nuisance parameters, regularized PQL employs an
#' iterative estimator based on maximizing the Laplace-approximated marginal
#' log-likelihood, assuming all other parameters are fixed, for estimating the
#' covariance matrix \eqn{\bm{D}_1, \bm{D}_2, \ldots}. This iterative estimator
#' was used in Hui et al., (2016) for independent clustered data specifically.
#' When they are multiple sets of random effects, each covariance matrix is
#' estimated conditionally on all others i.e., the random effect coefficients
#' corresponding to all other random effects are held constant. This can be
#' thought of as employing a series of conditional Laplace approximations to
#' obtain updates for \eqn{\bm{D}_1, \bm{D}_2, \ldots}.
#' 
#' %As nuisance parameters, regularized PQL currently offers three possible
#' estimates of the \eqn{\bm{D}_1, \bm{D}_2, \ldots}. The default is a simple
#' sample covariance estimator based on the estimated random effects
#' coefficients, "\code{bb}". That is, given estimates of \eqn{\bm{b}_{i1}} for
#' all \eqn{i = 1,\dots,n} based on regularized PQL, the estimate of
#' \eqn{\bm{D}_1} is given by
#' 
#' %\deqn{\hat{\bm{D}}_1 = \frac{1}{n} \sum\limits_{i=1}^n \hat{\bm{b}}_{i1}
#' \hat{\bm{b}}^T_{i1},}
#' 
#' %and similarly for \eqn{\bm{D}_2, \ldots}. The estimate of the random
#' effects covariance matrix is consistent, provided the number of observations
#' available to estimate each of the \eqn{\bm{b}_{i1}}'s i.e., the cluster
#' size, grows with \eqn{n}. The two other options for estimating the
#' covariance matrix are: 1) a weighted sample covariance estimator, where the
#' weights are based on the number of observations available to estimate each
#' random effect coefficient i.e., cluster size, such that more weight is given
#' to coefficients estimated with a larger cluster size. This estimator is
#' useful when you have unbalanced cluster sizes, which is almost always the
#' case in real life, and 2) an iterative estimator based on maximizing the
#' Laplace-approximated marginal log-likelihood assuming all other parameters
#' are fixed. This iterative estimator was used in Hui et al., (2016), but is
#' only applicable when you have independent clustered data e.g., longitudinal
#' studies.
#' 
#' \bold{A not so short discussion about information criterion}
#' 
#' How to choose the tuning parameters for penalized regression is an active
#' area of area of research in statistics (see for instance Zhang et al., 2010,
#' Hui et al., 2014), with the most popular solutions being cross validation
#' and information criteria. That is, a solution path is constructed and the
#' best submodel is then chosen by minimizing the value of the information
#' criterion. Anyway, \code{rpql} offers the following information criteria for
#' tuning parameter selection, as available in \code{ics} in the output. Please
#' note all of the criteria below use only the first part of the PQL function
#' as the loss function i.e., \eqn{IC = -2\sum\limits_{i=1}^n \log(f(y_i |
#' \bm{\beta}, \bm{b}_{i1}, \bm{b}_{i2}, \ldots)) +} model complexity terms.
#' 
#' \enumerate{ \item A AIC-type criterion that penalizes a values of 2 for
#' every non-zero fixed effect coefficient, and, for each set of random
#' effects, penalizes a value of 2 for every non-zero random effect coefficient
#' in that set.
#' 
#' \item A BIC-type criterion that penalizes a value of \eqn{\log(n)} for every
#' non-zero fixed effect coefficient, and, for each set of random effects,
#' penalizes a value of \eqn{\log(n_c)} for every non-zero, unique element in
#' covariance matrix for that set, where \code{n_c} denotes the number of
#' clusters corresponding to that random effect.
#' 
#' \item A BIC-type criterion that penalizes a value of \eqn{\log(n)} for every
#' non-zero fixed effect coefficient, and, for each set of random effects,
#' penalizes a value of \eqn{\log(n)} for every non-zero, unique element in
#' covariance matrix for that set. This combination of penalties is the one
#' used in the package \code{lme4}.
#' 
#' \item Three hybrid information criteria that penalizes a value \eqn{\log(n)}
#' for every non-zero fixed effect coefficient, and, for each set of random
#' effects, penalizes a value of 2/1/0.5 for every non-zero random effect
#' coefficient in that set. }
#' 
#' Selection consistency for all but the first AIC criteria have been
#' established, although empirically performance may differ. We generally
#' prefer the three hybrid criterion, although it is recommended that the user
#' tries several of them and see how results differ!
#' 
#' @aliases rpql rpql.default print.rpql
#' @param y A vector of responses
#' @param X A model matrix corresponding to the fixed effects. It should have
#' the same number of rows as the length of \code{y}. An intercept column must
#' be included if a fixed intercept is desired.
#' @param Z A list with each element being a model matrix for a set of random
#' effects. Each element of \code{Z} is referenced by a vector of IDs given by
#' the corresponding element in the list \code{id}. Each model matrix (element
#' of \code{Z}) should have the same number of rows as the length of \code{y}.
#' @param id A list with each element being a vector of IDs that reference the
#' model matrix in the corresponding element in the list \code{Z}. Each vector
#' of IDs \emph{must} be integers (but not factors).
#' @param x An object for class "rpql".
#' @param family The distribution for the responses in GLMM. The argument must
#' be applied as a object of class "family". Currently supported arguments
#' include: \code{gaussian()}, \code{poisson()}, \code{binomial()},
#' \code{Gamma()}, \code{nb2()} for negative binomial, \code{LOGNO()} for
#' log-normal, and \code{ZIP()} for zero-inflated Poisson.
#' @param trial.size The trial size if \code{family = binomial()}. Either takes
#' a single non-zero value or a vector of non-zero values with length the same
#' as the number of rows in \code{X}. The latter allows for differing trial
#' sizes across responses. Defaults to 1.
#' @param lambda A vector of length one or two specifying the tuning parameters
#' used in regularized PQL. If two elements are supplied, then first and second
#' elements are for the fixed and random effects penalty respectively. If one
#' element, then it is applied to both penalties.
#' @param pen.type A vector of one or two strings, specifying the penalty used
#' for variable selection. If two elements are supplied, then first and second
#' strings are the fixed and random effects penalty respectively. If one
#' element, the same type of penalty is used. Currently supported argument
#' include: "\code{lasso}" for standard lasso (Tibshirani, 1996), "\code{scad}"
#' for SCAD penalty with \eqn{a} controlled by \code{control$scad.a} (Fan and
#' Li, 2001), "\code{adl}" for adaptive lasso (Zou, 06), "\code{mcp}" for MC+
#' penalty with \eqn{\gamma} controlled by controlled by
#' \code{control$mcp.gamma} (Zhang, 2010). If the adaptive lasso is used, then
#' \code{pen.weights} must also be supplied. Defaults to standard lasso penalty
#' for both fixed and random effects.
#' @param start A list of starting values. It must contain the following
#' elements: \code{start$fixef} as starting values for the fixed effect
#' coefficients, \code{start$ranef} which is a list containing matrices of
#' starting values for the random effects coefficients. It may also contain
#' \code{start$D} which is a list of matrices to act as starting values for
#' random effects covariance matrices.
#' @param cov.groups A vector specifying if the columns of \code{X} (including
#' the intercept) should be regarded and therefore penalized in groups. For
#' example, if one or more of the fixed effect covariates are factors, then
#' \code{lme4} will automatically create dummy variables in the model matrix
#' and estimate coefficients for each level, using one level as the reference.
#' \code{cov.groups} is then used to identify all the coefficients that
#' corresponds to that factor, such that all of these coefficients are
#' penalized collectively as a group. Defaults to NULL, in which case it is
#' assumed all coefficients should be treated independently. Please see the
#' details and examples for more details.
#' @param pen.weights A list containing up to two elements for additional
#' (adaptive lasso) weights to be included for penalization. This must be
#' supplied if \code{pen.type} has one or both elements set to "adl", otherwise
#' it is optional. A weights equal to zero implies no penalization is applied
#' to the parameter. The two elements in the list are as follows: for fixed
#' effects, \code{pen.type$fixed} should be a vector with length equal to the
#' number of columns in \code{X}. For random effects, \code{pen.weights$ran}
#' should be a list of the same length as the list \code{Z}, where each element
#' in that list is a vector with length equal to the number of columns in the
#' corresponding element of the list \code{Z} (recall that each element of
#' \code{Z} is a model matrix). Defaults to NULL, in which case there are no
#' weights involved in the penalization.
#' @param hybrid.est Should a hybrid estimation approach be used? That is, once
#' model selection is performed using regularized PQL, should the submodel be
#' re-estimated using the \code{lme4} package, if possible? Defaults to FALSE.
#' @param offset This can be used to specify an \emph{a priori} known component
#' to be included in the linear predictor during fitting. It should be numeric
#' vector of length equal to \code{y}. Defaults to NULL.
#' @param intercept Is one of the columns of \code{X} an intercept term? This
#' is used to indicate the presence of a fixed intercept in the model, which
#' subsequently will NOT be penalized. Defaults to TRUE.
#' @param save.data Should \code{y, X}, and \code{Z}, be saved as part of the
#' output? Defaults to FALSE. The data is not saved by default in order to save
#' memory.
#' @param control A list controlling the finer details of the rPQL algorithm.
#' These include: \itemize{ \item\code{tol}: Tolerance value for convergence in
#' the regularized PQL to be declared, where convergence is measured as the
#' difference between the estimated parameters in successive iterations.
#' Defaults to a value of 1e-4.
#' 
#' \item\code{maxit}: The maximum number of update iterations for regularized
#' PQL. Defaults to 100.
#' 
#' \item\code{trace}: Should the update estimates of the fixed effect
#' coefficients and that random effect covariance matrices be printed at each
#' iteration? Defaults to FALSE.
#' 
#' \item\code{restarts}: The number of restarts to try in case the algorithm
#' diverges, i.e. the fixed effect coefficients and /or random effects
#' covariance matrices "blow up". Defaults to a value of 5. Divergence is
#' mostly likely to occur when you have count responses with some extremely
#' large counts in there, in which regularized PQL can throw a hissy fit.
#' 
#' \item\code{scad.a, mcp.gamma}: Controls the \eqn{a} and \eqn{\gamma}
#' parameters in the SCAD and MC+ penalty respectively. Defaults to \eqn{a =
#' 3.7} (Fan and Li, 2001) and \eqn{\gamma = 2} (Zhang, 2010) respectively.
#' Please note these parameters are only in use when \code{pen.type} involves
#' these penalties.
#' 
#' \item\code{seed}: A seed that can be used if results need to be replicated.
#' Defaults to NULL, in which case a random seed is used.  }
#' @param ... Not used.
#' @return An object of class "rpql" containing the following elements:
#' \item{call}{The matched call.}
#' 
#' \item{fixef}{A vector of estimated fixed effect coefficients, \eqn{\beta}.}
#' 
#' \item{ranef}{A list with each element being a matrix of estimated
#' (predicted) random effect coefficients, \eqn{\bm{b}_{i1}},
#' \eqn{\bm{b}_{i2}}, and so on.}
#' 
#' \item{ran.cov}{A list with each element being an estimated random effect
#' covariance matrices, \eqn{\bm{D}_1, \bm{D}_2, \ldots}.}
#' 
#' \item{logLik}{The (unpenalized) PQL likelihood value at convergence.}
#' 
#' \item{phi, shape, zeroprob}{Estimates of nuisance parameters (if
#' appropriate), including the variance and overdispersion parameter for
#' normal, lognormal and negative binomial families, the shape parameter for
#' the Gamma family, and the probability of a structural zero for zero-inflated
#' Poisson family.}
#' 
#' \item{family}{The family fitted.}
#' 
#' \item{n}{The length of \code{y}.}
#' 
#' \item{id}{The \code{id} argument.}
#' 
#' \item{lambda, pen.type}{The tuning parameters and penalties used.}
#' 
#' \item{ics}{A vector containing the number of estimated parameters in the
#' GLMM (note regularized PQL treats the random effects as "fixed"), and some
#' information criteria. Please see \code{details} above for more information.}
#' 
#' \item{nonzero.fixef}{A vector indexing which of the estimated fixed effect
#' coefficients are non-zero.}
#' 
#' \item{nonzero.ranef}{A list with each element being a vector indexing which
#' of the estimated random effects are non-zero, i.e. which of the diagonal
#' elements in the corresponding element of \code{ran.cov} are non-zero.}
#' 
#' \item{hybrid}{The estimated fit from \code{lme4}, if \code{hybrid.est =
#' TRUE}.}
#' 
#' \item{y,X,Z}{The data the GLMM is fitted to, if \code{save.data = TRUE}.}
#' @section Warnings: \itemize{ \item We strongly recommend you scale your
#' responses (if normally distributed) and any continuous covariates, otherwise
#' \code{rpql} like all penalized likelihood methods, may not make much sense!
#' 
#' \item Like its standard unpenalized counterpart, regularized PQL can produce
#' very bias parameter estimates in finite samples, especially if you do not
#' have a lot of data to estimate each random effect. We therefore envision
#' regularized PQL as a tool for fast model selection in GLMMs, and strongly
#' recommend you re-estimate the final submodel using more accurate estimation
#' methods i.e., use a hybrid estimation approach, in order to obtain better
#' final parameter estimates and predictions of the random effects.
#' 
#' \item If \code{save.data = TRUE}, the data you fitted the GLMM is also saved
#' as part of the output, and this can potentially take up a lot of memory.
#' 
#' \item If you are constantly suffering convergence issues with regularized
#' PQL, even after multiple restarts, consider increasing \code{lambda[2]} to
#' penalized the random effects more and stabilize the estimation algorithm.
#' You may also want to consider better starting values, in particular, smaller
#' values of \code{start$ranef}. Good luck!
#' 
#' }
#' @author
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_author(\"#1\")}",
#' "rpql")\Sexpr{tools:::Rd_package_author("rpql")}
#' 
#' Maintainer:
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_maintainer(\"#1\")}",
#' "rpql")\Sexpr{tools:::Rd_package_maintainer("rpql")}
#' @seealso \code{\link{rpqlseq}} for the wrapper function that runs
#' \code{rpql} multiple times on a sequence of tuning parameter values,
#' \code{\link{build.start.fit}} for building \code{start} lists from a GLMM
#' fitted using the \code{lme4} package, \code{\link{summary}} for a summary of
#' the regularized PQL fit. For alternative methods of fitting GLMMs, you may
#' also want be check out the packages \code{lme4, nlme, MCMCglmm} and
#' \code{glmmADMB}.
#' @references \itemize{ \item Breslow, N. E., and Clayton, D. G. (1993).
#' Approximate inference in generalized linear mixed models. Journal of the
#' American Statistical Association, 88, 9-25.
#' 
#' \item Fan, J., and Li, R. (2001). Variable selection via nonconcave
#' penalized likelihood and its oracle properties. Journal of the American
#' statistical Association, 96, 1348-1360.
#' 
#' \item Hui, F.K.C., Mueller, S., and Welsh, A.H. (2017). Joint Selection in
#' Mixed Models using Regularized PQL. Journal of the American Statistical
#' Association, 112, 1323-1333.
#' 
#' \item Hui, F.K.C., Mueller, S., and Welsh, A.H. (2017). Hierarchical
#' Selection of Fixed and Random Effects in Generalized Linear Mixed Models.
#' Statistica Sinica, 27, 501-518.
#' 
#' \item Hui, F. K., Warton, D. I., and Foster, S. D. (2014). Tuning parameter
#' selection for the adaptive lasso using ERIC. Journal of the American
#' Statistical Association, 110, 262-269.
#' 
#' \item Lin, X., and Breslow, N. E. (1996). Bias correction in generalized
#' linear mixed models with multiple components of dispersion. Journal of the
#' American Statistical Association, 91, 1007-1016.
#' 
#' \item Mueller, S., Scealy, J. L., and Welsh, A. H. (2013). Model selection
#' in linear mixed models. Statistical Science, 28, 135-167.
#' 
#' \item Tibshirani, R. (1996). Regression shrinkage and selection via the
#' lasso. Journal of the Royal Statistical Society. Series B (Methodological),
#' 58, 267-288.
#' 
#' \item Zhang, Y., Li, R., and Tsai, C. L. (2010). Regularization parameter
#' selections via generalized information criterion. Journal of the American
#' Statistical Association, 105, 312-323.
#' 
#' \item Zhang, C. H. (2010). Nearly unbiased variable selection under minimax
#' concave penalty. The Annals of Statistics, 38, 894-942.
#' 
#' \item Zou, H. (2006). The adaptive lasso and its oracle properties. Journal
#' of the American statistical association, 101, 1418-1429.  }
#' @examples
#' 
#' ## Please note all examples below use the \code{rpqlseq} wrapper function. 
#' 
#' library(lme4)
#' library(gamlss.dist)
#' 
#' ##################
#' ## Example 1: Poisson GLMM on simulated data 
#' ## Indepenent cluster model with 30 clusters and equal cluster sizes of 10
#' ## 9 fixed and random effect covariates including a fixed and random intercept
#' library(mvtnorm)
#' set.seed(1)
#' n <- 30; m <- 10; p <- 8; 
#' ## Generate rows of a model matrix from a multivariate normal distribution 
#' ## with AR1 covariance structure. 
#' 
#' H <- abs(outer(1:p, 1:p, "-")) 
#' X <- cbind(1,rmvnorm(n*m,rep(0,p),sigma=0.5^H)); 
#' Z <- X 
#' true_betas <- c(0.1,1,-1,-1,1,rep(0,p-4)) ## 5 truly important fixed effects
#' true_D <- matrix(0,ncol(Z),ncol(Z))
#' true_D[1:3,1:3] <- matrix(c(1,0.6,0.6,0.6,1,0.4,0.6,0.4,1),3,3,byrow=TRUE) 
#' ## 3 important random effects
#' 
#' simy <- gendat.glmm(id = list(cluster=rep(1:n,each=m)), X = X, beta = true_betas, 
#' 	Z = list(cluster=Z), D = list(cluster=true_D), family = poisson()) 
#' 
#' 	
#' \dontrun{
#' ## Construct a solution path using adaptive LASSO for selection 
#' dat <- data.frame(y = simy$y, simy$X, simy$Z$cluster, simy$id)
#' fit_satlme4 <- glmer(y ~ X - 1 + (Z - 1 | cluster), data = dat,
#' 	family = "poisson")
#' fit_sat <- build.start.fit(fit_satlme4, gamma = 2)
#' ## Please see example 3 for another way of constructing the adaptive weights
#' 
#' lambda_seq <- lseq(1e-6,1,length=100)
#' fit <- rpqlseq(y = simy$y, X = simy$X, Z = simy$Z, id = simy$id, 
#' 	family = poisson(), lambda = lambda_seq, pen.type = "adl", 
#' 	pen.weights = fit_sat$pen.weights, start = fit_sat)
#' 
#' summary(fit$best.fit[[3]])  
#' # apply(fit$collect.ics, 2, which.min) ## Look at best fit chosen by different ICs
#' 
#' ## Note, if you wanted to penalized the fixed effects only, this can achieved
#' ## by setting fit_sat$pen.weights$random$cluster <- rep(0,ncol(simy$Z$cluster))
#' 
#' ## An alternative way to construct the X and Z matrices for input into rpqlseq is as follows:
#' XMM <- unname(model.matrix(fit_satlme4)) 
#' ZMM <- getME(fit_satlme4,"mmList"); names(ZMM) <- "cluster"
#' lambda_seq <- lseq(1e-6,1,length=100)
#' fit <- rpqlseq(y = simy$y, X = XMM, Z = ZMM, id = simy$id, 
#'  	family = poisson(), lambda = lambda_seq, pen.type = "adl", 
#' 	pen.weights = fit_sat$pen.weights, start = fit_sat)
#' ## Big thanks for Andrew Olney for this suggestion!
#' }
#' 
#' 
#' ##################
#' ## Example 2: Similar to example 1 but with Bernoulli GLMMs 
#' ## 30 clusters, cluster size of 20
#' library(mvtnorm)
#' set.seed(1)
#' n <- 30; m <- 20; p <- 8; 
#' ## Generate rows of a model matrix from a multivariate normal distribution 
#' ## with AR1 covariance structure. 
#' 
#' H <- abs(outer(1:p, 1:p, "-")) 
#' X <- cbind(1,rmvnorm(n*m,rep(0,p),sigma=0.5^H)); 
#' Z <- X 
#' true_betas <- c(-0.1,1,-1,1,-1,rep(0,p-4)) ## 5 truly important fixed effects
#' true_D <- matrix(0,ncol(Z),ncol(Z))
#' true_D[1:3,1:3] <- diag(c(3,2,1), nrow = 3)
#' ## 3 important random effects
#' 
#' simy <- gendat.glmm(id = list(cluster=rep(1:n,each=m)), X = X, 
#'   beta = true_betas, Z = list(cluster=Z), D = list(cluster=true_D), family = binomial()) 
#' 
#' 	
#' \dontrun{
#' ## Construct a solution path using adaptive LASSO for selection 
#' dat <- data.frame(y = simy$y, simy$X, simy$Z$cluster, simy$id)
#' fit_satlme4 <- glmer(y ~ X - 1 + (Z - 1 | cluster), data = dat, 
#' 	family = "binomial")
#' fit_sat <- build.start.fit(fit_satlme4, gamma = 2)
#' 
#' lambda_seq <- lseq(1e-6,1,length=100)
#' best.fit <- list(ics = rep(Inf,6))
#' fit <- rpqlseq(y = simy$y, X = simy$X, Z = simy$Z, id = simy$id, 
#' 	family = binomial(), lambda = lambda_seq, pen.type = "adl", 
#' 	pen.weights = fit_sat$pen.weights, start = fit_sat)
#' 	
#' summary(fit$best.fit[[3]])  
#' # apply(fit$collect.ics, 2, which.min) ## Look at best fit chosen by different ICs
#' 
#' ## An alternative way to construct the X and Z matrices for input into rpqlseq is as follows:
#' XMM <- unname(model.matrix(fit_satlme4)) 
#' ZMM <- getME(fit_satlme4,"mmList"); names(ZMM) <- "cluster"
#' lambda_seq <- lseq(1e-6,1,length=100)
#' fit <- rpqlseq(y = simy$y, X = XMM, Z = ZMM, id = simy$id, 
#'  	family = binomial(), lambda = lambda_seq, pen.type = "adl", 
#' 	pen.weights = fit_sat$pen.weights, start = fit_sat)
#' }
#' 
#' 
#' ##################
#' ## Example 3: Bernoulli GLMMs on simulated data
#' ## Nested data with 200 observations in total: split into 10 creeks, 
#' ## 5 samples nested within each creek
#' ## Please see example in gendat.glmm for further details
#' mn <- 200; 
#' X <- matrix(1,mn,1); 
#' ids <- list(samples = rep(1:50,each=4), creek = rep(1:10,each=20)) 
#' ## We have two sets of random intercepts only, one for creek and one for 
#' ## samples nested within creek.
#' Zs <- list(samples = X, creek = X) 
#' 
#' true_betas <- -0.1 
#' true_D <- list(samples = as.matrix(0.001), creek = as.matrix(1)) 
#' ## Please ensure each element of true_D is a matrix
#' 
#' simy <- gendat.glmm(id = ids, X = X, beta = true_betas, Z = Zs, 
#' 	D = true_D, trial.size = 1, family = binomial())
#' 
#' \dontrun{
#' ## Construct a solution path use adaptive LASSO for selection
#' ## Here is another way of constructing the adaptive weights:
#' ## Use the fact that rpql can do a final fit based on maximum likelihood
#' ## to obtain a good saturated fit.
#' fit_sat <- rpql(y = simy$y, X = simy$X, Z = simy$Z, id = simy$id, 
#' 	family = binomial(), lambda = 0, hybrid = TRUE)
#' fit_sat <- build.start.fit(fit_sat$hybrid, gamma = 2)
#' 	
#' fit <- rpqlseq(y = simy$y, X = simy$X, Z = simy$Z, id = simy$id, 
#' 	family = binomial(), lambda = lambda_seq, pen.type = "adl", 
#' 	pen.weights = fit_sat$pen.weights, start = fit_sat)
#' 
#' summary(fit$best.fit[[3]])  
#' # apply(fit$collect.ics, 2, which.min) ## Look at best fit chosen by different ICs
#' }
#'    
#' 
#' ##################
#' ## Example 4: Linear mixed models on Alfalfa split-plot data
#' 
#' \dontrun{
#' 
#' library(nlme)
#' data(Alfalfa)
#' Alfalfa$Yield <- scale(Alfalfa$Yield)
#' X <- as.matrix(model.matrix(~ Date, data = Alfalfa)) 
#' ## Note Date is categorical variable!
#' colnames(X)[1] <- "x1"
#' Z <- list(BlockVariety = matrix(1,nrow(X),1), Block = matrix(1,nrow(X),1))
#' ## Four samples of each Block*Variety
#' ids <- list(BlockVariety = rep(1:(nrow(X)/4),each=4), 
#' 	Block = as.numeric(Alfalfa$Block)) 
#' 
#' ## How you would fit it in lme4
#' fit_satlme4 <- lmer(Yield ~ X - 1 + (1|Block/Variety), data = Alfalfa)
#' fit_sat <- build.start.fit(fit_satlme4, cov.groups = c(1,2,2,2), gamma = 2)
#' 
#' ## Construct a solution path using adaptive LASSO for selection
#' lambda_seq <- lseq(1e-5,2,length=100)
#' fit <- rpqlseq(y = Alfalfa$Yield, X = X, Z = Z, id = ids, 
#' 	lambda = lambda_seq, cov.groups = c(1,2,2,2), pen.type = "adl", 
#' 	pen.weights = fit_sat$pen.weights, start = fit_sat)
#' 
#' summary(fit$best.fit[[3]])  
#' # apply(fit$collect.ics, 2, which.min) ## Look at best fit chosen by different ICs
#' 
#' 
#' ## An alternative way to construct the X and Z matrices for input into rpqlseq is as follows:
#' X <- unname(model.matrix(fit_satlme4)) 
#' Z <- getME(fit_satlme4, "mmList"); names(Z) <- c("BlockVariety", "Block")
#' lambda_seq <- lseq(1e-6,1,length=100)
#' fit <- rpqlseq(y = Alfalfa$Yield, X = X, Z = Z, id = ids, 
#' 	lambda = lambda_seq, cov.groups = c(1,2,2,2), pen.type = "adl", 
#' 	pen.weights = fit_sat$pen.weights, start = fit_sat)
#' }
#' 
#' 
#' ##################
#' ## Example 5: Linear mixed models on sleep study dataset
#' 
#' \dontrun{
#' data(sleepstudy)
#' 
#' ## How you fit it in lme4
#' ## Response is scaled so as to avoid large variances and easier intepretation
#' sleepstudy$Reaction <- scale(sleepstudy$Reaction) 
#' sleepstudy$Days <- scale(sleepstudy$Days)
#' fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
#' 
#' ## How you fit it using rpql
#' ## Construct a solution path using adaptive LASSO for selection 
#' X <- cbind(1,sleepstudy$Days)
#' Z <- list(subject = X)
#' ids <- list(subject = as.numeric(sleepstudy$Subject))
#' fit_sat <- build.start.fit(fm1, gamma = 2)
#' 
#' lambda_seq <- lseq(1e-4,1,length=100)
#' fit <- rpqlseq(y = sleepstudy$Reaction, X = X, Z = Z, id = ids, 
#' 	lambda = lambda_seq, pen.type = "adl", 
#' 	pen.weights = fit_sat$pen.weights, start = fit_sat)
#' 
#' summary(fit$best.fit[[3]])  
#' # apply(fit$collect.ics, 2, which.min) ## Look at best fit chosen by different ICs
#' ## Best fit might well be the saturated fit! 
#' ## This is at least consistent with confint(fm1)
#' 
#' ## An alternative way to construct the X and Z matrices for input into rpqlseq is as follows:
#' X <- unname(model.matrix(fm1)) 
#' Z <- getME(fm1, "mmList"); names(Z) <- "subject"
#' lambda_seq <- lseq(1e-6,1,length=100)
#' fit <- rpqlseq(y = sleepstudy$Reaction, X = X, Z = Z, id = ids, 
#' 	lambda = lambda_seq, pen.type = "adl", 
#' 	pen.weights = fit_sat$pen.weights, start = fit_sat)
#' }
#' 
#' 
#' ##################
#' ## Example 6: GLMM with lognormal responses
#' ## Fixed effects selection only
#' 
#' \dontrun{
#' n <- 50; m <- 10; p <- 8; 
#' H <- abs(outer(1:p, 1:p, "-")) 
#' X <- cbind(1,rmvnorm(n*m,rep(0,p),sigma=0.5^H)); 
#' Z <- X[,1:3] ## 3 random effects all of which important
#' true_betas <- c(0.1,1,-1,-1,1,rep(0,p-4)) ## 5 important fixed effects
#' true_D <- matrix(0,ncol(Z),ncol(Z))
#' true_D[1:3,1:3] <- matrix(c(1,0.6,0.6,0.6,1,0.4,0.6,0.4,1),3,3,byrow=TRUE) 
#' 
#' simy <- gendat.glmm(id = list(cluster=rep(1:n,each=m)), X = X, 
#' 	beta = true_betas, Z = list(cluster=Z), D = list(cluster=true_D), 
#' 	family = LOGNO(), phi = 1) 
#' 
#' ## We will use the lasso penalty for fixed effects only with no weights
#' ## Note lognormal mixed models are usually hard to fit by maximum likelihood in R!
#' ## Hence adaptive weights are sightly hard to obtain
#' 
#' ## Note also that since random effects are not penalized, then generally 
#' ## the corresponding fixed effect covariates should not be penalized 
#' ## (at least in longitudinal studies), in keeping in line with the 
#' ## hierarchical principle of the effects.
#' ## To account for this in the above, we can use the pen.weights argument 
#' ## to prevent penalization of the first three fixed effect covariates
#' 
#' fit <- rpqlseq(y = simy$y, X = simy$X, Z = simy$Z, id = simy$id, 
#'   family = LOGNO(), lambda = lambda_seq, pen.type = "lasso", start = NULL, 
#' 	pen.weights = list(fixed = rep(c(0,1), c(3,ncol(X)-3))))
#' 
#' summary(fit$best.fit[[3]])  
#' # apply(fit$collect.ics, 2, which.min) ## Look at best fit chosen by different ICs
#' }
#' 
#' 
rpql <- function(y, ...) 
     UseMethod("rpql")




#' Wrapper function for joint effects selection in GLMMs using regularized PQL.
#' 
#' \code{rpql} offers fast joint selection of fixed and random effects in
#' Generalized Linear Mixed Model (GLMMs) via regularization. The penalized
#' quasi-likelihood (PQL) is used as a loss function, and penalties are added
#' on to perform fixed and random effects selection. This method of joint
#' selection in GLMMs, referred to regularized PQL, is fast compared to
#' information criterion and hypothesis testing (Hui et al., 2016).
#' 
#' \code{rpqlseq} is a wrapper function to permit a sequence of tuning
#' parameter values, which wraps around the code workhorse function
#' \code{rpql}.
#' 
#' Please see the help file for \code{rpql} for details on how regularized PQL
#' works. \code{rpqlseq} is simply a wrapper function to run the core
#' \code{rpql} function multiple times, on a sequence of tuning parameter
#' values, in order to construct a regularization path. The best models, based
#' on different information criteria for selecting the best tuning parameter
#' (degree of sparsity) are then returned.
#' 
#' @param y,X,Z,id,family,trial.size As per the \code{rpql} function. Please
#' see the help file for \code{rpql} for details on the arguments.
#' @param lambda Either a vector containing \bold{sequence} of tuning parameter
#' values, which is applied to both penalties, or two-column matrix containing
#' a \bold{sequence} of tuning parameter values for the fixed and random
#' effects penalty respectively.
#' @param
#' pen.type,start,cov.groups,pen.weights,offset,intercept,save.data,control As
#' per the \code{rpql} function. Please see the help file for \code{rpql} for
#' details on the arguments.
#' @param ... Not used.
#' @return An object of class "rpql" containing the following elements:
#' \item{best.fits}{A list containing the best fitted models as based on
#' different information criteria used to select the tuning parameter. Each
#' element in this list has the same structure as the output from the
#' \code{rpql} function. Please see the \code{rpql} function for details on the
#' information criteria available as well as the nature of the output.}
#' 
#' \item{collect.ics}{A matrix containing the values of various information
#' criteria calculated for the sequence of \code{lambda} values supplied. The
#' best fitted models found in \code{best.fits} is based off this matrix i.e.,
#' each element in \code{best.fits} corresponds to a model that was chosen
#' based on minimizing the corresponding information criterion in
#' \code{collect.ics}. Please see the \code{rpql} function for details on the
#' information criteria available.}
#' 
#' \item{lambda}{The sequence of tuning parameters considered.}
#' @author
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_author(\"#1\")}",
#' "rpql")\Sexpr{tools:::Rd_package_author("rpql")}
#' 
#' Maintainer:
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_maintainer(\"#1\")}",
#' "rpql")\Sexpr{tools:::Rd_package_maintainer("rpql")}
#' @seealso \code{\link{rpql}}, which is the core workhorse function that
#' performed regularized PQL for a single set of tuning parameter values.
#' @examples
#' 
#' ## Please see examples in help file for the \code{rpql} function for usage.
#' 
rpqlseq <- function(y, X, Z, id, family = gaussian(), trial.size = 1, lambda, pen.type = "lasso", start = NULL, cov.groups = NULL, 
     pen.weights = NULL, offset = NULL, intercept = TRUE, save.data = FALSE, 
     control = list(tol = 1e-4, maxit = 100, trace = FALSE, restarts = 5, scad.a = 3.7, mcp.gamma = 2, lasso.lambda.scale = TRUE, seed = NULL), ...) 
     {
	
	lambda <- as.matrix(lambda)
	if(!control$lasso.lambda.scale)
          warnings("Scaling factor for the second tuning parameter has been turned off for lasso and adaptive lassp penalties. 
               This may cause issues as the group coefficient penalty on the random effects is on a different scale to the single 
               coefficient penalty on the fixed effects.")
	if(ncol(lambda) == 1) 
          {
          lambda[,1] <- sort(lambda[,1], decreasing = FALSE)
          }
          
	
	control <- fillin.control(control) 
	
	## Do a starting fit for formatting reasons
	start_fit <- rpql(y = y, X = X, Z = Z, id = id, family = family, lambda = lambda[1,], pen.type = "lasso", start = start) 
	collect_ics <- matrix(0, nrow = nrow(lambda), ncol = length(start_fit$ics)+1) ## First column in start_fit$ics is DoF
	colnames(collect_ics) <- c("PQL Likelihood", names(start_fit$ics)) 

	best_fits <- vector("list", length(start_fit$ics)-1); 
	names(best_fits) <- names(start_fit$ics)[-1]
	rm(start_fit)
			
	for(l1 in 1:nrow(lambda)) 
          {
		message("Onto ", l1)
		if(l1 == 1) 
               prev_fit <- start
               
          new_fit <- rpql(y = y, X = X, Z = Z, id = id, family = family, lambda = lambda[l1,], pen.type = pen.type, pen.weights = pen.weights, start = prev_fit, 
               cov.groups = cov.groups, hybrid.est = FALSE, offset = offset, intercept = intercept, save.data = save.data, control = control)
		
		if(l1 == 1)
            {
            for(k in 1:length(best_fits)) 
                best_fits[[k]] <- new_fit
            }

		if(l1 > 1)
               {
               for(l2 in 1:length(best_fits)) 
                    {
                    if(new_fit$ics[l2+1] < best_fits[[l2]]$ics[l2+1]) 
                         best_fits[[l2]] <- new_fit
                    }
               }
		prev_fit <- new_fit
		collect_ics[l1,] <- c(new_fit$logLik1, new_fit$ics)
		}
	
	out <- list(best.fits = best_fits, collect.ics = collect_ics, lambda = lambda)	
	}



#y = simy$y
#X = simy$X
#Z = simy$Z
#id = simy$id
#family = binomial()
#lambda = 0
#pen.type = "lasso"
#trial.size = 10
#trace = TRUE
#start = gensatfit
#hybrid.est = FALSE; offset = NULL; 
#pen.weights = NULL; 
#cov.groups = NULL; 
#control = list(tol = 1e-4, maxit = 100, trace = TRUE, restarts = 5, scad.a = 3.7, mcp.gamma = 2, seed = NULL); 
#intercept = TRUE; 
#save.data = FALSE; 


rpql.default <- function(y, X, Z, id, family = gaussian(), trial.size = 1, lambda, pen.type = "lasso", start = NULL, cov.groups = NULL, 
     pen.weights = NULL, hybrid.est = FALSE, offset = NULL, intercept = TRUE, save.data = FALSE, 
     control = list(tol = 1e-4, maxit = 100, trace = FALSE, restarts = 5, scad.a = 3.7, mcp.gamma = 2, lasso.lambda.scale = TRUE, seed = NULL), ...) 
     {

     conv.eps <- 1e-3
     control <- fillin.control(control) 
	
     y <- as.vector(y)
     X <- as.matrix(X)
     num_fixed <- ncol(X)
     if(is.null(colnames(X))) 
          colnames(X) <- paste0("x",1:ncol(X)) 
     for(k in 1:length(Z)) 
          colnames(Z[[k]]) <- paste0("z",k,1:ncol(Z[[k]]))
     if(det(crossprod(X)) == 0) 
          warning("Is X rank-deficient? Please double check!", immediate. = TRUE) 
     if(!is.list(Z) || !is.list(id)) 
          stop("Please supply id and Z as lists. The element in the two lists correspond to a vector of IDs and the random effects model matrix for those IDs respectively.")
     if(sum(duplicated(names(id))) > 0 || sum(duplicated(names(Z))) > 0) 
          stop("Please ensure the names of the elements in the lists id (and Z) are unique.")
     if(length(Z) != length(id)) 
          stop("The number of elements in Z and id should be the same. Thanks.")
# 	if(any(apply(X,2,sd) != 1) || any(apply(X,2,sd) != 1))
# 		message("You may want to consider standardizing your X and Z matrix prior to penalization.")
	

          
     if(is.null(offset)) 
          offset <- rep(0,length(y))
     if(length(offset) != length(y)) 
          stop("Offset should have the same length as y. Thanks.")

     num_ran <- get_lengths_id <- get_nrows_Z <- length_uniids <- numeric(length(Z)); 
     for(k in 1:length(Z)) 
          { 
          Z[[k]] <- as.matrix(Z[[k]])
          num_ran[k] <- ncol(Z[[k]])
          get_nrows_Z[k] <- nrow(Z[[k]]) 
          }
     for(k in 1:length(id)) 
          { 
          id[[k]] <- as.integer(id[[k]])
          get_lengths_id[k] <- length(id[[k]])
          length_uniids[k] <- length(unique(id[[k]])) 
          }
     if((get_nrows_Z != get_lengths_id) || (get_nrows_Z != nrow(X)) || (nrow(X) != get_lengths_id)) 
          stop("The number of rows in X, the length of each element in list id, and the number of rows in each element in list Z should all be the same. Thanks.")
     rm(get_lengths_id, get_nrows_Z)

	
     if(!is.null(cov.groups)) 
          {
          actual_cov_groups <- cov.groups 
          if(control$trace) 
               message("Group penalization performed as cov.groups has been supplied.")
          if((ncol(X) != length(cov.groups))) 
               stop("If supplied, the length of cov.groups should be equal to the number of columns in X. Thanks") 
          }
     if(is.null(cov.groups)) 
          actual_cov_groups <- 1:ncol(X)
               
               
     if(length(pen.type) == 1) 
          pen.type <- rep(pen.type,2)
     if(sum(pen.type %in% c("lasso","scad","adl","mcp")) != 2) 
          stop("Current version does not permit specified penalty. Sorry!")		
     if(pen.type[1] == "adl" & !is.vector(pen.weights$fixed)) 
          stop("Please supply weights for the fixed effects in adaptive lasso. Specifically, pen.weights$fixed should be a vector.")
     if(pen.type[2] == "adl" & !is.list(pen.weights$ran)) 
          stop("Please supply weights for the random effects in adaptive lasso. Specifically, pen.weights$ran should be a list with the same number of elements as lists Z and id.")
     if(!is.null(pen.weights$fixed) & length(pen.weights$fixed) != ncol(X)) 
          stop("Weights provided for fixed effects must have the same length as the number of columns in X. Thanks.")
     if(!is.null(pen.weights$ran)) 
          { 
          get_lengths_weights <- sapply(pen.weights$ran, length)
          if(!all(get_lengths_weights == num_ran)) 
               stop("The length of each element in pen.weights$ran should equal the number of columns in the corresponding element of list Z. Thanks.") 
          rm(get_lengths_weights) 
          }
          
          
     # 	if(!(ran.cov.est %in% c("bb","wbb","laplace"))) 
     # 		stop("Current version does not permit specified estimator of random effects covariance matrix. Sorry!")
     # 	if(ran.cov.est == "laplace" & length(Z) > 1) {
     # 		message("Laplace type estimator of random effects covariance matrix is permitted only if there is one set of random effects, i.e. length(Z) == length(id) = 1.") 
     # 		message("Switching to ran.cov.est = wbb") 
     # 		}
          
          
     if(!(family$family[1] %in% c("gaussian","poisson","binomial","negative.binomial","Gamma","LOGNO","ZIP"))) 
          stop("Current version does not permit specified family. Sorry!")
     if(family$family[1] == "binomial" & !(length(trial.size) %in% c(1,nrow(X)))) 
          stop("Binomial family requires trial.size to either be a single non-zero value or a vector of non-zero values with length equal to the number of rows in X.")
     if(family$family[1] != "binomial" & length(trial.size) != 1) 
          {
          trial.size <- 1
          message("trial.size argument ignored for non-binomial families.")
          }
          
          
     if(length(lambda) == 1) 
          {
          lambda <- rep(lambda, 2)
          if(!control$lasso.lambda.scale)
               warnings("Scaling factor for the second tuning parameter has been turned off for lasso and adaptive lassp penalties. 
                    This may cause issues as the group coefficient penalty on the random effects is on a different scale to the single 
                    coefficient penalty on the fixed effects.")
          }
     lambda <- as.vector(unlist(lambda))
     if(!(length(lambda) %in% c(1,2))) 
          stop("lambda should be a vector of length 1 or 2.")
     runif(1)
     old <- .Random.seed
     on.exit( { .Random.seed <<- old } )
     if(!is.null(control$seed)) 
          set.seed(control$seed)
               
     init_fit_run <- 0
     if(is.null(start)) 
          {
          get_init <- start.fixed(y = y, X = X, Z = Z, id = id, family = family, offset = offset, trial.size = trial.size)
          start <- build.start.fit.rpql(fit = get_init, id = id, num.ran = num_ran, cov.groups = cov.groups)
          init_fit_run <- 1
          }

               
     new_beta <- beta <- start$fixef
     new_b <- b <- start$ranef
     if(!is.list(new_b)) 
          stop("start$ranef should be a list of random effects. Thanks.")
     for(k in 1:length(new_b)) 
          new_b[[k]] <- b[[k]] <- as.matrix(b[[k]])
     new_phi <- phi <- new_shape <- shape <- new_zeroprob <- zeroprob <- 1
     if(family$family[1] == "gaussian") 
          new_phi <- phi <- mean((y - X %*% beta - offset)^2) 
     if(family$family[1] == "LOGNO") 
          new_phi <- phi <- mean((log(y) - X %*% beta - offset)^2) 
     if(family$family[1] == "negative.binomial") 
          new_phi <- phi <- 0.01 
     if(family$family[1] == "Gamma") 
          new_shape <- shape <- 1
     if(family$family[1] == "ZIP") 
          new_zeroprob <- zeroprob <- 0.1

     if(is.null(start$ran.cov)) 
          { 
          new_ranef_covmat <- ranef_covmat <- vector("list",length(id))
          for(k in 1:length(ranef_covmat)) 
               new_ranef_covmat[[k]] <- ranef_covmat[[k]] <- cov(new_b[[k]]) + diag(x=rep(1e-4,num_ran[k])) 
          }
     if(!is.null(start$ran.cov)) 
          new_ranef_covmat <- ranef_covmat <- start$ran.cov
          

     ## Construcs big block diagonal matrix based on an a vector of id and a Z matrix. Basically converts Z to a block Z.
     make.bigZ <- function(Z, id) 
          {
          matlist <- split(as.data.frame(Z), id)
          matlist <- lapply(matlist, function(x) as.matrix(x))
          bigZ <- matrix(0, nrow = nrow(Z), ncol = ncol(Z)*length(matlist))
          for(k in 1:length(matlist)) 
               { 
               bigZ[as.numeric(rownames(matlist[[k]])), (k*ncol(Z)-ncol(Z)+1):(k*ncol(Z))] <- matlist[[k]] 
               }
          bigZ <- Matrix(bigZ, sparse = TRUE)
          return(bigZ)
          }

               
     diff <- 1e4
     cw_rpqllogLik <- 0
     counter <- true.counter <- restart_counter <- 0
     ranef_fullshrunk <- rep(0,length(Z)) ## Indicator if random effects has been shrunk to zero. If 1, then ignore updating it!
     if(family$family[1] %in% c("LOGNO","gaussian")) 
          XtWX <- crossprod(X) ## Define these outside as they do not change with iteration
     all_bigZ <- vector("list", length(id))
     for(k in 1:length(id)) 
          all_bigZ[[k]] <- make.bigZ(Z = Z[[k]], id = id[[k]]) ## Define these outside as they do not change with iteration, but save them as sparse matrices
	
    while(diff > control$tol & true.counter < control$maxit & restart_counter <= control$restarts) 
          {
          all_Zb <- numeric(length(y)) 
          for(k2 in 1:length(id)) 
               all_Zb <- all_Zb + as.numeric(all_bigZ[[k2]] %*% c(t(b[[k2]])))

            
          if(!(family$family[1] %in% c("gaussian","LOGNO"))) 
               {
               new_etas <- X %*% beta + all_Zb + offset
               new_etas[new_etas > 30] <- 30
               new_etas[new_etas < -30] <- -30

               if(family$family[1] == "negative.binomial") 
                    new_weights <- (family$mu.eta(new_etas))^2/(family$variance(family$linkinv(new_etas),phi=phi)) ## 1/(variance*g'(mu)^2)
               if(family$family[1] == "Gamma")
                    new_weights <- (family$mu.eta(new_etas))^2/(family$variance(family$linkinv(new_etas))/shape) 
               if(family$family[1] == "ZIP") 
                    {
                    new_weights <- (poisson()$mu.eta(new_etas))^2/(poisson()$variance(poisson()$linkinv(new_etas))) 
                    weightmulti <- (1-zeroprob)*dpois(y,lambda=family$mu.linkinv(new_etas))/dZIP(x=y,mu=family$mu.linkinv(new_etas),sigma=zeroprob)
                    weightmulti[!is.finite(weightmulti)] <- 0
                    new_weights <- new_weights*weightmulti
                    rm(weightmulti)
                    }
               if(!(family$family[1] %in% c("ZIP","negative.binomial","Gamma"))) 
                    new_weights <- trial.size*(family$mu.eta(new_etas))^2/(family$variance(family$linkinv(new_etas))) ## 1/(variance*g'(mu)^2)
               if(family$family[1] != "ZIP") 
                    new_workresp <- new_etas + (y/trial.size-family$linkinv(new_etas))/family$mu.eta(new_etas) ## z = eta + (y-mu)*g'(mu)
               if(family$family[1] == "ZIP") 
                    new_workresp <- new_etas + (y-family$mu.linkinv(new_etas))/poisson()$mu.eta(new_etas) 

               new_weights <- as.vector(new_weights)	
               new_weights[new_weights > 1e4] <- 1e4
               #XtWX <- crossprod(X*sqrt(new_weights)) ## REMOVE SINCE USING C++
               #XtWZ <- crossprod(X, (new_workresp-all_Zb-offset)*new_weights) ## REMOVE SINCE USING C++
               }

                         
          ## Update fixed effects conditional on random
          absbetas <- sapply(split(as.vector(beta),actual_cov_groups), l2.norm)[actual_cov_groups] + 1e-6
          penmat <- matrix(0,length(absbetas),length(absbetas))		
          if(pen.type[1] %in% c("adl","lasso")) 
               diag(penmat) <- length(y)*lambda[1]/as.vector(absbetas)
          if(pen.type[1] == "scad") 
               diag(penmat) <- length(y)*scad.deriv(absbetas, lambda=lambda[1], a = control$scad.a)/as.vector(absbetas)
          if(pen.type[1] == "mcp") 
               diag(penmat) <- length(y)*mcp.deriv(absbetas, lambda=lambda[1], gamma = control$mcp.gamma)/as.vector(absbetas)
          if(!is.null(pen.weights$fixed)) 
               diag(penmat) <- diag(penmat)*c(pen.weights$fixed)
          if(intercept) 
               diag(penmat)[1] <- 0 
          penmat[which(penmat > 1e8)] <- 1e8
        
          if(family$family[1] == "gaussian") 
               new_beta <- chol2inv(chol(XtWX + phi*penmat)) %*% crossprod(X, y-all_Zb-offset)
          if(family$family[1] == "LOGNO") 
                new_beta <- chol2inv(chol(XtWX + phi*penmat)) %*% crossprod(X, log(y)-all_Zb-offset)
          if(!(family$family[1] %in% c("gaussian","LOGNO"))) 
               {
               new_beta <- newcoeffn(X*sqrt(new_weights), penmat, (new_workresp-all_Zb-offset)*sqrt(new_weights))			
            #new_beta <- try(chol2inv(chol(XtWX + penmat)) %*% XtWZ, silent=TRUE) ## REMOVE SINCE USING C++
#  			if(inherits(new_beta,"try-error")) 
#                      new_beta <- ginv(XtWX + penmat) %*% XtWZ
               }
          if(lambda[1] > 0) 
               {
               if(intercept) 
                    new_betaint <- new_beta[1]
                    newbeta.l2 <- sapply(split(as.vector(new_beta),actual_cov_groups), l2.norm)
                    get.zeros <- which(newbeta.l2 < conv.eps)
                    new_beta[which(actual_cov_groups %in% get.zeros)] <- 0
               if(intercept) 
                    new_beta[1] <- new_betaint ## Do not allow intercept to be shrunk to zero
               rm(get.zeros, newbeta.l2)
               }
                    

          ## Update random effects conditional on fixed
          for(k in 1:length(id)) 
               {
               if(ranef_fullshrunk[k] == 1) 
                    next;
                
               all_Zb <- numeric(length(y)); 
               if(length(Z) > 1) 
                    { 
                    for(k2 in (1:length(id))[-k]) 
                         all_Zb <- all_Zb + as.numeric(all_bigZ[[k2]] %*% c(t(b[[k2]]))) 
                    } ## Conditioning on all linear predictors for other ids, if is more than one id
               bigZ <- as.matrix(all_bigZ[[k]])
        
               #do.inv.ranef_covmat <- try(solve(ranef_covmat[[k]] + diag(x=1e-5, nrow = num_ran[k])), silent=TRUE); 
               do.inv.ranef_covmat <- ginv(ranef_covmat[[k]]) # ## Need this when rows of ranef_covmat are zero
                
               if(pen.type[2] %in% c("adl","lasso")) 
                    {
                    if(!control$lasso.lambda.scale)
                         penmat2 <- diag(x=length(y)*lambda[2]/(sqrt(colSums(b[[k]]^2))+1e-6), nrow = num_ran[k]) 
                    if(control$lasso.lambda.scale)
                         penmat2 <- diag(x=length(y)*(lambda[2]/sqrt(length_uniids[k]))/(sqrt(colSums(b[[k]]^2))+1e-6), nrow = num_ran[k]) 
                    }
               if(pen.type[2] == "scad") 
                    penmat2 <- diag(x=length(y)*scad.deriv(sqrt(colSums(b[[k]]^2)),lambda=lambda[2],a=control$scad.a)/(sqrt(colSums(b[[k]]^2))+1e-6), nrow = num_ran[k]) 
               if(pen.type[2] == "mcp") 
                    penmat2 <- diag(x=length(y)*mcp.deriv(sqrt(colSums(b[[k]]^2)),lambda=lambda[2],gamma=control$mcp.gamma)/(sqrt(colSums(b[[k]]^2))+1e-6), nrow = num_ran[k]) 
               if(!is.null(pen.weights$ran)) 
                    diag(penmat2) <- diag(penmat2)*c(pen.weights$ran[[k]])
               penmat2 <- do.inv.ranef_covmat + penmat2
               penmat2[which(penmat2 > 1e8)] <- 1e8
                        
                        
               if(family$family[1] == "gaussian") 
                    {
                    #new_bvec <- chol2inv(chol(crossprod(bigZ) + kronecker(diag(x=length_uniids[k]),phi*penmat2))) %*% crossprod(bigZ, c(y - X%*%new_beta - all_Zb - offset)) 
                    new_bvec <- newcoeffn(bigZ, kronecker(diag(x=length_uniids[k]),phi*penmat2), c(y - X %*% new_beta - all_Zb - offset))
                    }
               if(family$family[1] == "LOGNO") 
                    {
                    #new_bvec <- chol2inv(chol(crossprod(bigZ) + kronecker(diag(x=length_uniids[k]),phi*penmat2))) %*% crossprod(bigZ, c(log(y) - X%*%new_beta -all_Zb-offset))
                    new_bvec <- newcoeffn(bigZ, kronecker(diag(x=length_uniids[k]),phi*penmat2), c(log(y) - X %*% new_beta - all_Zb - offset))
                    }
               if(!(family$family[1] %in% c("gaussian","LOGNO"))) 
                    {
                    #subworkresp <- crossprod(bigZ, diag(x=new_weights))%*%c(new_workresp - X%*%new_beta - all_Zb - offset) ## REMOVE SINCE USING C++
                    #do.bvec <- try(chol2inv(chol(as.matrix(crossprod(sqrt(new_weights)*bigZ)) + kronecker(diag(x=length_uniids[k]),penmat2))) %*% subworkresp, silent=TRUE) ## REMOVE SINCE USING C++
                    do.bvec <- newcoeffn(sqrt(new_weights)*bigZ, kronecker(diag(x=length_uniids[k]),penmat2), c(new_workresp - X %*% new_beta - all_Zb - offset)*sqrt(new_weights)) 

                    if(!inherits(do.bvec,"try-error")) 
                         new_bvec <- do.bvec
                    if(inherits(do.bvec,"try-error")) 
                         {
                         subworkresp <- crossprod(bigZ, diag(x=new_weights))%*%c(new_workresp - X %*% new_beta - all_Zb - offset)
                         new_bvec <- ginv(as.matrix(crossprod(sqrt(new_weights)*bigZ)) + kronecker(diag(x=length_uniids[k]),penmat2)) %*% subworkresp 
                         }
                    }
               new_b[[k]] <- matrix(new_bvec, nrow = length_uniids[k], ncol = num_ran[k], byrow = TRUE)		
               sel_nonzerob <- which(sqrt(colSums(new_b[[k]]^2)) > length_uniids[k]*conv.eps)
               if(length(sel_nonzerob) > 0) 
                    new_b[[k]][,-sel_nonzerob] <- 0
               if(length(sel_nonzerob) == 0) 
                    new_b[[k]][,] <- 0
                
                
                ## Update random effects covariance matrix 
                ## Iterative estimator based on maximizing the Laplace approximation 
                ## If multiple Z are applied, then the Laplace approximated estimator is done conditonally on all other random effects
                new_ranef_covmat[[k]] <- matrix(0, nrow = num_ran[k], ncol = num_ran[k])
                if(length(sel_nonzerob) == 0) 
                    next;

                if(family$family[1] %in% c("gaussian","LOGNO")) 
                    {
                    for(i in 1:length_uniids[k]) 
                         { 
                         sel.i <- which(id[[k]] == i)
                         new_ranef_covmat[[k]][sel_nonzerob,sel_nonzerob] <- new_ranef_covmat[[k]][sel_nonzerob,sel_nonzerob] + chol2inv(chol(crossprod(as.matrix(Z[[k]][sel.i,sel_nonzerob]))/phi + do.inv.ranef_covmat[sel_nonzerob,sel_nonzerob])) + tcrossprod(new_b[[k]][i,sel_nonzerob])
                         }
                    new_ranef_covmat[[k]] <- new_ranef_covmat[[k]]/length_uniids[k]
                    }

                if(!(family$family[1] %in% c("gaussian","LOGNO"))) 
                    {
                    for(i in 1:length_uniids[k]) 
                         { 
                         sel.i <- which(id[[k]] == i)
                         new_ranef_covmat[[k]][sel_nonzerob,sel_nonzerob] <- new_ranef_covmat[[k]][sel_nonzerob,sel_nonzerob] + chol2inv(chol(crossprod(sqrt(new_weights[sel.i])*as.matrix(Z[[k]][sel.i,sel_nonzerob])) + do.inv.ranef_covmat[sel_nonzerob,sel_nonzerob])) + tcrossprod(new_b[[k]][i,sel_nonzerob])
                         }
                    new_ranef_covmat[[k]] <- new_ranef_covmat[[k]]/length_uniids[k]
                    }			

# 			if(ran.cov.est == "bb") { ## Based on (1/n_k) sum_{i=1}^{n_k} b_{ik}%*%t(b_{ik})
# 				new_ranef_covmat[[k]] <- cov.wt(new_b[[k]], center = FALSE, method = "ML")$cov
# 				}
# 				
# 			if(ran.cov.est == "wbb") { ## Based on (1/n_k) sum_{i=1}^{n_k} w_{ik} b_{ik}%*%t(b_{ik}), where w_{ik} is based on cluster size
# 				clus.size <- unlist(as.vector(table(id[[k]])))
# 				new_ranef_covmat[[k]] <- cov.wt(new_b[[k]], center = FALSE, method = "ML", wt = clus.size)$cov
# 				}
                        
                if(length(sel_nonzerob) == 0) 
                    ranef_fullshrunk[k] <- 1
                rm(bigZ, new_bvec)
                }

                    
          ## Solve scale parameters if required
          all_Zb <- numeric(length(y)) 
          for(k2 in 1:length(id)) 
               {
               all_Zb <- all_Zb + as.numeric(all_bigZ[[k2]] %*% c(t(new_b[[k2]]))) 
               }
          if(family$family[1] == "gaussian") 
               new_phi <- mean((y - X %*% new_beta - all_Zb - offset)^2)
          if(family$family[1] == "LOGNO") 
               new_phi <- mean((log(y) - X %*% new_beta - all_Zb - offset)^2)

          if(family$family[1] == "negative.binomial") 
               {
               new_phi <- try(1/theta.ml(y = y, mu = family$linkinv(X %*% new_beta + all_Zb + offset), limit = 100), silent = TRUE)
               if(inherits(new_phi,"try-error")) 
                    new_phi <- phi
               }
                
          if(family$family[1] == "Gamma") 
               {
               prof.logl <- function(a, y, mu) 
                    {
                    out <- -y*a/mu + (a-1)*log(y) - lgamma(a) - a*log(mu) + a*log(a)
                    return(sum(out)) 
                    }
               update.shape <- try(optimize(f = prof.logl, interval = c(1e-3,1e3), y = y, mu = family$linkinv(X %*% new_beta+all_Zb+offset), maximum = TRUE), silent = TRUE)
               if(inherits(update.shape,"try-error"))    
                    new_shape <- shape 
               if(!inherits(update.shape,"try-error")) 
                    new_shape <- update.shape$max
               }
                
          if(family$family[1] == "ZIP") 
               {
               prof.logl <- function(a, y, mu) 
                    {
                    out <- sum(dZIP(x=y, mu=mu, sigma=a, log=TRUE)) 
                    return(out)
                    }
               update.zeroprob <- try(suppressWarnings(optimize(f = prof.logl, interval = c(1e-4,0.9999), y = y, mu = family$mu.linkinv(X %*% new_beta+all_Zb+offset), maximum = TRUE)), silent = TRUE)
               if(inherits(update.zeroprob,"try-error")) 
                    new_zeroprob <- zeroprob
               if(!inherits(update.zeroprob,"try-error")) 
                    new_zeroprob <- update.zeroprob$max
               }
                    

                    
          diff <- sum((new_beta-beta)^2) 
          for(k in 1:length(id)) 
               diff <- diff + sum((new_ranef_covmat[[k]][upper.tri(new_ranef_covmat[[k]],diag=TRUE)]-ranef_covmat[[k]][upper.tri(ranef_covmat[[k]],diag=TRUE)])^2)
          #for(k in 1:length(id)) { diff <- diff + sum((new_b[[k]]-b[[k]])^2) }
          if(control$trace) 
               { 
               cat("Iteration:", counter,"   Error:", round(diff,4), "\n")
                    
               new_beta2 <- new_beta
               names(new_beta2) <- colnames(X)
               message("Fixed effects:")
               print(c(new_beta2)); #print(c(newbeta.l2))
               message("Covariance Matrices:")
               print(new_ranef_covmat); 
               rm(new_beta2) 
               } 

                    
          beta <- new_beta
          b <- new_b; 
          ranef_covmat <- new_ranef_covmat; 
          shape <- new_shape
          zeroprob <- new_zeroprob
          phi <- new_phi
          counter <- counter + 1;
          true.counter <- true.counter + 1; 
          if(counter < 5) 
               true.counter <- 0

            
          if(any(unlist(lapply(ranef_covmat, function(x) any(abs(x) > 30)) == TRUE)) || sum(abs(new_beta)) < 1e-3) 
               {
               cat("Convergence issues. Restarting...\n")
               restart_counter <- restart_counter + 1
               if(init_fit_run == 0) 
                    {
                    get_init <- start.fixed(y = y, X = X, Z = Z, id = id, family = family, offset = offset, trial.size = trial.size)
                    start <- build.start.fit.rpql(fit = get_init, id = id, num.ran = num_ran, cov.groups = cov.groups)
                    init_fit_run <- 1
                    }
                         
               new_beta <- beta <- start$fixef
               new_b <- b <- vector("list",length(id))
               for(k in 1:length(b)) 
                    { 
                    new_b[[k]] <- b[[k]] <- matrix(rnorm(length_uniids[k]*num_ran[k], mean=0, sd = 0.2+0.2*(family$family!="poisson")),length_uniids[k],num_ran[k]) 
                    }
               new_ranef_covmat <- ranef_covmat <- vector("list",length(id))
               for(k in 1:length(ranef_covmat)) 
                    new_ranef_covmat[[k]] <- ranef_covmat[[k]] <- cov(new_b[[k]])

               new_phi <- phi <- new_shape <- shape <- new_zeroprob <- zeroprob <- 1
               if(family$family[1] == "gaussian") 
                    new_phi <- phi <- mean((y - X %*% beta - offset)^2) 
               if(family$family[1] == "LOGNO") 
                    new_phi <- phi <- mean((log(y) - X %*% beta - offset)^2) 
               if(family$family[1] == "negative.binomial") 
                    new_phi <- phi <- 0.01 
               if(family$family[1] == "Gamma") 
                    new_shape <- shape <- 1
               if(family$family[1] == "ZIP") 
                    new_zeroprob <- zeroprob <- 0.1
                         
               diff <- 1e4
               cw_rpqllogLik <- 0
               counter <- 0 
               }		

          if(restart_counter > control$restarts) 
               { 
               stop("Sorry, but convergence could not be reached within the allocated number of control$restarts. Consider increasing lambda[2].\n")
               return()
               }

        }
        ## DONE!

		
     sel_nonzerob <- sel.zero.b <- vector("list",length(id))
     for(k in 1:length(id)) 
          { 
          sel_nonzerob[[k]] <- which(sqrt(colSums(new_b[[k]]^2)) > length_uniids[k]*conv.eps) 
          if(length(sel_nonzerob[[k]]) > 0) 
               new_b[[k]][,-sel_nonzerob[[k]]] <- 0
          }
	
	
     ## The PQL likelihood
     new_etas <- X %*% new_beta + offset
     for(k2 in 1:length(id)) 
          new_etas <- new_etas + as.numeric(all_bigZ[[k2]] %*% c(t(new_b[[k2]])))

     if(family$family[1] == "gaussian") 
          pqllogLik <- sum(dnorm(y, mean=new_etas, sd=sqrt(new_phi), log=TRUE)) 
     if(family$family[1] == "Gamma") 
          pqllogLik <- sum(dgamma(y, shape=new_shape, scale=family$linkinv(new_etas)/new_shape, log=TRUE)) 
     if(family$family[1] == "binomial") 
          pqllogLik <- sum(dbinom(y, size=trial.size, prob=family$linkinv(new_etas), log=TRUE)) 
     if(family$family[1] == "poisson") 
          pqllogLik <- sum(dpois(y, lambda=family$linkinv(new_etas), log=TRUE)) 
     if(family$family[1] == "negative.binomial") 
          pqllogLik <- sum(dnbinom(y, mu=family$linkinv(new_etas), size=1/new_phi, log=TRUE))
     if(family$family[1] == "LOGNO") 
          pqllogLik <- sum(dlnorm(y, meanlog=new_etas, sdlog=sqrt(new_phi), log=TRUE))
     if(family$family[1] == "ZIP") 
          pqllogLik <- sum(dZIP(y, mu=family$mu.linkinv(new_etas), sigma=new_zeroprob, log=TRUE))
     loglik1 <- pqllogLik
    
     for(k in 1:length(id)) 
          { 
          if(length(sel_nonzerob[[k]]) > 0) 
               pqllogLik <- pqllogLik - 0.5*sum(rowSums(new_b[[k]][,sel_nonzerob[[k]]]*(new_b[[k]][,sel_nonzerob[[k]]]%*%solve(new_ranef_covmat[[k]][sel_nonzerob[[k]],sel_nonzerob[[k]]])))) 
          }
     names(sel_nonzerob) <- names(id)
               
               
     ## Information Criterion...argh!
     aic.v1 <- -2*loglik1 + 2*(sum(new_beta != 0) + sum(sapply(new_b, function(x) sum(x!=0))))
     bic.v2 <- -2*loglik1 + log(length(y))*sum(new_beta!=0) + sum(log(length_uniids)*sapply(new_ranef_covmat, function(x) sum(x[upper.tri(x,diag=T)]!=0))) ## Number of random effects to penalize is just number of non-zero elements in covariance matrix
     bic.v3 <- -2*loglik1 + log(length(y))*sum(new_beta!=0) + sum(log(length(y))*sapply(new_ranef_covmat, function(x) sum(x[upper.tri(x,diag=T)]!=0))) ## Number of random effects to penalize is just number of non-zero elements in covariance matrix
     hic.v1 <- -2*loglik1 + log(length(y))*sum(new_beta!=0) + 2*sum(sapply(new_b, function(x) sum(x!=0)))
     hic.v2 <- -2*loglik1 + log(length(y))*sum(new_beta!=0) + sum(sapply(new_b, function(x) sum(x!=0)))
     hic.v3 <- -2*loglik1 + log(length(y))*sum(new_beta!=0) + 0.5*sum(sapply(new_b, function(x) sum(x!=0)))
	
	
     dof <- sum(new_beta != 0) + sum(sapply(new_b, function(x) sum(x!=0)))
     ics <- c(dof, aic.v1, bic.v2, bic.v3, hic.v1, hic.v2, hic.v3) 
     names(ics) = c(
          "# of estimated parameters",
          "AIC: 2*sum(beta!=0) + 2*sum(b!=0)",
          "BIC: log(nrow(X))*sum(beta!=0) + log(number of clusters)*sum(ranef_covmat!=0)",
          "BIC: log(nrow(X))*sum(beta!=0) + log(nrow(X))*sum(ranef_covmat!=0)",
          "Hybrid IC: log(nrow(X))*sum(beta!=0) + 2*sum(b!=0)",
          "Hybrid IC: log(nrow(X))*sum(beta!=0) + sum(b!=0)",
          "Hybrid IC: log(nrow(X))*sum(beta!=0) + 0.5*sum(b!=0)")

               
     ## Garnishing
     out_list <- list(fixef = as.vector(new_beta), ranef = new_b, ran.cov = new_ranef_covmat, pql.logLik = pqllogLik, logLik1 = loglik1, 
          phi = new_phi, shape = new_shape, zeroprob = new_zeroprob, family = family, n = length(y), 
          trial.size = trial.size, id = id, lambda = lambda, pen.weights = pen.weights, pen.type = pen.type, ics = ics, 
          nonzero.fixef = which(as.vector(new_beta)!=0), nonzero.ranef = sel_nonzerob)

     if(save.data) 
          { 
          out_list$y <- y
          out_list$X <- X
          outlist$Z <- Z
          out_list$offset <- offset
          message("Please note the column names in Z have been overwritten for easier reference in the estimation algorithm. Apologies in advance") 
          }

		
     if(hybrid.est) 
          {
          if(family$family[1] %in% c("ZIP","LOGNO")) 
               { 
               message("Hybrid estimation not possible as lme4 does not do this. Sorry! Maybe try glmmTMB?")
               out_list$hybrid <- NULL 
               }
                    
          else {	
               if(control$trace) 
                    message("Performing hybrid estimation with lme4...might take a bit of time, and the labels for the parameters will be wrong. Apologies in advance!")

               make_dat <- data.frame(y, do.call(cbind,Z), X, do.call(cbind,id))
               restring1 <- NULL
               if(length(out_list$nonzero.fixef) > 0) 
                    restring1 <- paste(paste(colnames(X)[out_list$nonzero.fixef],collapse="+"),"-1")
               for(k in 1:length(id)) 
                    { 
                    if(length(out_list$nonzero.ranef[[k]])) restring1 <- paste(restring1, "+ (",paste0(colnames(Z[[k]])[out_list$nonzero.ranef[[k]]],collapse="+"),"-1|",names(id)[k],")")
                    }
               restring1 <- reformulate(restring1, response = "y")
    # 		print(restring1)
               if(family$family[1] == "gaussian") 
                    final_fit <- suppressWarnings(lmer(restring1, REML = TRUE, offset = offset, data = make_dat))
               if(family$family[1] != "gaussian" & family$family[1] != "negative.binomial") 
                    final_fit <- suppressWarnings(glmer(restring1, family = family, offset = offset, data = make_dat))
               if(family$family[1] == "negative.binomial") 
                    final_fit <- suppressWarnings(glmer.nb(restring1, family = family, offset = offset, data = make_dat))

               out_list$hybrid <- final_fit ## Access covariance matrix as VarCorr(final_fit)[[1]][1:nrow(VarCorr(final_fit)[[1]]),1:nrow(VarCorr(final_fit)[[1]])]
                    
               rm(make_dat)
               }
          }	
    

     names(out_list$fixef) <- colnames(X)
     names(out_list$ranef) <- names(out_list$ran.cov) <- names(id)
     for(k in 1:length(id)) 
          { 
          rownames(out_list$ranef[[k]]) <- 1:length_uniids[k]; colnames(out_list$ranef[[k]]) <- colnames(Z[[k]])
          rownames(out_list$ran.cov[[k]]) <- colnames(out_list$ran.cov[[k]]) <- colnames(Z[[k]]) 
          }
            
     class(out_list) <- "rpql"
     out_list$call <- match.call()
    
     return(out_list)
     }

	
	
print.rpql <- function(x, ...) 
     {
    message("Call:")
    print(x$call)
    message()

    cat("Family:", x$family$family[1], "\nPenalty type:", x$pen.type, "\nTuning parameters:", as.numeric(x$lambda), "\n") 
    cat("Total number of observations:", x$n, "\n IDs (groups):", paste(names(x$id),sapply(x$id,function(a) length(unique(a))), sep = ": ", collapse = ";"), "\n") 
    cat("Non-zero fixed effects:", x$nonzero.fixef, "\nNon-zero random effects by IDs (groups):", paste(names(x$id),x$nonzero.ranef, sep = ": ", collapse = ";"), "\n") 
    #message("Information Criterion:")
    #print(x$ics) 
    }
	
	
