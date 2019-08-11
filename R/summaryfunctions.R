print.summary.rpql <- function(x, ...) 
     {
    message("Call:")
    print(x$call)
    message()
    
    cat("Family:", x$family$family[1], "\nPenalty type:", x$pen.type, "\n Tuning parameters:", as.numeric(x$lambda), "\n") 
    cat("Value of PQL at convergence:", x$logLik, "\n\n")
    #cat("Estimated model\n\t Non-zero fixed effects --", x$nonzero.fixef, "\n\t Non-zero random effects --", paste(names(x$id),x$nonzero.ranef, sep = ": ", collapse = "; "), "\n\n") 
    message("Estimates of fixed effects:")
    print(x$fixef); 
    message()

    message("Estimates of variances for random effects (diagonal elements of the covariance matrices):")
    for(k in 1:length(x$ranef)) 
          { 
          message(names(x$id)[k], ":")
          print(diag(x$ran.cov[[k]]))
          message() 
          } 

    if(x$family$family[1] %in% c("gaussian","lognormal")) 
        message("Variance par ameter: ", x$phi, "\n")
    if(x$family$family[1] %in% c("negative.binomial")) 
        message("overdispersion parameter (V = mu + phi*mu^2): ", x$phi, "\n")
    if((x$family$family[1] == "Gamma")) 
        message("Shape parameter (V = mu^2/shape): ", x$shape, "\n")
    if((x$family$family[1] == "ZIP")) 
        message("Probability of structural zero: ", x$zeroprob, "\n") 
    }	
	
	


#' Summary of GLMM fitted using regularized PQL.
#' 
#' A summary of the results from applying \code{rpql}.
#' 
#' 
#' @aliases summary.rpql print.summary.rpql
#' @param object An object of class "rpql".
#' @param x An object of class "rpql".
#' @param ... Not used.
#' @return A list (some of which is printed) containing the following elements:
#' \item{Call}{The matched call.}
#' 
#' \item{fixed}{Estimated fixed effects coefficients.}
#' 
#' \item{ranef}{A list with each element being a matrix of estimated random
#' effects coefficients.}
#' 
#' \item{ran.cov}{A list with each element being a estimated random effects
#' covariance matrix.}
#' 
#' \item{logLik}{PQL log-likelihood value at convergence.}
#' 
#' \item{family}{The \code{family} argument, i.e. response type.}
#' 
#' \item{pen.type,lambda}{Penalties used for selection and the corresponding
#' tuning parameter values.}
#' 
#' \item{ics}{A vector containing the number of estimated, non-zero parameters,
#' and three information criterion. Please see the help file for \code{rpql}
#' for details on these criteria.}
#' 
#' \item{id}{The \code{id} argument, i.e. list of IDs.}
#' 
#' \item{nonzero.fixef}{A vector indexing which of the estimated fixed effect
#' coefficients are non-zero.}
#' 
#' \item{nonzero.ranef}{A list with each element being a vector indexing which
#' of the estimated random effects are non-zero, i.e. which of the diagonal
#' elements in the corresponding element of \code{ran.cov} are non-zero.}
#' @author
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_author(\"#1\")}",
#' "rpql")\Sexpr{tools:::Rd_package_author("rpql")}
#' 
#' Maintainer:
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_maintainer(\"#1\")}",
#' "rpql")\Sexpr{tools:::Rd_package_maintainer("rpql")}
#' @seealso \code{\link{rpql}} for fitting and performing model selection in
#' GLMMs using regularized PQL.
#' @examples
#' 
#' ## Please see examples in help file for the rpql function
#' 
summary.rpql <- function(object, ...) 
     {
     gather_output <- list(call = object$call, fixef = round(object$fixef,3), ranef = lapply(object$ranef,round,3), 
          ran.cov = lapply(object$ran.cov,round,3), logLik = round(object$logLik,3), family = object$family)
    
    if(object$family$family[1] %in% c("gaussian","lognormal","negative.binomial")) 
        gather_output$phi <- round(object$phi,3)
    if((object$family$family[1] == "Gamma")) 
        gather_output$shape <- round(object$shape,3)
    if((object$family$family[1] == "ZIP")) 
        gather_output$zeroprob <- round(object$zeroprob,3)

    gather_output$pen.type <- object$pen.type
    gather_output$lambda <- object$lambda
    gather_output$ics <- object$ics
    gather_output$id <- object$id 
    gather_output$nonzero.fixef <- object$nonzero.fixef
    gather_output$nonzero.ranef <- object$nonzero.ranef
    gather_output$ics <- object$ics 
    
    class(gather_output) <- "summary.rpql"
    gather_output 
    }	
	
	
	
