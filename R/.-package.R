

#' Joint effects selection in GLMMs using regularized PQL
#' 
#' \code{rpql} offers fast joint selection of fixed and random effects in
#' Generalized Linear Mixed Model (GLMMs) via regularization. Specifically the
#' penalized quasi-likelihood (PQL, Breslow and Clayton, 1993) is used as a
#' loss function, and penalties are added on to perform fixed and random
#' effects selection e.g., the lasso (Tibshirani, 1996) penalty. This method of
#' joint selection in GLMMs, referred to regularized PQL, is very fast compared
#' to information criterion and hypothesis testing, and has attractive large
#' sample properties (Hui et al., 2016). Its performance however may not be
#' great if the amount of data to estimate each random effect is not large,
#' i.e. the cluster size is not large.
#' 
#' \tabular{ll}{ Package: \tab rpql\cr Type: \tab Package\cr Version: \tab
#' 0.4\cr Date: \tab 2015-06-01\cr License: \tab GPL-2\cr }
#' 
#' @name rpql-package
#' @docType package
#' @author
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_author(\"#1\")}",
#' "rpql")\Sexpr{tools:::Rd_package_author("rpql")}
#' 
#' Maintainer:
#' c("\\Sexpr[results=rd,stage=build]{tools:::Rd_package_maintainer(\"#1\")}",
#' "rpql")\Sexpr{tools:::Rd_package_maintainer("rpql")}
#' @references \itemize{ \item Breslow, N. E., and Clayton, D. G. (1993).
#' Approximate inference in generalized linear mixed models. Journal of the
#' American Statistical Association, 88, 9-25.
#' 
#' \item Hui, F.K.C., Mueller, S., and Welsh, A.H. (2017). Joint Selection in
#' Mixed Models using Regularized PQL. Journal of the American Statistical
#' Association, 112, 1323-1333.
#' 
#' \item Tibshirani, R. (1996). Regression shrinkage and selection via the
#' lasso. Journal of the Royal Statistical Society. Series B (Methodological),
#' 58, 267-288.  }
#' @examples
#' 
#' ## Please see examples in help file for the rpql function
#' 
NULL



