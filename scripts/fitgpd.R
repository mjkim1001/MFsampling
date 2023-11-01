lgpd.weight = function (x,weight = NULL, u = 0, sigmau = 1, xi = 0, phiu = 1, log = TRUE) 
{
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.param(u)
  check.param(weight,allowvec = T,allownull = T)
  check.param(sigmau)
  check.param(xi)
  check.prob(phiu)
  check.logic(log)
  check.inputn(c(length(u), length(sigmau), length(xi), length(phiu)), 
               allowscalar = TRUE)
  if (any(!is.finite(x))) {
    warning("non-finite cases have been removed")
    x[!is.finite(x)] = NA
  }
  xu = x[which(x > u)]
  nu = length(xu)
  yu = (xu - u)/sigmau
  syu = 1 + xi * yu
  if ((min(syu) <= 0) | (sigmau <= 0) | (phiu <= 0) | (phiu > 
                                                       1)) {
    l = -Inf
  }
  else {
    if (abs(xi) < 1e-06) {
      if(is.null(weight)){
        l = -nu * log(sigmau) - sum(yu) + nu * log(phiu)
      }else{
        l = -log(sigmau) * sum(weight) - sum(yu*weight) + log(phiu) * sum(weight) 
      }
      
    }
    else {
      if(is.null(weight)){
        l = -nu * log(sigmau)  - (1/xi + 1) * sum(log(syu)) + 
          nu * log(phiu)
      }else{
        l = - log(sigmau) * sum(weight)  - (1/xi + 1) * sum(log(syu)*weight) + 
          log(phiu) * sum(weight) 
      }
      
    }
  }
  if (!log) 
    l = exp(l)
  l
}

nlgpd.weight<-function (pvector, x, weight=NULL, u = 0, phiu = 1, finitelik = FALSE) 
{
  np = 2
  check.nparam(pvector, nparam = np)
  check.param(u)
  check.param(weight,allowvec = T,allownull = T)
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.prob(phiu)
  check.logic(finitelik)
  sigmau = pvector[1]
  xi = pvector[2]
  nllh = -lgpd.weight(x,weight, u, sigmau, xi, phiu)
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e+06
  }
  nllh
}
fgpd.weight<-function (x, weight=1,u = 0, phiu = NULL, pvector = NULL, std.err = TRUE, 
          method = "BFGS", control = list(maxit = 10000), finitelik = TRUE, 
          ...) 
{
  call <- match.call()
  np = 2
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.param(u)
  check.param(weight,allowvec = T,allownull = T)
  check.prob(phiu, allownull = TRUE)
  check.nparam(pvector, nparam = np, allownull = TRUE)
  check.logic(std.err)
  check.optim(method)
  check.control(control)
  check.logic(finitelik)
  if (any(is.infinite(x))) 
    warning("infinite cases have been removed")
  x = x[!is.infinite(x)]
  if (any(is.na(x))) 
    warning("missing values are treated as below threshold when estimating tail fraction")
  check.quant(x, allowna = TRUE)
  n = length(x)
  if ((method == "L-BFGS-B") | (method == "BFGS")) 
    finitelik = TRUE
  xu = x[which(x > u)]
  if(length(weight)>0)  weights = weight[which(x > u)]
  nu = length(xu)
  if (nu < 1) 
    stop("no elements of x are above threshold")
  if (is.null(pvector)) {
    yu = xu - u
    pvector[1] = sqrt(6 * var(yu))/pi
    pvector[2] = 0.1
  }
  if (is.null(phiu)) {
    phiu = nu/n
    se.phiu = sqrt(phiu * (1 - phiu)/n)
  }
  else {
    if (phiu == 0) 
      stop("tail probability must be in (0, 1]")
    se.phiu = NA
  }
  nllh = nlgpd.weight(pvector, xu, weights, u, phiu)
  if (is.infinite(nllh)) {
    pvector[2] = 0.1
    nllh = nlgpd.weight(pvector, xu, weights, u, phiu)
  }
  if (is.infinite(nllh)) 
    stop("initial parameter values are invalid")
  fit = optim(par = as.vector(pvector), fn = nlgpd.weight, x = xu, weight=weights, 
              u = u, phiu = phiu, finitelik = finitelik, control = control, 
              method = method, hessian = TRUE, ...)
  conv = TRUE
  if ((fit$convergence != 0) | any(fit$par == pvector) | (abs(fit$value) >= 
                                                          1e+06)) {
    conv = FALSE
    warning("check convergence")
  }
  if (conv & std.err) {
    qrhess = qr(fit$hessian)
    if (qrhess$rank != ncol(qrhess$qr)) {
      warning("observed information matrix is singular")
      se = NULL
      invhess = NULL
    }
    else {
      invhess = solve(qrhess)
      vars = diag(invhess)
      if (any(vars <= 0)) {
        warning("observed information matrix is singular")
        invhess = NULL
        se = NULL
      }
      else {
        se = sqrt(vars)
      }
    }
  }
  else {
    invhess = NULL
    se = NULL
  }
  list(call = call, x = as.vector(x), init = as.vector(pvector), 
       optim = fit, conv = conv, cov = invhess, mle = fit$par, 
       se = se, rate = phiu, nllh = fit$value, n = n, u = u, 
       sigmau = fit$par[1], xi = fit$par[2], phiu = phiu, se.phiu = se.phiu)
}

