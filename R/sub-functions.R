Chao_Hill_abu = function(x,q){
  x = x[x>0]
  n = sum(x)
  f1 = sum(x==1)
  f2 = sum(x==2)
  p1 = ifelse(f2>0,2*f2/((n-1)*f1+2*f2),ifelse(f1>0,2/((n-1)*(f1-1)+2),1))
  r <- 0:(n-1)
  Sub <- function(q){
    if(q==0){
      sum(x>0) + (n-1)/n*ifelse(f2>0, f1^2/2/f2, f1*(f1-1)/2)
    }
    else if(q==1){
      A <- sum(x/n*(digamma(n)-digamma(x)))
      B <- ifelse(f1==0|p1==1,0,f1/n*(1-p1)^(1-n)*(-log(p1)-sum(sapply(1:(n-1), function(r)(1-p1)^r/r))))
      exp(A+B)
    }else if(abs(q-round(q))==0){
      A <- sum(exp(lchoose(x,q)-lchoose(n,q)))
      ifelse(A==0,NA,A^(1/(1-q)))
    }else {
      sort.data = sort(unique(x))
      tab = table(x)
      term = sapply(sort.data,function(z){
        k=0:(n-z)
        sum(choose(k-q,k)*exp(lchoose(n-k-1,z-1)-lchoose(n,z)))
      })
      A = sum(tab*term)
      B = ifelse(f1==0|p1==1,0,f1/n*(1-p1)^(1-n)*(p1^(q-1)-sum(choose(q-1,r)*(p1-1)^r)))
      (A+B)^(1/(1-q))
    }
  }
  sapply(q, Sub)
}
Hill <- function(x,q,datatype = c("abundance","incidence")){
  if(datatype=="incidence"){x = x[-1]}
  p <- x[x>0]/sum(x)
  Sub <- function(q){
    if(q==0) sum(p>0)
    else if(q==1) exp(-sum(p*log(p)))
    else exp(1/(1-q)*log(sum(p^q)))
  }
  sapply(q, Sub)
}
coverage_to_size = function(x, C, data_type='abundance'){

  if (data_type=='abundance'){

    n <- sum(x)
    refC <- iNEXT:::Chat.Ind(x, n)
    f <- function(m, C) abs(iNEXT:::Chat.Ind(x, m) - C)
    if (refC > C) {
      opt <- optimize(f, C = C, lower = 0, upper = sum(x))
      mm <- opt$minimum
      mm <- round(mm)
    }
    if (refC <= C) {
      f1 <- sum(x == 1)
      f2 <- sum(x == 2)
      if (f1 > 0 & f2 > 0) {
        A <- (n - 1) * f1/((n - 1) * f1 + 2 * f2)
      }
      if (f1 > 1 & f2 == 0) {
        A <- (n - 1) * (f1 - 1)/((n - 1) * (f1 - 1) + 2)
      }
      if (f1 == 1 & f2 == 0) {
        A <- 1
      }
      if (f1 == 0 & f2 == 0) {
        A <- 1
      }
      mm <- (log(n/f1) + log(1 - C))/log(A) - 1
      mm <- n + mm
      mm <- round(mm)
    }

  } else {

    m <- NULL
    n <- max(x)
    refC <- iNEXT:::Chat.Sam(x, n)
    f <- function(m, C) abs(iNEXT:::Chat.Sam(x, m) - C)
    if (refC > C) {
      opt <- optimize(f, C = C, lower = 0, upper = max(x))
      mm <- opt$minimum
      mm <- round(mm)
    }
    if (refC <= C) {
      f1 <- sum(x == 1)
      f2 <- sum(x == 2)
      U <- sum(x) - max(x)
      if (f1 > 0 & f2 > 0) {
        A <- (n - 1) * f1/((n - 1) * f1 + 2 * f2)
      }
      if (f1 > 1 & f2 == 0) {
        A <- (n - 1) * (f1 - 1)/((n - 1) * (f1 - 1) + 2)
      }
      if (f1 == 1 & f2 == 0) {
        A <- 1
      }
      if (f1 == 0 & f2 == 0) {
        A <- 1
      }
      mm <- (log(U/f1) + log(1 - C))/log(A) - 1
      mm <- n + mm
      mm <- round(mm)
    }

  }

  return(mm)
}
bootstrap_population_multiple_assemblage = function(data, data_gamma, data_type){

  if (data_type == 'abundance'){

    S_obs = sum(data_gamma > 0)
    n = sum(data_gamma)
    f1 = sum(data_gamma == 1)
    f2 = sum(data_gamma == 2)
    f0_hat = ifelse(f2 == 0, (n - 1)/n * f1 * (f1 - 1)/2, (n - 1)/n * f1^2/2/f2) %>% ceiling()

    output = apply(data, 2, function(x){

      p_i_hat = chaoUtility:::bootp_one_abu(Spec = x, zero = T)

      if(length(p_i_hat) != length(x)){

        p_i_hat_unobs = p_i_hat[length(p_i_hat)]
        p_i_hat_obs = p_i_hat[1:length(x)]
        p_i_hat = c(p_i_hat_obs, rep(0, f0_hat))
        candidate = which(p_i_hat==0)

        chosen = sample(x = candidate, size = min(length(p_i_hat) - length(x), length(candidate)), replace = F)
        p_i_hat[chosen] = (1-sum(p_i_hat))/length(chosen)

        p_i_hat

      } else {

        p_i_hat = c(p_i_hat, rep(0, f0_hat))
        p_i_hat

      }
    })

  }

  if (data_type == 'incidence'){

    S_obs = sum(data_gamma > 0)
    t = data_gamma[1]
    Q1 = sum(data_gamma == 1)
    Q2 = sum(data_gamma == 2)
    Q0_hat = if ( Q2 == 0 ){( (t-1)/t ) * ( Q1*(Q1-1)/2 )} else {( (t-1)/t ) * ( (Q1^2) / (2*Q2) )} %>% ceiling

    output = apply(data, 2, function(x){

      pi_i_hat = chaoUtility:::bootp_one_inc(Spec = x, zero = T)

      if(length(pi_i_hat) != length(x)){

        pi_i_hat_unobs = pi_i_hat[length(pi_i_hat)]
        pi_i_hat_obs = pi_i_hat[1:length(x)]
        pi_i_hat = c(pi_i_hat_obs, rep(0, Q0_hat))
        candidate = which(pi_i_hat==0)
        chosen = sample(x = candidate, size = min(length(pi_i_hat) - length(x), length(candidate)), replace = F)
        pi_i_hat[chosen] = pi_i_hat_unobs

        pi_i_hat

      } else {

        pi_i_hat = c(pi_i_hat, rep(0, Q0_hat))
        pi_i_hat

      }
    })

  }

  return(output)

}
raw2freq = function(x){

  output = lapply(x, function(data){

    nT = ncol(data[[1]])

    data_alpha_freq = lapply(data, rowSums) %>% do.call(c,.) %>% c(nT, .)

    gamma = Reduce('+', data)
    gamma[gamma>1] = 1
    data_gamma_freq = c(nT, rowSums(gamma))

    data_alpha_freq = data_alpha_freq[data_alpha_freq>0] %>% as.numeric
    data_gamma_freq = data_gamma_freq[data_gamma_freq>0] %>% as.numeric

    list(data_alpha_freq=data_alpha_freq, data_gamma_freq=data_gamma_freq)

  })

  return(output)

}
FD.m.est_0 = function (ai_vi, m, q, nT) {
  EFD = function(m, qs, obs, asy, beta, av) {
    m = m - nT
    out <- sapply(1:length(qs), function(i) {
      if (qs[i] != 2) {
        obs[i] + (asy[i] - obs[i]) * (1 - (1 - beta[i])^m)
      }
      else if (qs[i] == 2) {
        V_bar^2/sum((av[, 2]) * ((1/(nT + m)) * (av[,
                                                    1]/nT) + ((nT + m - 1)/(nT + m)) * (av[, 1] *
                                                                                          (av[, 1] - 1)/(nT * (nT - 1)))))
      }
    })
    return(out)
  }
  V_bar <- sum(ai_vi$ai[, 1] * ai_vi$vi[, 1])/nT
  asy <- FunD:::FD_est(ai_vi, q, nT)$est
  obs <- FunD:::FD_mle(ai_vi, q)
  out <- sapply(1:ncol(ai_vi$ai), function(i) {
    ai <- ai_vi$ai[, i]
    ai[ai < 1] <- 1
    av = cbind(ai = round(ai), vi = ai_vi$vi[, i])
    RFD_m = FunD:::RFD(av, nT, nT - 1, q, V_bar)
    beta <- rep(0, length(q))
    asy_i <- asy[, i]
    obs_i <- obs[, i]
    asy_i <- sapply(1:length(q), function(j) {
      max(asy_i[j], obs_i[j], RFD_m[j])
    })

    obs_i <- sapply(1:length(q), function(j) {
      max(RFD_m[j], obs_i[j])
    })

    beta0plus <- which(asy_i != obs_i)
    beta[beta0plus] <- (obs_i[beta0plus] - RFD_m[beta0plus])/(asy_i[beta0plus] - RFD_m[beta0plus])
    sapply(m, function(mm) {
      if (mm < nT) {
        FunD:::RFD(av, nT, mm, q, V_bar)
      }
      else if (mm == nT) {
        obs_i
      }
      else {
        EFD(m = mm, qs = q, obs = obs_i, asy = asy_i,
            beta = beta, av = av)
      }
    }) %>% t %>% as.numeric
  })
  matrix(out, ncol = ncol(ai_vi$ai))
}
Bootstrap_distance_matrix = function(data, distance_matrix, f0.hat, datatype){

  if (datatype=="incidence_freq") {
    n = data[1]
    X = data[-1]
    u = sum(data)
  } else if (datatype=="abundance") {
    n = sum(data)
    X = data
  }

  # n = sum(data)
  distance = as.matrix(distance_matrix)
  dij = distance
  # X = data

  F.1 <- sum(dij[, X==1]) ; F.2 <- sum(dij[, X==2])
  F11 <- sum(dij[X==1, X==1]) ; F22 <- sum(dij[X==2, X==2])

  if (datatype=="abundance") {
    F.0hat <- ifelse(F.2 > 0, ((n-1)/n) * (F.1^2/(2 * F.2)), ((n-1)/n)*(F.1*(F.1-0.01)/(2)))
    F00hat <- ifelse(F22 > 0, ((n-2)* (n-3)* (F11^2)/(4* n* (n-1)* F22)), ((n-2)* (n-3)* (F11*(F11-0.01))/(4 *n * (n-1))) )
  } else if (datatype=="incidence_freq") {
    F.0hat <- ifelse(F.2 > 0, ((n-1)/n) * (F.1^2/(2 * F.2)), ((n-1)/n)*(F.1*(F.1-0.01)/(2)))
    F00hat <- ifelse(F22 > 0, ((n-1)^2 * (F11^2)/(4* n* n* F22)), ((n-1)* (n-1)* (F11*(F11-0.01))/(4 *n * n)) )
  }

  # F.0hat <- ifelse(F.2 > 0, ((n-1)/n) * (F.1^2/(2 * F.2)), ((n-1)/n)*(F.1*(F.1-0.01)/(2)))
  # F00hat <- ifelse(F22 > 0, ((n-2)* (n-3)* (F11^2)/(4* n* (n-1)* F22)), ((n-2)* (n-3)* (F11*(F11-0.01))/(4 *n * (n-1))) )

  if (f0.hat==0) {
    d=dij
  } else if (f0.hat==1) {
    random_dij = as.vector(rmultinom(1, 1000, rep(1/(length(X)*f0.hat), length(X)*f0.hat) ) )/1000
    d.0bar <- matrix(random_dij*F.0hat, length(X), f0.hat, byrow = T)
    d00 = matrix(0, f0.hat, f0.hat)
    d <- cbind(dij, d.0bar )
    aa <- cbind(t(d.0bar), d00 )
    d <- rbind(d, aa)
    diag(d) = 0
  } else {
    random_dij = as.vector(rmultinom(1, 1000, rep(1/(length(X)*f0.hat), length(X)*f0.hat) ) )/1000
    d.0bar <- matrix(random_dij*F.0hat, length(X), f0.hat, byrow = T)

    fo.num = (f0.hat * (f0.hat-1) )/2
    random_d00 = as.vector(rmultinom(1, 1000, rep(1/fo.num, fo.num) ) )/1000
    d00 = matrix(0, f0.hat, f0.hat)
    d00[upper.tri(d00)] = (F00hat/2)*random_d00
    d00 <- pmax(d00, t(d00))###signmatrix
    d <- cbind(dij, d.0bar )
    aa <- cbind(t(d.0bar), d00 )
    d <- rbind(d, aa)
    diag(d) = 0
  }

  return(d)

}
