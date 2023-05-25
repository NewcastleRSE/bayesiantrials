model {

  for (i in 1:N){

    Y[i] ~ dnorm(mu[i],taub)

    mu[i] = alpha + X[i]*delta + c[Cl[i]]


  }

  for (j in 1:Nclust){

    c[j] ~ dnorm(0,tauw)

  }

  alpha ~ dnorm(1,0.001)
  delta ~ dnorm(0,0.001)

  taub ~ dgamma(eps,eps)
  tauw ~ dgamma(eps,eps)

}
