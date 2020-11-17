# Metropolis-Hastings Algorithm

```{r}

# Simulate data from mixture normal with $\mu_1=0$ and $\mu_2=2.5$
da = rbind(rnorm(10^3), 2.5+rnorm(3*10^3))
hist(da)

# Define the likelihood function
like = function(mu){
  sum(log((.25*dnorm(da-mu[1]) + .75*dnorm(da-mu[2]))))
}

```

# Objective 
- Use Metropolis-Hastings Alg to estimate the values of $\mu_1$ and $\mu_2$.
- Assume a uniform prior $U(-2, 5)$ for both $\mu_1$ and $\mu_2$. This determines the space state for the Markov Chain. 
- Use random walk as the transition kernel.

```{r}

# Define the state space of the Markov Chain
the = matrix(runif(2, -2, 5), ncol = 2)

# Standard deviation of the random walk
scale = 0.2

mu = c(0, 2.5)
curlike = hval = like(mu)

Niter = 10^4
for (iter in (1:Niter)){
  prop = the[iter,] + rnorm(2) * scale
  if ((max(-prop)>2) || (max(prop)>5) ||  # Proposed sample out of state space
      (log(runif(1)) > like(prop) - curlike)) # Accept with probability
    prop = the[iter,]
  curlike = like(prop)
  
  # Add new MH sample to simulations
  hval = c(hval, curlike)
  the = rbind(the, prop) 
}
```


# Plots

```{r}
hist(the[,1])
plot(the[,1])

hist(the[,2])
plot(the[,2])

hist(the)
```
# Add gradient drift to the code
```{r}
gradlike = function(mu){
  deno = .25 * dnorm(da-mu[1]) + .75 * dnorm(da - mu[2])
  gra = sum(.25 * (da - mu[1]) * dnorm(da - mu[1])/deno)
  grb = sum(.75 * (da - mu[2]) * dnorm(da - mu[2])/deno)
  return(c(gra, grb))
}
```

# Modified MH code with the gradient drift

```{r}

# Standard deviation of the random walk
scale = 0.2

mu = c(0, 2.5)
curlike = hval = like(mu)
curmean = c(0,0)

Niter = 10^4
for (iter in (1:Niter)){
  prop = rnorm(2) * scale
  meanprop = prop + .5 * scale^2 * gradlike(prop)
  if ((max(-prop)>2) || (max(prop)>5) ||  # Proposed sample out of state space
      (log(runif(1)) > like(prop) 
       - curlike - sum(dnorm(prop, curmean, lo=T)) +
       sum(dnorm(the[iter, ], meanprop, lo=T)))) {
          prop = the[iter,]
          meanprop = curmean
  }
  curlike = like(prop)
  curmean = meanprop

  # Add new MH sample to simulations
  hval = c(hval, curlike)
  the = rbind(the, prop) 
}

```

# Plots

```{r}
hist(the[,1])
plot(the[,1])

hist(the[,2])
plot(the[,2])

hist(the)
```
