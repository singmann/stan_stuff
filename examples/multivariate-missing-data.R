
library("rstan")
options(mc.cores = parallel::detectCores())


##----------------------------------------------------------------
##      Generate Multivariate Data With Missing Observations     -
##----------------------------------------------------------------


set.seed(567891234)
total_n <- 200
k <- 4
prop_missing <- 0.5
Sigma <- randcorr::randcorr(k)
mu <- rnorm(k)

dat_all <- MASS::mvrnorm(total_n, mu, Sigma)
dat_reduced <- dat_all  ## copy from which to remove observations

index <- rep(c(TRUE, FALSE), times = round(c(prop_missing, 1-prop_missing) * total_n))
select <- replicate(n = 3, sample(index))
dat_reduced[select] <- NA_real_

apply(dat_reduced, 2, function(x) mean(is.na(x)))

cor(dat_reduced, use = "pair")
cor(dat_reduced, use = "complete")

cor(dat_all)

## TRUE:
Sigma

##---------------------------------------------------------------
##                  Prepare Data and Stan Model                 -
##---------------------------------------------------------------

source("functions/functions-multi-missing.R")


dat_prep <- prep_data(as.data.frame(dat_reduced)) ## needs to be a data.frame
str(dat_prep, 1)

## make stan model:

zz <- file("models/multi-missing.stan", "w")

make_data_declaration(dat_prep, file = zz)
make_parameters_block(dat_prep, file = zz)
make_trans_parameters_block(dat_prep, file = zz)

make_model_block(dat_prep, file = zz, lkj_eta = 1, t_nu = 3, t_scale = 3.5)

close.connection(zz)

mod_miss <- stan(file = "models/multi-missing.stan", 
                 data = dat_prep)
print(mod_miss, pars = "Omega")
plot(mod_miss, pars = "Omega")

cor(dat_all)
cor(dat_reduced, use = "pair")
