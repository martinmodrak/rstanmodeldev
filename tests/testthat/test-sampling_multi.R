test_that("sampling multi runs", {
  simple_model_code = "
    data {
      real x[5];
    }

    parameters {
      real mu;
    }

    model {
      x ~ normal(mu, 1);
    }
    "
  simple_model = rstan::stan_model(model_code = simple_model_code)

  second_simple_model = rstan::stan_model(model_code = "
    data {
      real x[7];
    }

    parameters {
      real<lower=0> sigma;
    }

    model {
      x ~ normal(0, sigma);
    }
    ")

  num_data = 10
  data_list = list()
  for(i in 1:num_data) {
    data_list[[i]] = list(x = rnorm(5, i, 1))
  }

  res = sampling_multi(simple_model, data_list)

})
