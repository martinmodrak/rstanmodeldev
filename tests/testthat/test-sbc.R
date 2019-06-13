test_that("sbc runs", {
  simple_model = rstan::stan_model(model_code = "
    data {
      real x[5];
    }

    parameters {
      real mu;
    }

    model {
      x ~ normal(mu, 1);
      mu ~ normal(0, 1);
    }
    ")

  generator <- function() {
    mu = rnorm(1, 0, 1)

    list(
      true = list(
        mu = mu
      ),
      observed = list(
        x = rnorm(5, mu, 1)
      )
    )
  }

  sbc_res <- sbc(simple_model, generator, N_steps = 2)
})
