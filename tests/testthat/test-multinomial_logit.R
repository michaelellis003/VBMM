test_that("multinomial_logit function works", {

    probs <- list()
    probs[[1]] <- matrix(c(1/6, 1/6, 2/3,
                      1/6, 2/3, 1/6,
                      2/3, 1/6, 1/6),
                    nrow = 3, ncol = 3, byrow = TRUE)

    probs[[2]] <- matrix(c(0, 0, 0, 1,
                      0, 0, 1, 0,
                      0, 1, 0, 0,
                      1, 0, 0, 0),
                    nrow = 4, ncol = 4, byrow = TRUE)

    probs[[3]] <- matrix(c(1/2, 1/2, 0, 0,
                           1/2, 0, 1/2, 0,
                           1/2, 0, 0, 1/2,
                           0, 0, 1/2, 1/2),
                         nrow = 4, ncol = 4, byrow = TRUE)

    for(i in 1:3) {
        dens <- probs[[i]] * rnorm(1, mean = 1000, sd = 10)
        log_dens <- log(dens)

        output <- multinomial_logit(log_dens)

        expect_equal(probs[[i]], output)
    }

})
