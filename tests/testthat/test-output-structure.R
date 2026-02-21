# test-output-structure.R
# Tests for simulation output structure — canonical binary and genetic simulation cases

test_that("output has expected fields", {
  result <- run_simple_test_simulation()

  # Check main structure
  expect_named(result, c("summary", "survivors", "spatial", "events", "parameters"))

  # Check summary structure
  expect_named(result$summary, c("final_population_size", "total_events", "duration", "status"))

  # Check survivors structure (if any survived)
  if (result$summary$final_population_size > 0) {
    expect_true(is.data.frame(result$survivors))
    expect_true(all(c("id", "x", "y", "genotype", "phenotype", "width") %in% names(result$survivors)))
  }

  # Check events structure
  expect_true("times"          %in% names(result$events))
  expect_true("types"          %in% names(result$events))
  expect_true("individual_ids" %in% names(result$events))

  # Check parameters structure
  expect_true(all(c("landscape", "individual", "genetic", "simulation") %in% names(result$parameters)))
})

test_that("survivors data frame has correct structure", {
  # Genetic simulation for meaningful trait and coordinate checks
  result <- run_genetic_test_simulation()

  if (result$summary$final_population_size > 0) {
    # Check it's a data frame
    expect_s3_class(result$survivors, "data.frame")

    # Check number of rows matches population
    expect_equal(nrow(result$survivors), result$summary$final_population_size)

    # Check column types
    expect_type(result$survivors$id,        "integer")
    expect_type(result$survivors$x,         "double")
    expect_type(result$survivors$y,         "double")
    expect_type(result$survivors$genotype,  "double")
    expect_type(result$survivors$phenotype, "double")
    expect_type(result$survivors$width,     "double")

    # Check coordinate ranges against world geometry
    half_w <- result$spatial$world_width  / 2
    half_h <- result$spatial$world_height / 2
    expect_true(all(result$survivors$x >= -half_w & result$survivors$x <= half_w))
    expect_true(all(result$survivors$y >= -half_h & result$survivors$y <= half_h))

    # Genotype and phenotype must lie within the landscape value range [0, 100]
    expect_true(all(result$survivors$genotype  >= 0 & result$survivors$genotype  <= 100))
    expect_true(all(result$survivors$phenotype >= 0 & result$survivors$phenotype <= 100))
    expect_true(all(result$survivors$width > 0))
  }
})

test_that("events vectors have consistent lengths", {
  result <- run_simple_test_simulation()

  n_events <- length(result$events$times)

  # All event vectors should have same length
  expect_equal(length(result$events$types),          n_events)
  expect_equal(length(result$events$individual_ids), n_events)

  # For standard history detail
  if ("x_coordinates" %in% names(result$events)) {
    expect_equal(length(result$events$x_coordinates), n_events)
    expect_equal(length(result$events$y_coordinates), n_events)
    expect_equal(length(result$events$genotypes),     n_events)
  }

  # Check that event times are monotonically increasing
  if (n_events > 1) {
    expect_true(all(diff(result$events$times) >= 0))
  }
})

test_that("spatial history has correct structure", {
  landscape <- create_continuous_test_landscape()

  result <- twolife_simulation(
    landscape_params = landscape,
    individual_params = list(
      initial_population_size = 1000,
      neighbor_radius         = 2.0,
      vision_angle            = pi,
      step_length             = 5.0,
      base_dispersal_rate     = 0.4,
      base_birth_rate         = 0.6,
      base_mortality_rate     = 0.2,
      birth_density_slope     = 0.02,
      mortality_density_slope = 0.02
    ),
    simulation_params = list(max_events = 100),
    master_seed    = 45,
    history_detail = "full"
  )

  # result$spatial contains landscape geometry
  expect_true("world_width"  %in% names(result$spatial))
  expect_true("world_height" %in% names(result$spatial))
  expect_true("world_size"   %in% names(result$spatial))

  # World dimensions should be positive numbers
  expect_true(result$spatial$world_width  > 0)
  expect_true(result$spatial$world_height > 0)
  expect_true(result$spatial$world_size   > 0)
})

test_that("parameters are stored correctly", {
  landscape <- create_continuous_test_landscape()

  # Override matrix_mortality_multiplier to a distinct value for verification
  landscape$matrix_mortality_multiplier <- 2.5

  result <- twolife_simulation(
    landscape_params = landscape,
    individual_params = list(
      initial_population_size = 1000,
      neighbor_radius         = 2.0,
      vision_angle            = pi,
      step_length             = 5.0,
      base_dispersal_rate     = 0.4,
      base_birth_rate         = 0.6,
      base_mortality_rate     = 0.2,
      birth_density_slope     = 0.02,
      mortality_density_slope = 0.02
    ),
    genetic_params = list(
      genotype_means = rep(50, 1000),
      genotype_sds   = 0.4,
      mutation_rates = 0.01,
      plasticities   = 0.05
    ),
    simulation_params = list(
      max_events  = 100,
      neutral_mode = FALSE
    ),
    master_seed = 49
  )

  # Check individual parameters
  expect_equal(result$parameters$individual$initial_population_size, 1000)
  expect_equal(result$parameters$individual$base_birth_rate,         0.6)
  expect_equal(result$parameters$individual$base_mortality_rate,     0.2)

  # matrix_mortality_multiplier is a landscape param, stored under parameters$landscape
  expect_equal(result$parameters$landscape$matrix_mortality_multiplier, 2.5)

  # Genetic parameters are expanded to full population-length vectors internally
  expect_equal(length(result$parameters$genetic$genotype_means), 1000)
  expect_equal(length(result$parameters$genetic$genotype_sds),   1000)
  expect_true(all(result$parameters$genetic$genotype_sds  == 0.4))
  expect_true(all(result$parameters$genetic$mutation_rates == 0.01))
  expect_true(all(result$parameters$genetic$plasticities   == 0.05))

  # Check simulation parameters
  expect_equal(result$parameters$simulation$max_events,   100)
  expect_equal(result$parameters$simulation$neutral_mode, FALSE)

  # Check landscape matrix is 5x5
  expect_true(is.matrix(result$parameters$landscape$habitat))
  expect_equal(nrow(result$parameters$landscape$habitat), 5)
  expect_equal(ncol(result$parameters$landscape$habitat), 5)
})

test_that("summary statistics are reasonable", {
  result <- run_simple_test_simulation()

  # Check summary values are non-negative
  expect_true(result$summary$final_population_size >= 0)
  expect_true(result$summary$total_events          >= 0)
  expect_true(result$summary$duration              >= 0)

  # Check final population is integer
  expect_equal(result$summary$final_population_size,
               as.integer(result$summary$final_population_size))

  # Check status is valid
  expect_true(result$summary$status %in% c("surviving", "extinct"))

  # If population survived, check survivors match summary
  if (result$summary$final_population_size > 0) {
    expect_equal(nrow(result$survivors), result$summary$final_population_size)
  }
})
