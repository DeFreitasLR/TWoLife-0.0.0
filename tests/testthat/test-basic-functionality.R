# test-basic-functionality.R
# Basic functionality tests for TWoLife package

test_that("Package loads without error", {
  expect_true(require(TWoLife, quietly = TRUE))
})

test_that("Basic simulation runs without error", {
  # Create a simple landscape
  test_landscape <- create_test_landscape(size = 10, habitat_prop = 0.5)
  
  # Run a short simulation
  expect_error({
    result <- twolife_simulation(
      landscape_params = list(habitat = test_landscape),
      individual_params = list(
        initial_population_size = 10,
        base_birth_rate = 0.35,
        base_mortality_rate = 0.25
      ),
      simulation_params = list(max_events = 100)
    )
  }, NA)  # NA means "no error expected"
})

test_that("Simulation returns expected structure", {
  # Use helper function for cleaner test
  result <- run_simple_test_simulation(steps = 100, n = 10)
  
  # Check that result is the right class
  expect_s3_class(result, "twolife_result")
  
  # Check for required components
  expect_true("summary" %in% names(result))
  expect_true("survivors" %in% names(result))
  expect_true("spatial" %in% names(result))
  expect_true("events" %in% names(result))
  expect_true("parameters" %in% names(result))
  
  # Check summary fields
  expect_true("final_population_size" %in% names(result$summary))
  expect_true("total_events" %in% names(result$summary))
  expect_true("duration" %in% names(result$summary))
  expect_true("status" %in% names(result$summary))
})

test_that("Simulation completes successfully with max_events limit", {
  # max_events prevents infinite loops but includes initialization events
  result <- run_simple_test_simulation(steps = 50, n = 5)
  
  # Simulation should complete successfully
  expect_s3_class(result, "twolife_result")
  expect_true("summary" %in% names(result))
  
  # Should have some events recorded (at minimum, initialization)
  expect_gt(result$summary$total_events, 0)
  
  # Note: total_events includes initialization events (placing initial individuals)
  # so it may exceed max_events parameter, which limits simulation events only
})

test_that("Population size can be computed", {
  result <- run_simple_test_simulation(steps = 100, n = 10)
  
  pop_trajectory <- compute_population_size(result)
  
  expect_true(is.data.frame(pop_trajectory))
  expect_true("time" %in% names(pop_trajectory))
  expect_true("population_size" %in% names(pop_trajectory))
  expect_true(nrow(pop_trajectory) > 0)
})