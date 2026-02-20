# test-basic-functionality.R
# Basic functionality tests using consistent parameters and seeding

test_that("twolife_simulation runs and returns expected structure", {
  # Use standard continuous landscape (seed 200)
  landscape <- create_continuous_test_landscape()
  
  # Run simple simulation (seed 300)
  set.seed(300)
  result <- twolife_simulation(
    landscape_params = landscape,
    individual_params = list(
      initial_population_size = 50,
      base_birth_rate = 0.5,
      base_mortality_rate = 0.2
    ),
    simulation_params = list(max_events = 100)
  )
  
  # Test structure
  expect_s3_class(result, "twolife_result")
  expect_type(result, "list")
  
  # Test main components
  expect_true("summary" %in% names(result))
  expect_true("survivors" %in% names(result))
  expect_true("events" %in% names(result))
  expect_true("parameters" %in% names(result))
  expect_true("spatial" %in% names(result))
  
  # Test summary contents
  expect_true("final_population_size" %in% names(result$summary))
  expect_true("total_events" %in% names(result$summary))
  expect_true("status" %in% names(result$summary))
  expect_true("duration" %in% names(result$summary))
  
  # Test survivors structure
  expect_type(result$survivors, "list")
  expect_true("x" %in% names(result$survivors))
  expect_true("y" %in% names(result$survivors))
  expect_true("id" %in% names(result$survivors))
})

test_that("population_size extracts trajectory correctly", {
  # Run simple test simulation
  result <- run_simple_test_simulation(
    population = 50,
    events = 100,
    seed = 300
  )
  
  # Extract population trajectory
  pop_trajectory <- population_size(result)
  
  # Test trajectory structure
  expect_s3_class(pop_trajectory, "data.frame")
  expect_true("time" %in% names(pop_trajectory))
  expect_true("population" %in% names(pop_trajectory))
  expect_true(nrow(pop_trajectory) > 0)
  expect_true(all(pop_trajectory$population >= 0))
  
  # Test that population values are integers
  expect_type(pop_trajectory$population, "integer")
})

test_that("simulation handles genetic parameters", {
  # Run genetic test simulation
  result <- run_genetic_test_simulation(
    population = 50,
    events = 100,
    seed = 400
  )
  
  # Test genetic fields are present
  expect_true("genotype" %in% names(result$survivors))
  expect_true("phenotype" %in% names(result$survivors))
  expect_true("width" %in% names(result$survivors))
  
  # Test field lengths match final population
  expect_equal(length(result$survivors$genotype), result$summary$final_population_size)
  expect_equal(length(result$survivors$phenotype), result$summary$final_population_size)
  expect_equal(length(result$survivors$width), result$summary$final_population_size)
  
  # Test that genetic values are numeric and reasonable
  if (result$summary$final_population_size > 0) {
    expect_type(result$survivors$genotype, "double")
    expect_type(result$survivors$phenotype, "double")
    expect_true(all(result$survivors$width > 0))
  }
})

test_that("different history detail levels work", {
  landscape <- create_continuous_test_landscape()
  
  # Test minimal history
  set.seed(500)
  result_minimal <- twolife_simulation(
    landscape_params = landscape,
    individual_params = list(initial_population_size = 30),
    simulation_params = list(max_events = 50),
    history_detail = "minimal"
  )
  
  # Test standard history
  set.seed(500)
  result_standard <- twolife_simulation(
    landscape_params = landscape,
    individual_params = list(initial_population_size = 30),
    simulation_params = list(max_events = 50),
    history_detail = "standard"
  )
  
  # Test full history
  set.seed(500)
  result_full <- twolife_simulation(
    landscape_params = landscape,
    individual_params = list(initial_population_size = 30),
    simulation_params = list(max_events = 50),
    history_detail = "full"
  )
  
  # All should complete successfully
  expect_s3_class(result_minimal, "twolife_result")
  expect_s3_class(result_standard, "twolife_result")
  expect_s3_class(result_full, "twolife_result")
  
  # Minimal should have fewer event details
  expect_true(length(names(result_minimal$events)) <= length(names(result_full$events)))
})

test_that("snapshot_at_time extracts snapshots correctly", {
  # Run simulation
  result <- run_simple_test_simulation(
    population = 50,
    events = 100,
    seed = 300
  )
  
  # Get mid-simulation snapshot
  mid_time <- result$summary$duration / 2
  snapshot_mid <- snapshot_at_time(result, target_time = mid_time)
  
  # Get final snapshot
  snapshot_final <- snapshot_at_time(result, target_time = Inf)
  
  # Test snapshot structure
  expect_type(snapshot_mid, "list")
  expect_true("time" %in% names(snapshot_mid))
  expect_true("individuals" %in% names(snapshot_mid))
  
  # Test final snapshot matches survivors
  expect_equal(nrow(snapshot_final$individuals), result$summary$final_population_size)
})

test_that("binary and continuous landscapes both work", {
  # Binary landscape
  landscape_binary <- create_binary_test_landscape()
  set.seed(600)
  result_binary <- twolife_simulation(
    landscape_params = landscape_binary,
    individual_params = list(
      initial_population_size = 30,
      base_birth_rate = 0.5,
      base_mortality_rate = 0.2,
      matrix_mortality_multiplier = 3.0
    ),
    simulation_params = list(max_events = 50)
  )
  
  # Continuous landscape
  landscape_continuous <- create_continuous_test_landscape()
  set.seed(600)
  result_continuous <- twolife_simulation(
    landscape_params = landscape_continuous,
    individual_params = list(
      initial_population_size = 30,
      base_birth_rate = 0.5,
      base_mortality_rate = 0.2
    ),
    simulation_params = list(max_events = 50)
  )
  
  # Both should complete
  expect_s3_class(result_binary, "twolife_result")
  expect_s3_class(result_continuous, "twolife_result")
})

test_that("master_seed ensures reproducibility", {
  landscape <- create_continuous_test_landscape()
  
  # Run same simulation twice with same seed
  result1 <- twolife_simulation(
    landscape_params = landscape,
    individual_params = list(initial_population_size = 30),
    simulation_params = list(max_events = 100),
    master_seed = 12345
  )
  
  result2 <- twolife_simulation(
    landscape_params = landscape,
    individual_params = list(initial_population_size = 30),
    simulation_params = list(max_events = 100),
    master_seed = 12345
  )
  
  # Results should be identical
  expect_equal(result1$summary$final_population_size, result2$summary$final_population_size)
  expect_equal(result1$summary$total_events, result2$summary$total_events)
  expect_equal(result1$summary$duration, result2$summary$duration)
})
