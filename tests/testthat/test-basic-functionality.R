# test-basic-functionality.R
# Basic functionality tests — canonical binary and genetic simulation cases

test_that("twolife_simulation runs and returns expected structure", {
  result <- run_simple_test_simulation()

  # Test class and type
  expect_s3_class(result, "twolife_result")
  expect_type(result, "list")

  # Test main components
  expect_true("summary"    %in% names(result))
  expect_true("survivors"  %in% names(result))
  expect_true("events"     %in% names(result))
  expect_true("parameters" %in% names(result))
  expect_true("spatial"    %in% names(result))

  # Test summary contents
  expect_true("final_population_size" %in% names(result$summary))
  expect_true("total_events"          %in% names(result$summary))
  expect_true("status"                %in% names(result$summary))
  expect_true("duration"              %in% names(result$summary))

  # Test survivors structure
  expect_type(result$survivors, "list")
  expect_true("x"  %in% names(result$survivors))
  expect_true("y"  %in% names(result$survivors))
  expect_true("id" %in% names(result$survivors))
})

test_that("population_size extracts trajectory correctly", {
  result <- run_simple_test_simulation()

  pop_trajectory <- population_size(result)

  # Test trajectory structure
  expect_s3_class(pop_trajectory, "data.frame")
  expect_true("time"            %in% names(pop_trajectory))
  expect_true("population_size" %in% names(pop_trajectory))
  expect_true(nrow(pop_trajectory) > 0)
  expect_true(all(pop_trajectory$population_size >= 0))

  # Test that population values are numeric
  expect_type(pop_trajectory$population_size, "double")
})

test_that("simulation handles genetic parameters", {
  result <- run_genetic_test_simulation()

  # Test genetic fields are present
  expect_true("genotype"  %in% names(result$survivors))
  expect_true("phenotype" %in% names(result$survivors))
  expect_true("width"     %in% names(result$survivors))

  # Test field lengths match final population
  expect_equal(length(result$survivors$genotype),  result$summary$final_population_size)
  expect_equal(length(result$survivors$phenotype), result$summary$final_population_size)
  expect_equal(length(result$survivors$width),     result$summary$final_population_size)

  # Test that genetic values are numeric and reasonable
  if (result$summary$final_population_size > 0) {
    expect_type(result$survivors$genotype,  "double")
    expect_type(result$survivors$phenotype, "double")
    expect_true(all(result$survivors$width > 0))
  }
})

test_that("different history detail levels work", {
  landscape <- create_continuous_test_landscape()

  result_minimal <- twolife_simulation(
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
    history_detail = "minimal"
  )

  result_standard <- twolife_simulation(
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
    history_detail = "standard"
  )

  result_full <- twolife_simulation(
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

  # All should complete successfully
  expect_s3_class(result_minimal,  "twolife_result")
  expect_s3_class(result_standard, "twolife_result")
  expect_s3_class(result_full,     "twolife_result")

  # Minimal should have fewer event details than full
  expect_true(length(names(result_minimal$events)) <= length(names(result_full$events)))
})

test_that("snapshot_at_time extracts snapshots correctly", {
  result <- run_genetic_test_simulation()

  mid_time       <- result$summary$duration / 2
  snapshot_mid   <- snapshot_at_time(result, target_time = mid_time, show_plot = FALSE)
  snapshot_final <- snapshot_at_time(result, target_time = Inf,      show_plot = FALSE)

  # Test snapshot structure
  expect_type(snapshot_mid, "list")
  expect_true("time"       %in% names(snapshot_mid))
  expect_true("population" %in% names(snapshot_mid))
  expect_true("n_alive"    %in% names(snapshot_mid))

  # n_alive must be internally consistent with the returned population data frame.
  # Note: snapshot_at_time reconstructs state from the event log, which may not
  # exactly match result$summary$final_population_size (that comes directly from
  # the C++ output). The reconstruction is a best-effort replay, not a guarantee.
  expect_equal(snapshot_final$n_alive, nrow(snapshot_final$population))
  expect_true(snapshot_final$n_alive >= 0)
})

test_that("binary and continuous landscapes both work", {
  result_binary  <- run_simple_test_simulation()
  result_genetic <- run_genetic_test_simulation()

  expect_s3_class(result_binary,  "twolife_result")
  expect_s3_class(result_genetic, "twolife_result")
})

test_that("master_seed ensures reproducibility", {
  landscape <- create_continuous_test_landscape()

  run_once <- function() {
    twolife_simulation(
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
      master_seed = 12345
    )
  }

  result1 <- run_once()
  result2 <- run_once()

  expect_equal(result1$summary$final_population_size, result2$summary$final_population_size)
  expect_equal(result1$summary$total_events,          result2$summary$total_events)
  expect_equal(result1$summary$duration,              result2$summary$duration)
})
