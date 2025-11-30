#' Example Landscape: Multiple Small Habitat Patches
#'
#' A 15x15 binary landscape matrix representing fragmented habitat with
#' several small patches of suitable habitat (value = 1) scattered across
#' a matrix landscape (value = 0).
#'
#' @format A 15x15 numeric matrix with binary values:
#' \describe{
#'   \item{0}{Matrix (non-habitat)}
#'   \item{1}{Suitable habitat}
#' }
#'
#' @details
#' This landscape demonstrates habitat fragmentation with multiple small,
#' isolated patches. It is useful for testing:
#' \itemize{
#'   \item Population persistence in fragmented landscapes
#'   \item Dispersal between isolated patches
#'   \item Metapopulation dynamics
#' }
#'
#' The landscape has approximately 22% habitat coverage distributed across
#' several distinct patches.
#'
#' @examples
#' # Load and examine the dataset
#' data(several_small)
#' dim(several_small)
#' table(several_small)
#'
#' # Visualize the landscape
#' plot_landscape(several_small,
#'                main = "Fragmented Habitat",
#'                colors = "habitat")
#'
#' # Run a simulation on this landscape
#' result <- twolife_simulation(
#'   landscape_params = list(habitat = several_small),
#'   individual_params = list(initial_population_size = 15),
#'   simulation_params = list(max_events = 200),
#'   master_seed = 111
#' )
#'
#' # Check if population persisted
#' result$summary$final_population_size
#'
#' @keywords datasets
"several_small"


#' Example Landscape: Single Large Habitat Patch
#'
#' A 15x15 binary landscape matrix representing a single large, contiguous
#' patch of suitable habitat (value = 1) surrounded by matrix (value = 0).
#'
#' @format A 15x15 numeric matrix with binary values:
#' \describe{
#'   \item{0}{Matrix (non-habitat)}
#'   \item{1}{Suitable habitat}
#' }
#'
#' @details
#' This landscape demonstrates a single large habitat patch configuration.
#' It is useful for testing:
#' \itemize{
#'   \item Population dynamics in continuous habitat
#'   \item Comparison with fragmented landscapes (SLOSS debate)
#'   \item Baseline population persistence
#' }
#'
#' The landscape has approximately 30% habitat coverage in a single
#' contiguous patch.
#'
#' @examples
#' # Load and examine the dataset
#' data(single_large)
#' dim(single_large)
#' table(single_large)
#'
#' # Visualize the landscape
#' plot_landscape(single_large,
#'                main = "Single Large Patch",
#'                colors = "habitat")
#'
#' # Run a simulation on this landscape
#' result <- twolife_simulation(
#'   landscape_params = list(habitat = single_large),
#'   individual_params = list(initial_population_size = 20),
#'   simulation_params = list(max_events = 300),
#'   master_seed = 222
#' )
#'
#' # Check if population persisted
#' result$summary$final_population_size
#'
#' @keywords datasets
"single_large"
