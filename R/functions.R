# functions.R - Complete TWoLife R interface
# UPDATED: Added all missing @param documentation and proper imports
# UPDATED: Added show_legend parameter to check_habitat_match functions
# FIXED: Replaced non-ASCII characters with ASCII equivalents

#' @importFrom grDevices colorRampPalette heat.colors terrain.colors
#' @importFrom graphics abline grid image legend mtext par points
#' @importFrom stats cor median quantile rnorm runif sd
#' @importFrom utils modifyList
NULL

# ============================================================================
# CORE SIMULATION FUNCTIONS
# ============================================================================

#' Run TWoLife Individual-Based Simulation
#'
#' Runs a spatially-explicit individual-based simulation with habitat selection,
#' genetic variation, phenotypic plasticity, and demographic processes. Supports
#' rectangular landscapes and customizable habitat selection parameters.
#'
#' @param landscape_params List containing landscape parameters:
#'   \describe{
#'     \item{habitat}{Matrix (integer, logical or numeric). World Habitat Pixel Matrix. Matrix values should be binary (0/1) or continuous habitat values (usually from 0 to 1). Each value represents an environmental value of a pixel. Individuals experience fitness effects based on how well their phenotype matches the habitat value of the pixel they are located.  Matrix dimensions define the landscape size and grid resolution. Required.}
#'     \item{cell_size}{Numeric. Length of the side of a landscape cell in world units. World units are the arbitrary spatial measurement scale for the simulation (e.g., meters, kilometers, or abstract distance units) - all spatial parameters must use the same unit for consistency. Defines the spatial resolution of the landscape. If habitat is a 100x100 matrix and cell_size = 1.0, the simulated world dimensions are 100x100 world units (world_width = 100, world_height = 100). All spatial parameters (step_length, neighbor_radius, coordinates) are expressed in these same world units. Individual positions are continuous (x,y) coordinates and are not restricted to pixel centers.}
#'     \item{boundary_condition}{Integer. Defines what happens to individuals that reach the edges of world borders:
#'       \itemize{
#'         \item 1 = reflective: agent keeps its speed, but its path "bounces" off the border at the same angle it hit it (angle of incidence equals angle of reflection);
#'         \item 2 = absorbing: individuals that cross borders permanently exit the simulation and are removed from the population, generating emigration events (recorded as event type = 3)
#'         \item 3 = periodic: individuals moving beyond one edge wrap to the opposite edge, creating torus topology (infinite landscape approximation)
#'       }}
#'     \item{density_type}{Integer. Determines how population density is calculated for density-dependent demographic processes to be applied to each individual:
#'       \itemize{
#'         \item 1 = local: density calculated as number of individuals within a neighbor_radius distance of the focal individual. Represents local, spatially-explicit local competition (territoriality, local resource depletion).
#'         \item 2 = global: density calculated as total population size divided by total landscape area. Represents population-wide resource limitation.
#'       }
#'       
#'       This density value is then used in birth and mortality rate calculations via birth_density_slope and mortality_density_slope parameters. See Details for mathematical formulas.}
#'     \item{matrix_mortality_multiplier}{Numeric. Mortality rate multiplier applied based on habitat suitability. Controls how mortality scales from optimal to unsuitable habitat. Values > 1 increase mortality in poor-quality habitat, creating "hostile matrix" effects. For perfect specialists (genotype_sds = 0), the multiplier is applied directly in non-optimal habitat. For generalists (genotype_sds > 0), mortality is interpolated based on fitness, creating a smooth gradient. See Details for mathematical formulas.}
#'     \item{matrix_dispersal_multiplier}{Numeric. Dispersal rate multiplier applied based on habitat suitability. Controls the frequency of dispersal events in unsuitable habitat. Values < 1 reduce dispersal frequency in matrix, while values > 1 increase it. For perfect specialists (genotype_sds = 0), the multiplier is applied in non-optimal habitat. For generalists (genotype_sds > 0), dispersal rate remains at base_dispersal_rate regardless of habitat. Important: This affects the rate of dispersal events, NOT the distance (step_length). See Details for mathematical formulas.}
#'   }
#' @param individual_params List containing individual-level parameters:
#'   \describe{
#'     \item{initial_population_size}{Integer. Number of individuals at simulation start.}
#'     \item{neighbor_radius}{Numeric. Distance (in world units) within which other individuals are counted as neighbors for local density calculations. Only used when density_type = 1 (local). See Details section 'Density Calculations' for formula.}
#'     \item{vision_angle}{Numeric. Angular range (in radians) within which an individual can change direction during random walk dispersal (sampling_points = 0). See Details section 'Habitat Selection' for movement formulas.}
#'     \item{step_length}{Numeric. Maximum distance (in world units) an individual moves during each dispersal event. For random walk (sampling_points = 0), individual moves exactly this distance in chosen direction. For habitat selection (sampling_points > 0), defines the radius for sampling candidate locations. See Details section 'Habitat Selection' for formulas.}
#'     \item{base_dispersal_rate}{Numeric. Baseline probability per time unit that a dispersal event occurs (usual range: 0-1). Modified by habitat suitability for perfect specialists (genotype_sds = 0) via matrix_dispersal_multiplier.}
#'     \item{base_birth_rate}{Numeric. Baseline probability per time unit that a birth event occurs (usual range: 0-1). Modified by density-dependence (see birth_density_slope) and by habitat suitability for perfect specialists (specialists cannot reproduce in non-optimal habitat). See Details section 'Birth Rate Calculations' for formulas.}
#'     \item{base_mortality_rate}{Numeric. Baseline probability per time unit that a mortality event occurs (usual range: 0-1). Mathematical application shown in matrix_mortality_multiplier parameter.}
#'     \item{birth_density_slope}{Numeric. Controls the strength of negative density-dependence on birth rate. Higher values cause birth rate to decrease more rapidly as local density increases. See Details section 'Birth Rate Calculations' for formula.}
#'     \item{mortality_density_slope}{Numeric. Controls the strength of positive density-dependence on mortality rate. Higher values cause mortality rate to increase more rapidly as local density increases. Applied only to the baseline model (genotype_sds = 0). See Details sections 'Matrix Mortality Multiplier' and 'Density Calculations' for formulas.}
#'     \item{generalism_cost_scale}{Numeric. Half-saturation parameter (k) controlling the strength of the specialist-generalist mortality trade-off by modulating baseline mortality costs at the environmental optimum. See Details for mathematical formulas.}
#'     \item{initial_placement_mode}{Integer. Determines how individuals are positioned at simulation start:
#'       \itemize{
#'         \item 1 = uniform random placement 
#'         \item 2 = random placement following bivariate normal distribution centered on landscape center, with \eqn{\sigma_{placement} = L \sqrt{d_0}}
#'         \item 3 = custom coordinates (requires initial_x_coordinates and initial_y_coordinates)
#'       }
#'       
#'       Where \eqn{L} = step_length and \eqn{d_0} = base_dispersal_rate
#'}
#'     \item{initial_x_coordinates}{Numeric vector. Custom x-coordinates for initial individual positions. Required if initial_placement_mode = 3. Length must equal initial_population_size.}
#'     \item{initial_y_coordinates}{Numeric vector. Custom y-coordinates for initial individual positions. Required if initial_placement_mode = 3. Length must equal initial_population_size.}
#'   }
#' @param genetic_params List containing genetic parameters. Each parameter can be either a single value (applied to all individuals) or a vector with length equal to initial_population_size (one value per individual). If a vector is provided with length different from initial_population_size, an error is raised. R's automatic vector recycling is not used to avoid unintended parameter assignments:
#'   \describe{
#'     \item{genotype_means}{Numeric. The genetic environmental optimum value(s) for each individual. Represents the habitat value at which fitness is maximized for each genotype. Should be on the same scale as habitat values (e.g., 0-1). Can be single value (all individuals same genotype) or vector (different genotypes per individual). See Details section 'Fitness Function' for how genotype determines fitness.}
#'     \item{genotype_sds}{Numeric. Niche width parameter controlling tolerance to habitat mismatch. When = 0, individual is perfect specialist (fitness = 1 only at exact optimum, enables matrix_dispersal_multiplier effects). When > 0, individual is generalist (fitness decreases gradually with habitat deviation, larger values = broader tolerance). See Details section 'Fitness Function' for mathematical formula.}
#'     \item{mutation_rates}{Numeric. Standard deviation of mutations added to offspring genotype at each birth event. When = 0, offspring genotype identical to parent. When > 0, offspring genotype varies from parent enabling genetic evolution. Mutations are heritable (passed to future generations). See Details section 'Genetic Mechanisms' for formula.}
#'     \item{plasticities}{Numeric. Controls non-heritable phenotypic variation. Standard deviation of phenotypic noise added to genotype at birth (once per individual, fixed for lifetime). When = 0, phenotype equals genotype. When > 0, phenotype varies around genotype. Phenotype (not genotype) determines fitness and habitat selection. Plasticity is NOT heritable (unlike mutation). See Details section 'Genetic Mechanisms' for formula and comparison with mutation.}
#'     \item{sampling_points}{Integer. Number of candidate locations sampled within step_length distance during each dispersal event. When sampling_points = 0, individuals perform vision_angle-constrained random walk. When sampling_points > 0, individuals perform habitat selection by sampling candidate locations and choosing the best according to the habitat_selection_temperatures parameter. Larger values enable more informed habitat choice but slower computation. See Details for mathematical formulas.}
#'     \item{habitat_selection_temperatures}{Positive real. Temperature parameter for softmax function in habitat selection (sampling_points > 0). Controls strength of preference for high-quality habitat. Lower values -> stronger selection for best habitat (nearly deterministic). Higher values -> more random exploration. T = 1 gives balanced selection proportional to relative fitness. See Details section 'Habitat Selection' for mathematical formula.}
#'   }
#' @param simulation_params List containing simulation control parameters:
#'   \describe{
#'     \item{max_events}{Numeric. Maximum simulation time to run before stopping. Time advances continuously in the Gillespie algorithm based on exponentially distributed waiting times drawn from individual event rates (birth, death, dispersal). Higher demographic rates lead to faster time progression. The simulation stops when world_time reaches max_events, population goes extinct, or population exceeds 1,000,000. See Details for information on the Gillespie algorithm and time advancement.}
#'     \item{neutral_mode}{Logical. If TRUE, creates a null model without genetic variation or habitat selection. Effects: (1) all individuals are assigned the mean genotype of the initial population with no variation. Demographic processes (birth, death, dispersal rates) and movement mechanics remain unchanged. Useful for null model comparisons to isolate effects of genetic variation and habitat selection.}
#'   }
#' @param history_detail Character. Level of detail recorded in event history outputed by the function. Higher detail levels enable more analyses but use more memory.
#'
#'   Options:
#'   \itemize{
#'     \item minimal: Records only time, event type, and individual ID. Fastest, smallest memory footprint. Sufficient for population size trajectories.
#'     \item standard: Adds spatial coordinates (x, y), pixel ID (the landscape grid cell where the individual is located), and genotype for each event. Enables spatial and genetic analyses. Recommended for most uses.
#'     \item full: Adds phenotype and niche width for each event. The niche width is the phenotypic standard deviation (genotype_sds) of each individual at each moment, representing the breadth of habitats where the individual maintains high fitness. Enables complete historical reconstruction. Use when analyzing plasticity or detailed evolutionary dynamics.
#'   }
#'
#' @param master_seed Integer. Random seed for reproducible simulations. If NULL, results are stochastic.
#'
#' @section Default Values:
#' When parameters are not specified, the following defaults are used:
#'
#' \preformatted{
#' landscape_params = list(
#'   habitat = NULL,               # Required - must provide matrix
#'   cell_size = 1.0,
#'   boundary_condition = 1,
#'   density_type = 1,
#'   matrix_mortality_multiplier = 2.0,
#'   matrix_dispersal_multiplier = 0.5
#' )
#'
#' individual_params = list(
#'   initial_population_size = 200,
#'   neighbor_radius = 2.0,
#'   vision_angle = pi,
#'   step_length = 1.0,
#'   base_dispersal_rate = 0.1,
#'   base_birth_rate = 0.3,
#'   base_mortality_rate = 0.20,
#'   birth_density_slope = 0.02,
#'   mortality_density_slope = 0.02,
#'   generalism_cost_scale = 1.0,
#'   initial_placement_mode = 1,
#'   initial_x_coordinates = NULL,
#'   initial_y_coordinates = NULL
#' )
#'
#' genetic_params = list(
#'   genotype_means = 1,
#'   genotype_sds = 0,
#'   mutation_rates = 0,
#'   plasticities = 0,
#'   sampling_points = 0,
#'   habitat_selection_temperatures = 1.0
#' )
#'
#' simulation_params = list(
#'   max_events = 50 * initial_population_size,
#'   neutral_mode = FALSE
#' )
#'
#' history_detail = "standard"
#' master_seed = NULL
#' }
#'
#' @return A list of class 'twolife_result' with components:
#'   \describe{
#'     \item{summary}{List with status (surviving or extinct), final_population_size (integer),
#'       total_events (integer), and duration (numeric time units)}
#'     \item{survivors}{Data frame with columns: id (integer), x (numeric), y (numeric),
#'       genotype (numeric), phenotype (numeric), width (numeric). Empty data frame if extinct.}
#'     \item{spatial}{List containing world_width (numeric), world_height (numeric),
#'       world_size (numeric, maximum dimension: max(world_width, world_height)), num_patches (integer)}
#'     \item{events}{List with event history. Content depends on history_detail:
#'       \itemize{
#'         \item Always included: times (numeric vector), types (integer vector: -1=initial, 0=death,
#'           1=birth, 2=movement, 3=emigration), individual_ids (integer vector)
#'         \item If history_detail is "standard" or "full": patch_ids, x_coordinates, y_coordinates, genotypes
#'         \item If history_detail == "full": phenotypes, widths
#'       }}
#'     \item{parameters}{Nested list preserving all input parameters (landscape, individual, genetic, simulation)}
#'   }
#'
#' @details
#' Simulation Algorithm:
#'   The simulation uses a Gillespie algorithm (stochastic simulation algorithm) where:
#'   \enumerate{
#'     \item Each individual has rates for three possible events: birth, death, and dispersal
#'     \item Time advances exponentially:
#'       \deqn{\Delta t \sim \text{Exp}(\lambda)}
#'       Where:
#'       \deqn{\lambda = \sum_{i=1}^{N} (d_i + b_i + \mu_i)}
#'       And for each individual i:
#'       \itemize{
#'         \item \eqn{d_i} = dispersal rate
#'         \item \eqn{b_i} = birth rate
#'         \item \eqn{\mu_i} = mortality rate
#'         \item \eqn{N} = current population size
#'       }
#'     \item One event occurs per time step, chosen proportionally to rates
#'     \item Individual rates update after density or location changes
#'   }
#'
#' Stopping Conditions:
#'   The simulation stops when:
#'   \itemize{
#'     \item Population reaches 0 (extinction)
#'     \item Population exceeds 1,000,000 individuals (overflow protection)
#'     \item Simulation time reaches max_events
#'   }
#'   User interruption (Ctrl+C or Esc) is checked periodically during execution.
#'
#' Event Types:
#'   Events are recorded with integer type codes:
#'   \itemize{
#'     \item -1: Initial placement (simulation start)
#'     \item 0: Death
#'     \item 1: Birth
#'     \item 2: Movement (dispersal within landscape)
#'     \item 3: Emigration (exit through absorbing boundary)
#'   }
#'
#' Fitness Function:
#'   Individual fitness determines survival, reproduction, and habitat selection success. Fitness is calculated using a Gaussian function:
#'   
#'   \deqn{W = \exp\left(-\frac{(h - p)^2}{2\sigma_g^2}\right)}
#'   
#'   Where:
#'   \itemize{
#'     \item \eqn{W} = fitness (ranges 0 to 1)
#'     \item \eqn{h} = habitat value at current location
#'     \item \eqn{p} = phenotype (individual's expressed environmental optimum)
#'     \item \eqn{\sigma_g} = genotype_sds (niche width parameter)
#'   }
#'   
#'   Specialist vs Generalist Behavior:
#'   \itemize{
#'     \item \eqn{\sigma_g = 0} (Perfect Specialist): 
#'       \itemize{
#'         \item Fitness = 1 when \eqn{h = p} (exact match)
#'         \item Fitness = 0 otherwise (any mismatch)
#'         \item Enables matrix_dispersal_multiplier and strict habitat-dependent reproduction
#'         \item Strong selection pressure for matching habitat
#'       }
#'     \item \eqn{\sigma_g > 0} (Generalist):
#'       \itemize{
#'         \item Fitness = 1 when \eqn{h = p} (optimal habitat)
#'         \item Fitness decreases gradually as \eqn{|h - p|} increases
#'         \item Larger \eqn{\sigma_g} = broader tolerance (flatter fitness curve)
#'         \item Can survive and reproduce at constant rates across habitat gradients (birth rate not habitat-dependent)
#'       }
#'   }
#'   
#'   Fitness affects multiple processes:
#'   \itemize{
#'     \item Mortality rates (via matrix_mortality_multiplier for generalists)
#'     \item Birth rates (perfect specialists cannot reproduce in non-optimal habitat)
#'     \item Habitat selection (via habitat_selection_temperatures during dispersal)
#'   }
#'
#' Density Calculations:
#'   Population density is calculated differently depending on density_type:
#'   
#'   For local density (density_type = 1):
#'   \deqn{\rho_{local} = \frac{N_{neighbors}}{\pi r^2}}
#'   
#'   For global density (density_type = 2):
#'   \deqn{\rho_{global} = \frac{N_{total}}{\text{width} \times \text{height}}}
#'   
#'   Where:
#'   \itemize{
#'     \item \eqn{\rho} = density
#'     \item \eqn{N_{neighbors}} = number of neighbors within radius
#'     \item \eqn{N_{total}} = total population size
#'     \item \eqn{r} = neighbor_radius
#'     \item width and height are the landscape dimensions (world_width, world_height)
#'   }
#'
#' Birth Rate Calculations:
#'   Birth rates vary with density and habitat suitability:
#'   
#'   For perfect specialists (genotype_sds = 0):
#'   \itemize{
#'     \item In optimal habitat: \eqn{b = b_0 - \beta_b \rho}
#'     \item In non-optimal habitat: \eqn{b = 0} (cannot reproduce)
#'   }
#'   
#'   For generalists (genotype_sds > 0):
#'   \deqn{b = \max(0, b_0 - \beta_b \rho)}
#'   
#'   Where:
#'   \itemize{
#'     \item \eqn{b} = birth_rate
#'     \item \eqn{b_0} = base_birth_rate
#'     \item \eqn{\beta_b} = birth_density_slope
#'     \item \eqn{\rho} = density (from density_type calculation)
#'   }
#'   
#'   Birth rate is constrained to non-negative values. Perfect specialists cannot reproduce in non-optimal habitat, while generalists' reproduction depends only on density, not habitat quality.
#'
#' Mortality Rate Calculations:
#'   The mortality framework implements a tradeoff mechanism where maximum fitness 
#'   is modulated by the specialist-generalist continuum. Individuals experience 
#'   different baseline mortality rates at their environmental optimum depending on 
#'   their niche width, with mortality then interpolating to a universal maximum 
#'   based on habitat mismatch.
#'   
#'   For perfect specialists (genotype_sds = 0):
#'   \itemize{
#'     \item In optimal habitat: \eqn{\mu = \mu_0 + \beta_\mu \rho}
#'     \item In non-optimal habitat: \eqn{\mu = m \mu_0 + \beta_\mu \rho}
#'   }
#'   
#'   For generalists (genotype_sds > 0):
#'   \deqn{\mu = \mu_{realized} + (1 - W)(\mu_{max} - \mu_{realized})}
#'   
#'   Where the realized minimum (baseline mortality at optimum) is:
#'   \deqn{\mu_{realized} = \mu_0 (1 + \text{penalty} \times (m - 1))}
#'   
#'   Where:
#'   \itemize{
#'     \item \eqn{\mu} = mortality rate
#'     \item \eqn{\mu_0} = base_mortality_rate
#'     \item \eqn{m} = matrix_mortality_multiplier
#'     \item \eqn{W} = relative fitness (0 to 1, from Gaussian function)
#'     \item \eqn{\beta_\mu} = mortality_density_slope
#'     \item \eqn{\rho} = density
#'     \item \eqn{\mu_{max} = m \mu_0} (universal maximum mortality)
#'     \item \eqn{\mu_{realized}} = individual-specific baseline at optimum
#'   }
#'   
#'   The generalism penalty \eqn{\text{penalty} = \sigma_g/(\sigma_g + k)}, where 
#'   \eqn{\sigma_g} = genotype_sds (niche width) and \eqn{k} = generalism_cost_scale, 
#'   implements Michaelis-Menten saturation kinetics, ensuring the penalty ranges 
#'   from 0 (specialists, \eqn{\sigma_g = 0}) to approaching 1 (extreme generalists, 
#'   \eqn{\sigma_g >> k}). At \eqn{\sigma_g = k}, the penalty reaches 50\% of its 
#'   maximum. All individuals converge to the same maximum mortality \eqn{\mu_{max}} 
#'   in completely unsuitable habitat (\eqn{W = 0}), but specialists start from 
#'   lower baseline values at their optimum, creating the specialist-generalist 
#'   tradeoff: specialists achieve low mortality in optimal conditions but suffer 
#'   steep penalties when displaced, while generalists tolerate environmental 
#'   variation at the cost of elevated baseline mortality.
#'
#' Matrix Dispersal Multiplier:
#'   The matrix_dispersal_multiplier controls the frequency of dispersal events:
#'   
#'   For perfect specialists (genotype_sds = 0):
#'   \itemize{
#'     \item In optimal habitat: \eqn{d = d_0}
#'     \item In non-optimal habitat: \eqn{d = m_d d_0}
#'   }
#'   
#'   For generalists (genotype_sds > 0), dispersal rate remains \eqn{d_0} regardless of habitat.
#'   
#'   Where:
#'   \itemize{
#'     \item \eqn{d} = dispersal_rate
#'     \item \eqn{d_0} = base_dispersal_rate
#'     \item \eqn{m_d} = matrix_dispersal_multiplier
#'   }
#'
#' Habitat Selection:
#'   When sampling_points = 0, individuals perform random walk:
#'   \deqn{\theta_{new} = \theta_{current} + U\left(-\frac{\alpha}{2}, \frac{\alpha}{2}\right)}
#'   \deqn{\Delta x = L \cos(\theta_{new}), \quad \Delta y = L \sin(\theta_{new})}
#'   
#'   When sampling_points > 0, individuals perform habitat selection:
#'   \enumerate{
#'     \item Sample n = sampling_points candidate locations uniformly within step_length radius:
#'       \deqn{r \sim U(0, L), \quad \theta \sim U(0, 2\pi)}
#'       \deqn{x_{candidate} = x_{current} + r \cos(\theta)}
#'       \deqn{y_{candidate} = y_{current} + r \sin(\theta)}
#'     \item Calculate mortality rate at each candidate location
#'     \item Choose location probabilistically using softmax based on inverse mortality:
#'       \deqn{P(location_i) = \frac{\exp(-\mu_i / T)}{\sum_{j=1}^{n} \exp(-\mu_j / T)}}
#'   }
#'   
#'   Where:
#'   \itemize{
#'     \item \eqn{L} = step_length
#'     \item \eqn{\alpha} = vision_angle
#'     \item \eqn{T} = habitat_selection_temperatures
#'     \item \eqn{\mu_i} = mortality rate at location i (from Mortality Rate Calculations)
#'   }
#'   
#'   Lower mortality locations have higher selection probability. When T is small, 
#'   selection strongly favors low-mortality sites. When T is large, selection becomes 
#'   more exploratory with less discrimination between locations.
#'
#' Genetic Mechanisms:
#'   TWoLife implements two sources of variation:
#'   
#'   Mutation (Heritable Genetic Variation):
#'   \deqn{g_{offspring} = g_{parent} + \varepsilon_{mutation}}
#'   
#'   Where \eqn{\varepsilon_{mutation} \sim N(0, \mu_r)} and \eqn{\mu_r} = mutation_rate
#'   
#'   \itemize{
#'     \item Occurs at birth event
#'     \item Modifies offspring's genotype
#'     \item Passed to future generations (HERITABLE)
#'     \item Enables evolutionary adaptation
#'     \item When \eqn{\mu_r = 0}, offspring genotype = parent genotype exactly
#'   }
#'   
#'   Plasticity (Non-heritable Phenotypic Variation):
#'   \deqn{p = g + \varepsilon_{plasticity}}
#'   
#'   Where \eqn{\varepsilon_{plasticity} \sim N(0, \psi)} and \eqn{\psi} = plasticity parameter
#'   
#'   \itemize{
#'     \item Occurs once at birth, fixed for individual's lifetime
#'     \item Determines phenotype from genotype
#'     \item NOT passed to offspring (offspring get parent's genotype, not phenotype)
#'     \item Represents environmental/developmental noise
#'     \item When \eqn{\psi = 0}, phenotype = genotype exactly
#'   }
#'   
#'   Key Distinction:
#'   \itemize{
#'     \item Genotype: Inherited genetic value (evolves via mutation)
#'     \item Phenotype: Expressed trait value (genotype + plastic noise)
#'     \item Fitness calculation uses phenotype
#'     \item Inheritance passes genotype
#'     \item Both mutation and plasticity add variation, but only mutation is evolutionary
#'   }
#'
#' @examples
#' set.seed(100)
#' landscape <- create_fractal_landscape(
#'   cells_per_row = 5,
#'   fractality = 0.5,
#'   habitat_proportion = 0.6,
#'   return_as_landscape_params = TRUE
#' )
#' result <- twolife_simulation(
#'   landscape_params = landscape,
#'   individual_params = list(
#'     initial_population_size = 15,
#'     base_birth_rate = 0.4,
#'     base_mortality_rate = 0.15
#'   ),
#'   simulation_params = list(max_events = 5),
#'   master_seed = 123
#' )
#' print(result)
#' summary(result)
#' head(result$survivors)
#' landscape_genetic <- create_fractal_landscape(
#'   cells_per_row = 5,
#'   fractality = 0.5,
#'   min_value = 0.35,
#'   max_value = 0.64,
#'   return_as_landscape_params = TRUE
#' )
#' result_genetic <- twolife_simulation(
#'   landscape_params = landscape_genetic,
#'   individual_params = list(
#'     initial_population_size = 15,
#'     base_birth_rate = 0.4,
#'     base_mortality_rate = 0.15
#'   ),
#'   genetic_params = list(
#'     genotype_means = rnorm(15, mean = 0.5, sd = 0.15),
#'     genotype_sds = 0.15
#'   ),
#'   simulation_params = list(max_events = 5),
#'   master_seed = 456
#' )
#' summary(result_genetic)
#' head(result_genetic$survivors)
#' hist(result$survivors$genotype, main = "No Genetic Variation")
#' hist(result_genetic$survivors$genotype, main = "With Genetic Variation")
#' @export
twolife_simulation <- function(landscape_params = list(),
                               individual_params = list(),
                               genetic_params = list(),
                               simulation_params = list(),
                               history_detail = "standard",
                               master_seed = NULL) {
  
  # Validate history_detail parameter
  valid_levels <- c("minimal", "standard", "full")
  if (!history_detail %in% valid_levels) {
    stop("history_detail must be one of: ",
         paste(valid_levels, collapse = ", "),
         call. = FALSE)
  }
  
  if (!is.list(landscape_params)) {
    stop("landscape_params must be a list", call. = FALSE)
  }
  
  if (is.null(landscape_params$habitat)) {
    stop("habitat matrix is required in landscape_params$habitat", call. = FALSE)
  }
  
  habitat_grid <- landscape_params$habitat
  
  if (!is.matrix(habitat_grid)) {
    stop("landscape_params$habitat must be a matrix", call. = FALSE)
  }
  
  if (nrow(habitat_grid) <= 0 || ncol(habitat_grid) <= 0) {
    stop("landscape_params$habitat must have positive dimensions", call. = FALSE)
  }
  
  if (any(!is.finite(habitat_grid))) {
    stop("landscape_params$habitat contains non-finite values", call. = FALSE)
  }
  
  if (!is.null(master_seed)) {
    if (!is.numeric(master_seed) || length(master_seed) != 1 || master_seed != round(master_seed)) {
      stop("master_seed must be a single integer", call. = FALSE)
    }
  }
  
  defaults <- list(
    cell_size = 1.0,
    boundary_condition = 1,
    density_type = 1,
    matrix_mortality_multiplier = 2.0,
    matrix_dispersal_multiplier = 0.5,
    neighbor_radius = 2.0,
    vision_angle = pi,
    step_length = 1.0,
    base_dispersal_rate = 0.1,
    base_birth_rate = 0.3,
    base_mortality_rate = 0.20,
    birth_density_slope = 0.02,
    mortality_density_slope = 0.02,
    generalism_cost_scale = 1.0,
    initial_population_size = 200,
    initial_placement_mode = 1,
    neutral_mode = FALSE
  )
  
  all_params <- modifyList(defaults, c(
    landscape_params[names(landscape_params) != "habitat"],
    individual_params,
    simulation_params
  ))
  
  pop_size <- all_params$initial_population_size
  if (is.null(simulation_params$max_events)) {
    all_params$max_events <- (50 * pop_size)
  } else if (is.null(all_params$max_events)) {
    all_params$max_events <- (50 * pop_size)
  }
  
  if (all_params$base_birth_rate <= all_params$base_mortality_rate) {
    stop("Invalid parameters: base_birth_rate must be > base_mortality_rate", call. = FALSE)
  }
  
  if (all_params$step_length > all_params$neighbor_radius) {
    warning("step_length should be <= neighbor_radius", call. = FALSE)
  }
  
  process_genetic_parameter <- function(param_value, param_name, pop_size, default_value) {
    if (is.null(param_value)) {
      return(rep(default_value, pop_size))
    } else if (length(param_value) == 1) {
      return(rep(param_value, pop_size))
    } else if (length(param_value) == pop_size) {
      return(param_value)
    } else {
      stop("Parameter '", param_name, "' must have length 1 or ", pop_size,
           ", got length ", length(param_value), call. = FALSE)
    }
  }
  
  genotype_means <- process_genetic_parameter(
    genetic_params$genotype_means, "genotype_means", pop_size, 1
  )
  
  genotype_sds <- process_genetic_parameter(
    genetic_params$genotype_sds, "genotype_sds", pop_size, 0
  )
  
  mutation_rates <- process_genetic_parameter(
    genetic_params$mutation_rates, "mutation_rates", pop_size, 0
  )
  
  plasticities <- process_genetic_parameter(
    genetic_params$plasticities, "plasticities", pop_size, 0
  )
  
  sampling_points_raw <- genetic_params$sampling_points
  if (is.null(sampling_points_raw)) {
    sampling_points <- rep(0L, pop_size)
  } else if (length(sampling_points_raw) == 1) {
    sampling_points <- rep(as.integer(sampling_points_raw), pop_size)
  } else if (length(sampling_points_raw) == pop_size) {
    sampling_points <- as.integer(sampling_points_raw)
  } else {
    stop("Parameter 'sampling_points' must have length 1 or ", pop_size,
         ", got length ", length(sampling_points_raw), call. = FALSE)
  }
  
  habitat_selection_temperatures <- process_genetic_parameter(
    genetic_params$habitat_selection_temperatures, "habitat_selection_temperatures", pop_size, 1.0
  )
  
  if (any(habitat_selection_temperatures <= 0)) {
    stop("All habitat_selection_temperatures must be positive", call. = FALSE)
  }
  
  if (all_params$initial_placement_mode == 3) {
    if (is.null(individual_params$initial_x_coordinates) ||
        is.null(individual_params$initial_y_coordinates)) {
      stop("initial_x_coordinates and initial_y_coordinates required for placement mode 3", call. = FALSE)
    }
    initial_x <- individual_params$initial_x_coordinates
    initial_y <- individual_params$initial_y_coordinates
  } else {
    initial_x <- rep(0, pop_size)
    initial_y <- rep(0, pop_size)
  }
  
  result <- run_twolife_simulation(
    neighbor_radius = all_params$neighbor_radius,
    initial_population_size = as.integer(all_params$initial_population_size),
    vision_angle = all_params$vision_angle,
    step_length = all_params$step_length,
    base_dispersal_rate = all_params$base_dispersal_rate,
    base_birth_rate = all_params$base_birth_rate,
    base_mortality_rate = all_params$base_mortality_rate,
    birth_density_slope = all_params$birth_density_slope,
    mortality_density_slope = all_params$mortality_density_slope,
    habitat = habitat_grid,
    cell_size = all_params$cell_size,
    density_type = as.integer(all_params$density_type),
    matrix_mortality_multiplier = all_params$matrix_mortality_multiplier,
    matrix_dispersal_multiplier = all_params$matrix_dispersal_multiplier,
    generalism_cost_scale = all_params$generalism_cost_scale,
    initial_placement_mode = as.integer(all_params$initial_placement_mode),
    boundary_condition = as.integer(all_params$boundary_condition),
    max_events = all_params$max_events,
    initial_x_coordinates = initial_x,
    initial_y_coordinates = initial_y,
    genotype_means = genotype_means,
    genotype_sds = genotype_sds,
    mutation_rates = mutation_rates,
    plasticities = plasticities,
    sampling_points = sampling_points,
    habitat_selection_temperatures = habitat_selection_temperatures,
    neutral_mode = all_params$neutral_mode,
    history_detail = history_detail,
    master_seed = master_seed
  )
  
  final_pop <- as.integer(length(result$survivor_x))
  total_events <- length(result$event_times)
  duration <- if(total_events > 0) max(result$event_times) else 0
  
  # Build events list conditionally based on history_detail
  events_list <- list(
    times = result$event_times,
    types = result$event_types,
    individual_ids = result$individual_ids
  )
  
  if (history_detail != "minimal") {
    events_list$patch_ids <- result$patch_ids
    events_list$x_coordinates <- result$x_coordinates
    events_list$y_coordinates <- result$y_coordinates
    events_list$genotypes <- result$genotypes
  }
  
  if (history_detail == "full") {
    events_list$phenotypes <- result$phenotypes
    events_list$widths <- result$widths
  }
  
  lean_result <- list(
    summary = list(
      final_population_size = final_pop,
      total_events = total_events,
      duration = duration,
      status = if(final_pop > 0) "surviving" else "extinct"
    ),
    survivors = if(final_pop > 0) {
      data.frame(
        id = result$survivor_ids,
        x = result$survivor_x,
        y = result$survivor_y,
        genotype = result$survivor_genotypes,
        phenotype = result$survivor_phenotypes,
        width = result$survivor_widths,
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(id = integer(0), x = numeric(0), y = numeric(0),
                 genotype = numeric(0), phenotype = numeric(0),
                 width = numeric(0), stringsAsFactors = FALSE)
    },
    spatial = list(
      world_width = result$world_width,
      world_height = result$world_height,
      world_size = max(result$world_width, result$world_height),
      num_patches = result$num_patches
    ),
    events = events_list,
    parameters = list(
      landscape = list(
        habitat = habitat_grid,
        cell_size = all_params$cell_size,
        boundary_condition = all_params$boundary_condition,
        density_type = all_params$density_type,
        matrix_mortality_multiplier = all_params$matrix_mortality_multiplier,
        matrix_dispersal_multiplier = all_params$matrix_dispersal_multiplier
      ),
      individual = list(
        initial_population_size = all_params$initial_population_size,
        neighbor_radius = all_params$neighbor_radius,
        vision_angle = all_params$vision_angle,
        step_length = all_params$step_length,
        base_dispersal_rate = all_params$base_dispersal_rate,
        base_birth_rate = all_params$base_birth_rate,
        base_mortality_rate = all_params$base_mortality_rate,
        birth_density_slope = all_params$birth_density_slope,
        mortality_density_slope = all_params$mortality_density_slope,
        initial_placement_mode = all_params$initial_placement_mode,
        initial_x_coordinates = if(all_params$initial_placement_mode == 3) initial_x else NULL,
        initial_y_coordinates = if(all_params$initial_placement_mode == 3) initial_y else NULL
      ),
      genetic = list(
        genotype_means = genotype_means,
        genotype_sds = genotype_sds,
        mutation_rates = mutation_rates,
        plasticities = plasticities,
        sampling_points = sampling_points,
        habitat_selection_temperatures = habitat_selection_temperatures
      ),
      simulation = list(
        max_events = all_params$max_events,
        neutral_mode = all_params$neutral_mode,
        history_detail = history_detail,
        master_seed = master_seed
      )
    )
  )
  
  class(lean_result) <- c("twolife_result", "list")
  return(lean_result)
}

#' Calculate Population Trajectory Over Time
#'
#' Computes population size at each event time from simulation results by 
#' tracking the cumulative effects of demographic events. Essential for 
#' analyzing population dynamics, growth rates, extinction risk, and 
#' identifying critical demographic transitions. Compatible with all 
#' history_detail levels.
#'
#' @param result A 'twolife_result' object returned by \code{\link{twolife_simulation}}.
#'   Must contain \code{$events$times} and \code{$events$types} components.
#'
#' @return A data frame with two columns:
#'   \describe{
#'     \item{time}{Numeric. Event times in chronological order (in simulation time units)}
#'     \item{population_size}{Integer. Population size at each event time, calculated
#'       as cumulative sum of demographic changes}
#'   }
#'
#' @details
#' ## Algorithm
#' 
#' Population size at time \eqn{t} is calculated using:
#' 
#' \deqn{N(t) = N_0 + \sum_{i=1}^{t} \Delta N_i}
#' 
#' Where:
#' \itemize{
#'   \item \eqn{N(t)} = population size at time t
#'   \item \eqn{N_0} = initial population size
#'   \item \eqn{\Delta N_i} = population change from event i
#' }
#' 
#' ## Event Types and Population Change
#' 
#' Each event type contributes differently to \eqn{\Delta N}:
#' 
#' \itemize{
#'   \item **-1 (initial)**: \eqn{\Delta N = +1}. Records each founding individual
#'   \item **0 (death)**: \eqn{\Delta N = -1}. Mortality event
#'   \item **1 (birth)**: \eqn{\Delta N = +1}. Successful reproduction
#'   \item **2 (movement/dispersal)**: \eqn{\Delta N = 0}. Spatial redistribution only
#'   \item **3 (emigration)**: \eqn{\Delta N = -1}. Individual crosses absorbing boundary
#' }
#' 
#' ## Boundary Condition Effects
#' 
#' Emigration events (type 3) only occur with **boundary_condition = 2** (absorbing):
#' \itemize{
#'   \item boundary_condition = 1 (reflective): Individuals bounce back, no emigration
#'   \item boundary_condition = 2 (absorbing): Crossing boundary causes emigration (\eqn{\Delta N = -1})
#'   \item boundary_condition = 3 (periodic): Individuals wrap around, no emigration
#' }
#' 
#' ## Time Scaling
#' 
#' Event times follow the Gillespie algorithm (see \code{\link{twolife_simulation}} Details).
#' Time intervals between events are exponentially distributed with rate proportional 
#' to total demographic rates. Faster dynamics (higher birth/death rates) compress 
#' more events into less time.
#' 
#' ## Analysis Applications
#' 
#' This function enables:
#' \itemize{
#'   \item **Population viability analysis**: Detect extinction times and probabilities
#'   \item **Growth rate estimation**: Calculate \eqn{r = \frac{1}{t}\log\frac{N(t)}{N_0}}
#'   \item **Demographic transitions**: Identify when population stabilizes or crashes
#'   \item **Density-dependent effects**: Compare trajectories with different density parameters
#'   \item **Animation frames**: Provides time points for \code{\link{snapshot_at_time}}
#' }
#' 
#' @examples
#' set.seed(100)
#' landscape <- create_fractal_landscape(
#'   cells_per_row = 5,
#'   fractality = 0.5,
#'   habitat_proportion = 0.6,
#'   return_as_landscape_params = TRUE
#' )
#' landscape_genetic <- create_fractal_landscape(
#'   cells_per_row = 5,
#'   fractality = 0.5,
#'   min_value = 0.35,
#'   max_value = 0.64,
#'   return_as_landscape_params = TRUE
#' )
#' result <- twolife_simulation(
#'   landscape_params = landscape,
#'   individual_params = list(
#'     initial_population_size = 15,
#'     base_birth_rate = 0.4,
#'     base_mortality_rate = 0.15
#'   ),
#'   simulation_params = list(max_events = 5),
#'   master_seed = 123
#' )
#' result_genetic <- twolife_simulation(
#'   landscape_params = landscape_genetic,
#'   individual_params = list(
#'     initial_population_size = 15,
#'     base_birth_rate = 0.4,
#'     base_mortality_rate = 0.15
#'   ),
#'   genetic_params = list(
#'     genotype_means = rnorm(15, mean = 0.5, sd = 0.15),
#'     genotype_sds = 0.15
#'   ),
#'   simulation_params = list(max_events = 5),
#'   master_seed = 456
#' )
#' trajectory <- population_size(result)
#' trajectory_genetic <- population_size(result_genetic)
#' ylim_max <- max(trajectory$population_size, trajectory_genetic$population_size)
#' plot(trajectory, main = "No Genetic Variation", ylim = c(0, ylim_max))
#' plot(trajectory_genetic, main = "With Genetic Variation", ylim = c(0, ylim_max))
#' cat("Final size (no variation):", tail(trajectory$population_size, 1), "\n")
#' cat("Final size (with variation):", tail(trajectory_genetic$population_size, 1), "\n")
#' @seealso \code{\link{twolife_simulation}} for running simulations,
#'   \code{\link{snapshot_at_time}} for reconstructing population state at specific times
#'
#' @export
population_size <- function(result) {
  # Validate input
  if (!inherits(result, "twolife_result")) {
    stop("result must be a twolife_result object from twolife_simulation()", call. = FALSE)
  }
  
  if (is.null(result$events) || is.null(result$events$times) || is.null(result$events$types)) {
    stop("result$events must contain 'times' and 'types' vectors", call. = FALSE)
  }
  
  events_df <- data.frame(
    time = result$events$times,
    event_type = result$events$types,
    stringsAsFactors = FALSE
  )
  
  events_df <- events_df[order(events_df$time), ]
  
  events_df$pop_change <- ifelse(events_df$event_type == -1, 1,
                                 ifelse(events_df$event_type == 0, -1,
                                        ifelse(events_df$event_type == 1, 1,
                                               ifelse(events_df$event_type == 2, 0,
                                                      ifelse(events_df$event_type == 3, -1, 0)))))
  
  events_df$population_size <- cumsum(events_df$pop_change)
  
  return(events_df[, c("time", "population_size")])
}

#' Reconstruct Population State at Specific Time
#'
#' Replays the event log from a simulation to reconstruct the exact spatial distribution,
#' genotypes, and phenotypes of all living individuals at a specified time point. 
#' Essential for temporal analysis, creating animations, studying critical demographic 
#' transitions, and visualizing spatial-genetic patterns through time.
#'
#' @param simulation_result A 'twolife_result' object from \code{\link{twolife_simulation}}.
#'   Must have been created with \code{history_detail = "standard"} or \code{"full"}.
#'   The "minimal" history level lacks spatial coordinates and cannot be used.
#' @param target_time Numeric. Time point to reconstruct population state (in simulation time units).
#'   Must be \eqn{\geq 0} and \eqn{\leq} maximum simulation time. If larger than max time,
#'   automatically uses final time with a warning.
#' @param color_by Character. Trait used to color points in visualization:
#'   \itemize{
#'     \item "genotype" - Color by genetic environmental optimum (\eqn{\mu_g})
#'     \item "phenotype" - Color by expressed phenotype (\eqn{\mu_p = \mu_g + \epsilon}). 
#'       Requires history_detail = "full"
#'     \item "none" - Single color (red) for all individuals
#'   }
#' @param show_plot Logical. If TRUE, displays spatial visualization of reconstructed population
#'   with landscape background.
#' @param point_size Numeric. Size of points representing individuals in the plot.
#'
#' @return Invisibly returns a list with components:
#'   \describe{
#'     \item{time}{Numeric. The actual time used (may differ from input if capped at max time)}
#'     \item{n_alive}{Integer. Number of living individuals at target_time}
#'     \item{population}{Data frame with columns:
#'       \itemize{
#'         \item id: Integer individual ID
#'         \item x, y: Numeric world coordinates
#'         \item genotype: Numeric genetic optimum value
#'         \item phenotype: Numeric phenotype (NA if history_detail != "full")
#'         \item width: Numeric niche width (\eqn{\sigma_g}, NA if history_detail != "full")
#'         \item patch_id: Integer landscape cell ID
#'       }}
#'     \item{events_processed}{Integer. Number of events processed up to target_time}
#'     \item{history_level}{Character. The history_detail level: "standard" or "full"}
#'   }
#'
#' @details
#' ## Reconstruction Algorithm
#' 
#' The function reconstructs population state using:
#' 
#' \deqn{P(t) = \{i : t_{\text{birth},i} \leq t < t_{\text{death},i}\}}
#' 
#' Where:
#' \itemize{
#'   \item \eqn{P(t)} = set of individuals alive at time \eqn{t}
#'   \item \eqn{t_{\text{birth},i}} = birth time of individual \eqn{i}
#'   \item \eqn{t_{\text{death},i}} = death/emigration time of individual \eqn{i} 
#'     (or \eqn{+\infty} if still alive)
#' }
#' 
#' Pseudocode:
#' ```
#' 1. Initialize empty population tracker
#' 2. Filter events: keep only events with time <= target_time
#' 3. For each event in chronological order:
#'    - If birth (type = 1) or initial (type = -1): Add individual to alive set
#'    - If death (type = 0) or emigration (type = 3): Remove from alive set
#'    - If movement (type = 2): Update position of individual
#' 4. For each alive individual, retrieve most recent state (position, genotype, phenotype)
#' 5. Return reconstructed population
#' ```
#' 
#' ## Temporal Resolution
#' 
#' Reconstruction accuracy depends on event density:
#' \itemize{
#'   \item **High temporal resolution**: Many events per time unit -> precise state reconstruction
#'   \item **Low temporal resolution**: Few events per time unit -> states may span long intervals
#'   \item Between events, individual states remain constant (positions frozen until next movement)
#' }
#' 
#' This reflects the event-based nature of individual-based models: changes only occur 
#' at discrete event times.
#' 
#' ## History Detail Requirements
#' 
#' \tabular{lll}{
#'   \strong{history_detail} \tab \strong{Available Data} \tab \strong{color_by Options} \cr
#'   "minimal" \tab Event times, types only \tab Cannot reconstruct \cr
#'   "standard" \tab + positions, genotypes \tab "genotype", "none" \cr
#'   "full" \tab + phenotypes, niche widths \tab All options \cr
#' }
#' 
#' ## Visualization Components
#' 
#' When show_plot = TRUE:
#' \itemize{
#'   \item **Background**: Habitat values from landscape matrix
#'   \item **Points**: Individual positions colored by selected trait
#'   \item **Title**: Shows target_time and number of individuals alive
#'   \item **Legend**: Maps colors to trait values (if color_by != "none")
#' }
#' 
#' ## Animation Workflow
#' 
#' To create animations:
#' 1. Get time points: \code{times <- seq(0, max_time, length.out = 50)}
#' 2. For each time: \code{snapshot_at_time(result, t)}
#' 3. Save plots or compile frames into video
#' 
#' See Examples for complete animation code.
#' 
#' @examples
#' set.seed(100)
#' landscape <- create_fractal_landscape(
#'   cells_per_row = 5,
#'   fractality = 0.5,
#'   habitat_proportion = 0.6,
#'   return_as_landscape_params = TRUE
#' )
#' landscape_genetic <- create_fractal_landscape(
#'   cells_per_row = 5,
#'   fractality = 0.5,
#'   min_value = 0.35,
#'   max_value = 0.64,
#'   return_as_landscape_params = TRUE
#' )
#' result <- twolife_simulation(
#'   landscape_params = landscape,
#'   individual_params = list(
#'     initial_population_size = 15,
#'     base_birth_rate = 0.4,
#'     base_mortality_rate = 0.15
#'   ),
#'   simulation_params = list(max_events = 5),
#'   master_seed = 123
#' )
#' result_genetic <- twolife_simulation(
#'   landscape_params = landscape_genetic,
#'   individual_params = list(
#'     initial_population_size = 15,
#'     base_birth_rate = 0.4,
#'     base_mortality_rate = 0.15
#'   ),
#'   genetic_params = list(
#'     genotype_means = rnorm(15, mean = 0.5, sd = 0.15),
#'     genotype_sds = 0.15
#'   ),
#'   simulation_params = list(max_events = 5),
#'   master_seed = 456
#' )
#' state_5 <- snapshot_at_time(result, 5, show_plot = TRUE)
#' state_genetic_5 <- snapshot_at_time(result_genetic, 5, show_plot = TRUE)
#' cat("At time 5 (no variation):", state_5$n_alive, "individuals\n")
#' cat("At time 5 (with variation):", state_genetic_5$n_alive, "individuals\n")
#' if (state_genetic_5$n_alive > 0) {
#'   cat("Mean genotype:", mean(state_genetic_5$population$genotype), "\n")
#'   cat("SD genotype:", sd(state_genetic_5$population$genotype), "\n")
#' }
#' @seealso \code{\link{population_size}} for population trajectories,
#'   \code{\link{twolife_simulation}} for running simulations with appropriate history_detail,
#'   \code{\link{check_habitat_match}} for analyzing trait-habitat matching
#'
#' @export
snapshot_at_time <- function(simulation_result,
                             target_time,
                             color_by = "genotype",
                             show_plot = TRUE,
                             point_size = 2) {
  
  if (!inherits(simulation_result, "twolife_result")) {
    stop("simulation_result must be a twolife_result object", call. = FALSE)
  }
  
  # Check history detail level
  history_level <- simulation_result$parameters$simulation$history_detail
  
  if (is.null(history_level)) {
    # Backward compatibility - infer from available fields
    if ("phenotypes" %in% names(simulation_result$events)) {
      history_level <- "full"
    } else if ("x_coordinates" %in% names(simulation_result$events)) {
      history_level <- "standard"
    } else {
      history_level <- "minimal"
    }
  }
  
  # Validate capability
  if (history_level == "minimal") {
    stop("Cannot reconstruct spatial positions with history_detail='minimal'.\n",
         "Re-run simulation with history_detail='standard' or 'full'.",
         call. = FALSE)
  }
  
  # Validate color_by parameter
  if (!color_by %in% c("genotype", "phenotype", "none")) {
    stop("color_by must be one of: 'genotype', 'phenotype', or 'none'", call. = FALSE)
  }
  
  # Check if phenotype coloring is requested but not available
  if (color_by == "phenotype" && history_level != "full") {
    stop("color_by='phenotype' requires history_detail='full'.\n",
         "Either re-run simulation with history_detail='full' or use color_by='genotype'.",
         call. = FALSE)
  }
  
  # Extract available fields based on history level
  events <- data.frame(
    time = simulation_result$events$times,
    type = simulation_result$events$types,
    id = simulation_result$events$individual_ids,
    stringsAsFactors = FALSE
  )
  
  if (history_level != "minimal") {
    events$patch_id <- simulation_result$events$patch_ids
    events$x <- simulation_result$events$x_coordinates
    events$y <- simulation_result$events$y_coordinates
    events$genotype <- simulation_result$events$genotypes
  }
  
  use_exact_phenotypes <- FALSE
  if (history_level == "full") {
    events$phenotype <- simulation_result$events$phenotypes
    events$width <- simulation_result$events$widths
    use_exact_phenotypes <- TRUE
  }
  
  # Validate target_time
  if (target_time < 0) {
    stop("target_time must be non-negative", call. = FALSE)
  }
  
  if (target_time > max(events$time)) {
    warning("target_time exceeds simulation duration. Using final time.", call. = FALSE)
    target_time <- max(events$time)
  }
  
  # Filter events up to target time
  events_up_to <- events[events$time <= target_time, ]
  
  if (nrow(events_up_to) == 0) {
    stop("No events occurred before target_time", call. = FALSE)
  }
  
  # Track population state by replaying events
  alive_ids <- c()
  alive_data <- list()
  
  for (i in seq_len(nrow(events_up_to))) {
    event <- events_up_to[i, ]
    
    if (event$type == -1) {
      # Initial state
      alive_ids <- c(alive_ids, event$id)
      alive_data[[as.character(event$id)]] <- event
      
    } else if (event$type == 1) {
      # Birth
      alive_ids <- c(alive_ids, event$id)
      alive_data[[as.character(event$id)]] <- event
      
    } else if (event$type == 0 || event$type == 3) {
      # Death or emigration
      alive_ids <- setdiff(alive_ids, event$id)
      alive_data[[as.character(event$id)]] <- NULL
      
    } else if (event$type == 2) {
      # Movement - update position
      if (as.character(event$id) %in% names(alive_data)) {
        alive_data[[as.character(event$id)]]$x <- event$x
        alive_data[[as.character(event$id)]]$y <- event$y
        alive_data[[as.character(event$id)]]$patch_id <- event$patch_id
      }
    }
  }
  
  # Convert to data frame - return EXACTLY what's in the event history (no approximations)
  if (length(alive_ids) == 0) {
    alive_pop <- data.frame(
      id = integer(0),
      x = numeric(0),
      y = numeric(0),
      genotype = numeric(0),
      phenotype = numeric(0),
      width = numeric(0),
      patch_id = integer(0)
    )
  } else {
    alive_pop <- do.call(rbind, lapply(alive_data, function(row) {
      data.frame(
        id = row$id,
        x = row$x,
        y = row$y,
        genotype = row$genotype,
        phenotype = if(use_exact_phenotypes) row$phenotype else NA_real_,
        width = if(use_exact_phenotypes) row$width else NA_real_,
        patch_id = row$patch_id,
        stringsAsFactors = FALSE
      )
    }))
    rownames(alive_pop) <- NULL
  }
  
  # Visualization
  if (show_plot && nrow(alive_pop) > 0) {
    habitat_grid <- simulation_result$parameters$landscape$habitat
    cell_size <- simulation_result$parameters$landscape$cell_size
    
    n_rows <- nrow(habitat_grid)
    n_cols <- ncol(habitat_grid)
    world_width <- n_cols * cell_size
    world_height <- n_rows * cell_size
    
    x_coords <- seq(-world_width/2, world_width/2, length.out = n_cols)
    y_coords <- seq(-world_height/2, world_height/2, length.out = n_rows)
    z <- t(apply(habitat_grid, 2, rev))
    
    # Determine if binary and set colors accordingly
    is_binary <- all(habitat_grid %in% c(0, 1))
    z_range <- range(habitat_grid)
    
    if (is_binary) {
      # Binary: 0 = white (matrix), 1 = green (habitat)
      landscape_cols <- c("white", "#228B22")
    } else {
      # Continuous: terrain colors
      landscape_cols <- terrain.colors(100)
    }
    
    par(mfrow = c(1, 2))
    
    # Left panel: Landscape
    image(x = x_coords, y = y_coords, z = z,
          col = landscape_cols,
          main = "Landscape",
          xlab = "X", ylab = "Y", asp = 1)
    grid(col = "gray80", lty = 1, lwd = 0.5)
    
    # Right panel: Population with color mapping
    trait_values <- if (color_by == "phenotype") alive_pop$phenotype else alive_pop$genotype
    
    if (color_by != "none" && z_range[1] != z_range[2]) {
      # Normalize trait values to [0, 1] based on landscape range
      trait_normalized <- (trait_values - z_range[1]) / (z_range[2] - z_range[1])
      trait_normalized <- pmax(0, pmin(1, trait_normalized))
      
      # Map to the SAME color palette as the landscape
      if (is_binary) {
        # For binary: traits < 0.5 -> white, >= 0.5 -> green
        color_indices <- ifelse(trait_normalized < 0.5, 1, 2)
        point_colors <- landscape_cols[color_indices]
      } else {
        # For continuous: map to terrain.colors consistently
        color_indices <- pmax(1, pmin(100, round(trait_normalized * 99) + 1))
        point_colors <- landscape_cols[color_indices]
      }
    } else {
      point_colors <- rep("red", nrow(alive_pop))
    }
    
    trait_name <- if (color_by == "phenotype") "Phenotype" else if (color_by == "genotype") "Genotype" else "Trait"
    
    plot(alive_pop$x, alive_pop$y,
         col = point_colors, pch = 16, cex = point_size,
         main = paste("Population at t =", round(target_time, 2),
                      if (color_by != "none") paste("\n(colored by", trait_name, ")") else ""),
         xlab = "X", ylab = "Y", asp = 1,
         xlim = c(-world_width/2, world_width/2),
         ylim = c(-world_height/2, world_height/2))
    grid(col = "gray80", lty = 1, lwd = 0.5)
    
    par(mfrow = c(1, 1))
  }
  
  return(invisible(list(
    time = target_time,
    n_alive = nrow(alive_pop),
    population = alive_pop,
    events_processed = nrow(events_up_to),
    history_level = history_level
  )))
}

# ============================================================================
# LANDSCAPE GENERATION FUNCTIONS
# ============================================================================

#' Create Fractal Landscape with Spatial Autocorrelation
#'
#' Generates spatially autocorrelated landscapes using iterative neighbor averaging,
#' a simplified fractal-like algorithm. Creates realistic habitat patterns with
#' controllable fragmentation and clustering. Supports both square and rectangular
#' landscapes, and can generate either continuous environmental gradients or
#' binary habitat/matrix patterns.
#'
#' @param cells_per_row Integer. Number of cells per row (matrix rows). Must be positive.
#'   Typical values: 10-100 for most simulations. Larger values increase resolution
#'   but also memory requirements and computation time.
#' @param cells_per_col Integer. Number of cells per column (matrix columns). If NULL,
#'   creates square landscape (cells_per_col = cells_per_row). Use different values
#'   to create rectangular landscapes matching specific study regions.
#' @param fractality Numeric between 0 and 1. Controls degree of spatial autocorrelation
#'   and clumping:
#'   \itemize{
#'     \item **0.0-0.2**: Nearly random, little spatial structure. Isolated habitat cells
#'       scattered throughout matrix.
#'     \item **0.3-0.5**: Moderate clumping. Small habitat patches with some aggregation.
#'       Good for testing intermediate fragmentation scenarios.
#'     \item **0.6-0.7**: Strong clumping. Larger, more contiguous patches with clear
#'       spatial structure. Typical of many real landscapes.
#'     \item **0.8-1.0**: Very strong clumping. Large, highly aggregated patches or gradients.
#'       Useful for testing extreme fragmentation or strong environmental gradients.
#'   }
#'   
#'   Higher fractality -> more smoothing iterations -> stronger spatial autocorrelation.
#'   
#' @param min_value Numeric. Minimum habitat value for continuous landscapes.
#'   Only used when habitat_proportion = NULL. Typically 0.0.
#' @param max_value Numeric. Maximum habitat value for continuous landscapes.
#'   Only used when habitat_proportion = NULL. Must be > min_value. Typically 1.0.
#' @param habitat_proportion Numeric between 0 and 1, or NULL. Controls landscape type:
#'   \itemize{
#'     \item **If provided**: Creates **binary** landscape (0/1, matrix/habitat).
#'       Value specifies proportion of cells designated as habitat (value = 1).
#'       Remaining cells become matrix (value = 0). Fractality determines spatial
#'       clustering of habitat cells.
#'     \item **If NULL**: Creates **continuous** landscape with values scaled to
#'       [min_value, max_value]. Represents smooth environmental gradients.
#'   }
#' @param return_as_landscape_params Logical. Controls return format:
#'   \itemize{
#'     \item **TRUE**: Returns \code{list(habitat = matrix)}, ready for direct use
#'       as landscape_params in \code{\link{twolife_simulation}}
#'     \item **FALSE**: Returns matrix only, for custom processing or visualization
#'   }
#'
#' @return If \code{return_as_landscape_params = FALSE}:
#'   \itemize{
#'     \item Numeric matrix with dimensions \code{cells_per_row} x \code{cells_per_col}
#'     \item Values are either continuous ([min_value, max_value]) or binary (0, 1)
#'   }
#'
#'   If \code{return_as_landscape_params = TRUE}:
#'   \itemize{
#'     \item List with component: \code{$habitat} (the generated landscape matrix)
#'     \item Compatible with \code{twolife_simulation(landscape_params = ...)}
#'   }
#'
#' @details
#' ## Fractal Generation Algorithm
#' 
#' The function uses iterative neighbor averaging to create spatial autocorrelation:
#' 
#' **Step 1: Initialization**
#' 
#' Each cell \eqn{(i,j)} starts with random value:
#' \deqn{v_0(i,j) \sim \text{Uniform}(0, 1)}
#' 
#' **Step 2: Iterative Smoothing**
#' 
#' The algorithm performs \eqn{n_{\text{iter}}} smoothing iterations:
#' 
#' \deqn{n_{\text{iter}} = \text{round}(10 \times \text{fractality})}
#' 
#' At each iteration \eqn{t}, cell values are updated:
#' 
#' \deqn{v_{t+1}(i,j) = \alpha \cdot \bar{v}_t(i,j) + (1-\alpha) \cdot \varepsilon_{i,j}}
#' 
#' Where:
#' \itemize{
#'   \item \eqn{v_t(i,j)} = value at cell \eqn{(i,j)} at iteration \eqn{t}
#'   \item \eqn{\bar{v}_t(i,j)} = mean of 8-connected neighbors (Moore neighborhood):
#'     \deqn{\bar{v}_t(i,j) = \frac{1}{N_{\text{neighbors}}}\sum_{(i',j') \in N(i,j)} v_t(i',j')}
#'   \item \eqn{\alpha = 0.9} = smoothing weight (high value preserves spatial structure)
#'   \item \eqn{\varepsilon_{i,j} \sim \text{Uniform}(0,1)} = random noise (maintains variability)
#'   \item Edge cells use only available neighbors (\eqn{N_{\text{neighbors}} < 8})
#' }
#' 
#' **Step 3: Finalization**
#' 
#' After all iterations:
#' 
#' For **continuous landscapes** (habitat_proportion = NULL):
#' \deqn{v_{\text{final}}(i,j) = \text{min\_value} + v_n(i,j) \times (\text{max\_value} - \text{min\_value})}
#' 
#' For **binary landscapes** (habitat_proportion specified):
#' \itemize{
#'   \item Sort all cell values
#'   \item Set top \code{habitat_proportion} quantile to 1 (habitat)
#'   \item Set remaining cells to 0 (matrix)
#'   \item Result: spatially clustered binary pattern
#' }
#' 
#' ## Spatial Autocorrelation
#' 
#' The iterative averaging creates spatial autocorrelation described by:
#' 
#' \deqn{\text{Cor}(v(i,j), v(i',j')) \approx \exp\left(-\frac{d}{\lambda}\right)}
#' 
#' Where:
#' \itemize{
#'   \item \eqn{d = \sqrt{(i-i')^2 + (j-j')^2}} = distance between cells
#'   \item \eqn{\lambda \propto \text{fractality}} = correlation length scale
#'   \item Higher fractality -> larger \eqn{\lambda} -> longer-range correlations -> bigger patches
#' }
#' 
#' ## Comparison to Other Methods
#' 
#' \tabular{llll}{
#'   \strong{Method} \tab \strong{Spatial Structure} \tab \strong{Computation} \tab \strong{Parameters} \cr
#'   This function \tab Moderate \tab Fast \tab Simple (fractality) \cr
#'   Midpoint displacement \tab Strong fractal \tab Medium \tab Hurst exponent \cr
#'   Random fields (NLMR) \tab Tunable \tab Slow \tab Multiple parameters \cr
#'   Percolation \tab Binary only \tab Fast \tab Threshold only \cr
#' }
#' 
#' Advantages of this method:
#' \itemize{
#'   \item Fast computation (suitable for parameter sweeps)
#'   \item Intuitive fractality parameter (0-1 scale)
#'   \item Supports both continuous and binary landscapes
#'   \item No external dependencies
#' }
#' 
#' ## Parameter Selection Guidelines
#' 
#' **For habitat fragmentation studies:**
#' \itemize{
#'   \item Low fragmentation: fractality = 0.7-0.9, habitat_proportion = 0.5-0.7
#'   \item Medium fragmentation: fractality = 0.4-0.6, habitat_proportion = 0.3-0.5
#'   \item High fragmentation: fractality = 0.1-0.3, habitat_proportion = 0.2-0.4
#' }
#' 
#' **For environmental gradient studies:**
#' \itemize{
#'   \item Smooth gradients: fractality = 0.8-1.0, habitat_proportion = NULL
#'   \item Patchy gradients: fractality = 0.4-0.6, habitat_proportion = NULL
#' }
#' 
#' **For realistic landscapes:**
#' \itemize{
#'   \item Forest fragments: fractality = 0.5-0.7, habitat_proportion = 0.3-0.5
#'   \item Agricultural landscapes: fractality = 0.4-0.6, habitat_proportion = 0.4-0.6
#'   \item Urban-rural gradients: fractality = 0.6-0.8, continuous landscape
#' }
#' 
#' ## Edge Effects
#' 
#' Cells at landscape edges have fewer neighbors, which can create slight edge artifacts:
#' \itemize{
#'   \item Effect is minor for most fractality values
#'   \item Becomes more noticeable at very high fractality (> 0.9)
#'   \item Use larger landscapes if edge effects are a concern
#'   \item Consider periodic boundary conditions in simulations if modeling infinite landscapes
#' }
#'
#' @examples
#' set.seed(100)
#' landscape1 <- create_fractal_landscape(
#'   cells_per_row = 5,
#'   fractality = 0.5,
#'   habitat_proportion = 0.6,
#'   return_as_landscape_params = TRUE
#' )
#' plot_landscape(landscape1, main = "Binary Habitat Landscape")
#' landscape2 <- create_fractal_landscape(
#'   cells_per_row = 5,
#'   fractality = 0.5,
#'   min_value = 0.35,
#'   max_value = 0.64,
#'   return_as_landscape_params = TRUE
#' )
#' plot_landscape(landscape2, main = "Continuous Gradient for Genetics")
#' print(landscape1)  # Binary habitat values
#' print(landscape2)  # Continuous range
#' @seealso \code{\link{plot_landscape}} for visualization,
#'   \code{\link{twolife_simulation}} for using landscapes in simulations,
#'   \code{\link{plot_simulation_on_landscape}} for overlaying simulation results
#'
#' @export
create_fractal_landscape <- function(cells_per_row,
                                     cells_per_col = NULL,
                                     fractality,
                                     min_value = 0.0,
                                     max_value = 1.0,
                                     habitat_proportion = NULL,
                                     return_as_landscape_params = FALSE) {
  
  if (!is.numeric(cells_per_row) || cells_per_row <= 0 || cells_per_row != round(cells_per_row)) {
    stop("cells_per_row must be a positive integer", call. = FALSE)
  }
  
  if (is.null(cells_per_col)) {
    cells_per_col <- cells_per_row
  }
  
  if (!is.numeric(cells_per_col) || cells_per_col <= 0 || cells_per_col != round(cells_per_col)) {
    stop("cells_per_col must be a positive integer", call. = FALSE)
  }
  
  if (!is.numeric(fractality) || fractality < 0 || fractality > 1) {
    stop("fractality must be between 0 and 1", call. = FALSE)
  }
  
  if (!is.numeric(min_value) || !is.numeric(max_value)) {
    stop("min_value and max_value must be numeric", call. = FALSE)
  }
  
  if (min_value >= max_value) {
    stop("min_value must be less than max_value", call. = FALSE)
  }
  
  if (!is.null(habitat_proportion)) {
    if (!is.numeric(habitat_proportion) || habitat_proportion < 0 || habitat_proportion > 1) {
      stop("habitat_proportion must be between 0 and 1", call. = FALSE)
    }
  }
  
  range01 <- function(x) {
    (x - min(x)) / (max(x) - min(x))
  }
  
  # Generate fractal pattern using iterative smoothing
  fractal_matrix <- matrix(runif(cells_per_row * cells_per_col), nrow = cells_per_row, ncol = cells_per_col)
  n_iterations <- max(1, round(fractality * 10))
  
  for (iter in 1:n_iterations) {
    smoothed <- matrix(0, nrow = cells_per_row, ncol = cells_per_col)
    
    for (i in 1:cells_per_row) {
      for (j in 1:cells_per_col) {
        neighbors <- c()
        if (i > 1) neighbors <- c(neighbors, fractal_matrix[i-1, j])
        if (i < cells_per_row) neighbors <- c(neighbors, fractal_matrix[i+1, j])
        if (j > 1) neighbors <- c(neighbors, fractal_matrix[i, j-1])
        if (j < cells_per_col) neighbors <- c(neighbors, fractal_matrix[i, j+1])
        if (i > 1 && j > 1) neighbors <- c(neighbors, fractal_matrix[i-1, j-1])
        if (i > 1 && j < cells_per_col) neighbors <- c(neighbors, fractal_matrix[i-1, j+1])
        if (i < cells_per_row && j > 1) neighbors <- c(neighbors, fractal_matrix[i+1, j-1])
        if (i < cells_per_row && j < cells_per_col) neighbors <- c(neighbors, fractal_matrix[i+1, j+1])
        
        weight_original <- 1 - fractality
        weight_neighbors <- fractality / length(neighbors)
        
        smoothed[i, j] <- weight_original * fractal_matrix[i, j] + sum(neighbors * weight_neighbors)
      }
    }
    fractal_matrix <- smoothed
  }
  
  normalized_fractal <- range01(fractal_matrix)
  
  if (!is.null(habitat_proportion)) {
    if (habitat_proportion >= 1.0) {
      binary_landscape <- matrix(1, nrow = cells_per_row, ncol = cells_per_col)
    } else if (habitat_proportion <= 0.0) {
      binary_landscape <- matrix(0, nrow = cells_per_row, ncol = cells_per_col)
    } else {
      threshold <- quantile(normalized_fractal, 1 - habitat_proportion, na.rm = TRUE)
      binary_landscape <- matrix(0, nrow = cells_per_row, ncol = cells_per_col)
      binary_landscape[normalized_fractal >= threshold] <- 1
    }
    
    if (return_as_landscape_params) {
      return(list(habitat = binary_landscape))
    } else {
      return(binary_landscape)
    }
    
  } else {
    scaled_fractal <- min_value + (normalized_fractal * (max_value - min_value))
    
    habitat_grid <- matrix(nrow = cells_per_row, ncol = cells_per_col)
    count <- 0
    for (i in 1:cells_per_row) {
      for (j in 1:cells_per_col) {
        count <- count + 1
        habitat_grid[i, j] <- scaled_fractal[count]
      }
    }
    
    if (return_as_landscape_params) {
      return(list(habitat = habitat_grid))
    } else {
      return(habitat_grid)
    }
  }
}


# ============================================================================
# VISUALIZATION FUNCTIONS
# ============================================================================

#' Plot Landscape with World Coordinates
#'
#' Visualizes landscape habitat values in the world coordinate system used by 
#' TWoLife simulations. Displays proper axis scaling, centered origin, and 
#' optional grid overlay. Essential for understanding spatial context of 
#' simulation results and verifying landscape structure.
#'
#' @param landscape_data Either:
#'   \itemize{
#'     \item A numeric matrix of habitat values (dimensions: rows x columns)
#'     \item A list with a \code{$habitat} component (output from \code{\link{create_fractal_landscape}}
#'       or landscape creation functions)
#'   }
#' @param cell_size Numeric. Size of each landscape cell in world units. Should match
#'   the cell_size parameter used in \code{\link{twolife_simulation}}. Together with
#'   matrix dimensions, defines world_width and world_height.
#' @param filename Character. Optional file path to save plot (e.g., "landscape.png").
#'   If NULL, displays interactively in the active graphics device.
#' @param main Character. Plot title displayed at the top of the figure.
#' @param colors Character. Color scheme for habitat visualization:
#'   \itemize{
#'     \item "habitat" - Custom green gradient (white = 0 = matrix/unsuitable, 
#'       dark green = 1 = optimal habitat). Emphasizes habitat vs. matrix distinction.
#'     \item "terrain" - R's built-in terrain.colors() palette (green-brown-white gradient).
#'       General-purpose landscape visualization.
#'     \item "viridis" - Perceptually uniform viridis color scale (requires viridisLite package).
#'       Best for continuous habitat value gradients and colorblind-friendly displays.
#'   }
#' @param show_legend Logical. If TRUE, displays color legend showing mapping between
#'   habitat values and colors.
#' @param add_grid Logical. If TRUE, adds:
#'   \itemize{
#'     \item Gray grid lines at cell boundaries
#'     \item Red crosshairs marking origin (0, 0)
#'   }
#'
#' @return NULL, invisibly. Function is called for its side effect (creating a plot).
#'
#' @details
#' ## World Coordinate System
#' 
#' TWoLife uses a centered coordinate system where:
#' 
#' \deqn{x \in \left[-\frac{W}{2}, +\frac{W}{2}\right], \quad y \in \left[-\frac{H}{2}, +\frac{H}{2}\right]}
#' 
#' Where:
#' \itemize{
#'   \item \eqn{W = \text{ncol}(\text{habitat}) \times \text{cell\_size}} = world_width
#'   \item \eqn{H = \text{nrow}(\text{habitat}) \times \text{cell\_size}} = world_height
#'   \item Origin (0, 0) is at the landscape center
#'   \item Individual positions (x, y) are continuous within these bounds
#' }
#' 
#' ## Coordinate Mapping
#' 
#' Matrix indices map to world coordinates via:
#' 
#' \deqn{x_{\text{center}} = \left(\text{col} - \frac{\text{ncol}}{2}\right) \times \text{cell\_size}}
#' \deqn{y_{\text{center}} = \left(\text{row} - \frac{\text{nrow}}{2}\right) \times \text{cell\_size}}
#' 
#' Important notes:
#' \itemize{
#'   \item Matrix row 1 (top) maps to positive y values
#'   \item Matrix column 1 (left) maps to negative x values
#'   \item Each cell spans +/-cell_size/2 around its center
#' }
#' 
#' ## Color Palette Details
#' 
#' \tabular{llll}{
#'   \strong{Palette} \tab \strong{Range} \tab \strong{Best For} \tab \strong{Properties} \cr
#'   "habitat" \tab White -> Dark Green \tab Binary landscapes \tab Clear matrix/habitat contrast \cr
#'   "terrain" \tab Green -> Brown -> White \tab Natural landscapes \tab Traditional terrain visualization \cr
#'   "viridis" \tab Purple -> Green -> Yellow \tab Continuous gradients \tab Perceptually uniform, colorblind-safe \cr
#' }
#' 
#' ## Rectangular Landscapes
#' 
#' The function automatically handles non-square landscapes:
#' \itemize{
#'   \item Maintains correct aspect ratio
#'   \item Adjusts axis limits independently for x and y
#'   \item Preserves cell_size scaling in both dimensions
#' }
#' 
#' ## Grid Overlay
#' 
#' When add_grid = TRUE:
#' \itemize{
#'   \item **Gray lines**: Mark boundaries between landscape cells (aid in counting cells,
#'     measuring distances)
#'   \item **Red crosshairs**: Show (0, 0) origin (reference point for individual positions)
#' }
#' 
#' ## Visualization Best Practices
#' 
#' \itemize{
#'   \item **Match simulation parameters**: Always use the same cell_size as in your simulation
#'   \item **Verify scale**: Check that world dimensions match simulation output spatial range
#'   \item **Color choice**: Use "habitat" for binary, "viridis" for continuous gradients
#'   \item **Grid for precision**: Enable grid when analyzing specific spatial positions
#' }
#'
#' @examples
#' set.seed(100)
#' landscape <- create_fractal_landscape(
#'   cells_per_row = 5,
#'   fractality = 0.5,
#'   habitat_proportion = 0.6,
#'   return_as_landscape_params = TRUE
#' )
#' plot_landscape(landscape, main = "Binary Habitat Landscape")
#' landscape_genetic <- create_fractal_landscape(
#'   cells_per_row = 5,
#'   fractality = 0.5,
#'   min_value = 0.35,
#'   max_value = 0.64,
#'   return_as_landscape_params = TRUE
#' )
#' plot_landscape(
#'   landscape_genetic,
#'   main = "Continuous Gradient",
#'   colors = "viridis"
#' )
#' @seealso \code{\link{plot_simulation_on_landscape}} to overlay simulation results,
#'   \code{\link{create_fractal_landscape}} for landscape generation,
#'   \code{\link{twolife_simulation}} for running simulations with landscapes
#'
#' @export
plot_landscape <- function(landscape_data,
                           cell_size = 1.0,
                           filename = NULL,
                           main = "Landscape (World Coordinates)",
                           colors = "habitat",
                           show_legend = TRUE,
                           add_grid = TRUE) {
  
  if (is.list(landscape_data) && "habitat" %in% names(landscape_data)) {
    habitat_grid <- landscape_data$habitat
  } else if (is.matrix(landscape_data)) {
    habitat_grid <- landscape_data
  } else {
    stop("landscape_data must be a matrix or list with $habitat component", call. = FALSE)
  }
  
  if (!is.matrix(habitat_grid)) {
    stop("habitat_grid must be a matrix", call. = FALSE)
  }
  
  n_rows <- nrow(habitat_grid)
  n_cols <- ncol(habitat_grid)
  world_width <- n_cols * cell_size
  world_height <- n_rows * cell_size
  x_coords <- seq(-world_width/2, world_width/2, length.out = n_cols)
  y_coords <- seq(-world_height/2, world_height/2, length.out = n_rows)
  
  z <- t(apply(habitat_grid, 2, rev))
  
  is_binary <- all(habitat_grid %in% c(0, 1))
  z_range <- range(z, na.rm = TRUE)
  is_uniform <- z_range[1] == z_range[2]
  
  if (colors == "habitat") {
    if (is_binary) {
      # Binary landscapes: 0 = white (matrix), 1 = green (habitat)
      color_palette <- c("white", "#228B22")
      if (is_uniform) {
        color_palette <- if(z_range[1] == 0) "white" else "#228B22"
      }
    } else {
      color_palette <- colorRampPalette(c("#F5F5DC", "#90EE90", "#228B22", "#006400"))(100)
    }
  } else {
    n_colors <- if (is_binary && !is_uniform) 2 else if (is_uniform) 1 else 100
    color_palette <- switch(colors,
                            "terrain" = terrain.colors(n_colors),
                            "viridis" = {
                              if (requireNamespace("viridisLite", quietly = TRUE)) {
                                viridisLite::viridis(n_colors)
                              } else {
                                terrain.colors(n_colors)
                              }
                            },
                            terrain.colors(n_colors)
    )
  }
  
  image(x = x_coords, y = y_coords, z = z,
        col = color_palette,
        main = main,
        xlab = "X (World Coordinates)",
        ylab = "Y (World Coordinates)",
        asp = 1,
        useRaster = TRUE)
  
  if (add_grid) {
    grid(col = "gray80", lty = 1, lwd = 0.5)
    abline(h = 0, v = 0, col = "red", lty = 2, lwd = 2)
  }
  
  invisible(NULL)
}

#' Visualize Simulation Results Overlaid on Landscape
#'
#' Creates a comprehensive visualization showing the spatial distribution of
#' surviving individuals on the habitat landscape. Displays habitat values
#' as background colors with survivors as colored points, enabling assessment
#' of habitat selection, spatial clustering, and genotype-environment matching.
#'
#' @param simulation_result A 'twolife_result' object from \code{\link{twolife_simulation}}.
#' @param point_size Numeric. Size of points representing surviving individuals in the plot.
#' @param point_color Character. Color for points when color_by="none". Any valid R color name or hex code.
#' @param point_shape Integer. Point shape following R's pch codes (16 = filled circle, 1 = open circle, 17 = filled triangle).
#' @param color_by Character. Trait used to color survivor points:
#'   \itemize{
#'     \item "none" - Single color specified by point_color
#'     \item "genotype" - Color gradient based on genotype values (heat colors)
#'     \item "phenotype" - Color gradient based on phenotype values (heat colors)
#'   }
#' @param landscape_colors Character. Color scheme for landscape background:
#'   \itemize{
#'     \item "habitat" - Green gradient (white = matrix, dark green = high-quality habitat). Best for binary landscapes.
#'     \item "terrain" - Built-in terrain.colors palette (green-brown-white). General-purpose visualization.
#'     \item "viridis" - Viridis color scale (requires viridisLite package). Perceptually uniform and colorblind-friendly.
#'   }
#' @param main Character. Plot title. If NULL, automatically generates title showing final time and population size.
#' @param filename Character. Optional file path to save plot. If NULL, displays interactively.
#' @param add_stats Logical. If TRUE, adds statistics text to plot (currently not implemented).
#'
#' @return The input simulation_result object, invisibly. Allows for piping/chaining operations.
#'
#' @details
#' ## Visualization Components
#'
#' 1. **Background Layer**: Landscape habitat values
#'    \itemize{
#'      \item Color intensity represents habitat values
#'      \item Binary: white (matrix) vs green (habitat)
#'      \item Continuous: gradient from low to high values
#'    }
#'
#' 2. **Foreground Layer**: Surviving individuals
#'    \itemize{
#'      \item Each point represents one individual at final position
#'      \item Point color represents trait value (genotype or phenotype)
#'      \item Point size controlled by point_size parameter
#'    }
#'
#' 3. **Coordinate System**:
#'    \itemize{
#'      \item X-axis: Horizontal world coordinates (centered at 0)
#'      \item Y-axis: Vertical world coordinates (centered at 0)
#'      \item Range: \eqn{\pm \text{world\_dimension}/2}
#'      \item Grid overlay marks cell boundaries
#'    }
#'
#' ## Interpretation Guidelines
#'
#' **Spatial Pattern Analysis:**
#' \itemize{
#'   \item **Clustered survivors**: Suggests habitat selection or density-dependent processes
#'   \item **Even distribution**: May indicate weak habitat selection or high dispersal
#'   \item **Edge concentration**: Possible boundary condition effects or edge habitat preference
#'   \item **Patch occupancy**: Compare survivor locations to habitat patches in binary landscapes
#' }
#'
#' **Genotype-Environment Matching:**
#'
#' When \code{color_by = "genotype"}:
#' \itemize{
#'   \item **Matching colors**: Point color similar to background -> good habitat matching
#'   \item **Contrasting colors**: Point color differs from background -> potential mismatch
#'   \item **Gradient tracking**: Points follow habitat gradient -> successful adaptation
#'   \item **Spatial sorting**: Different genotypes in different patches -> local adaptation
#' }
#'
#' Visual pattern interpretation depends on:
#' \deqn{\text{Matching Quality} \propto \frac{\text{sampling\_points} \times (1 - \text{genotype\_sds})}{\text{habitat\_selection\_temperatures}}}
#'
#' Strong matching expected with:
#' \itemize{
#'   \item High sampling_points (informed habitat choice)
#'   \item Low genotype_sds (narrow niches, specialists)
#'   \item Low habitat_selection_temperatures (deterministic selection)
#'   \item High mutation_rates (genetic variation for adaptation)
#' }
#'
#' **Phenotype vs Genotype:**
#'
#' When plasticity > 0:
#' \itemize{
#'   \item \code{color_by = "genotype"}: Shows inherited genetic values
#'   \item \code{color_by = "phenotype"}: Shows expressed traits (genotype + plastic response)
#'   \item Phenotypes typically show stronger habitat matching than genotypes
#'   \item Difference reveals plastic response magnitude
#' }
#'
#' ## Color Palette Selection
#'
#' \tabular{lll}{
#'   \strong{Palette} \tab \strong{Best For} \tab \strong{Interpretation} \cr
#'   "habitat" \tab Binary landscapes \tab Clear habitat/matrix distinction \cr
#'   "terrain" \tab General use \tab Traditional landscape colors \cr
#'   "viridis" \tab Continuous gradients \tab Accurate gradient perception \cr
#' }
#'
#' For publications:
#' \itemize{
#'   \item Use "viridis" for colorblind accessibility
#'   \item Match color scheme to manuscript style guidelines
#'   \item Consider grayscale conversion for print journals
#' }
#'
#' @examples
#' set.seed(100)
#' landscape <- create_fractal_landscape(
#'   cells_per_row = 5,
#'   fractality = 0.5,
#'   habitat_proportion = 0.6,
#'   return_as_landscape_params = TRUE
#' )
#' landscape_genetic <- create_fractal_landscape(
#'   cells_per_row = 5,
#'   fractality = 0.5,
#'   min_value = 0.35,
#'   max_value = 0.64,
#'   return_as_landscape_params = TRUE
#' )
#' result <- twolife_simulation(
#'   landscape_params = landscape,
#'   individual_params = list(
#'     initial_population_size = 15,
#'     base_birth_rate = 0.4,
#'     base_mortality_rate = 0.15
#'   ),
#'   simulation_params = list(max_events = 5),
#'   master_seed = 123
#' )
#' result_genetic <- twolife_simulation(
#'   landscape_params = landscape_genetic,
#'   individual_params = list(
#'     initial_population_size = 15,
#'     base_birth_rate = 0.4,
#'     base_mortality_rate = 0.15
#'   ),
#'   genetic_params = list(
#'     genotype_means = rnorm(15, mean = 0.5, sd = 0.15),
#'     genotype_sds = 0.15
#'   ),
#'   simulation_params = list(max_events = 5),
#'   master_seed = 456
#' )
#' plot_simulation_on_landscape(result, main = "No Genetic Variation")
#' plot_simulation_on_landscape(result_genetic, main = "With Genetic Variation")
#' @seealso \code{\link{twolife_simulation}} for running simulations,
#'   \code{\link{plot_landscape}} for landscape visualization,
#'   \code{\link{snapshot_at_time}} for temporal snapshots
#'
#' @export
plot_simulation_on_landscape <- function(simulation_result,
                                         point_size = 2,
                                         point_color = "red",
                                         point_shape = 16,
                                         color_by = "none",
                                         landscape_colors = "habitat",
                                         main = NULL,
                                         filename = NULL,
                                         add_stats = TRUE) {
  if (!inherits(simulation_result, "twolife_result")) {
    stop("simulation_result must be a twolife_result object", call. = FALSE)
  }
  
  if (!color_by %in% c("none", "genotype", "phenotype")) {
    stop("color_by must be one of: 'none', 'genotype', or 'phenotype'", call. = FALSE)
  }
  
  habitat_grid <- simulation_result$parameters$landscape$habitat
  cell_size <- simulation_result$parameters$landscape$cell_size
  survivors <- simulation_result$survivors
  n_survivors <- if (is.null(survivors)) 0L else nrow(survivors)
  if (is.null(main)) main <- paste("Simulation Results:", n_survivors, "Survivors")
  
  n_rows <- nrow(habitat_grid)
  n_cols <- ncol(habitat_grid)
  world_width <- n_cols * cell_size
  world_height <- n_rows * cell_size
  x_coords <- seq(-world_width/2, world_width/2, length.out = n_cols)
  y_coords <- seq(-world_height/2, world_height/2, length.out = n_rows)
  
  z <- t(apply(habitat_grid, 2, rev))
  
  is_binary <- all(habitat_grid %in% c(0, 1))
  z_range <- range(z, na.rm = TRUE)
  is_uniform <- z_range[1] == z_range[2]
  
  if (landscape_colors == "habitat") {
    if (is_binary) {
      # Binary landscapes: 0 = white (matrix), 1 = green (habitat)
      color_palette <- c("white", "#228B22")
      if (is_uniform) {
        color_palette <- if (z_range[1] == 0) "white" else "#228B22"
      }
    } else {
      color_palette <- colorRampPalette(c("#F5F5DC", "#90EE90", "#228B22", "#006400"))(100)
    }
  } else {
    n_colors <- if (is_binary && !is_uniform) 2 else if (is_uniform) 1 else 100
    color_palette <- switch(landscape_colors,
                            "terrain" = terrain.colors(n_colors),
                            "viridis" = {
                              if (requireNamespace("viridisLite", quietly = TRUE)) {
                                viridisLite::viridis(n_colors)
                              } else {
                                terrain.colors(n_colors)
                              }
                            },
                            terrain.colors(n_colors)
    )
  }
  
  image(x = x_coords, y = y_coords, z = z,
        col = color_palette,
        main = main,
        xlab = "X (World Coordinates)",
        ylab = "Y (World Coordinates)",
        asp = 1,
        useRaster = TRUE)
  
  grid(col = "gray80", lty = 1, lwd = 0.5)
  abline(h = 0, v = 0, col = "red", lty = 2, lwd = 1)
  
  if (n_survivors > 0) {
    if (color_by == "genotype" && !is.null(survivors$genotype)) {
      gpal <- heat.colors(100)
      gr <- range(survivors$genotype, na.rm = TRUE)
      if (gr[1] != gr[2]) {
        gscaled <- (survivors$genotype - gr[1])/(gr[2]-gr[1])
        idx <- pmax(1, pmin(100, round(gscaled*99)+1))
        ptcols <- gpal[idx]
      } else {
        ptcols <- rep(gpal[50], n_survivors)
      }
      points(survivors$x, survivors$y, col = ptcols, pch = point_shape, cex = point_size)
    } else if (color_by == "phenotype" && !is.null(survivors$phenotype)) {
      gpal <- heat.colors(100)
      gr <- range(survivors$phenotype, na.rm = TRUE)
      if (gr[1] != gr[2]) {
        gscaled <- (survivors$phenotype - gr[1])/(gr[2]-gr[1])
        idx <- pmax(1, pmin(100, round(gscaled*99)+1))
        ptcols <- gpal[idx]
      } else {
        ptcols <- rep(gpal[50], n_survivors)
      }
      points(survivors$x, survivors$y, col = ptcols, pch = point_shape, cex = point_size)
    } else {
      points(survivors$x, survivors$y, col = point_color, pch = point_shape, cex = point_size)
    }
  }
  
  invisible(simulation_result)
}

# ============================================================================
# S3 METHODS FOR TWOLIFE_RESULT CLASS
# ============================================================================

#' Print Method for TWoLife Results
#'
#' Provides a clean summary when a twolife_result object is printed to console.
#' Shows key simulation outcomes without overwhelming detail.
#'
#' @param x A 'twolife_result' object from \code{\link{twolife_simulation}}
#' @param ... Additional arguments (currently unused)
#'
#' @return The input object \code{x}, invisibly
#'
#' @examples
#' set.seed(500)
#' landscape <- create_fractal_landscape(
#'   cells_per_row = 5,
#'   fractality = 0.5,
#'   habitat_proportion = 0.6,
#'   return_as_landscape_params = TRUE
#' )
#'
#' result <- twolife_simulation(
#'   landscape_params = landscape,
#'   individual_params = list(
#'     initial_population_size = 15,
#'     base_birth_rate = 0.4,
#'     base_mortality_rate = 0.15
#'   ),
#'   simulation_params = list(max_events = 150),
#'   master_seed = 123
#' )
#'
#' # Just type the object name to see summary
#' result
#'
#' # Or explicitly call print
#' print(result)
#'
#' @export
print.twolife_result <- function(x, ...) {
  cat("TWoLife Simulation Result\n")
  cat("==========================\n\n")
  cat("Status:", x$summary$status, "\n")
  cat("Final population:", x$summary$final_population_size, "\n")
  cat("Total events:", x$summary$total_events, "\n")
  cat("Duration:", round(x$summary$duration, 2), "time units\n")
  
  if (x$summary$status == "surviving" && !is.null(x$survivors)) {
    cat("\nSurvivors:", nrow(x$survivors), "individuals\n")
    if ("genotype" %in% names(x$survivors)) {
      cat("Genotype range:",
          round(range(x$survivors$genotype), 3), "\n")
    }
  }
  
  cat("\nLandscape:",
      nrow(x$parameters$landscape$habitat), "x",
      ncol(x$parameters$landscape$habitat), "cells\n")
  
  cat("\nUse summary() for more details\n")
  cat("Use plot() to visualize results\n")
  
  invisible(x)
}

#' Summary Method for TWoLife Results
#'
#' Provides detailed statistical summary of simulation outcomes including
#' population dynamics, spatial patterns, and genetic characteristics.
#'
#' @param object A 'twolife_result' object from \code{\link{twolife_simulation}}
#' @param ... Additional arguments (currently unused)
#'
#' @return A list of class 'summary.twolife_result' containing:
#'   \describe{
#'     \item{status}{Simulation outcome ("surviving" or "extinct")}
#'     \item{n_survivors}{Final population size}
#'     \item{initial_pop}{Starting population size}
#'     \item{events}{Total count of demographic/movement events executed (births, deaths, movements, emigrations)}
#'     \item{duration}{Simulation duration in time units}
#'     \item{extinction_time}{Time of extinction (NA if surviving)}
#'     \item{landscape_size}{Dimensions of landscape (rows x columns)}
#'     \item{history_detail}{Level of event history recorded}
#'   }
#'
#' @examples
#' landscape <- create_fractal_landscape(
#'   cells_per_row = 5,
#'   fractality = 0.7,
#'   habitat_proportion = 0.4,
#'   return_as_landscape_params = TRUE
#' )
#'
#' result <- twolife_simulation(
#'   landscape_params = landscape,
#'   individual_params = list(initial_population_size = 10),
#'   simulation_params = list(max_events = 150),
#'   master_seed = 456
#' )
#'
#' # Get detailed summary
#' summary(result)
#'
#' # Store summary for further use
#' sim_summary <- summary(result)
#' sim_summary$n_survivors
#'
#' @export
summary.twolife_result <- function(object, ...) {
  structure(
    list(
      status = object$summary$status,
      n_survivors = object$summary$final_population_size,
      initial_pop = object$parameters$individual$initial_population_size,
      events = object$summary$total_events,
      duration = object$summary$duration,
      extinction_time = if (object$summary$status == "extinct")
        max(object$events$times) else NA,
      landscape_size = dim(object$parameters$landscape$habitat),
      history_detail = object$parameters$simulation$history_detail
    ),
    class = "summary.twolife_result"
  )
}

#' Print Method for TWoLife Summary
#'
#' Prints the summary of a twolife_result object in a readable format.
#'
#' @param x A 'summary.twolife_result' object from \code{summary.twolife_result()}
#' @param ... Additional arguments (currently unused)
#'
#' @return The input object \code{x}, invisibly
#'
#' @export
print.summary.twolife_result <- function(x, ...) {
  cat("TWoLife Simulation Summary\n")
  cat("===========================\n\n")
  
  cat("Population Dynamics:\n")
  cat("  Initial population:", x$initial_pop, "\n")
  cat("  Final population:", x$n_survivors, "\n")
  cat("  Status:", x$status, "\n")
  
  if (!is.na(x$extinction_time)) {
    cat("  Extinction time:", round(x$extinction_time, 2), "time units\n")
  }
  
  cat("\nSimulation Details:\n")
  cat("  Total events:", x$events, "\n")
  cat("  Duration:", round(x$duration, 2), "time units\n")
  cat("  Events per time unit:", round(x$events / x$duration, 2), "\n")
  
  cat("\nLandscape:\n")
  cat("  Dimensions:", x$landscape_size[1], "x", x$landscape_size[2], "cells\n")
  cat("  Total cells:", prod(x$landscape_size), "\n")
  
  cat("\nHistory Detail:", x$history_detail, "\n")
  
  invisible(x)
}

#' Plot Method for TWoLife Results
#'
#' Generic plot method for twolife_result objects. Calls
#' \code{\link{plot_simulation_on_landscape}} with convenient syntax.
#'
#' @param x A 'twolife_result' object from \code{\link{twolife_simulation}}
#' @param ... Additional arguments passed to \code{\link{plot_simulation_on_landscape}}.
#'   Common options include: point_size, point_color, color_by, landscape_colors, main
#'
#' @return The input object \code{x}, invisibly
#'
#' @examples
#' set.seed(200)
#' landscape <- create_fractal_landscape(
#'   cells_per_row = 5,
#'   fractality = 0.5,
#'   habitat_proportion = 0.6,
#'   return_as_landscape_params = TRUE
#' )
#'
#' result <- twolife_simulation(
#'   landscape_params = landscape,
#'   individual_params = list(
#'     initial_population_size = 15,
#'     base_birth_rate = 0.4,
#'     base_mortality_rate = 0.15
#'   ),
#'   simulation_params = list(max_events = 150),
#'   master_seed = 789
#' )
#'
#' # Simple plot using S3 method
#' plot(result)
#'
#' # With options
#' plot(result, point_size = 3, color_by = "genotype")
#'
#' # Equivalent to:
#' plot_simulation_on_landscape(result, point_size = 3, color_by = "genotype")
#'
#' @seealso \code{\link{plot_simulation_on_landscape}} for full documentation of plotting options
#'
#' @export
plot.twolife_result <- function(x, ...) {
  plot_simulation_on_landscape(x, ...)
  invisible(x)
}

# ============================================================================
# VALIDATION AND ANALYSIS FUNCTIONS
# ============================================================================

#' Validate Habitat Matching for Simulation Results
#'
#' Creates side-by-side visualizations to assess whether surviving individuals are positioned
#' in habitat matching their genotype or phenotype values. The left plot displays the
#' landscape with habitat values, while the right plot shows survivor locations
#' colored by their trait values. Calculates correlation between individual trait values
#' and habitat environmental values to quantify effectiveness of habitat selection or adaptation.
#'
#' @param simulation_result A 'twolife_result' object from \code{\link{twolife_simulation}}
#' @param main Character. Plot title. If NULL, automatically generated.
#' @param point_size Numeric. Size of points in scatter plot.
#' @param landscape_colors Character. Color scheme for landscape background.
#'   Options: "habitat", "terrain", "viridis". Note: Binary landscapes (habitat
#'   values of 0 and 1 only) automatically use the "habitat" color scheme
#'   (0 = white matrix, 1 = green habitat) regardless of this parameter, ensuring
#'   consistent visual interpretation.
#' @param color_by Character. Which trait to analyze:
#'   \itemize{
#'     \item "genotype" - Analyze genotype-habitat matching
#'     \item "phenotype" - Analyze phenotype-habitat matching (recommended for plasticity)
#'   }
#' @param show_stats Logical. If TRUE, prints statistical summary to console.
#' @param show_legend Logical. If TRUE, displays color legend on plot.
#'
#' @return Invisibly returns a data frame with columns:
#'   \describe{
#'     \item{id}{Integer. Individual ID}
#'     \item{x}{Numeric. X coordinate in world space}
#'     \item{y}{Numeric. Y coordinate in world space}
#'     \item{genotype}{Numeric. Genotype value}
#'     \item{phenotype}{Numeric. Phenotype value}
#'     \item{habitat_value}{Numeric. Habitat value at individual's location}
#'     \item{row_index}{Integer. Landscape row index (for debugging)}
#'     \item{col_index}{Integer. Landscape column index (for debugging)}
#'   }
#'
#' @details
#' What This Function Measures:
#'   This function tests whether individuals are found in habitats matching their trait
#'   values. Perfect matching would show:
#'   \itemize{
#'     \item Individuals with trait value 0.8 in habitats with values ~0.8
#'     \item Individuals with trait value 0.2 in habitats with values ~0.2
#'     \item Strong positive correlation (r close to 1.0)
#'   }
#'
#' Mathematical Basis:
#'   Calculates Pearson correlation coefficient:
#'
#'   \deqn{r = \frac{\text{Cov}(X, H)}{\sigma_X \sigma_H}}
#'
#'   Where:
#'   \itemize{
#'     \item \eqn{X} = trait values (genotype or phenotype) of survivors
#'     \item \eqn{H} = habitat values at survivor locations
#'     \item \eqn{\sigma_X}, \eqn{\sigma_H} = standard deviations
#'     \item \eqn{r} ranges from -1 (perfect negative) to +1 (perfect positive correlation)
#'   }
#'
#'   Also provides Spearman rank correlation (robust to outliers and non-linear relationships).
#'
#' Visualization Components:
#'   Creates a two-panel plot:
#'
#'   Left panel - Landscape visualization:
#'   \itemize{
#'     \item Displays the spatial landscape with habitat values as background colors
#'     \item X and Y axes show world coordinates
#'     \item Color intensity represents habitat values (darker/greener = higher values)
#'     \item Binary landscapes (0/1 values): white = matrix, green = habitat
#'     \item Continuous landscapes: color gradient shows habitat values variation
#'   }
#'
#'   Right panel - Survivor distribution:
#'   \itemize{
#'     \item Shows coordinates of all surviving individuals as points
#'     \item Each point represents one survivor's position
#'     \item Point color represents the individual's trait value (genotype or phenotype)
#'     \item Color mapping matches landscape color scheme for visual comparison
#'     \item Good matching: point colors should match the background colors in left panel
#'     \item Poor matching: points show colors different from their corresponding landscape locations
#'   }
#'
#' Connection to Simulation Parameters:
#'   Correlation strength depends on:
#'   \itemize{
#'     \item sampling_points: Higher values enable better habitat selection (stronger r)
#'     \item habitat_selection_temperatures: Lower values create stronger selection (stronger r)
#'     \item mutation_rates: Higher values create more trait variation to match habitat variation
#'     \item genotype_sds: Niche width - narrower niches (lower values) favor stronger matching
#'     \item simulation duration: Longer simulations allow more adaptation (stronger r over time)
#'   }
#'
#' Statistical Output:
#'   When show_stats = TRUE, prints to console:
#'   \itemize{
#'     \item Sample size (number of survivors)
#'     \item Pearson correlation with p-value (tests H0: r = 0)
#'     \item Spearman rank correlation (non-parametric alternative)
#'     \item Mean and range of trait values
#'     \item Mean and range of habitat values at survivor locations
#'   }
#'
#' @examples
#' set.seed(100)
#' landscape <- create_fractal_landscape(
#'   cells_per_row = 5,
#'   fractality = 0.5,
#'   habitat_proportion = 0.6,
#'   return_as_landscape_params = TRUE
#' )
#' landscape_genetic <- create_fractal_landscape(
#'   cells_per_row = 5,
#'   fractality = 0.5,
#'   min_value = 0.35,
#'   max_value = 0.64,
#'   return_as_landscape_params = TRUE
#' )
#' result <- twolife_simulation(
#'   landscape_params = landscape,
#'   individual_params = list(
#'     initial_population_size = 15,
#'     base_birth_rate = 0.4,
#'     base_mortality_rate = 0.15
#'   ),
#'   simulation_params = list(max_events = 5),
#'   master_seed = 123
#' )
#' result_genetic <- twolife_simulation(
#'   landscape_params = landscape_genetic,
#'   individual_params = list(
#'     initial_population_size = 15,
#'     base_birth_rate = 0.4,
#'     base_mortality_rate = 0.15
#'   ),
#'   genetic_params = list(
#'     genotype_means = rnorm(15, mean = 0.5, sd = 0.15),
#'     genotype_sds = 0.15
#'   ),
#'   simulation_params = list(max_events = 5),
#'   master_seed = 456
#' )
#' validation <- check_habitat_match(result)
#' validation_genetic <- check_habitat_match(result_genetic, color_by = "genotype")
#' head(validation)
#' head(validation_genetic)
#' if (!is.null(validation_genetic) && nrow(validation_genetic) > 0) {
#'   cor_value <- cor(validation_genetic$genotype, validation_genetic$habitat_value)
#'   cat("Genotype-habitat correlation:", round(cor_value, 3), "\n")
#' }
#' @seealso \code{\link{habitat_mismatch}} for fitness-based analysis
#'
#' @export
check_habitat_match <- function(simulation_result,
                                main = NULL,
                                point_size = 2,
                                landscape_colors = "terrain",
                                color_by = "phenotype",
                                show_stats = TRUE,
                                show_legend = TRUE) {
  
  if (!inherits(simulation_result, "twolife_result")) {
    stop("simulation_result must be a twolife_result object", call. = FALSE)
  }
  
  if (!color_by %in% c("genotype", "phenotype")) {
    stop("color_by must be either 'genotype' or 'phenotype'", call. = FALSE)
  }
  
  habitat_grid <- simulation_result$parameters$landscape$habitat
  cell_size <- simulation_result$parameters$landscape$cell_size
  survivors <- simulation_result$survivors
  n_survivors <- nrow(survivors)
  
  if (is.null(main)) {
    main <- paste("Habitat Matching Validation:", n_survivors, "Survivors")
  }
  
  if (n_survivors == 0) {
    cat("No survivors to validate\n")
    return(invisible(NULL))
  }
  
  n_rows <- nrow(habitat_grid)
  n_cols <- ncol(habitat_grid)
  world_width <- n_cols * cell_size
  world_height <- n_rows * cell_size
  
  col_indices <- round((survivors$x + world_width/2) / cell_size) + 1
  row_indices <- n_rows - round((survivors$y + world_height/2) / cell_size)
  
  col_indices <- pmax(1, pmin(n_cols, col_indices))
  row_indices <- pmax(1, pmin(n_rows, row_indices))
  
  habitat_at_survivors <- numeric(n_survivors)
  for (i in 1:n_survivors) {
    habitat_at_survivors[i] <- habitat_grid[row_indices[i], col_indices[i]]
  }
  
  validation_data <- data.frame(
    id = survivors$id,
    x = survivors$x,
    y = survivors$y,
    genotype = survivors$genotype,
    phenotype = survivors$phenotype,
    habitat_value = habitat_at_survivors,
    row_index = row_indices,
    col_index = col_indices
  )
  
  trait_values <- if (color_by == "phenotype") survivors$phenotype else survivors$genotype
  trait_name <- if (color_by == "phenotype") "Phenotype" else "Genotype"
  
  if (show_stats) {
    cat("\nHabitat Matching Statistics (colored by", trait_name, "):\n")
    cat("----------------------------\n")
    cat("Number of survivors:", n_survivors, "\n")
    cat(trait_name, "distribution:\n")
    cat("  Mean:", round(mean(trait_values), 3), "\n")
    cat("  Range:", round(min(trait_values), 3), "to",
        round(max(trait_values), 3), "\n")
    cat("Habitat value at survivor locations:\n")
    cat("  Mean:", round(mean(habitat_at_survivors), 3), "\n")
    cat("  Range:", round(min(habitat_at_survivors), 3), "to",
        round(max(habitat_at_survivors), 3), "\n")
    
    if (n_survivors > 1) {
      correlation <- cor(trait_values, habitat_at_survivors)
      cat("\n", trait_name, "-Habitat correlation:", round(correlation, 3), "\n", sep="")
      if (!is.na(correlation) && correlation > 0.3) {
        cat("  -> Positive match: individuals in habitat matching their", tolower(trait_name), "\n")
      }
    }
  }
  
  z_range <- range(habitat_grid, na.rm = TRUE)
  is_uniform <- z_range[1] == z_range[2]
  is_binary <- all(habitat_grid %in% c(0, 1))
  
  # FIXED: Force habitat colors for binary landscapes
  if (is_binary) {
    # Binary landscapes: always use habitat colors (0 = white, 1 = green)
    color_palette <- c("white", "#228B22")
    if (is_uniform) {
      color_palette <- if (z_range[1] == 0) "white" else "#228B22"
    }
  } else if (landscape_colors == "habitat") {
    color_palette <- colorRampPalette(c("#F5F5DC", "#90EE90", "#228B22", "#006400"))(100)
  } else {
    n_colors <- if (is_uniform) 1 else 100
    color_palette <- switch(landscape_colors,
                            "terrain" = terrain.colors(n_colors),
                            "viridis" = {
                              if (requireNamespace("viridisLite", quietly = TRUE)) {
                                viridisLite::viridis(n_colors)
                              } else {
                                terrain.colors(n_colors)
                              }
                            },
                            terrain.colors(n_colors))
  }
  
  if (!is_uniform) {
    trait_normalized <- (trait_values - z_range[1]) / (z_range[2] - z_range[1])
    trait_normalized <- pmax(0, pmin(1, trait_normalized))
    
    if (is_binary) {
      color_indices <- ifelse(trait_normalized < 0.5, 1, 2)  # FIXED: was (2, 1)
    } else {
      color_indices <- pmax(1, pmin(100, round(trait_normalized * 99) + 1))
    }
    survivor_colors <- color_palette[color_indices]
  } else {
    survivor_colors <- rep(color_palette, n_survivors)
  }
  
  par(mfrow = c(1, 2))
  
  x_coords <- seq(-world_width/2, world_width/2, length.out = n_cols)
  y_coords <- seq(-world_height/2, world_height/2, length.out = n_rows)
  z <- t(apply(habitat_grid, 2, rev))
  
  image(x = x_coords, y = y_coords, z = z,
        col = color_palette,
        main = "Landscape (Habitat Values)",
        xlab = "X (World Coordinates)",
        ylab = "Y (World Coordinates)",
        asp = 1,
        useRaster = TRUE)
  
  grid(col = "gray80", lty = 1, lwd = 0.5)
  abline(h = 0, v = 0, col = "red", lty = 2, lwd = 1)
  
  plot(survivors$x, survivors$y,
       col = survivor_colors,
       pch = 16,
       cex = point_size,
       main = paste("Survivors (Colored by", trait_name, ")"),
       xlab = "X (World Coordinates)",
       ylab = "Y (World Coordinates)",
       asp = 1,
       xlim = c(-world_width/2, world_width/2),
       ylim = c(-world_height/2, world_height/2))
  
  grid(col = "gray80", lty = 1, lwd = 0.5)
  abline(h = 0, v = 0, col = "red", lty = 2, lwd = 1)
  
  if (!is_uniform && !is_binary && show_legend) {
    legend_vals <- seq(z_range[1], z_range[2], length.out = 5)
    legend("topright",
           legend = round(legend_vals, 2),
           col = color_palette[seq(1, 100, length.out = 5)],
           pch = 16,
           title = paste(trait_name, "Value"),
           cex = 0.8,
           bg = "white")
  }
  
  mtext("If habitat selection works: dot colors should match background colors",
        side = 1, line = -1.5, outer = TRUE, cex = 0.9, col = "blue")
  
  mtext(main, side = 3, line = -1.5, outer = TRUE, cex = 1.2, font = 2)
  
  par(mfrow = c(1, 1))
  
  return(invisible(validation_data))
}

#' Calculate Genotype-Habitat Fitness Statistics
#'
#' Computes fitness-based metrics quantifying how well surviving individuals match their
#' habitat. Calculates fitness using the same Gaussian function used during the simulation,
#' providing comprehensive statistics on habitat matching quality and population fitness.
#'
#' @param simulation_result A 'twolife_result' object from \code{\link{twolife_simulation}}
#' @param return_individuals Logical. If TRUE, includes individual-level data in output.
#'
#' @return A list of class 'genotype_habitat_mismatch' with components:
#'   \describe{
#'     \item{n_survivors}{Integer. Number of surviving individuals analyzed}
#'     \item{mean_fitness}{Numeric. Mean fitness across all survivors (0-1 scale)}
#'     \item{median_fitness}{Numeric. Median fitness}
#'     \item{sd_fitness}{Numeric. Standard deviation of fitness}
#'     \item{min_fitness}{Numeric. Minimum fitness value}
#'     \item{max_fitness}{Numeric. Maximum fitness value}
#'     \item{high_fitness_count}{Integer. Number with fitness > 0.8}
#'     \item{medium_fitness_count}{Integer. Number with fitness 0.5-0.8}
#'     \item{low_fitness_count}{Integer. Number with fitness < 0.5}
#'     \item{percent_high_fitness}{Numeric. Percentage with fitness > 0.8}
#'     \item{percent_medium_fitness}{Numeric. Percentage with fitness 0.5-0.8}
#'     \item{percent_low_fitness}{Numeric. Percentage with fitness < 0.5}
#'     \item{correlation}{Numeric. Pearson correlation between phenotype and habitat (-1 to 1)}
#'     \item{mean_genotype}{Numeric. Mean genotype value of survivors}
#'     \item{sd_genotype}{Numeric. Standard deviation of genotype}
#'     \item{mean_phenotype}{Numeric. Mean phenotype value of survivors}
#'     \item{sd_phenotype}{Numeric. Standard deviation of phenotype}
#'     \item{mean_habitat}{Numeric. Mean habitat value at survivor locations}
#'     \item{sd_habitat}{Numeric. Standard deviation of habitat at survivor locations}
#'     \item{mean_niche_width}{Numeric. Mean niche width (genotype_sds) of survivors}
#'     \item{sd_niche_width}{Numeric. Standard deviation of niche width}
#'     \item{mean_absolute_mismatch}{Numeric. Mean |phenotype - habitat|}
#'     \item{median_absolute_mismatch}{Numeric. Median |phenotype - habitat|}
#'     \item{individuals}{Data frame with per-individual data (only if return_individuals=TRUE).
#'       Contains columns: id, x, y, genotype, phenotype, width, habitat_value, fitness,
#'       absolute_mismatch, raw_mismatch, row_index, col_index}
#'   }
#'
#'   The object has a custom print method that displays formatted statistics.
#'
#' @details
#' Fitness Calculation:
#'   This function uses the same fitness formula that determines demographic rates during
#'   the simulation. Fitness quantifies how well an individual's phenotype matches the
#'   habitat value at its location.
#'
#'   For generalists (niche_width > 0):
#'
#'   \deqn{f = \exp\left(-\frac{(p - h)^2}{2\sigma^2}\right)}
#'
#'   Where:
#'   \itemize{
#'     \item \eqn{f} = fitness (ranges 0 to 1)
#'     \item \eqn{p} = phenotype (individual's expressed trait value)
#'     \item \eqn{h} = habitat value at individual's location (0 to 1)
#'     \item \eqn{\sigma} = niche_width (genotype_sds parameter, tolerance to mismatch)
#'   }
#'
#'   For perfect specialists (niche_width = 0):
#'
#'   \deqn{f = \begin{cases} 
#'   1 & \text{if } |p - h| < \epsilon \\
#'   0 & \text{otherwise}
#'   \end{cases}}
#'
#'   Where \eqn{\epsilon = 0.001} is the tolerance threshold.
#'
#'   This binary fitness function means perfect specialists only survive in exactly matching habitat.
#'
#' Fitness Interpretation:
#'   \itemize{
#'     \item fitness = 1.0: Perfect match. Phenotype exactly equals habitat value.
#'     \item fitness > 0.8: High fitness. Within ~1 niche width of optimum.
#'     \item fitness = 0.5-0.8: Medium fitness. 1-2 niche widths from optimum.
#'     \item fitness < 0.5: Low fitness. More than 2 niche widths from optimum.
#'     \item fitness = 0.368 (1/e): Exactly 1 niche width from optimum (Gaussian inflection point).
#'   }
#'
#' Connection to Simulation Demography:
#'   During simulation, this fitness value affects demographic rates. For generalists
#'   (genotype_sds > 0), the mortality rate is interpolated based on fitness:
#'
#'   \deqn{\mu = \mu_{max} - f_{rel}(\mu_{max} - \mu_{min})}
#'
#'   Where:
#'   \itemize{
#'     \item \eqn{\mu_{max} = m \mu_0} (matrix_mortality_multiplier x base_mortality_rate)
#'     \item \eqn{\mu_{min} = \mu_0} (base_mortality_rate)
#'     \item \eqn{f_{rel} = f/f_{max}} (relative fitness: fitness at current habitat / fitness at optimal habitat)
#'   }
#'
#'   Higher fitness leads to lower mortality rate, higher birth rate (for generalists),
#'   and greater probability of leaving offspring.
#'
#' @examples
#' set.seed(100)
#' landscape <- create_fractal_landscape(
#'   cells_per_row = 5,
#'   fractality = 0.5,
#'   habitat_proportion = 0.6,
#'   return_as_landscape_params = TRUE
#' )
#' landscape_genetic <- create_fractal_landscape(
#'   cells_per_row = 5,
#'   fractality = 0.5,
#'   min_value = 0.35,
#'   max_value = 0.64,
#'   return_as_landscape_params = TRUE
#' )
#' result <- twolife_simulation(
#'   landscape_params = landscape,
#'   individual_params = list(
#'     initial_population_size = 15,
#'     base_birth_rate = 0.4,
#'     base_mortality_rate = 0.15
#'   ),
#'   simulation_params = list(max_events = 5),
#'   master_seed = 123
#' )
#' result_genetic <- twolife_simulation(
#'   landscape_params = landscape_genetic,
#'   individual_params = list(
#'     initial_population_size = 15,
#'     base_birth_rate = 0.4,
#'     base_mortality_rate = 0.15
#'   ),
#'   genetic_params = list(
#'     genotype_means = rnorm(15, mean = 0.5, sd = 0.15),
#'     genotype_sds = 0.15
#'   ),
#'   simulation_params = list(max_events = 5),
#'   master_seed = 456
#' )
#' mismatch <- habitat_mismatch(result)
#' mismatch_genetic <- habitat_mismatch(result_genetic)
#' print(mismatch)
#' print(mismatch_genetic)
#' cat("\nNo genetic variation:\n")
#' cat("  Mean fitness:", mismatch$mean_fitness, "\n")
#' cat("  Correlation:", mismatch$correlation, "\n")
#' cat("\nWith genetic variation:\n")
#' cat("  Mean fitness:", mismatch_genetic$mean_fitness, "\n")
#' cat("  Correlation:", mismatch_genetic$correlation, "\n")
#' cat("  % high fitness:", mismatch_genetic$percent_high_fitness, "\n")
#' @seealso \code{\link{check_habitat_match}} for visual validation
#'
#' @export
habitat_mismatch <- function(simulation_result,
                             return_individuals = FALSE) {
  
  if (!inherits(simulation_result, "twolife_result")) {
    stop("simulation_result must be a twolife_result object", call. = FALSE)
  }
  
  habitat_grid <- simulation_result$parameters$landscape$habitat
  cell_size <- simulation_result$parameters$landscape$cell_size
  survivors <- simulation_result$survivors
  n_survivors <- nrow(survivors)
  
  if (n_survivors == 0) {
    warning("No survivors to analyze", call. = FALSE)
    return(list(
      n_survivors = 0,
      mean_fitness = NA,
      median_fitness = NA,
      sd_fitness = NA,
      min_fitness = NA,
      high_fitness_count = 0,
      medium_fitness_count = 0,
      low_fitness_count = 0,
      mean_absolute_mismatch = NA,
      correlation = NA,
      mean_genotype = NA,
      mean_phenotype = NA,
      mean_habitat = NA
    ))
  }
  
  n_rows <- nrow(habitat_grid)
  n_cols <- ncol(habitat_grid)
  world_width <- n_cols * cell_size
  world_height <- n_rows * cell_size
  
  col_indices <- round((survivors$x + world_width/2) / cell_size) + 1
  row_indices <- n_rows - round((survivors$y + world_height/2) / cell_size)
  
  col_indices <- pmax(1, pmin(n_cols, col_indices))
  row_indices <- pmax(1, pmin(n_rows, row_indices))
  
  habitat_at_survivors <- numeric(n_survivors)
  for (i in 1:n_survivors) {
    habitat_at_survivors[i] <- habitat_grid[row_indices[i], col_indices[i]]
  }
  
  survivor_phenotypes <- survivors$phenotype
  survivor_widths <- survivors$width
  
  fitness_values <- numeric(n_survivors)
  for (i in 1:n_survivors) {
    if (survivor_widths[i] > 0) {
      deviation <- habitat_at_survivors[i] - survivor_phenotypes[i]
      variance <- survivor_widths[i]^2
      fitness_values[i] <- exp(-(deviation^2) / (2 * variance))
    } else {
      fitness_values[i] <- ifelse(
        abs(habitat_at_survivors[i] - survivor_phenotypes[i]) < 0.001,
        1.0,
        0.0
      )
    }
  }
  
  absolute_mismatch <- abs(survivor_phenotypes - habitat_at_survivors)
  raw_mismatch <- survivor_phenotypes - habitat_at_survivors
  
  correlation <- if (n_survivors > 1) {
    cor(survivor_phenotypes, habitat_at_survivors)
  } else {
    NA
  }
  
  high_fitness <- sum(fitness_values > 0.8)
  medium_fitness <- sum(fitness_values >= 0.5 & fitness_values <= 0.8)
  low_fitness <- sum(fitness_values < 0.5)
  
  results <- list(
    n_survivors = n_survivors,
    mean_fitness = mean(fitness_values),
    median_fitness = median(fitness_values),
    sd_fitness = sd(fitness_values),
    min_fitness = min(fitness_values),
    max_fitness = max(fitness_values),
    high_fitness_count = high_fitness,
    medium_fitness_count = medium_fitness,
    low_fitness_count = low_fitness,
    percent_high_fitness = (high_fitness / n_survivors) * 100,
    percent_medium_fitness = (medium_fitness / n_survivors) * 100,
    percent_low_fitness = (low_fitness / n_survivors) * 100,
    mean_absolute_mismatch = mean(absolute_mismatch),
    median_absolute_mismatch = median(absolute_mismatch),
    correlation = correlation,
    mean_genotype = mean(survivors$genotype),
    sd_genotype = sd(survivors$genotype),
    mean_phenotype = mean(survivor_phenotypes),
    sd_phenotype = sd(survivor_phenotypes),
    mean_habitat = mean(habitat_at_survivors),
    sd_habitat = sd(habitat_at_survivors),
    mean_niche_width = mean(survivor_widths),
    sd_niche_width = sd(survivor_widths)
  )
  
  if (return_individuals) {
    results$individuals <- data.frame(
      id = survivors$id,
      x = survivors$x,
      y = survivors$y,
      genotype = survivors$genotype,
      phenotype = survivor_phenotypes,
      width = survivor_widths,
      habitat_value = habitat_at_survivors,
      fitness = fitness_values,
      absolute_mismatch = absolute_mismatch,
      raw_mismatch = raw_mismatch,
      row_index = row_indices,
      col_index = col_indices
    )
  }
  
  class(results) <- c("genotype_habitat_mismatch", "list")
  return(results)
}


#' Print Method for Mismatch Statistics
#'
#' @param x A genotype_habitat_mismatch object
#' @param ... Additional arguments (ignored)
#'
#' @return x (invisibly)
#'
#' @export
print.genotype_habitat_mismatch <- function(x, ...) {
  cat("Genotype-Habitat Fitness Analysis\n")
  cat("==================================\n\n")
  
  if (x$n_survivors == 0) {
    cat("No survivors to analyze\n")
    return(invisible(x))
  }
  
  cat("Sample size:", x$n_survivors, "individuals\n\n")
  
  cat("Fitness Metrics (0-1 scale, based on PHENOTYPE-environment match):\n")
  cat("------------------------------------------------------------------\n")
  cat("Mean fitness:", round(x$mean_fitness, 4), "\n")
  cat("Median fitness:", round(x$median_fitness, 4), "\n")
  cat("SD of fitness:", round(x$sd_fitness, 4), "\n")
  cat("Range:", round(x$min_fitness, 4), "to", round(x$max_fitness, 4), "\n\n")
  
  cat("Fitness Quality Distribution:\n")
  cat("-----------------------------\n")
  cat("High fitness (>0.8):", x$high_fitness_count,
      paste0("(", round(x$percent_high_fitness, 1), "%)"), "\n")
  cat("Medium fitness (0.5-0.8):", x$medium_fitness_count,
      paste0("(", round(x$percent_medium_fitness, 1), "%)"), "\n")
  cat("Low fitness (<0.5):", x$low_fitness_count,
      paste0("(", round(x$percent_low_fitness, 1), "%)"), "\n\n")
  
  cat("Phenotype-Habitat Relationship:\n")
  cat("-------------------------------\n")
  cat("Correlation:", round(x$correlation, 4), "\n")
  cat("Mean genotype:", round(x$mean_genotype, 4),
      "+/- SD:", round(x$sd_genotype, 4), "\n")
  cat("Mean phenotype:", round(x$mean_phenotype, 4),
      "+/- SD:", round(x$sd_phenotype, 4), "\n")
  cat("Mean habitat at survivors:", round(x$mean_habitat, 4),
      "+/- SD:", round(x$sd_habitat, 4), "\n")
  cat("Mean niche width:", round(x$mean_niche_width, 4),
      "+/- SD:", round(x$sd_niche_width, 4), "\n")
  cat("Mean absolute mismatch:", round(x$mean_absolute_mismatch, 4),
      "(for reference)\n\n")
  
  cat("Interpretation:\n")
  cat("---------------\n")
  
  if (x$mean_fitness > 0.8) {
    cat("Excellent fitness: Individuals well-matched to their habitat\n")
  } else if (x$mean_fitness > 0.6) {
    cat("Good fitness: Most individuals reasonably well-matched\n")
  } else if (x$mean_fitness > 0.4) {
    cat("Moderate fitness: Individuals somewhat matched to habitat\n")
  } else {
    cat("Poor fitness: Individuals poorly matched to habitat\n")
  }
  
  if (!is.na(x$correlation)) {
    if (x$correlation > 0.5) {
      cat("Strong positive correlation: Effective habitat selection\n")
    } else if (x$correlation > 0.3) {
      cat("Moderate positive correlation: Some habitat selection\n")
    } else if (x$correlation > 0) {
      cat("Weak positive correlation: Limited habitat selection\n")
    } else {
      cat("No/negative correlation: Poor or absent habitat selection\n")
    }
  } else {
    cat("Correlation: NA (insufficient data variation)\n")
  }
  
  if (x$percent_high_fitness > 50) {
    cat("Majority in high-fitness locations (>50% above 0.8 fitness)\n")
  } else if (x$percent_low_fitness > 30) {
    cat("Many individuals in low-fitness locations (>30% below 0.5 fitness)\n")
  }
  
  invisible(x)
}