# TWoLife

> Individual-based spatial population simulations with genetic evolution and habitat selection

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

TWoLife is an R package for running spatially-explicit, individual-based population models. It simulates organisms that move, reproduce, and die across heterogeneous landscapes, with genetic evolution and behavioral habitat selection.

## Features

- **Spatial Population Dynamics**: Individuals move across realistic landscapes
- **Genetic Evolution**: Mutation, selection, and adaptation to local conditions  
- **Habitat Selection**: Active behavioral choice of optimal environments
- **Demographic Processes**: Density-dependent birth and death rates
- **Fast Simulation**: C++ backend via Rcpp for computational efficiency
- **Flexible Landscapes**: Binary, continuous, or fractal habitat patterns
- **Comprehensive Visualization**: Built-in plotting functions for results

## Installation

### From CRAN (Coming Soon)

TWoLife is currently under review for CRAN. Once accepted, install with:

```r
install.packages("TWoLife")
```

### Development Version from GitHub

```r
# Install development tools if needed
install.packages("devtools")

# Install TWoLife from GitHub
devtools::install_github("DeFreitasLR/TWoLife-0.0.0")
```

## Dependencies

The package requires:
- R >= 4.0.0
- Rcpp >= 1.0.0
- mathjaxr (for documentation)

Optional for enhanced visualization:
- viridisLite

Note: Fractal landscape generation is implemented internally using base R functions.

## Quick Start

```r
library(TWoLife)

# 1. Create a test landscape
landscape <- create_fractal_landscape(
  cells_per_row = 20,
  fractality = 0.6,
  habitat_proportion = 0.4
)

# 2. Run a basic simulation
result <- twolife_simulation(
  landscape_params = list(habitat = landscape),
  individual_params = list(initial_population_size = 50),
  simulation_params = list(max_events = 1000)
)

# 3. View results
print(result)
plot_simulation_on_landscape(result)

# 4. Analyze population trajectory
trajectory <- compute_population_size(result)
plot(trajectory$time, trajectory$population_size,
     type = "l", main = "Population Over Time",
     xlab = "Time", ylab = "Population Size")
```

## Advanced Usage

### Simulation with Genetic Diversity

```r
# Create landscape
landscape <- create_fractal_landscape(
  cells_per_row = 20,
  fractality = 0.6,
  habitat_proportion = 0.5
)

# Run simulation with genetic variation
result <- twolife_simulation(
  landscape_params = list(habitat = landscape),
  individual_params = list(
    initial_population_size = 40,
    base_birth_rate = 0.4,
    base_mortality_rate = 0.15
  ),
  genetic_params = list(
    genotype_means = runif(40, 0.2, 0.8),  # Diverse starting genotypes
    genotype_sds = rep(0.15, 40),           # Niche width tolerance
    mutation_rates = rep(0.03, 40),         # Evolution rate
    plasticities = rep(0.02, 40),           # Phenotypic flexibility
    sampling_points = rep(15, 40)           # Habitat selection intensity
  ),
  simulation_params = list(max_events = 2000),
  master_seed = 123
)

# Visualize results
plot_simulation_on_landscape(result, color_by = "genotype")
```

### Compare Different Landscapes

```r
# Generate different landscape types
continuous <- create_fractal_landscape(
  cells_per_row = 15, 
  fractality = 0.7,
  habitat_proportion = 0.6
)

fragmented <- create_fractal_landscape(
  cells_per_row = 15, 
  fractality = 0.3, 
  habitat_proportion = 0.3
)

# Run simulations
results <- lapply(list(continuous = continuous, fragmented = fragmented), 
  function(landscape) {
    twolife_simulation(
      landscape_params = list(habitat = landscape),
      individual_params = list(initial_population_size = 30),
      simulation_params = list(max_events = 500),
      master_seed = 456
    )
  }
)

# Compare final population sizes
sapply(results, function(r) r$summary$final_population_size)
```

## Documentation

For detailed examples and biological interpretation:

```r
# View vignette
vignette("introduction", package = "TWoLife")

# Function help
?twolife_simulation        # Main simulation function
?create_fractal_landscape  # Landscape generation
?plot_simulation_on_landscape  # Result visualization
?compute_population_size   # Extract population trajectories

# Complete function list
help(package = "TWoLife")
```

## Main Functions

| Function | Purpose |
|----------|---------|
| `twolife_simulation()` | Run complete individual-based simulation |
| `create_fractal_landscape()` | Generate realistic fractal landscapes |
| `create_corner_landscape()` | Create simple test landscapes |
| `plot_simulation_on_landscape()` | Visualize simulation results on landscape |
| `plot_landscape_image()` | Display landscape patterns |
| `plot_landscape_world_coords()` | Plot landscape with world coordinates |
| `compute_population_size()` | Extract population trajectories over time |
| `quick_plot_result()` | Quick visualization of results |

## Key Parameters

### Landscape Parameters
- `habitat`: Habitat quality matrix (required)
- `cell_size`: Spatial resolution (default: 1.0)
- `boundary_condition`: Edge behavior ("periodic", "absorbing", "reflective")
- `density_type`: Density dependence scale ("local" or "global")

### Individual Parameters
- `initial_population_size`: Starting number of individuals (default: 200)
- `neighbor_radius`: Density dependence spatial scale (default: 2.0)
- `step_length`: Movement distance per dispersal event (default: 1.0)
- `base_birth_rate`, `base_mortality_rate`: Demographic rates
- `base_dispersal_rate`: Movement rate

### Genetic Parameters
- `genotype_means`: Optimal habitat values for each individual
- `genotype_sds`: Niche width (habitat tolerance)
- `mutation_rates`: Evolution rate per generation
- `plasticities`: Phenotypic plasticity
- `sampling_points`: Habitat selection behavior intensity

### Simulation Parameters
- `max_events`: Maximum number of events to simulate (auto-calculated if not specified)
- `neutral_mode`: Disable habitat selection (default: FALSE)
- `output_file`: Optional file path for detailed output

## Examples and Tutorials

See the package vignette for:
- Basic workflow examples
- Genetic evolution demonstrations
- Landscape effect comparisons
- Parameter sensitivity analysis
- Biological interpretation of results

Access with: `vignette("introduction", package = "TWoLife")`

## Citation

If you use TWoLife in your research, please cite:

```r
citation("TWoLife")
```

Or in text:

Freitas, L. (2025). TWoLife: Individual-Based Spatial Population Simulations with Genetic Evolution. R package version 0.1.0. https://github.com/DeFreitasLR/TWoLife-0.0.0

## License

GPL-3

See the LICENSE file for details.

## Bug Reports and Feature Requests

Please report bugs or request features via [GitHub Issues](https://github.com/DeFreitasLR/TWoLife-0.0.0/issues).

## Contact

**Maintainer**: Lucas Freitas  
**Email**: rodriguesdefreitas@gmail.com  
**ORCID**: [0000-0002-2773-0981](https://orcid.org/0000-0002-2773-0981)

## Acknowledgments

Contributors: piLaboratory

## Development

For developers working on the package:

```r
# Generate documentation
devtools::document()

# Run tests
devtools::test()

# Check package
devtools::check()

# Install locally
devtools::install()
```

## Contributing

Contributions are welcome! To contribute:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Make changes with tests
4. Commit your changes (`git commit -m 'Add amazing feature'`)
5. Push to the branch (`git push origin feature/amazing-feature`)
6. Open a Pull Request

Please ensure all tests pass and documentation is updated.

---

**Note**: TWoLife is designed for ecological and evolutionary research. The package implements spatially-explicit individual-based models with continuous genetics, making it suitable for studying eco-evolutionary dynamics, range expansion, and adaptive landscapes in heterogeneous environments.
