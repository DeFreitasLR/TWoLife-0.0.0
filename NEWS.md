# TWoLife 0.1.0

## Initial CRAN Release

### Core Features

* Individual-based spatial population simulation engine with C++ backend for computational efficiency
* Full genetic evolution system with mutation, selection, and quantitative traits
* Habitat selection behavior with sampling-based decision making
* Density-dependent demographic rates (birth, death, dispersal) with local or global density effects
* Multiple boundary conditions: absorbing, periodic, and reflective edges
* Fractal landscape generation with configurable spatial autocorrelation
* Comprehensive visualization system with multiple color schemes

### Main Functions

* `twolife_simulation()` - Primary simulation function with flexible parameter system
* `create_fractal_landscape()` - Generate realistic heterogeneous landscapes
* `create_corner_landscape()` - Create simple test landscapes
* `plot_simulation_on_landscape()` - Visualize population distributions on landscapes
* `plot_landscape_image()` - Display landscape patterns with custom coloring
* `compute_population_size()` - Extract and analyze population trajectories
* `quick_plot_result()` - Simplified result visualization

### Simulation Capabilities

* Spatial dynamics in continuous 2D space with flexible landscape sizes
* Continuous genetic traits with normal distributions and phenotypic plasticity
* Active habitat choice through multi-point sampling
* Directional movement with configurable vision angles and step lengths
* Differential mortality and dispersal rates in habitat versus matrix
* Stochastic event timing using exponential distributions
* Population tracking and demographic analysis

### Documentation

* Comprehensive vignette with step-by-step tutorials and biological interpretation
* Complete help documentation for all 16 exported functions
* Working examples demonstrating progressive complexity
* Detailed parameter descriptions with biological context
* Mathematical formulations for demographic and genetic processes

### Technical Implementation

* Fast C++ core via Rcpp for performance-critical operations
* Memory-safe implementation with smart pointers and bounds checking
* Proper random number generation with seed control for reproducibility
* Input validation with informative error messages
* Cross-platform compatibility (Windows, macOS, Linux)
* Clean namespace with selective imports

### Package Structure

* 5 R source files with clear organization
* 9 C++ source files with modular design
* 16 documentation files (.Rd format)
* Test suite using testthat framework
* Example datasets: `several_small` and `single_large`
* Example scripts in inst/examples directory

### Dependencies

* R >= 4.0.0
* Rcpp >= 1.0.0 for C++ integration
* mathjaxr for mathematical notation in documentation
* Standard R packages: stats, graphics, grDevices, utils

### Suggested Packages

* testthat >= 3.0.0 for testing
* knitr for vignette building
* rmarkdown for documentation
* viridisLite for enhanced color palettes

### System Requirements

* C++11 compiler (automatically handled by R)
* Windows: Rtools for compilation
* macOS: Xcode command line tools
* Linux: g++ compiler

### License and Citation

* GPL-3 license
* CITATION file included for proper academic attribution
* ORCID identifier for author verification

### Quality Assurance

* Comprehensive input validation
* Memory-safe C++ implementation
* Thread-safe atomic operations
* Numerical stability in edge cases
* Complete roxygen2 documentation
* Reproducible examples in all help files
