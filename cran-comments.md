# CRAN Submission Comments — TWoLife 0.0.0

## Test Environments

- Windows 11 x64 (build 26200), R 4.4.2 (2024-10-31 ucrt), GCC 13.3.0
- R-hub (to be run before submission)

## R CMD check Results

0 errors | 0 warnings | 2 notes

### Note 1: C++14 specification

```
checking C++ specification ... NOTE
    Specified C++14: please drop specification unless essential
```

C++14 is required and cannot be dropped. The package uses `std::make_unique`
(introduced in C++14) throughout the C++ backend for memory-safe management
of `Individual` and `Landscape` objects in the Gillespie simulation engine
(see `src/landscape.cpp` and `src/twolife_rcpp.cpp`). Replacing `std::make_unique`
with C++11-compatible alternatives would require reimplementing smart-pointer
ownership semantics that are central to the simulation's correctness.

### Note 2: Unable to verify current time

```
checking for future file timestamps ... NOTE
    unable to verify current time
```

This is an infrastructure issue unrelated to the package itself. It occurs
when the check environment cannot reach a time server to verify file timestamps.

## Downstream Dependencies

This is a new package with no downstream dependencies.
