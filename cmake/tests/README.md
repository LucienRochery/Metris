Consumer tests for Metris

Purpose
-------
This folder contains two small consumer CMake projects used to exercise the installed
Metris package and its exported targets. They check that the package config generated
by Metris exposes any imported targets (e.g., `Eigen3::Eigen`, `NLopt::nlopt`) so that
consumer projects can link to `Metris::libMetris` without extra manual steps.

Tests
-----
- `no_deps/` — a minimal consumer that only does `find_package(Metris CONFIG REQUIRED)` and
  links to `Metris::libMetris`. This simulates a consumer that expects Metris to bring
  its own dependency information.

- `with_deps/` — a consumer that finds `Eigen3` and `NLopt` first, then `find_package(Metris)`.
  This verifies Metris works correctly when the consumer already provides/imports
  these dependencies.

How the test harness works
--------------------------
The top-level `cmake/tests/CMakeLists.txt` defines a custom target `consumer_tests` that:
1. Configures each consumer in a separate build directory.
2. Builds each consumer executable.

By default the harness assumes you have installed Metris into the project's install prefix
(i.e. you ran `make install` in your build dir). The harness passes `-DMETRIS_DIR` and
`-DCMAKE_PREFIX_PATH` pointing to the install prefix when configuring the consumers so
`find_package(Metris CONFIG REQUIRED)` finds the installed package.

Running the tests
-----------------
From the root of the repository, using the same build dir you normally use:

```bash
# configure and build Metris (done already in your workflow)
cmake -S . -B build/release_clang -DCMAKE_BUILD_TYPE=Release
cmake --build build/release_clang --target install

# configure and build the consumer tests
cmake --build build/release_clang --target consumer_tests
```

Or run via ctest (runs the `consumer_tests` custom target):

```bash
cd build/release_clang
ctest -R consumer_tests -V
```

Notes and troubleshooting
-------------------------
- If you prefer to consume Metris from the build-tree (without installing), edit `cmake/tests/CMakeLists.txt`
  and set `-DCMAKE_PREFIX_PATH=${CMAKE_BINARY_DIR}` for consumer configuration, or pass `-DMETRIS_DIR` explicitly
  to point to the build `CMAKE_CURRENT_BINARY_DIR` where `MetrisConfig.cmake` is generated.

- `with_deps/` requires system-available `Eigen3` and `NLopt` (or CMake to fetch/build them). If `find_package` for
  these fails on your machine, either install those dependencies or change the consumer CMakeLists to skip finding
  them.

- The goal of these consumers is to surface issues like the missing `Eigen3::Eigen` target in the exported package.
  If you see an error like "The link interface of target Metris::metris_deps contains Eigen3::Eigen but the target
  was not found", then the package config does not declare a `find_dependency(Eigen3)` and the fix needs to ensure
  the package config injects that find_dependency. The repository includes changes addressing this; re-run configure
  and reinstall Metris after pulling those changes.

Want improvements?
------------------
- I can add execution (running the consumer binaries) as part of the tests and add separate ctest entries for each consumer.
- I can add a CI job or a GitHub Actions workflow to run these consumers automatically.
- I can make the harness support both install-tree and build-tree consumer scenarios automatically.

If you want any of the above, tell me which and I'll add it.
