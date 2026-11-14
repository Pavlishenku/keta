cmake_minimum_required(VERSION 3.20)
project(keta VERSION 1.0.0 LANGUAGES CXX)

# keta (桁) — high-performance bar-structure analysis kernel.
# SI units (m, N, Pa, kg, rad), right-handed global XYZ with Z up.

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(KETA_NATIVE "Build with -march=native" OFF)
# AddressSanitizer + UndefinedBehaviorSanitizer on the library and the programs:
# the runtime net behind the compile-time lifetime contract of keta/param
# (keta/param/datatree.h). Use a separate build directory:
#   cmake -B build-asan -G Ninja -DKETA_SANITIZE=ON && ninja -C build-asan
option(KETA_SANITIZE "Build with -fsanitize=address,undefined" OFF)
if(KETA_SANITIZE)
  add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer)
  add_link_options(-fsanitize=address,undefined)
endif()

add_library(keta STATIC
  src/keta/json.cpp
  src/keta/model.cpp
  src/keta/sections.cpp
  src/keta/element.cpp
  src/keta/assembly.cpp
  src/keta/linsolve.cpp
  src/keta/ordering.cpp
  src/keta/supernodal.cpp
  src/keta/drying.cpp
  src/keta/statics.cpp
  src/keta/section_cut.cpp
  src/keta/plastic_hinge.cpp
  src/keta/rc_beam_nl.cpp
  src/keta/shell_plastic.cpp
  src/keta/profiles.cpp
  src/keta/pushover.cpp
  src/keta/nldyn.cpp
  src/keta/thermal_profile.cpp
  src/keta/ec3.cpp
  src/keta/ec4.cpp
  src/keta/composite.cpp
  src/keta/ec2_section.cpp
  src/keta/contour_integrals.cpp
  src/keta/section_losses.cpp
  src/keta/section_service.cpp
  src/keta/section_sls.cpp
  src/keta/resal.cpp
  src/keta/wood_armer.cpp
  src/keta/soil_springs.cpp
  src/keta/link_history.cpp
  src/keta/track_bridge.cpp
  src/keta/qp.cpp
  src/keta/cable_tuning.cpp
  src/keta/rail_dynamics.cpp
  src/keta/nldyn_model.cpp
  src/keta/traffic_surface.cpp
  src/keta/link_nl.cpp
  src/keta/area_load.cpp
  src/keta/shell_tendon.cpp
  src/keta/edge_springs.cpp
  src/keta/corot_warping.cpp
  src/keta/ec2_design.cpp
  src/keta/warping.cpp
  src/keta/eigsolve.cpp
  src/keta/dynamics.cpp
  src/keta/nonlinear.cpp
  src/keta/cable.cpp
  src/keta/influence.cpp
  src/keta/launching.cpp
  src/keta/footbridge.cpp
  src/keta/tangent.cpp
  src/keta/accidental.cpp
  src/keta/traffic.cpp
  src/keta/prestress.cpp
  src/keta/prestress_tools.cpp
  src/keta/phasing.cpp
  src/keta/phasing_check.cpp
  src/keta/traveler.cpp
  src/keta/creep.cpp
  src/keta/envelope.cpp
  src/keta/etude.cpp
  src/keta/combos_ec0.cpp
  src/keta/fiber.cpp
  src/keta/fire.cpp
  src/keta/wind.cpp
  src/keta/diagnostics.cpp
  src/keta/st1.cpp
  src/keta/st1_lang.cpp
  src/keta/st1_run.cpp
  src/keta/st1_domains.cpp
  src/keta/st1_edit.cpp
  src/keta/tapered.cpp
  src/keta/bef.cpp
  src/keta/matlaws.cpp
  src/keta/surch.cpp
  src/keta/draw.cpp
  src/keta/soilzones.cpp
  src/keta/eig_ext.cpp
  src/keta/heat2d.cpp
  src/keta/substructure.cpp
  src/keta/rotating.cpp
  src/keta/hinge_law.cpp
  src/keta/continuum.cpp
  src/keta/traffic_groups.cpp
  src/keta/traffic_fatigue.cpp
  src/keta/rail_fatigue.cpp
  src/keta/traffic_multi.cpp
  src/keta/traffic_transverse.cpp
  src/keta/rail_actions.cpp
  src/keta/intensity_law.cpp
  src/keta/influence_scatter.cpp
  src/keta/wind_turbulence.cpp
  src/keta/wind_buffeting.cpp
  src/keta/wind_equivalent.cpp
  src/keta/wind_flutter.cpp
  src/keta/wind_synthesis.cpp
  src/keta/combination.cpp
  src/keta/firechain.cpp
  src/keta/springmass.cpp
  src/keta/io.cpp
  src/keta/mesh2d.cpp
  src/keta/suspension.cpp
  src/keta/influence_surface.cpp
  src/keta/section_exact.cpp
  src/keta/section_stress.cpp
  src/keta/thinwall.cpp
  src/keta/thinwall_morphology.cpp
  src/keta/shell.cpp
  src/keta/hereditary.cpp
  src/keta/fiber_frame.cpp
  src/keta/rc.cpp
  src/keta/slabmesh.cpp
  src/keta/fiber_hyst.cpp
  src/keta/catalog.cpp
  src/keta/parallel.cpp
  src/keta/script.cpp
  src/keta/param/value.cpp
  src/keta/param/helpers.cpp
  src/keta/param/components_model.cpp
  src/keta/param/components_analyze.cpp
  src/keta/param/document.cpp
  src/keta/param/expr.cpp
  src/keta/param/components_math.cpp
  src/keta/param/components_lists.cpp
  src/keta/param/components_sets.cpp
  src/keta/param/components_grid.cpp
  src/keta/param/components_stability.cpp
  src/keta/param/components_results.cpp
  src/keta/param/components_section.cpp
  src/keta/param/components_prestress.cpp
  src/keta/param/components_import.cpp
  src/keta/param/components_baked.cpp
  src/keta/param/components_cluster.cpp
  src/keta/param/units.cpp
  src/keta/param/diff.cpp
  src/keta/param/components_moving.cpp
)
target_include_directories(keta PUBLIC src third_party/eigen)
# std::thread: before glibc 2.34 the thread functions live in libpthread
find_package(Threads REQUIRED)
target_link_libraries(keta PUBLIC Threads::Threads)
if(MSVC)
  # Windows / Visual Studio: big translation units (the Python bindings), M_PI,
  # UTF-8 sources, standard conformance
  target_compile_options(keta PUBLIC /bigobj /utf-8 /permissive- /W3)
  target_compile_definitions(keta PUBLIC _USE_MATH_DEFINES NOMINMAX)
else()
  target_compile_options(keta PRIVATE -Wall -Wextra)
endif()
target_compile_definitions(keta PRIVATE KETA_VERSION="${PROJECT_VERSION}")
if(KETA_NATIVE)
  # -ffp-contract=off: no automatic FMA contraction of scalar expressions
  # (the rounding of the portable build is kept there; Eigen kernels still
  # use their explicit FMA packets)
  target_compile_options(keta PUBLIC -march=native -ffp-contract=off)
endif()

add_executable(keta_cli cli/main.cpp)
target_link_libraries(keta_cli PRIVATE keta)
add_executable(keta_bench EXCLUDE_FROM_ALL bench/keta_bench.cpp)
target_link_libraries(keta_bench PRIVATE keta)
set_target_properties(keta_cli PROPERTIES OUTPUT_NAME keta)
if(SKBUILD)
  # `pip install` also puts the command-line program `keta` on the PATH
  install(TARGETS keta_cli RUNTIME DESTINATION ${SKBUILD_SCRIPTS_DIR})
endif()


# The Python module (`import keta`), built when
# Python's development files are found. pybind11 is vendored (headers only,
# third_party/pybind11). The package lands in <build>/python/keta:
#   PYTHONPATH=build/python python3 -c "import keta"
# Not built under the sanitizers (an instrumented module needs a preloaded
# runtime in the interpreter).
option(KETA_PYTHON "Build the keta Python module" ON)
if(KETA_PYTHON AND NOT KETA_SANITIZE)
  find_package(Python3 COMPONENTS Interpreter Development.Module)
  if(Python3_Development.Module_FOUND)
    set_target_properties(keta PROPERTIES POSITION_INDEPENDENT_CODE ON)
    Python3_add_library(_keta MODULE WITH_SOABI src/keta/py/keta_py.cpp)
    target_include_directories(_keta PRIVATE third_party/pybind11/include)
    target_link_libraries(_keta PRIVATE keta)
    target_compile_definitions(_keta PRIVATE KETA_VERSION="${PROJECT_VERSION}")
    set_target_properties(_keta PROPERTIES
      LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/python/keta
      CXX_VISIBILITY_PRESET hidden)
    # the pure-Python part of the package and its typing stubs (tools/gen_pyi.py)
    set(KETA_PY_FILES __init__.py __init__.pyi _keta.pyi py.typed mcp.py docs_corpus.json.gz)
    foreach(f ${KETA_PY_FILES})
      add_custom_command(TARGET _keta POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                ${CMAKE_SOURCE_DIR}/python/keta/${f} ${CMAKE_BINARY_DIR}/python/keta/${f})
    endforeach()
    if(SKBUILD)
      # wheel layout: site-packages/keta/{_keta.*, __init__.py, stubs}
      install(TARGETS _keta LIBRARY DESTINATION keta)
      foreach(f ${KETA_PY_FILES})
        install(FILES ${CMAKE_SOURCE_DIR}/python/keta/${f} DESTINATION keta)
      endforeach()
    endif()
  else()
    message(STATUS "keta: Python development files not found — the Python module is not built")
  endif()
endif()
