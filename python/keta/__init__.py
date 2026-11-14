"""keta (桁) — bar-structure analysis kernel, scriptable from Python.

The C++ kernel behind a Python scripting surface: a `Model` built call by
call (keyword arguments with defaults, auto ids) and analyses returning
result objects::

    import keta as k

    m = k.Model()
    a = m.add_node(0, 0, 0)
    b = m.add_node(6, 0, 0)
    steel = m.steel()
    ipe = m.add_section("IPE300")
    beam = m.add_beam(a, b, section=ipe, material=steel)
    m.fix(a)
    m.add_point_load(b, fz=-10_000)
    r = m.analyze()
    print(r.case(1).displacement(b))

Analyses: linear statics with cases and combinations, linear buckling and
imperfections, modal, second-order (P-Delta), geometrically nonlinear
(co-rotational), arc-length continuation and ruin load, elastic-catenary
cables, harmonic and PSD forced vibration, response spectra. SI units
throughout (m, N, Pa, kg, rad); results are plain tuples and lists; a missing
id raises KeyError, a bad argument or a failed analysis ValueError.

Conventions:
  * default gravity is the standard 9.80665 m/s^2;
  * catalog profiles carry their root fillets;
  * the shear areas of the parametric sections are the energy-consistent
    ones of Saint-Venant flexure for the material's Poisson ratio —
    `keta.set_shear_areas("shape_factor")` selects the classical shape
    factors 5/6 (rectangle), 0.9 (circle) for the beams added afterwards;
  * Model.to_json / Model.from_json use keta's model JSON;
  * validate() reports the first problem keta's validation finds.

Build: the CMake target `_keta` (Python development files required) writes
the package to <build>/python/keta; put that directory on PYTHONPATH.
"""

try:
    from keta._keta import (  # type: ignore[import-not-found]
        ArcTrace,
        Buckling,
        CableStatic,
        CaseResult,
        Harmonic,
        Modal,
        Model,
        Nonlinear,
        Psd,
        Results,
        SecondOrder,
        Spectrum,
        LinkLawResult,
        ec4_deflection,
        ec4_eta_min,
        ec4_hogging,
        ec4_sagging,
        rc_design_section,
        rc_shell_reinforcement,
        section_props,
        set_shear_areas,
        shear_areas,
        wood_armer,
        wood_armer_assess,
        clark_nielsen,
        response_spectrum_of,
        uniaxial_path,
        annex_a_vehicle,
        lm2_wheels,
        combine_extreme,
        set_solver_threads,
        solver_threads,
        hinge_law_at,
        track_bridge_uic774,
        soil_ka_rankine,
        soil_kp_rankine,
        soil_k0_jaky,
        soil_ka_coulomb,
        soil_kp_coulomb,
        soil_py_soft_clay,
        soil_py_stiff_clay,
        soil_py_sand,
        soil_tz_clay,
        soil_tz_sand,
        soil_qz,
        soil_earth_law,
        __version__,
    )
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        "the compiled keta._keta extension is missing — build the CMake target `_keta` "
        "and put <build>/python on PYTHONPATH"
    ) from exc

__all__ = [
    "Model",
    "Results",
    "CaseResult",
    "Buckling",
    "SecondOrder",
    "Nonlinear",
    "ArcTrace",
    "Modal",
    "CableStatic",
    "Harmonic",
    "Psd",
    "Spectrum",
    "LinkLawResult",
    "rc_design_section",
    "rc_shell_reinforcement",
    "ec4_sagging",
    "ec4_hogging",
    "ec4_deflection",
    "ec4_eta_min",
    "section_props",
    "set_shear_areas",
    "shear_areas",
    "wood_armer",
    "wood_armer_assess",
    "clark_nielsen",
    "response_spectrum_of",
    "uniaxial_path",
    "annex_a_vehicle",
    "lm2_wheels",
    "combine_extreme",
    "set_solver_threads",
    "solver_threads",
    "hinge_law_at",
    "track_bridge_uic774",
    "soil_ka_rankine",
    "soil_kp_rankine",
    "soil_k0_jaky",
    "soil_ka_coulomb",
    "soil_kp_coulomb",
    "soil_py_soft_clay",
    "soil_py_stiff_clay",
    "soil_py_sand",
    "soil_tz_clay",
    "soil_tz_sand",
    "soil_qz",
    "soil_earth_law",
    "__version__",
]

# Every other public name of the extension (the bindings added after the list
# above: phasing, influence lines, load models, fatigue, railway, wind ...),
# re-exported so `keta.<name>` reaches them without maintaining two lists.
from keta import _keta as _ext  # noqa: E402

for _name in dir(_ext):
    if not _name.startswith("_") and _name not in globals():
        globals()[_name] = getattr(_ext, _name)
        __all__.append(_name)
del _name, _ext
