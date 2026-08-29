from __future__ import annotations

import numpy as np

from .core import C_A_TO_M, C_M_TO_A

ENGINE_VERSION = "2.0.0"
BUILT_IN_TRANSFORMATION = "B2 austenite -> B19' martensite"


def method_provenance(cc1_tol: float, cc2_tol: float) -> dict[str, object]:
    """Machine-readable method metadata for result exports."""
    return {
        "engine_version": ENGINE_VERSION,
        "built_in_transformation": BUILT_IN_TRANSFORMATION,
        "correspondence_A_to_M": np.asarray(C_A_TO_M, float).tolist(),
        "correspondence_M_to_A": np.asarray(C_M_TO_A, float).tolist(),
        "numerical_tolerances": {
            "cc1_abs_lambda2_minus_1": float(cc1_tol),
            "cc2_normalized": float(cc2_tol),
        },
        "scientific_separation": {
            "austenite_martensite": "lambda2 and CMC degeneracy",
            "classical_supercompatibility": "CC1-CC3 for classified Type-I, Type-II and Compound domains",
            "metric_intercompatibility": "CMC habit plane + SMC shear + twin shear/shear epsilon",
        },
    }
