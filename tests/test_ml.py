import numpy as np
import pandas as pd

from src.ml import TARGET_COLUMNS, physics_screen_predictions, predict_lattice, train_lattice_model


def test_linear_ml_pipeline_and_physics_screen():
    rng = np.random.default_rng(3)
    n = 30
    cu = rng.uniform(0,10,n)
    fe = rng.uniform(0,4,n)
    df = pd.DataFrame({
        "Cu_atpct": cu,
        "Fe_atpct": fe,
        "a_B2_A": 3.01 + 0.001*cu - 0.0004*fe,
        "a_B19p_A": 2.90 + 0.002*cu,
        "b_B19p_A": 4.10 + 0.01*cu,
        "c_B19p_A": 4.64 + 0.004*fe,
        "beta_deg": 97.8 - 0.03*cu + 0.02*fe,
    })
    res = train_lattice_model(df, ["Cu_atpct","Fe_atpct"], model_name="Linear regression", folds=5)
    assert (res.cv_metrics["R2"] > 0.99).all()
    cand = pd.DataFrame({"Cu_atpct":[3.0,6.0],"Fe_atpct":[1.0,2.0]})
    pred = predict_lattice(res,cand)
    for t in TARGET_COLUMNS:
        assert f"pred_{t}" in pred
    screened = physics_screen_predictions(pred)
    assert screened["physics_valid"].all()
