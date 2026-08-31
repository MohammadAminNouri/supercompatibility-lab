# V4.1 deterministic bundle fix

No scientific formula, numerical engine, Workspace 3 content, or V4 UX calculation was changed.

V4.1 fixes only deployment integration:

- the repository requires `app.py`, `streamlit_app.py`, `supercompatibility_r7.py`, and
  `supercompatibility_final.py` to be byte-identical;
- V4 installed only the last file, so `verify_embedded_app.py` correctly failed;
- V4.1 backs up all four, installs the exact same verified V4 file to all four,
  verifies their SHA-256 hashes, then runs the repository's own audit suite;
- any failure restores all four automatically.
