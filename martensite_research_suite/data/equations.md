# Equation ledger

## Variant/groupoid engine
Common parent subgroup:
H^P = G^P ∩ R G^D R^{-1}

Variant count:
N_V = |G^P| / |H^P|

Operator classes:
O_k = H^P g_k H^P

Operator composition (CT-reference-author convention):
(O_m,O_n) -> O_m^{-1} O_n
The reduced operator composition is generally multivalued.

## Interaction work
For a prescribed stress tensor sigma and transformation/deformation strain epsilon:
IW = sigma : epsilon

For uniaxial loading:
sigma = sigma0 (n ⊗ n)

The app reports IW in MJ/m^3 when sigma is entered in MPa.

## Compatibility diagnostics
Right stretch:
U = sqrt(F^T F)

Middle stretch residual:
r_lambda2 = |lambda_2(U) - 1|

The current public build treats this as a diagnostic, NOT as proof of full cofactor/supercompatibility.

## Reconstruction validation
Known-truth partition recovery uses ARI, NMI, homogeneity, completeness and V-measure.
Method-vs-method agreement is never labeled accuracy.
