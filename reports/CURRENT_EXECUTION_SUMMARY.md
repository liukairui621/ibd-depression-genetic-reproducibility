# Current corrected execution summary

Date: 2026-08-06

## Input correction

The Howard depression input was verified as the public PGC-UK Biobank
genome-wide release excluding 23andMe (University of Edinburgh DataShare,
doi:10.7488/ds/2458). The corrected sample size is 500,199: 170,756 cases and
329,443 controls. The raw file contains 8,483,301 variants and has SHA256
`c422cbcdd4296958273fb38e4a29e98a498022acd21f3e4b13d93232511f17e3`.

## Corrected global analysis

- Howard observed-scale h2: 0.0480 (SE 0.0021; z 22.86).
- Howard depression x FinnGen strict IBD: rg 0.1205 (SE 0.0336; P 0.0003).
- FinnGen depression x de Lange IBD: rg 0.0829 (SE 0.0351; P 0.0181).
- The five Howard-involving rg estimates were numerically unchanged after the
  N correction.
- The reciprocal estimates are reported separately. Their standardized
  sampling correlation is 0.120, so fixed-effect pooling is not used.

## Corrected local analysis

- Eligible bivariate tests: 47 for Howard depression x FinnGen IBD and 63 for
  FinnGen depression x de Lange IBD; 110 jointly corrected tests.
- Four pair-specific signals pass BH q<0.05.
- Three blocks are estimable in both reciprocal pairs.
- Zero blocks meet the same-region, same-direction, q<0.05-in-both rule.
- Independent R verification reproduced all BH q values with maximum absolute
  difference 4.94e-13.
- Approximate median absolute rho detectable with 80% power at 0.05/110 is
  0.715 and 0.741 for the two pairs. This is a descriptive normal-theory
  diagnostic, not a formal LAVA power simulation.

## Extension sensitivity

After excluding six cells with absolute genetic-covariance-intercept z at
least 2, all 18 retained estimates remained positive and nine passed BH-FDR.
The direction of the global result therefore did not depend on the flagged
cells. No stable depression-definition or IBD-subtype mechanism is inferred.

## Evidence boundary

The corrected data support small positive genome-wide shared susceptibility.
They do not establish individual prediction, causal direction, a stable local
component, a candidate gene, a shared molecular mechanism, or a common
therapeutic target.
