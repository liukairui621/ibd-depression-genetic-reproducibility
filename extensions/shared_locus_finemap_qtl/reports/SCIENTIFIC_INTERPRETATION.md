# Fine-mapping, candidate-gene mapping and molecular colocalisation report

## Scope

The analysis was restricted a priori to blocks 154, 464 and 1671 identified in the reciprocal PLACO screen. Four GWAS were fine-mapped with SuSiE-RSS using the same 1000 Genomes EUR LD reference. Candidate genes were mapped only from credible-set/lead-variant overlap or a 100-kb window, followed by local eQTL/sQTL colocalisation and a separate Open Targets pQTL annotation.

## Fine-mapping result

- Block 154: the IBD component replicated across deLange and FinnGen (top-variant r2=0.858; six exact credible-set variants). Neither depression GWAS produced a credible set.
- Block 464: the IBD component replicated across deLange and FinnGen (top-variant r2=1.000; twelve exact credible-set variants). FinnGen depression colocalised with the IBD component, but Howard depression had no estimable credible set, so cross-depression replication was not established.
- Block 1671: Howard and FinnGen depression produced different credible sets (top-variant r2=0.108; no exact overlap); neither IBD GWAS produced a credible set.

## Candidate-gene and QTL result

- GPR25, MST1 met the primary H4>=0.8 criterion with both IBD cohorts. GPR25 was supported by whole-blood eQTL; MST1 was supported by amygdala eQTL and repeated sQTL signals in amygdala, sigmoid/transverse colon, terminal ileum and blood. Under the stricter p12=1e-6 sensitivity, the repeated MST1 sQTL findings remained above H4=0.8 in both IBD cohorts, whereas GPR25 remained strong for deLange IBD (H4=0.903) but fell below the strong threshold for FinnGen IBD (H4=0.768).
- No gene met H4>=0.8 with both depression cohorts. Single-cohort molecular candidates were FADS1, TMEM258.
- No candidate met the same molecular-QTL criterion across all four GWAS. Therefore no replicated cross-disease shared causal gene was identified.

## Protein layer

Open Targets returned 77 pQTL credible sets for MST1 and none for GPR25, FADS1 or TMEM258. The MST1 pQTL led by rs11130213 colocalised with multiple IBD GWAS (maximum H4=0.987). Different MST1 pQTL signals colocalised with depression GWAS (maximum H4=0.995). No single pQTL credible set had H4>=0.8 for both disease families (n=0). In 1000 Genomes EUR, the principal IBD pQTL lead rs11130213 was weakly correlated with depression-linked pQTL leads (r2=0.022-0.098), indicating distinct protein-regulatory signals rather than one shared causal component.

## Boundary of inference

The three PLACO blocks reproduced at the variant or broad LD-block level, with uneven molecular resolution across traits. In block 464, FinnGen depression colocalised with the IBD component, whereas Howard depression had no credible set. GPR25 and MST1 had molecular support in IBD analyses; FADS1 and TMEM258 had support in the FinnGen depression analysis. These observations did not establish a component supported by all four GWAS. Missing credible sets indicate unresolved evidence, not disease specificity.

## Technical coverage

Seventy-eight molecular traits passed the frozen local-QTL eligibility gate. GTEx and BLUEPRINT were analyzable. Nine CEDAR locus-dataset requests returned HTTP 500 from the official API and are recorded as technically unavailable, not negative. All primary ABF comparisons were retained, and SuSiE was run only after the frozen ABF H3/H4 gate.
