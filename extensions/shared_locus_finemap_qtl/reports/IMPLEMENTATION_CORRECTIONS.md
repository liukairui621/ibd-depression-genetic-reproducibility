# Implementation corrections and non-result-affecting events

The following events occurred during pipeline implementation. Each was corrected before the final result set was generated.

1. The first R data.table pass used an invalid `..trait` scoping expression. The error occurred before any SuSiE or colocalisation model was fitted. The expression was corrected and the complete GWAS stage was rerun.
2. The initial GENCODE parser assumed chromosome labels without the `chr` prefix and therefore produced no mapped genes. Chromosome normalization was corrected before candidate-mapping outputs were generated.
3. The first local-QTL preparation pass handled duplicate variants globally rather than within each molecular trait. Duplicate handling was changed to operate separately for every molecular trait, and all QTL preparation and colocalisation outputs were regenerated.
4. The system tabix build lacked HTTPS support. A project-specific Python environment was created with pysam and an explicit CA bundle; final remote tabix retrievals used this environment.
5. Repeated remote FTP connections were consolidated after an inefficient early retrieval attempt. The orphaned process was stopped. Only the complete final retrieval is represented in the manifests and checksums.
6. Direct tabix access to the approximately 5-GB CEDAR resources stalled. The workflow was switched to the official gene-and-region REST endpoint. All nine requested CEDAR locus-dataset combinations returned HTTP 500 and are retained as technical failures, not biological negatives.
7. The summary plotting stage encountered a pandas/MKL compatibility issue in the base conda environment. It was rerun in the project environment with pinned pandas, numpy and matplotlib versions. Numerical source tables were unchanged.

Final validation covers model convergence, result-grid completeness, prespecified evidence gates, disease-family replication and the separation of MST1 pQTL signals. All 20 checks passed; see `reports/FINAL_VALIDATION.md`.

