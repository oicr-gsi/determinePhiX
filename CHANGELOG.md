## 1.1.0 - 2022-06-14
- Estimate PhiX content using bbalign instead of bbduk.
- Scatter alignment over read 1 and read 2.
- Add fastQC as subworkflow for QC.
- Add generatePhixFastqs task that runs when phiXindices argument is used.

## 1.0.1 - 2022-04-22
- [GP-3244](https://jira.oicr.on.ca/browse/GP-3244) - Workflow for PhiX determination
    - A workflow similar to bcl2barcode that accepts BCL data and report PhiX percentages.
    - Designed to get accurate contamination reading on failed NextSeq runs.
