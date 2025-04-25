# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased] - 2025-04-23
### Added
- [GRD-795](https://jira.oicr.on.ca/browse/GRD-795) - Expanded built-in documentation (metadata changes only).

## [1.2.0] - 2024-06-25
### Added
- [GRD-797](https://jira.oicr.on.ca/browse/GRD-797) - add vidarr labels to outputs (changes to medata only)

## [1.1.1] - 2022-06-21
### Changed
- Change input lane (string) to lanes (array) to allow the olive to run the workflow with no-lane-splitting.

## [1.1.0] - 2022-06-14
### Changed
- Estimate PhiX content using bbalign instead of bbduk.
- Scatter alignment over read 1 and read 2.
### Added
- Add fastQC as subworkflow for QC.
- Add generatePhixFastqs task that runs when phiXindices argument is used.

## [1.0.1] - 2022-04-22
### Added
- [GP-3244](https://jira.oicr.on.ca/browse/GP-3244) - Workflow for PhiX determination
- A workflow similar to bcl2barcode that accepts BCL data and report PhiX percentages.
- Designed to get accurate contamination reading on failed NextSeq runs.
