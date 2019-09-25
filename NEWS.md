# immunarch 0.3.2

Remove MonetDBLite from dependencies because it got removed from CRAN.

## Bug fixes
* Fix a bug in MiXCR parser.


# immunarch 0.3.1

## Features
* Boxplots and barplots are now support statistical tests via the `.test` argument.

* Add parsers for old VDJtools formats.

## Documentation
* Update docs and vignettes with statistical tests information.

* Add a note for list names to vignettes.

* Documentation for clustering.

* Documentation for dimension reduction.

* Minor fix for `repOverlap` documentation.

## Bug fixes
* Fix a grouping bug in visualisations.

* Fix statistical tests from `ggpubr`.

* Fix for `geneUsage` with `.type="family"`.


# immunarch 0.3.0

## Features
* `fixVis` now supports the following legends: size, shape, color, fill, linetype.

* `fixVis` can plot figures to R console / RStudio "Plots" tab.

* `fixVis` now supports the number of columns in legends.

* Support for the AIRR file format.

* Experimental support for the 10xGenomics format.

* Save and load `immunarch` format via `repSave` and `repLoad`.

* Save and load VDJtools format via `repSave` and `repLoad`.

## Bug fixes
* `.a` and `.b` didn't passed to Tversky index.

* `fixVis` - fix a bug when users apply X/Y settings to the other axis.
