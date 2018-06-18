# ahnr 0.3.1

* Removing generics.
* plot function name changed to visualize.

# ahnr 0.3.0

* Generic functions updated: summary, predict, plot.
* AHNnD function name changed to fit.
* SimAHNnD function name changed to predict.
* Fix in `predict` function: ncol(H[[1]]) --> max(unlist(sapply(H, ncol)))
* Fix in `predict` function to handle empty molecules: if (is.null(H[[i]])) {next}

# ahnr 0.2.1

* Correction to the labels of the second plot in the vignette.

# ahnr 0.2.0

* Added legends in the plots of the vignette file.
* Correction due to the changes to the R development sources: "arithmetic with zero-column data.frames".

# ahnr 0.1.0

* Added a `NEWS.md` file to track changes to the package.

