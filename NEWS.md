# ggmanh 1.10

## New Feature

* Added a visualization that is an extension of manhattan plot called binned manhattan plot.
Instead of plotting each individual variant, the variants are binned based on position (x-axis) and 
log10 p-value (y-axis). More information can be found in the vignette.

* Preprocessing function is `binned_manhattan_preprocess` and plotting function is `binned_manhattan_plot`.

## Internal Changes

* Change the way positions of the variants are calculated. `MPdata` now saves
the unscaled positions and the scaling is done in `manhattan_plot` function.
(this is for the x-axis)

* To get the actual relative positions used in the manhattan plot, use `calc_new_pos()`
on the `MPdata` object.

* Gap between chromosomes can be rescaled by using `chr.gap.scaling=` argument
inside `manhattan_data_preprocess` and `manhattan_plot` functions.

* Now accounts for negative positions.

# ggmanh 0.99

* Added a `NEWS.md` file to track changes to the package.
