# ggmanh 1.10

* Change the way positions of the variants are calculated. `MPdata` now saves
the unscaled positions and the scaling is done in `manhattan_plot` function.
(this is for the x-axis in)

* To get the actual positions used in the manhattan plot, use `calc_new_pos()`
on the `MPdata` object.

* Gap between chromosomes can be rescaled by using `chr.gap.scaling=` argument
inside `manhattan_data_preprocess` and `manhattan_plot` functions.

* Now accounts for negative positions.

# ggmanh 0.99

* Added a `NEWS.md` file to track changes to the package.
