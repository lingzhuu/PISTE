# =============================================================================
# mesh.R
#
# Triangular mesh for the SPDE spatial fields of screen_joint_fit().
# =============================================================================


#' Build the mesh for the spatial fields
#'
#' Builds a 2D triangular mesh over the spots with [fmesher::fm_mesh_2d()].
#' Build it once for a section and pass it to [screen_joint_fit()] for every
#' gene set. Coordinates are assumed rescaled to the unit square (shifted to
#' start at zero and divided by their largest value).
#'
#' The mesh must resolve the spatial field at its practical range `h`: a
#' Matern field cannot be represented on triangles that are not much smaller
#' than its range. Unless given explicitly, the mesh parameters are derived
#' from `h` with the usual SPDE rules of thumb,
#' * `max_edge = c(h / 3, h)`: inner edge at most a third of the range, outer
#'   extension edge equal to the range;
#' * `cutoff = h / 9`: spots closer than this share a node (also prevents
#'   tiny triangles);
#' * `offset = c(h / 2, 1.5 * h)`: inner and outer extension, the outer at
#'   least the range so that boundary effects stay outside the tissue.
#' For Visium (spot spacing 100 um, about 0.015 of a 6-7 mm section) this
#' gives roughly 250 nodes at h = 0.3, 1,000 at h = 0.1 and every spot its
#' own node at h <= 0.06. Computational cost grows with the node count times
#' the number of genes.
#'
#' @param xy N x 2 matrix of spot coordinates.
#' @param h Practical range of the spatial field in the same units as `xy`;
#'   used only to derive the defaults below.
#' @param max_edge Largest triangle edge; one value or `c(inner, outer)`.
#' @param cutoff Minimum distance between mesh nodes; closer spots share a
#'   node.
#' @param offset Inner and outer extension of the mesh around the spots.
#' @return An `fm_mesh_2d` mesh; `mesh$n` is its number of nodes.
#' @export
build_mesh <- function(xy, h = 0.3, max_edge = NULL, cutoff = NULL,
                       offset = NULL) {
  xy <- as.matrix(xy)
  stopifnot(is.numeric(h), length(h) == 1, is.finite(h), h > 0)
  if (is.null(max_edge)) max_edge <- c(h / 3, h)
  if (is.null(cutoff)) cutoff <- h / 9
  if (is.null(offset)) offset <- c(h / 2, 1.5 * h)
  stopifnot(all(max_edge > 0), length(max_edge) %in% 1:2,
            cutoff >= 0, length(offset) == 2)
  mesh <- fmesher::fm_mesh_2d(loc = xy, max.edge = max_edge, cutoff = cutoff,
                              offset = offset)
  attr(mesh, "spapath_mesh_args") <- list(h = h, max_edge = max_edge,
                                          cutoff = cutoff, offset = offset)
  mesh
}
