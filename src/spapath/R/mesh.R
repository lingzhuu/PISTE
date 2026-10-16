# =============================================================================
# mesh.R
#
# Triangular mesh for the SPDE spatial fields of screen_joint_fit().
# =============================================================================


#' Build the mesh for the spatial fields
#'
#' Builds a 2D triangular mesh over the spots with [fmesher::fm_mesh_2d()].
#' Build it once for a section and pass it to [screen_joint_fit()] for every
#' gene set. The defaults assume coordinates rescaled to the unit square
#' (shifted to start at zero and divided by their largest value).
#'
#' @param xy N x 2 matrix of spot coordinates.
#' @param max_edge Largest triangle edge.
#' @param cutoff Minimum distance between mesh nodes; closer spots share a
#'   node.
#' @param offset Inner and outer extension of the mesh around the spots.
#' @return An `fm_mesh_2d` mesh; `mesh$n` is its number of nodes.
#' @export
build_mesh <- function(xy, max_edge = 0.25, cutoff = 0.08,
                       offset = c(0.1, 0.2)) {
  xy <- as.matrix(xy)
  stopifnot(max_edge > 0, cutoff >= 0, length(offset) == 2)
  fmesher::fm_mesh_2d(loc = xy, max.edge = max_edge, cutoff = cutoff,
                      offset = offset)
}
