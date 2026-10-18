# =============================================================================
# 02 - Voxel-wise subject-prevalence lower bounds (closed testing)
# -----------------------------------------------------------------------------
# For each voxel i, computes TD(S_i): a simultaneous (1 - alpha) lower
# confidence bound on the number of subjects active at that voxel.
# Input:    output/pvals.RData (from 01_download_and_preprocess.R)
# Output:   output/prevalence_bounds.RData
#           output/prevalence_lower_bound.nii   (used by the Shiny app)
#           output/prevalence_summary.csv       (Table 4.1 in the thesis)
# =============================================================================

library(RNifti)
library(hommel)
library(ARIbrain)

out_dir <- "output"
alpha   <- 0.05

load(file.path(out_dir, "pvals.RData"))       # p, mask, copes

mask_v <- as.vector(mask)
m      <- sum(mask_v)                          # voxels in the mask
n      <- length(p) / m                        # subjects


# -----------------------------------------------------------------------------
# 1. Closed testing with Hommel / Simes on all m * N hypotheses
# -----------------------------------------------------------------------------

hom <- hommel(p)

td <- numeric(m)
for (j in seq_len(m)) {
  if (j %% 5000 == 0) message(j, " / ", m)
  ix    <- (1:n) + (j - 1) * n                 # the N subjects for voxel j
  td[j] <- discoveries(hom, ix = ix, alpha = alpha)
}

save(td, file = file.path(out_dir, "prevalence_bounds.RData"))


# -----------------------------------------------------------------------------
# 2. Write the prevalence-bound map back to NIfTI
# -----------------------------------------------------------------------------

td_img <- mask
td_img[] <- NA
td_img[mask == 1] <- td
writeNifti(td_img, file.path(out_dir, "prevalence_lower_bound.nii"))

# 3D array version for plotting / clustering (NA outside the mask)
td_map <- array(NA, dim = dim(mask))
td_map[mask == 1] <- td


# -----------------------------------------------------------------------------
# 3. Summary table across prevalence thresholds u
# -----------------------------------------------------------------------------

exact_counts    <- table(factor(td, levels = 0:n))
at_least_counts <- rev(cumsum(rev(as.vector(exact_counts))))

prevalence_summary <- data.frame(
  u               = 0:n,
  voxels_exact    = as.vector(exact_counts),
  voxels_at_least = at_least_counts,
  pct_of_mask     = round(at_least_counts / m * 100, 2)
)

print(prevalence_summary)
write.csv(prevalence_summary, file.path(out_dir, "prevalence_summary.csv"), row.names = FALSE)


# -----------------------------------------------------------------------------
# 4. u-out-of-N cluster maps
# -----------------------------------------------------------------------------
# For each u, keeps connected clusters (>= min_size voxels) of voxels with
# TD(S_i) >= u. The final map stores, for each voxel, the highest u for which
# it belongs to such a cluster.

make_cl <- function(u, min_size = 10) {
  cond <- td_map >= u
  cond[is.na(cond)] <- FALSE
  if (sum(cond) < min_size) return(array(0, dim = dim(td_map)))

  cl    <- cluster_threshold(cond)
  sizes <- table(cl)
  keep  <- setdiff(names(sizes[sizes >= min_size]), "0")

  cl[!cl %in% keep] <- 0
  (cl > 0) * u
}

cl_final <- make_cl(1)
for (u in 2:n) cl_final <- pmax(cl_final, make_cl(u), na.rm = TRUE)

table(as.vector(cl_final))

# Plotting (Figures 4.5 and 4.6) uses vis_brain() from R/helpers.R:
# vis_brain(td_map,   mask, x = 31, y = 6, z = 20, type = "abs")
# vis_brain(cl_final, mask, x = 31, y = 6, z = 20)
