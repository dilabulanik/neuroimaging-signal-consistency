# =============================================================================
# 03 - Comparison with All-Resolutions Inference (ARI)
# -----------------------------------------------------------------------------
# Group-level one-sample t-tests, supra-threshold clusters (tau = 3.2,
# >= 20 voxels) and ARI lower bounds on the true discovery proportion (TDP)
# within each cluster (Figure 4.7 in the thesis).
# Input:    output/pvals.RData (from 01_download_and_preprocess.R)
# =============================================================================

library(ARIbrain)
library(pARI)

out_dir <- "output"
load(file.path(out_dir, "pvals.RData"))       # p, mask, copes

mask_v <- as.vector(mask)


# -----------------------------------------------------------------------------
# 1. Group-level test statistics and p-values
# -----------------------------------------------------------------------------

tmp      <- pARI:::oneSamplePar(t(copes), alternative = "two.sided")
ts_group <- tmp$Test
p_group  <- tmp$pv

# Map a vector over masked voxels back to a 3D array
from_v_to_map <- function(X, pvalue = FALSE) {
  out <- mask_v
  out[mask_v == 1] <- X
  out[mask_v == 0] <- ifelse(pvalue, 1, 0)
  array(out, dim = dim(mask))
}

tmap_group <- from_v_to_map(ts_group)
pmap_group <- from_v_to_map(p_group)


# -----------------------------------------------------------------------------
# 2. Supra-threshold clusters
# -----------------------------------------------------------------------------

get_cl <- function(thr = 3.2, min_size = 20) {
  cl <- cluster_threshold(abs(tmap_group) > thr)
  tb <- table(cl)
  small <- as.numeric(names(tb[tb < min_size]))
  cl[array(cl %in% small, dim = dim(cl))] <- 0
  cl
}

clusters <- get_cl(thr = 3.2, min_size = 20)


# -----------------------------------------------------------------------------
# 3. ARI: lower bound on TDP within each cluster
# -----------------------------------------------------------------------------

ari <- ARI(pmap_group, clusters = clusters, mask = mask)
ari <- ari[-nrow(ari), ]                      # drop the "outside clusters" row
print(ari)

tdp_map <- array(0, dim = dim(clusters))
cl_ids  <- as.numeric(sub("^cl", "", rownames(ari)))
for (i in seq_along(cl_ids)) tdp_map[clusters == cl_ids[i]] <- ari[i, 4]
tdp_map[mask == 0] <- NA

# Plotting (Figure 4.7) uses vis_brain() from R/helpers.R:
# vis_brain(tdp_map, mask, x = 31, y = 6, z = 20, type = "tdp",
#           limits = c(0, max(tdp_map, na.rm = TRUE)))
