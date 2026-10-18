# =============================================================================
# 01 - Download data, build the analysis mask and compute p-values
# -----------------------------------------------------------------------------
# Data:     NeuroVault collection #2447, contrast "noun/verb judgement > baseline"
# Input:    subject-level t-statistic maps (downloaded via the NeuroVault API)
# Output:   data/neurovault_2447/   tmap_XX.nii.gz, pmap_XX.nii.gz, mask.nii.gz
#           output/pvals.RData      p (vector of m*N p-values), mask, copes
# =============================================================================

library(httr)
library(jsonlite)
library(RNifti)

collection_id <- 2447
contrast      <- "noun/verb judgement > baseline"
data_dir      <- file.path("data", paste0("neurovault_", collection_id))
out_dir       <- "output"
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, showWarnings = FALSE)


# -----------------------------------------------------------------------------
# 1. Download subject-level t-maps from NeuroVault
# -----------------------------------------------------------------------------

api_url <- paste0("https://neurovault.org/api/collections/", collection_id, "/images/")
res  <- GET(api_url)
meta <- fromJSON(content(res, "text", encoding = "UTF-8"))

sel <- meta$results$contrast_definition == contrast
sel[is.na(sel)] <- FALSE
image_urls <- meta$results$file[sel]

n    <- length(image_urls)                    # number of subjects (N = 10)
sids <- sprintf("%02d", seq_len(n))

for (i in seq_len(n)) {
  fname <- file.path(data_dir, paste0("tmap_", sids[i], ".nii.gz"))
  if (!file.exists(fname)) download.file(image_urls[i], destfile = fname, mode = "wb")
}


# -----------------------------------------------------------------------------
# 2. Build the analysis mask
# -----------------------------------------------------------------------------
# No mask is provided with the collection, so the mask is the union of voxels
# that are non-zero in at least one subject, excluding voxels with infinite
# values in any subject.

any_nonzero  <- NULL
any_infinite <- NULL

for (i in seq_len(n)) {
  nii <- readNifti(file.path(data_dir, paste0("tmap_", sids[i], ".nii.gz")))
  arr <- as.array(nii)

  if (is.null(any_nonzero)) {
    any_nonzero  <- arr != 0
    any_infinite <- is.infinite(arr)
  } else {
    any_nonzero  <- any_nonzero  | (arr != 0)
    any_infinite <- any_infinite | is.infinite(arr)
  }
}

mask_arr <- any_nonzero & !any_infinite
writeNifti(asNifti(mask_arr * 1, reference = nii), file.path(data_dir, "mask.nii.gz"))

mask   <- readNifti(file.path(data_dir, "mask.nii.gz"))
mask_v <- as.vector(mask)
m      <- sum(mask_v)                         # number of voxels in the mask


# -----------------------------------------------------------------------------
# 3. Voxel-wise p-values
# -----------------------------------------------------------------------------
# t-statistics are stored in an N x m matrix (rows = subjects, cols = voxels).
# Flattening column-wise gives a vector where x[(1:n) + (j-1)*n] are the N
# subjects for voxel j, which is the input format expected by hommel.
# Degrees of freedom are not provided, so two-sided p-values use the
# standard normal approximation.

copes <- matrix(0, nrow = n, ncol = m)
for (i in seq_len(n)) {
  fname <- file.path(data_dir, paste0("tmap_", sids[i], ".nii.gz"))
  copes[i, ] <- as.vector(readNifti(fname))[mask_v == 1]
}

x <- as.vector(copes)
p <- 2 * pnorm(abs(x), lower.tail = FALSE)

save(p, mask, copes, file = file.path(out_dir, "pvals.RData"))


# -----------------------------------------------------------------------------
# 4. Subject-level p-value maps (for visualisation)
# -----------------------------------------------------------------------------

for (i in seq_len(n)) {
  tmap <- readNifti(file.path(data_dir, paste0("tmap_", sids[i], ".nii.gz")))
  parr <- 2 * pnorm(abs(as.array(tmap)), lower.tail = FALSE)
  parr[!mask_arr] <- NA
  writeNifti(asNifti(parr, reference = tmap),
             file.path(data_dir, paste0("pmap_", sids[i], ".nii.gz")))
}
