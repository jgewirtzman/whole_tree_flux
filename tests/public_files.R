#!/usr/bin/env Rscript
# Check the Git index (including staged changes) before sharing: only project code, data, results,
# figures and documentation are allowed; no Word or Quarto files, Office lock files or Mac metadata.
inventory <- tempfile()
status <- system2("git", c("ls-files", "-z"), stdout = inventory)
if (status != 0) stop("Run this check from the Git repository.")
bytes <- readBin(inventory, "raw", n = file.info(inventory)$size)
unlink(inventory)
paths <- strsplit(rawToChar(replace(bytes, bytes == as.raw(0), charToRaw("\n"))), "\n", fixed = TRUE)[[1]]
allowed <- "^(\\.gitignore|\\.githooks/|README\\.md$|LICENSE$|CITATION\\.cff$|\\.zenodo\\.json$|run_analysis\\.R$|renv\\.lock$|whole_tree_flux\\.Rproj$|data/|flux/|scaling/|figures/|docs/|scripts/|tests/)"
blocked <- "(\\.(docx|qmd|doc|pptx)$|(^|/)~\\$|(^|/)\\.DS_Store$|manuscript|(^|/)notes/|(^|/)archive/|(^|/)flux/cache/)"
bad <- paths[!grepl(allowed, paths) | grepl(blocked, paths, ignore.case = TRUE)]
if (length(bad)) stop("Files outside the shared project layout are staged:\n", paste(bad, collapse = "\n"))
cat("PASS: only project code, data, results, figures and documentation are in the Git index\n")
