#!/usr/bin/env Rscript
# Audit the current Git index, including staged additions and removals.
# Read a NUL-delimited inventory without locale-dependent filename quoting.
inventory <- tempfile()
status <- system2("git", c("ls-files", "-z"), stdout = inventory)
if (status != 0) stop("Run this check from the public Git repository.")
bytes <- readBin(inventory, "raw", n = file.info(inventory)$size)
unlink(inventory)
paths <- strsplit(rawToChar(replace(bytes, bytes == as.raw(0), charToRaw("\n"))), "\n", fixed = TRUE)[[1]]
private <- paste(c("(^|/)(manuscript[^/]*|DRAFT_[^/]*|WORKLOG_[^/]*|STATUS_TRACKER[^/]*|CLAUDE[^/]*|comment_responses[^/]*|prompt_for_[^/]*)(/|$)",
 "(^|/)(synthesis_merge|grl_ge2m|\\.git-manuscript|\\.claude|\\.codex)(/|$)",
 "(^|/)(build_manuscript\\.sh|ms-git\\.sh|comment_response_stats\\.R)$", "\\.(docx|qmd)$", "_writeup\\.", "(^|/)\\.DS_Store$"), collapse = "|")
bad <- paths[grepl(private, paths, ignore.case = TRUE)]
if (length(bad)) stop("Private/editorial files in public index:\n", paste(bad, collapse = "\n"))
cat("PASS: no manuscript, editorial, internal-planning or Mac metadata files in the public index\n")
