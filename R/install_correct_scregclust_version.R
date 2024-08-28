repo   <- "scmethods/scregclust"
#
# #commit with outputs of expected format as of 8JUN
# commit <- "e0516fc"#dd89a6a549d7906f6502f41586c3ed47f

# commit with initial variance estimate that works in n < p case
#commit <- "e490910"#a1398294c1164c5ddd0965c7b9b42e881

# commit with modules terminology and some updates
commit <- "62dcd25"#42d5784c3760ed09e0aa630a60d5b626b



# important: build rtools manually.
# https://cran.r-project.org/bin/windows/Rtools/rtools44/rtools.html

# pkgbuild::has_build_tools()

# Install the package from the specific commit
devtools::install_github(paste0(repo, "@", commit), force = TRUE)


# note:  if updating, might need to remove (and above code will reinstall) screclust, cli, rcpp
