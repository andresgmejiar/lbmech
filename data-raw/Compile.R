# Make sure all function documentation files have been compiled
roxygen2::roxygenize()

# Build the PDF manual
devtools::build_manual()

# Build the README
devtools::build_readme()

# Run check to make sure that everything is compatible with systems and
# all of the tests within examples are running OK
devtools::check()

# Build the website
pkgdown::build_site()
