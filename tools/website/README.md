# Website maintenance

The package website uses pkgdown with the Bootstrap 5 Cosmo theme. Install
pkgdown and roxygen2 in your usual shared R user library. Use the shared R
installation to install the current source package before building the site.

User-facing tutorial prose and code live in `inst/examples/getting-started.R`.
`tools/build-getting-started.R` runs that script and regenerates its Markdown
and figure. Function documentation is maintained in the roxygen comments in
`R/`; regenerate the corresponding Rd files with `roxygen2::roxygenise()`.
README.md is maintained directly.

From the repository root, run `Rscript tools/build-site.R` to build a clean
site in `docs/`. This also updates the guide and copies the downloadable R
scripts. After regenerating the guide, reinstall the source package before
running `Rscript tests/examples.R`, which reads the installed example files.
Confirm that the installed scripts match the workspace; a matching version
number alone does not establish that they are current.

GitHub Actions validates the examples and builds the site on pushes and pull
requests. Pushes to `master` deploy to GitHub Pages; pull requests do not publish.

Keep maintainer notes and planning documents in this directory. pkgdown can
publish Markdown files from the repository root and `.github`, even when they
are excluded from the R package by `.Rbuildignore`.
