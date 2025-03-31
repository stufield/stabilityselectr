# h/t to @jimhester and @yihui for this parse block:
# https://github.com/yihui/knitr/blob/dc5ead7bcfc0ebd2789fe99c527c7d91afb3de4a/Makefile#L1-L4
# Note the portability change as suggested in the manual:
# https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Writing-portable-packages
PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename `pwd`)
RM = rm -rf
RCMD = R --vanilla CMD
RSCRIPT = Rscript --vanilla


all: check clean
roxygen: docs

docs:
	@ $(RSCRIPT) -e "roxygen2::roxygenise(roclets = c('collate', 'namespace', 'rd'))"

readme:
	@ echo "Rendering README.Rmd"
	@ $(RSCRIPT) \
	-e "Sys.setenv(RSTUDIO_PANDOC='/usr/bin/pandoc/')" \
	-e "options(cli.width = 80L)" \
	-e "rmarkdown::render('README.Rmd', quiet = TRUE)"
	@ $(RM) README.html

test:
	@ $(RSCRIPT) \
	-e "Sys.setenv(TZ = 'America/Denver')" \
	-e "devtools::test(reporter = 'summary', stop_on_failure = TRUE)"

test_file:
	@ $(RSCRIPT) \
	-e "Sys.setenv(TZ = 'America/Denver', NOT_CRAN = 'true')" \
	-e "devtools::load_all()" \
	-e "testthat::test_file('$(FILE)', reporter = 'progress', stop_on_failure = TRUE)"

accept_snapshots:
	@ Rscript -e "testthat::snapshot_accept()"

build: docs
	@ cd ..;\
	$(RCMD) build --resave-data $(PKGSRC)

check: build
	@ cd ..;\
	$(RCMD) check --no-manual $(PKGNAME)_$(PKGVERS).tar.gz

# saving data/cluster.rda from progeny clustering
# paper data set example: `clust_data` & `progeny_data`
# Comes from the Progeny Clustering paper: pg 3.
objects:
	@ echo "Creating 'data/cluster.rda' ..."
	@ $(RSCRIPT) \
	-e "clust_data <- data.frame(" \
	-e "  F1 = c(0.4, 0.98, 0.35, 0.62, 0.48, 0.63, 0.22," \
	-e "         0.92, 0.38, 0.05, 0.91, 0.57, 0.77, 0.41," \
	-e "         0.57, 0.02, 0.19, 0.63, 0.95, 0.16)," \
	-e "  F2 = c(0.85, 0.24, 0.98, 0.08, 0.79, 0.38, 0.88," \
	-e "         0.14, 0.75, 0.56, 0.33, 0.32, 0.22, 0.75," \
	-e "         0.35, 0.61, 0.63, 0.34, 0.09, 0.67)" \
	-e ")" \
	-e "set.seed(101)" \
	-e "Sigma        <- matrix(c(1, 0.25, 0.25, 1), 2, 2)" \
	-e "progeny_data <- rbind(MASS::mvrnorm(50, mu = c(-1, 2), Sigma)," \
	-e "                      MASS::mvrnorm(50, mu = c(2, 0), Sigma)," \
	-e "                      MASS::mvrnorm(50, mu = c(-1, -2), Sigma)" \
	-e ")" \
	-e "colnames(progeny_data) <- c('F1', 'F2')" \
	-e "save(clust_data, progeny_data, file = 'data/cluster.rda', compress = 'xz')"
	@ echo "Saving 'data/cluster.rda' ..."
	@ echo "Creating 'inst/prgeny-plots/*.pdf' ..."
	@ $(RSCRIPT) inst/progeny-plots/create-plots.R
	@ echo "Saving 'inst/prgeny-plots/*.pdf' ..."

install:
	@ R CMD INSTALL --use-vanilla --preclean --resave-data .

increment:
	@ echo "Adding Version '$(ver)' to DESCRIPTION"
	@ $(shell sed -i 's/^Version: .*/Version: $(ver)/' DESCRIPTION)
	@ echo "Adding new heading to 'NEWS.md'"
	@ $(shell sed -i '1s/^/# $(PKGNAME) $(ver)\n\n/' NEWS.md)

release:
	@ echo "Adding release commit"
	@ git add -u
	@ git commit -m "Increment version number"
	@ git push origin main
	@ git tag -a v$(PKGVERS) -m "Release of $(PKGVERS)"
	@ git push origin v$(PKGVERS)
	@ echo "Remember to bump the DESCRIPTION file with bump_to_dev()"

clean:
	@ cd ..;\
	$(RM) $(PKGNAME)_$(PKGVERS).tar.gz $(PKGNAME).Rcheck
