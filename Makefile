# ~
#
# Makefile for managing development/deployment
#
# ------------------------------------------------------


# config
# -----
PROJECT     = `grep 'name *= *' setup.cfg | sed 's/.*=//g' | tr -d '[:space:]'`
VERSION     = `git show -s --format="%cd-%h" --date=format:'%Y.%m.%d'`


# targets
# -------
.PHONY: help docs info clean init build

help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'


info: ## list info about package
	@echo "$(PROJECT), version $(VERSION)"
	@echo last updated: `git log | grep --color=never 'Date:' | head -1 | sed 's/Date:   //g'`


clean: ## clean unnecessary files from repository
	rm -rf build dist
	rm -f ubuntu*.log
	rm -rf tests/__pycache__ .pytest_cache .ipynb_*
	rm -f roles/*.retry


init: ## intialize repository for development
	git submodule init
	git submodule update


update: ## update submodules with most recent versions
	git submodule update --remote
	git add packages
	git commit -m "Updated submodules to most recent remote ref."


requirements: init ## install all requirements for using repository
	python -m pip install ansible
	python -m pip install -r requirements.txt

test: ## run testing suite for module using py.test
	py.test

