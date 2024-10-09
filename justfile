import 'files/overleaf.just'

default:
   just --list

ensure-env: clone-overleaf
   pixi install --frozen

publish:
	nbqa isort nbs/*.ipynb
	nbqa black nbs/*.ipynb
	quarto publish gh-pages --no-prompt

copy-figure:
	#only copy the figures that are used in the manuscript
	cp -r figures/fig* overleaf/draftfigures/

publish-poster:
  Rscript -e 'pagedown::chrome_print("manuscripts/.AGU23_poster.rmd")'

publish-qrcode:
  segno "https://beforerr.github.io/EMIC_multi-satellites_study/" -o=images/qrcode.png --light transparent --scale 10

check_os:
	@if [ "$(shell uname)" != "Darwin" ]; then echo "This makefile is intended for MacOS only." && exit 1; fi

idl:
	idl -e ".run ./src/elf_diffusion_coefficient.pro"

poes:
	cd src && idl -e ".run ./poes_process.pro"

# As long as numpy and C/Fortran compilers are present everything should build and spacepy.irbempy should be available.
install_spacepy: check_os install_cdf
	@echo "Installing spacepy..."
	@SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk pip install spacepy
	@echo "spacepy installed successfully."

install_cdf:
	@if [ ! -f CDF3_9_0-binary-signed.pkg ]; then \
		echo "CDF3_9_0-binary-signed.pkg not found. Downloading..."; \
		curl -O https://spdf.gsfc.nasa.gov/pub/software/cdf/dist/latest/macosx/CDF3_9_0-binary-signed.pkg; \
	fi
	@echo "Installing CDF3_9_0-binary-signed.pkg"
	# Must be run as root to install this package.
	@sudo installer -pkg CDF3_9_0-binary-signed.pkg -target /
	@echo "CDF V3.9.0 installed successfully."