# SFREEMAP

SFREEMAP is an analytical approach to obtain accurate, per-branch
expectations of numbers of state transitions and dwelling times. It also
provides an intuitive way of visualizing the results by integrating over
the posterior and summarizing the parameters onto a target reference
topology (such as a consensus or MAP tree) provided by the user.

Paper: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1554-7

## Installation

First install R and dependencies according to your operating system (see below), then install
the package from CRAN, github or source, as described in the following sections.


##### Debian/Ubuntu

If you are using a Debian/Ubuntu based distribution, install the following packages.

```bash
sudo apt-get install r-base-core libblas-dev liblapack-dev libmagick++-dev
```

##### Windows

Install the latest version of Rtools for Windows following
http://cran.r-project.org/bin/windows/Rtools/.

##### MAC OS

You might need to install install Xcode command line tools:
```bash
xcode-select --install
```

#### From Github using devtools

The most up to date version from Github can be installed using `devtools`. Devtools itself has a different set of
dependencies depending on the OS you are on, for some `install.packages('devtools')` on R console might be enough, 
for Debian/Ubuntu the easiest way is with `sudo apt install r-cran-devtools`.

After installing `devtools` you can run the lines below on R to install `sfreemap`.
Use `build_vignettes = TRUE` parameter to `install_github` if you want to build usage examples as well.
```R
require(devtools)
install_github('dpasqualin/sfreemap')
```

#### From source

If you don't want to install `devtools` you can also clone the repository and install it using R command line.
You will need `git` installed on your system.
```
git clone https://github.com/dpasqualin/sfreemap.git
R CMD build sfreemap && R CMD INSTALL sfreemap
```

## Documentation

If you want to look at the reference manual, load R and type the following
commands. A new tab will open in your web browser with the description of
all available functions.

```
require(sfreemap)
help(package = "sfreemap", help_type = "html")
```

The following command will open a new tab in your browser with a link to the vignettes, if you chose
to install it.

```
browseVignettes("sfreemap")
```

## Development

For development you are gonna need a few extra packages:
```
# On Fedora libssl-dev might be called libssh2-devel, and on OS X might be just libssh2
sudo apt-get install texlive-full zlib1g-dev libcurl4-openssl-dev libssl-dev
```
