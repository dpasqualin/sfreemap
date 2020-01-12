# SFREEMAP

SFREEMAP is an analytical approach to obtain accurate, per-branch
expectations of numbers of state transitions and dwelling times. It also
provides an intuitive way of visualizing the results by integrating over
the posterior and summarizing the parameters onto a target reference
topology (such as a consensus or MAP tree) provided by the user.

## Installation

You need to have R installed on your system, plus some dependencies.
If you are using a Debian/Ubuntu based distribution, just type the following command in a terminal.

`sudo apt-get install r-base-core libblas-dev liblapack-dev libmagick++-dev`

You can install the stable version of the package from CRAN by opening R console and typing:
```
install.packages(sfreemap)
```

The most up to date version from github can be installed using devtools.
If you want to build vignettes (usage examples), provide `build_vignettes = FALSE`
parameter to `install_github` command (installation will take a few minutes more):
```
install.packages('devtools')
require(devtools)
install_github('dpasqualin/sfreemap')
```

You can also clone the repository and install it like this:
```
git clone https://github.com/dpasqualin/sfreemap.git
R CMD build sfreemap && R CMD INSTALL sfreemap
```

If you want to look at the reference manual, load R and type the following
commands. A new tab will open in your web browser with the description of
all available functions.

```
require(sfreemap)
help(package = "sfreemap", help_type = "html")
```

The following command will open a new tab in your web browser with a link to
the vignettes (HTML).

```
browseVignettes("sfreemap")
```


## Development

For development you are gonna need a few extra packages:
```
# On Fedora libssl-dev might be called libssh2-devel, and on OS X might be just libssh2
sudo apt-get install texlive-full zlib1g-dev libcurl4-openssl-dev libssl-dev
```