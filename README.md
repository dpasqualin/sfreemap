### WARNING

This package is not production ready. There are still some improvements and tests that need to be done.

### Package Requirements for Production

You need to have R installed on your system, plus some dependencies. If you are using a debian/ubuntu based distribution, just type the following command in a terminal.

`sudo apt-get install r-base-core libblas-dev liblapack-dev` 

### Package Requirements for Development

`sudo apt-get install r-base-core texlive-full zlib1g-dev libcurl4-openssl-dev`

### Install

```
install.packages('devtools')
install_github('dpasqualin/sfreemap')
```

If you have troubles installing the `devtools` package, try downloading
`sfreemap` and then building and installing it using the following commands:

```
git clone https://github.com/dpasqualin/sfreemap.git
R CMD build sfreemap && R CMD INSTALL sfreemap
```

If you choose to install using the command above, the documentation will be
available in the directory `sfreemap.Rcheck`.

Some people might have problems with package `Briostrings` as well, which is
a dependency of `phangorn`, which is a dependency of `sfreemap`. If you do,
the official Biostrings website suggest the following commands to install
it:

```
## try http:// if https:// URLs are not supported
source('https://bioconductor.org/biocLite.R')
biocLite('Biostrings')
```
