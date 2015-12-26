### WARNING

This package is not production ready. There are still some improvements and tests that need to be done.

### Package Requirements for Production

You need to have R installed on your system. If you are using a debian/ubuntu based distribution, just type the following command in a terminal.

`sudo apt-get install r-base-core`

### Package Requirements for Development

`sudo apt-get install r-base-core texlive-full`

###  R Dependencies

### Install

Before installing *sfreemap* make sure you have the package *phytools*:

```
install.packages('phytools')
install.packages('devtools')
install_github('dpasqualin/sfreemap')
```

If you have troubles installing the `devtools` package, try downloading
`sfreemap` and then building and installing it using the following commands:

```
git clone https://github.com/dpasqualin/sfreemap.git
R CMD check sfreemap && R CMD build sfreemap && R CMD INSTALL sfreemap
```

If you choose to install using the command above, the documentation will be
available in the directory `sfreemap.Rcheck`.
