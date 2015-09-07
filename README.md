### WARNING

This package is not production ready. There are still some improvements and tests that need to be done.

### Package Requirements for Production

You need to have R installed on your system. If you are using a debian/ubuntu based distribution, just type the following command in a terminal.

`sudo apt-get install r-base-core`

### Package Requirements for Development

`sudo apt-get install r-base-core texlive-full`

###  R Dependencies

### Install

Before installing *sfreemapc* make sure you have the package *phytools*:

```
install.packages('phytools')
install.packages('devtools')
install_github('dpasqualin/sfreemapc')
```

If you have troubles installing the `devtools` package, try downloading
`sfreemapc` and then building and installing it using the following commands:

```
git clone https://github.com/dpasqualin/sfreemapc.git
R CMD check sfreemapc && R CMD build sfreemapc && R CMD INSTALL sfreemapc
```

If you choose to install using the command above, the documentation will be
available in the directory `sfreemapc.Rcheck`.


### Example

```
require(sfreemapc) # load package
tree <- pbtree(n=100,scale=1) # create a tree with 100 taxa
Q <- matrix(c(-1,1,1,-1),2,2) # create a transition rate matrix
rownames(Q)<-colnames(Q)<-letters[1:nrow(Q)] # give name to the states
tree <- sim.history(tree,Q,anc="a") # simulate a history for the tree

# estimate the history
sm <- sfreemap.map(tree, tree$states, Q='empirical')
```
