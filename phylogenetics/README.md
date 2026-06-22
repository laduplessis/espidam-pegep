This is a short tutorial to explore some R functions to build and display phylogenetic trees. Over the course of the tutorial you'll build neighbour-joining, maximum parsimony and maximum-likelihood trees. 

# Required packages

The functions we'll be using are part of the `ape` and `phangorn` packages. If you haven't installed them install them now.

```r
install.packages("ape")
install.packages("phangorn")
```

Now load the packages!

```r
library(ape)
library(phangorn)
```
# Data

For this tutorial we'll be using an HIV dataset that was used in a criminal case. The dataset comes from [Metzker et al. (2002) PNAS](http://www.pnas.org/content/99/22/14292). From the abstract of the paper:

> _A gastroenterologist was convicted of attempted second-degree 
> murder by injecting his former girlfriend with blood or blood-products obtained from an HIV type 1 (HIV-1)-infected patient under
> his care. Phylogenetic analyses of HIV-1 sequences were admitted and
> used as evidence in this case, representing the first use of phylogenetic analyses in a criminal court case in the United States._

The dataset consists of two alignments, from the reverse transcriptase (pol) and from the envelope (env). Each alignment contains multiple sequences from the victim, the patient and the community. 

In the original court case (which dates from the late 1990s) maximum parsimony and least evolution trees were used (this is a different type of distance based tree). The publication, which came out a few years later additionally used Bayesian phylogenetics. 

---
# Pol
## Loading the data

The first step is to get the data into R. You can use the `read.dna()` or `read.FASTA()` functions in `ape`. If you're using `read.dna()` remember to supply an argument for the file format! Alternatively, you can use the `read.phyDat()` function in `phangorn`.

**Hint:** You can always use `?` to get quick help for a function in R, e.g. type in `?read.dna` to display the help for it. 

At first we'll just load the pol alignment (`hiv_pol.fasta`). Explore the alignment and find out how many sequences and sites it contains. Note that the sequence ids are of format  `<Genbank id>|<Patient id>`. The Genbank id can be used to find the sequence records on NCBI Genbank. The patient ids start with either V (sequences from the victim), P (sequences from the HIV positive patient), or LA (sequences from the community). 

**Before doing any phylogenetics, what is your hypothesis for what we'll observe if the victim really was injected with blood from the HIV positive patient?** 

You can extract the sequence ids with the `labels()` function. We'll need to find out which sequences are from the victim, the patient or the community, so that we can annotate them on the tree later. Use the `strsplit()` and  `grep()` functions to do that. You may also need to use a for loop or `sapply()` to get there. (There are probably many other ways to do that, and you should use whatever way is easiest for you to implement and understand). 

---
## Distance-based trees 

Use the `nj()` function in `phangorn` to build a neighbour-joining tree using the TN93 (Tamura Nei 1993) model. This model is similar to an HKY model, but allows the two transition rates to be different. 

Before you can do that you first need to construct a pairwise distance matrix. The easiest way is using the `ape` function `dist.dna()`. (There is also the `dist.ml()` function in `phangorn`, but it has a very limited range of nucleotide substitution models). 

**Does neighbour-joining produce a rooted tree?** Use the `is.rooted()` function to check if the tree is rooted or not. 

We can plot the tree using `plot.phylo()` (or just `plot()`). If the tree is unrooted we have to also supply an argument  `type = "unrooted"` otherwise an arbitrary root will be used. Unrooted trees are also more legible when using the `lab4ut = "axial"` argument. In general you may also want to scale down the labels to make it more readable. 

Since the neighbour-joining does in fact not produce a rooted tree it has to be rooted. We will root the tree on the sequence of the community patient `LA02`. Which community patient we root on is not especially important for deciding whether or not the doctor is guilty - it is simply important to use _any_ outgroup. 

To root the tree use the function `root.phylo()`. Now plot the rooted tree. You can colour the labels if you supply a vector of colours to the `color` argument in `plot()`. To do that you'll need to match the tip labels of the tree to the patient, victim, and community sequence names you extracted earlier. The `match()` function could be useful. Tip labels of the tree can be extracted using, e.g. `poltree$tip.label` if your tree is called `poltree`. 

**What are the taxonomic relationships of the patient and victim sequences?Do they form monopyletic, paraphyletic, or polyphyletic groups? Does the phylogeny look like your hypothesis?** 

Finally, let's do the same using a Jukes-Cantor model and see how the two trees compare to each other. You can plot a tanglegram connecting the two trees using the `cophyloplot()` function. It'll take some adjusting to get it to look good - in particular you'll need to increase the `space` argument quite a bit or else you won't be able to see the "tangle". 

**Do both models lead to the same conclusion? How different are the trees?** 

_**(optional)**_
Also use the `upgma()` function to look at a UPGMA tree. **Which tree do you trust more between the NJ and UPGMA trees?** 

---
## Maximum parsimony

To use the parsimony methods implemented in the `phangorn` package we have to convert our alignment to a `phyDat` object using the `as.phyDat()` function. Alternatively we could have just used `read.phyDat()` to read it in as a `phyDat` object. 

We can use the `parsimony()` function to calculate the parsimony score of a tree using Fitch's algoritm. **What is the parsimony score of the neighbour joining tree?** 

To actually find the maximum parsimony tree we can use the function `pratchet()` (parsimony ratchet). This function returns the topology of the most parsimonious tree, but it doesn't have branch lengths, i.e. it only does the post-order traversal. (plot it!). To add branch-lengths we also have to do the pre-order traversal, which we can do using the `acctran()` function. 

**Does the maximum parsimony tree support the same conclusion as the neighbour-joining trees?** 

Now plot the tanglegram between the TN93 neighbour-joining tree and the maximum parsimony tree. Note that while the maximum parsimony tree is rooted, it's not rooted on the same sequence (the root seems arbitrary), so you'll need to first reroot it to the same outgroup as the neighbour-joining tree. 

---
## Maximum-likelihood

We can use the `pml_bb()` function to find the maximum-likelihood tree. Use a GTR substitution model. 

Once again, the tree is unrooted. Root it and plot it, and plot the tanglegram with the TN93 neighbour-joining tree. 

**Does the maximum-likelihood tree lead to any different conclusions?** 

---
# Env 

Now repeat all of the above with the env alignment. The alignment is bigger, so it will take considerably longer for the maximum parsimony and maximum-likelihood trees! 

**Does env support the same conclusions as pol?** 

---
# Extra stuff

## Prettier trees
The tree plotting functions in `ape` are quite bare-bones and it requires a lot of effort to make them output pretty trees. `ggtree` uses the same grammar for graphics as `ggplot2` to plot pretty trees, but the learning curve is quite steep if you're not used to `ggplot2`. Fortunately it's quite widely used, and it's easy to get help online. You can also look at the `ggplot` book by Guangchang Yu: https://yulab-smu.top/treedata-book/index.html, especially chapters 4 and 5 (it's quite short). 

## IQ-TREE
IQ-TREE is one of the most widely used programs for building maximum-likelihood trees. It also automatically uses ModelFinder to find the best-fitting nucleotide substitution model. Does it lead to the same conclusions? You can build the maximum-likelihood tree using default settings by simply typing in `iqtree3 -s pol.fasta` (if IQ-TREE version 3 is installed). 