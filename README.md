# EpiStat

EpiStat is perl pipeline for searching epistatic interactions on phylogenetic trees. The input data is a phylogenetic tree with a ist of substitutions on branches in XPARR format:

`Header`\
`RootID`\
`BranchLine+`

`Header="child\tparent\tlength"`

`RootID="\w+"` - any string containing letters or digits\
`BranchLine = "NodeID\tParentNodeID\tBranchLength\tSynSiteList\tNonSynMutList\tNonSynMutList"`\
`site=\d+` - site in the reference sequence, positive number\
`SynSiteList = "site(;site)*"` - a list of sites, in which synonymous substitutions has occured on this edge\
`NonSynMutList = "mutation(;mutation)*"` - alist of amino acid substitutions on this edge\
`mutation = "AsiteB"`\
`A` and `B` is one letter amino acid code

You can find an example of XPARR tree in /epistat.7.1/prm_templates/intragene/nophen/

XPARR formatted trees can be generated from the output of pylognetic reconstruction program which can reconstruct both the tree and the ancestor states. The *mk_xparr.pl* script can be used to generate XPARR from tree and alignment:\
`mk_xparr.pl -t <tree_NWK> [options] <alignment_FASTA>`\
Here <tree_NWK> is a newick formated tree, where the bootstrap-slot is filled with internal nodes' identifiers, and <alignment_FASTA> is an FASTA alignment of both terminal and ancestor sequences.

## EpiStat workflow

EpiStat is implemented in perl 5 and use R libraries. The workflow itself consists of scenarios whic shuold be executed subsequently:

*estimate_tau.pl* - Estimates the scaling parameter *tau*, which reflects how fast or slow do the sequential substitutions happen one after another for an average pair of sites.\
*run_epistat.pl - Runs the main analysis.\
*[estimate_fdr.pl]* - Estimates FDR values for nominal p-values. This script is not required for obtaining of final result, but it is necessary to evaluate the presence or absence of signal in the data.\
*mk_coevolmtx.pl* - Forms a covariance matrix of substitutions in pairs of sites\
*minvert/cor2pcor.R* - Inverts the covariance matrix, normalizes the result.

For each script, the description of options and values of input parameters passed in files is described in the comments at the beginning of \*.pl files. The end result generated with *minvert/cor2pcor.R* is a matrix containing the *wi_score* statistics that is used to order the pairs. Larger statistic values correspond to pairs of sites with a higher probability of being close on the 3D structure of the protein. When using scripts from the package, the parameter file naming convention is adopted (optional). Let the source file be named example.xparr, then the parameter files are named: *example.<script_name>.prm*. The parameted values may depend on the data being processed.\
Examples of parameter files that can be used as templates can be found in: /epistat.7.1/prm_templates/intragene/nophen/.

## Installation

- Install perl 5 (you can use https://perlbrew.pl/ if you like)
- Install R (https://www.r-project.org/)
- Install *cpanminus* and use it to install a couple of Perl libraries from CPAN:
  - `cpan App::cpanminus`
  - `cpanm Bio::Phylo`
  - `cpanm Time::Progress`
  - `cpanm List::BinarySearch`
  - `cpanm Math::Gradient`
  - `cpanm Color::Rgb`
- Install R packages:
  - BiRewire package for R:

    To install this package, start R (version "4.0") and enter:
    ```
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("BiRewire")
    ```
    For older versions of R, please refer to the appropriate Bioconductor release.

  - `BiocManager::install("corpcor")`
  - `BiocManager::install("getopt")`

- Clone the repository or download and unzip the release archive. The epistat.7.1 contains the scripts, while three other directories contain libraries.
- Create two system variables: *EPISTAT_HOME* pointing to epistat.7.1 directory and *EPISTAT_LIB* pointing to the directory with pdbUtils, TreeUtils and DnaUtilities. On linux you can add `export EPISTAT_HOME=your_epistat_home` and `export EPISTAT_LIB=your_epistat_lib` to the .bash_profile, for example:
  ```
  export EPISTAT_HOME=$HOME/epistat.7.1/
  export EPISTAT_LIB=$HOME/epistat_lib/
  ```
- On linux: make all perl scripts executable (chmod)

## Quick start

- Copy the example files to a separate directory:
```
  mkdir example
  cp ~$HOME/epistat.7.1/prm_templates/intragene/nophen/* ./example
  cd example
```
- Run these commands from the shell:
  - `~$HOME/epistat.7.1/estimate_tau.pl -b 0.95 example.xparr > example.tau`\
  This script estimates 'tau' parameter and prints it to example.tau. Copy the value of "TAU=" from example.tau to the field "TAU=" in example.epistat.prm.
  - `~$HOME/epistat.7.1/run_epistat.pl -x example.xparr -m 3 -p 10000 example.epistat.prm run_epistat.prm`
  - `~$HOME/epistat.7.1/estimate_fdr.pl example.upper.pvalue.sites.fdr.prm > example.upper.pvalue.sites.fdr`
  

## Advanced configaration

### estimate_fdr.pl

`estimate_fdr.pl example.upper.pvalue.(sites|pairs).fdr.prm> example.upper.pvalue.(sites|pairs).fdr`

The parameter files *example.upper.pvalue.(sites|pairs).fdr.prm* requires the appropriate names of the output files generated by run_epistat.pl. All output files have the "example" prefix to replace the prefixes of the corresponding filenames in the template file. The script allows you to apply filters to pairs of sites when calculating the FDR, which allows you to detect a signal in a subset of sites with certain properties. For some types of data, it may be required for signal detection to limit the minimum number of substitutions in each pair, which are specified in the parameter file:
```
MutationNumbersFilter = "sites"
BGR_SiteMinMutations = "2"
FGR_SiteMinMutations = "2"

MutationNumbersFilter="pairs"
BGR_SiteMinMutations="2"
FGR_SiteMinMutations="2"
```
The "sites" filter requires at least BGR_SiteMinMutations substitution and FGR_SiteMinMutations substitution in the first (background) and the second (foreground) sites in the pair. The "pairs" filter is more restrictive, it requires that consecutive pairs of mutations in pairs of sites contain at least BGR_SiteMinMutations and FGR_SiteMinMutations of *different* mutations. By default, you should use the filter "pairs" and 2 mutations in sites, if N/L >= 5, where N is the number of sequences in the alignment, L is the length of the alignment. Disabling the filter for the minimum number of pairs is necessary if the initial alignment contains few sequences (N / L <5), for this you should set:
```
MutationNumbersFilter="sites"
BGR_SiteMinMutations="1"
FGR_SiteMinMutations="1"
```
