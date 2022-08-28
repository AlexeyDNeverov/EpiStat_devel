# EpiStat

EpiStat is a pipeline for searching epistatic interactions on phylogenetic trees. 

The method is based on the following publications:

1. Kryazhimsky S., et. al. Prevalence of Epistasis in the Evolution of Influenza A Surface Proteins. 2011 PLoS Genetics 7(2):e1001301
2. Neverov A.D., et al. Coordinated Evolution of Influenza A Surface Proteins. 2015 PLoS Genetics 11(8):e1005404
3. Neverov A.D., et al. Episodic evolution of coadapted sets of amino acid sites in mitochondrial proteins. 2021 PLoS Genetics 17(1):e1008711.


The input data is a phylogenetic tree with a list of substitutions on branches in XPARR format:

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
Here <tree_NWK> is a newick formated tree, where the bootstrap-slot is filled with internal nodes' identifiers, and <alignment_FASTA> is a FASTA alignment of both terminal and ancestor sequences.

## EpiStat workflow

EpiStat is implemented in perl 5, R and python 2. The workflow itself consists of scenarios which should be executed subsequently:

- *estimate_tau.pl* - Estimates the scaling parameter *tau*, which reflects how fast or slow do the sequential substitutions happen one after another for an average pair of sites.
- *run_epistat.pl* - Runs the main analysis.
- *estimate_fdr.pl* - Estimates FDR values for nominal p-values. This script is not required for obtaining of final result, but it is necessary to evaluate the presence or absence of signal in the data.
- *mk_coevolmtx.pl* - Forms a covariance matrix of substitutions in pairs of sites
- *minvert/cor2pcor.R* - Inverts the covariance matrix, normalizes the result.
- *graph/mk_graph.pl* - Cerates coevolution graphs in GRAPHML format based on lists of weighted edges.
- *graph/find_partition.py* - Searches for clusters of graph vertices with a high density of positive connections within and negative connections between clusters.
- *stat/xparr_site_groups_test.pl* - Calculates on which edges of the tree the relative evolutionary rates have changed in the given site groups.
- *stat/summ_site_groups_tests.pl* - Determine the coordinated changes in the rates of substitutions in different genes.

For each script, the description of options and values of input parameters passed in files is described in the comments at the beginning of \*.pl files. The end result generated with *minvert/cor2pcor.R* is a matrix containing the *wi_score* statistics that is used to order the pairs. Larger statistic values correspond to pairs of sites with a higher probability of being close on the 3D structure of the protein. When using scripts from the package, the parameter file naming convention is adopted (optional). Let the source file be named example.xparr, then the parameter files are named: *example.<script_name>.prm*. The parameted values may depend on the data being processed.\
Examples of parameter files that can be used as templates can be found in: /epistat.7.1/prm_templates/intragene/nophen/.

## Installation

- Install perl 5 (you can use [perlbrew](https://perlbrew.pl/) if you like)
- Install [R](https://www.r-project.org/)
- Install python 2 (typically system python for linux, you can also use [pyenv](https://github.com/pyenv/pyenv) with [miniconda](https://docs.conda.io/en/latest/miniconda.html))
- Install [GNU Parallel](https://www.gnu.org/software/parallel/)
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

These instructions are for linux.
- Put *epistat.7.1* to place you like (for example your home directory) and *pdbUtils*, *TreeUtils* and *DnaUtilities* to *epistat_lib* dir somewhere (for example your home directory), than add EPISTAT_HOME and EPISTAT_LIB to your .bash_profile, for example:
  ```
  export EPISTAT_HOME=$HOME/epistat.7.1/
  export EPISTAT_LIB=$HOME/epistat_lib/
  ```
- Copy the example files to a separate directory:
```
  mkdir example
  cp $EPISTAT_HOME/prm_templates/intragene/nophen/* ./example
  cd example
```
- Run these commands from the shell:
  - `$EPISTAT_HOME/estimate_tau.pl -b 0.95 example.xparr > example.tau`\
  This script estimates 'tau' parameter and prints it to example.tau. Copy the value of "TAU=" from example.tau to the field "TAU=" in example.epistat.prm.
  - `$EPISTAT_HOME/run_epistat.pl -x example.xparr -m 3 -p 10000 example.epistat.prm run_epistat.prm`
  - `$EPISTAT_HOME/estimate_fdr.pl example.upper.pvalue.sites.fdr.prm > example.upper.pvalue.sites.fdr`
  - `$EPISTAT_HOME/mk_coevolmtx.pl –m Z-score example.mk_coevolmtx.prm`
  - `$EPISTAT_HOME/minvert/cor2pcor.R -f example.block.mtx -l 0.9 -n 0`
  - `$EPISTAT_HOME/graph/mk_graph.pl –u –p pos.edges –n neg.edges –o example –f graphml`
  - `$EPISTAT_HOME/graph/find_partition.py –n example.negative_edges.graphml  –p example.positive_edges.graphml –r 10000`
  - `$EPISTAT_HOME/stat/xparr_site_groups_test.pl –x example.xparr –a crosstab –b site2group.tab –с $EPISTAT_HOME/stat/figtree.colors > example.site_groups_test.out`
  - `$EPISTAT_HOME/stat/summ_site_groups_tests.pl –a parent –x example.xparr site_groups_test.out.list`
  
## Advanced configaration

### Estimation of the presence of signal in the data: estimate_fdr.pl

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
The "sites" filter requires at least BGR_SiteMinMutations substitution and FGR_SiteMinMutations substitution in the first (background) and the second (foreground) sites in the pair. The "pairs" filter is more restrictive, it requires that consecutive pairs of mutations in pairs of sites contain at least BGR_SiteMinMutations and FGR_SiteMinMutations of *different* mutations. By default, you should use the filter "pairs" and 2 mutations in sites, if N/L >= 5, where N is the number of sequences in the alignment, L is the length of the alignment. Disabling the filter for the minimum number of pairs is necessary if the initial alignment contains few sequences (N/L < 5), for this you should set:
```
MutationNumbersFilter="sites"
BGR_SiteMinMutations="1"
FGR_SiteMinMutations="1"
```
*example.upper.pvalue.(sites|pairs).fdr* contains the table:\
`pvalue	#obs	#exp	P`\
Here pvalue is the nominal p-value; #obs is the number of observations in real data corresponding to a given or smaller nominal p-value; #exp is the average number of observations in the null model, in which sites do not interact; P is the proportion of random realizations of the null model, where #obs <= #exp. If P <0.05, then for the corresponding nominal pvalue we have FDR < 1. We can assume that a signal is present in the data if in the range 0 <= pvalue <= 0.05 a significant proportion of pairs have nominal p-values, such that FDR < 1.

### Coevolution matrix creation: mk_coevolmtx.pl

`mk_coevolmtx.pl –m Z-score example.mk_coevolmtx.prm`

In the parameter file *example.mk_coevolmtx.prm*, set the filter to the minimum mutation value in the sites corresponding to the strongest signal in the data (see above). For example:
```
MutationNumbersFilter="pairs"
BGR_SiteMinMutations="2"
FGR_SiteMinMutations="2"
```
or
```
MutationNumbersFilter="sites"
BGR_SiteMinMutations="1"
FGR_SiteMinMutations="1"
```
The resulting covariance matrix will be in files: *example.block.mtx* and *example.mtx*

### Site interaction estimation: cor2pcor.R

`cor2pcor.R -f example.block.mtx -l 0.9 -n 0`

The output file *example.cor2pcor.R.out* will contain the result of the analysis: the values of the ordering statistics *wi_score*, reflecting the strength of the interaction between sites in a pair.

### Statistical analysis of results

The result of the work of EpiStat is a list of ordered or unordered pairs of sites, the evolution of which is coordinated (positive epistasis) or discoordinated (negative epistasis) in time. We measure the strength of a positive or negative association by the elements of the covariance matrix or *wi_score* statistics, which is proportional to the values of the inverse covariance matrix taken in modulus. The p-value is a measure of the statistical significance of associations. These statistics can be used to order the list of site pairs from the most strongly associated sites to pairs of independently evolving sites.

### Building a coevolution graph: mk_graph.pl 

`mk_graph.pl –u –p pos.edges –n neg.edges –o example –f graphml`

The coevolution graph is used to find groups of sites with a similar substitution patterns in a phylogenetic tree. The nodes of the graph are sites, and the edges with positive (negative) weights determine the presence of coordinated (discordant) evolution in pairs of sites. To build a graph, you need to define lists of edges and their weights. Edges with positive and negative weights are specified in separate files.

As a result, two files in GRAPHML format will be generated containing graph slices with positive and negative edges, respectively: *example.positive_edges.graphml* and *example.negative_edges.graphml*

### Coevolving site groups detection: find_partition.py 

`find_partition.py –n example.negative_edges.graphml –p example.positive_edges.graphml –r 10000`

The find_partition.py script implements the search for clusters of graph vertices with a high density of positive connections within and negative connections between clusters.

### Changing of substitution rates in site groups: xparr_site_groups_test.pl

`xparr_site_groups_test.pl –x example.xparr –a crosstab –b site2group.tab [–с $EPISTAT_HOME/stat/figtree.colors] > example.site_groups_test.out`

The *xparr_site_groups_test.pl* script calculates on which edges of the tree the relative evolutionary rates have changed in the given site groups. The *site2group.tab* file maps the site ID to the group ID. The optional parameter '-c' sets the correspondence between the group number and the color for displaying the results on the tree in FigTree format. The main result is printed as a table to standard output.

### Determination of coordinated changes in the rates of substitutions in different genes: summ_site_groups_tests.pl

`summ_site_groups_tests.pl –a parent –x example.xparr site_groups_test.out.list`

In case when there are several genes whose evolution corresponds to one phylogenetic tree, for example, example.xparr, it is possible to match the edges, in which in each gene there was a change in the relative rates of evolution. The summ_site_groups_tests.pl script determines the probability of a random match. The *site_groups_test.out.list* file should contain a list of file paths with the results of the *xparr_site_groups_test.pl* script for each gene.

## Reproducing publication results

To reproduce publication results for mitochondrian proteins you can use configuration files from the [data](https://github.com/gFedonin/EpiStat/tree/master/data) folder.
