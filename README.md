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

XPARR formatted trees can be generated from the output of pylognetic reconstruction program which can reconstruct both the tree and the ancestor states. The *mk_xparr.pl* script can be used:\
`mk_xparr.pl -t <tree_NWK> [options] <alignment_FASTA>`\
Here <tree_NWK> is a newick formated tree, where the bootstrap-slot is filled with internal nodes' identifiers, and <alignment_FASTA> is an FASTA alignment of both terminal and ancestor sequences.

# EpiStat workflow

EpiStat is implemented in perl 5 and use R libraries. The workflow itself consists of scenarios whic shuold be executed subsequently:

*estimate_tau.pl* - Estimates the scaling parameter *tau*, which reflects how fast or slow do the sequential substitutions happen one after another for an average pair of sites.\
*run_epistat.pl - Runs the main analysis.\
*[estimate_fdr.pl]* - Estimates FDR values for nominal p-values. This script is not required for obtaining of final result, but it is necessary to evaluate the presence or absence of signal in the data.\
*mk_coevolmtx.pl* - Forms a covariance matrix of substitutions in pairs of sites\
*minvert/cor2pcor.R* - Inverts the covariance matrix, normalizes the result.

For each script, the description of options and values of input parameters passed in files is described in the comments at the beginning of \*.pl files. The end result generated with *minvert/cor2pcor.R* is a matrix containing the *wi_score* statistics that is used to order the pairs. Larger statistic values correspond to pairs of sites with a higher probability of being close on the 3D structure of the protein. When using scripts from the package, the parameter file naming convention is adopted (optional). Let the source file be named example.xparr, then the parameter files are named: *example.<script_name>.prm*. The parameted values may depend on the data being processed.\
Examples of parameter files that can be used as templates can be found in: /epistat.7.1/prm_templates/intragene/nophen/.


