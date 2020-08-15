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
