use strict;
use warnings;
use Module::Build;
 
use 5.008;

system("chmod -R u+x script");

my %stat_files;
open INPF, "<script_stat_list.txt" or die "\nUnable to open input file script_stat_list.txt!";
while(<INPF>){
	chomp;
	s/^\s+//;
	s/\s+$//;
	if(/\S+/){
		$stat_files{$_}=$_;
	}
}
close INPF;
my %minvert_files;
open INPF, "<script_minvert_list.txt" or die "\nUnable to open input file script_minvert_list.txt!";
while(<INPF>){
	chomp;
	s/^\s+//;
	s/\s+$//;
	if(/\S+/){
		$minvert_files{$_}=$_;
	}
}
close INPF;
my %utils_files;
open INPF, "<script_utils_list.txt" or die "\nUnable to open input file script_utils_list.txt!";
while(<INPF>){
	chomp;
	s/^\s+//;
	s/\s+$//;
	if(/\S+/){
		$utils_files{$_}=$_;
	}
}
close INPF;
my %graph_files;
open INPF, "<script_graph_list.txt" or die "\nUnable to open input file script_graph_list.txt!";
while(<INPF>){
	chomp;
	s/^\s+//;
	s/\s+$//;
	if(/\S+/){
		$graph_files{$_}=$_;
	}
}
close INPF;
my %examples_files;
open INPF, "<examples_list.txt" or die "\nUnable to open input file examples_list.txt!";
while(<INPF>){
	chomp;
	s/^\s+//;
	s/\s+$//;
	if(/\S+/){
		$examples_files{$_}=$_;
	}
}
close INPF;

my $class = Module::Build->subclass(
    class => 'My::Builder',
    code => q{
        sub ACTION_docs {
            use Pod::Text;
            use Pod::Select;

            my $self = shift;

            podselect({ -output => 'README.pod' },
                'lib/Epistat.pm');

            my $parser = Pod::Text->new();
            $parser->parse_from_file('README.pod', 'README');

            return $self->SUPER::ACTION_docs;
        }
    }
);

my $builder = $class->new(
    module_name         => 'Epistat',
    license             => 'gpl_3',
    dist_author         => 'Alexey D. Neverov <neva_2000(at)mail.ru>',
    dist_abstract       => 'Epistat - the set of utilities for the searching for associations on a phylogenetic tree',
    create_makefile_pl  => 0,     #'traditional',
    script_files        => 'script',
	dist_version_from   => 'lib/Epistat.pm',
    create_readme       => 1,
    requires            => {
		'Bio::Phylo'       => 0,
        'Time::Progress'    => 0,
        'List::BinarySearch'  => 0,
		'Color::Rgb'  => 0,
        'Math::Gradient' => 0,
        'File::Basename'      => 0,
        'File::Path'          => 0,
		'Getopt::Std'          => 0,
		'Statistics::Multtest' => 0,
		'Data::Float' => 0
    },
	stat_files => \%stat_files,
	graph_files => \%graph_files,
	minvert_files => \%minvert_files,
	utils_files => \%utils_files,
	examples_files => \%examples_files,
	configure_requires => {
		'Pod::Text' => 0,
		'Pod::Select' => 0,
	},
    build_requires      => {
        'Test::More'          => '0.47',
		'Pod::Text' => 0,
		'Pod::Select' => 0,
    },
    #meta_merge => {
    #    resources => {
    #        repository => 'https://github.com/szabgab/App',
    #        bugtracker => 'https://github.com/szabgab/App/issues'
    #    }
    #},
);
$builder->add_build_element('stat');
$builder->add_build_element('graph');
$builder->add_build_element('minvert');
$builder->add_build_element('utils');
$builder->add_build_element('examples');
$builder->install_base_relpaths('examples' => 'epistat_examples');
$builder->create_build_script();
