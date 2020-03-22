#!/usr/bin/perl -w
# Bruno Contreras-Moreira, Pablo Vinuesa, ...
# 2005-6 UNAM, Mexico

# program that translates nucleotide sequences, aligns them, builds a distance-based tree
# and finally produces clusters of sequences under a defined distance cutoff

## (old and finally runs branchclust in order to define sequence clusters)

use strict;
use FindBin '$Bin';
use lib "$Bin";
use lib "$Bin/bioperl-1.5.2_102/";
use Getopt::Std;
use phyTools;
use File::Basename;
use Bio::TreeIO; 

set_phyTools_env();

my $progname = $0;
my $MINSEQUENCES = 4;
my $MAXPAIRDIST = $ENV{"MAX_PAIR_DIST"};

########################################################################


my (%opts,$FASTA_DNA_FILE,$gencode,$maxdist,$filter_redundant,$print_matrix,$skip_muscle);
my ($userMATRIX,$user_cluster_boundaries,$user_excluded,$server_output,$useDNAdistances);

getopts('hDSsmfi:c:d:M:b:e:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0)) 
{ 
	print   "usage: $progname [options]\n";
	print   "-h \t this message\n";
	print   "-i \t file with nucleotide sequences in FASTA format\n";
	print   "-c \t genetic translation code                            (optional, default universal [2-16])\n";
	print   "-d \t maximum pairwise distance                           (optional, default = $MAXPAIRDIST)\n";
	print   "-D \t use DNA (HKY85+G) distances                         (optional, default = WAG+G)\n";
	print   "-M \t take user-provided distance matrix in Phylip format (optional)\n";
	print   "-f \t filter redundant sequences                          (optional)\n";
	print   "-m \t print distance matrix                               (optional)\n";
	print   "-s \t format output for the web server                    (optional)\n";
	print   "-S \t skip multiple alignment if sequences are aligned    (optional)\n"; 
	print   "-b \t comma-separated pair of cluster boundaries          (optional, requires -M)\n";
	print   "-e \t list of comma-separated sequences to be excluded    (optional, requires -M and -b)\n";
	print   "   \t within a -b delimited clade\n";
	exit; 
}

if(defined($opts{'i'})){ $FASTA_DNA_FILE = $opts{'i'}; }
else{ die "# $0 : need a nucleotide file in FASTA format\n" }
if(defined($opts{'c'})){ $gencode = $opts{'c'}; }
else{ $gencode = 'universal'; }
if(defined($opts{'d'})){ $maxdist = $opts{'d'}; }
else{ $maxdist = $MAXPAIRDIST; }
if(defined($opts{'D'})){ $useDNAdistances = 1; }
else{ $useDNAdistances = 0; }
if(defined($opts{'S'})){ $skip_muscle = 1; }
else{ $skip_muscle = 0; }
if(defined($opts{'f'})){ $filter_redundant = 1; }
else{ $filter_redundant = 0; }
if(defined($opts{'m'})){ $print_matrix = 1; }
else{ $print_matrix = 0; }
if(defined($opts{'s'})){ $server_output = 1; }
else{ $server_output = 0; }

$user_cluster_boundaries=$user_excluded='';
if(defined($opts{'M'}))
{ 
	if(-e $opts{'M'}){ $userMATRIX = $opts{'M'}; }
	else{ die "# $0 : need a valid distance matrix file in Phylip format\n" } 

	if(defined($opts{'b'}))
	{ 
		$user_cluster_boundaries = $opts{'b'}; 
		if(defined($opts{'e'}))
		{	
			$user_excluded = $opts{'e'};
		}
	}
}
else{ $userMATRIX = ''; }

if(!$server_output)
{
	print "# $progname -i $FASTA_DNA_FILE -c $gencode -d $maxdist -f $filter_redundant -m $print_matrix ".
	"-s $server_output -M $userMATRIX -b $user_cluster_boundaries -e $user_excluded -D $useDNAdistances -S $skip_muscle\n\n";
}

############################################################################

my (%FASTAseqs,%TAXA,$n_of_taxa,$n_of_sequences,%FASTA_AA_seqs,$dist_matrix,$label);
my ($FASTA_PROT_FILE,$seq,$sequence,$FASTA_PROT_ALN_FILE,$PHY_PROT_ALN_FILE,$PROT_DIST_FILE);
my ($alpha,$NEWICK_TREE_FILE,%CLUSTERS,$cluster,$clusterfile,$clusterfile_AA,%user_cluster);
my ($display_labelled_tree,$TREE_GRAPH_FILE,$display_graph_tree,$scale_string,$alignment_stats_string,$mean_distance_string);

# 0) parse FASTA file, detect and count present taxa
%FASTAseqs = read_FASTA_sequence( $FASTA_DNA_FILE, 0, $gencode, 1 ); # check CDS quality and removes identical sequences
%TAXA = find_taxa_FASTA_headers(\%FASTAseqs);

$n_of_taxa = scalar(keys(%TAXA));
$n_of_sequences = scalar(keys(%FASTAseqs));

print  "# number of sequences read = $n_of_sequences\n";
print  "# number of recognised taxa = $n_of_taxa\n";

if($n_of_sequences < $MINSEQUENCES)
{
	die "# $progname : abort, need at least $MINSEQUENCES sequences to run this program\n";
}

if($FASTAseqs{(keys(%FASTAseqs))[0]}{'SEQ'} =~ /RNDQEHILKMFPSWYVBZ/) 
{
	die "# $progname : abort, need nucleotide sequences to run this program\n";
} 

# 1) translate dna to protein sequence using defined $gencode
%FASTA_AA_seqs = convert_FASTA_dna2prot(\%FASTAseqs, $gencode);

$FASTA_PROT_FILE = (split(/\./,basename($FASTA_DNA_FILE)))[0] . ".faa";
open(FAA, ">$FASTA_PROT_FILE") || die "# $progname: can't write to file $FASTA_PROT_FILE $!\n";
foreach $seq (sort {$a<=>$b} (keys(%FASTA_AA_seqs)))
{
	$sequence = $FASTA_AA_seqs{$seq}{'SEQ'};
	$sequence =~ s/\*(-*)$/$1/g; # prealigned input might have trailing gaps! 
	if($sequence =~ /\*/)
	{
		print "\n# ERROR: Translated protein sequences contain internal STOP codon(s).\n";
		print "# Make sure that the first nucleotide in your DNA frames corresponds\n";
		print "# to the first codon position or/and that an appropriate translation table\n";
		print "# is applied.\n\n";
		print "# translation: $sequence\n\n";

		close(FAA);
		exit(-1);
	}

	print FAA ">$seq ".substr($FASTA_AA_seqs{$seq}{'NAME'},1)."$sequence\n";
}	
close(FAA);

if($skip_muscle)
{ 
	if(!check_aligned_FASTA_sequences(%FASTA_AA_seqs)){ $skip_muscle = 0; print "# input sequences are not aligned\n"; } 
	else{ print "# input sequences seem to be aligned\n"; }
}
 
if($userMATRIX eq '')
{
	# 2) align protein sequences (if required)
	if($skip_muscle){ $FASTA_PROT_ALN_FILE = $FASTA_PROT_FILE }
	else
	{
		print "# aligning translated sequences...\n";
		$FASTA_PROT_ALN_FILE = run_muscle($FASTA_PROT_FILE,1); 
		if($FASTA_PROT_ALN_FILE eq $ERROR)
		{
			die "# $progname : sorry, run_muscle failed\n";
		}
	}

	# 3) calculate distance matrix (if required) and NJ tree
	print "# computing distance matrix...\n";
	if($useDNAdistances)
	{
		my $myDNAfile = "my_" . basename($FASTA_DNA_FILE);
		open(MYDNA,">$myDNAfile") || die "# $progname : cannot write to $myDNAfile\n";
		foreach my $seq (sort {$a<=>$b} (keys(%FASTAseqs)))
		{
	        	$FASTAseqs{$seq}{'NAME'} =~ s/>//;   
		        print MYDNA ">$seq " . $FASTAseqs{$seq}{'NAME'} . $FASTAseqs{$seq}{'SEQ'},"\n";
		}
		close(MYDNA);
		
		my $DNA_aln_file = convert_AAaln2NTaln($myDNAfile,$FASTA_PROT_ALN_FILE);
		if($DNA_aln_file eq $ERROR){ die "# $progname : sorry, convert_AAaln2NTaln failed\n"; }
		
		my $DNA_phylip_file = convert_FAS2PHY($DNA_aln_file,1);
		if($DNA_phylip_file eq $ERROR){ die "# $progname : sorry, FAS2PHY failed ($DNA_aln_file)\n"; }
		
		($PROT_DIST_FILE,$alpha,$alignment_stats_string,$mean_distance_string) = run_PUZZLE_DIST_HKYG($DNA_phylip_file,1);
                if($PROT_DIST_FILE eq $ERROR)
                {
                        die "# $progname : sorry, run_PUZZLE_DIST_HKYG failed: $!\n";
                }
		else
                { 
                        my $correct_dist_filename = $PROT_DIST_FILE;
			$correct_dist_filename =~ s/my_//;
			system("mv -f $PROT_DIST_FILE $correct_dist_filename"); 
                        $PROT_DIST_FILE = $correct_dist_filename;
                }
                if(!$server_output)
		{ 
			print "# run_PUZZLE_DIST_HKYG alpha = $alpha\n"; 
			unlink($myDNAfile,$DNA_aln_file,$DNA_phylip_file); #print "$myDNAfile,$DNA_aln_file,$DNA_phylip_file\n";
		}
		else
		{
			print "# multiple alignment file = $DNA_aln_file\n";
                        unlink($myDNAfile,$DNA_phylip_file);
		}
	}
	else
	{
		$PHY_PROT_ALN_FILE = convert_FAS2PHY($FASTA_PROT_ALN_FILE,1);
        	if($PHY_PROT_ALN_FILE eq $ERROR){ die "# $progname : sorry, convert_FAS2PHY failed\n"; }

		($PROT_DIST_FILE,$alpha,$alignment_stats_string,$mean_distance_string) = run_PUZZLE_DIST($PHY_PROT_ALN_FILE,1);
		if($PROT_DIST_FILE eq $ERROR)
		{
			die "# $progname : sorry, run_PUZZLE_DIST failed: $!\n";
		}
		if(!$server_output){ print "# run_PUZZLE_DIST alpha = $alpha\n"; }
		else{ print "# multiple alignment file = $FASTA_PROT_ALN_FILE\n"; }
	}

	($NEWICK_TREE_FILE,$TREE_GRAPH_FILE) = run_NEIGHBOR($PROT_DIST_FILE,1);
	print "$alignment_stats_string\n";
}
else
{ 
	$PROT_DIST_FILE = $userMATRIX;
	$NEWICK_TREE_FILE = dirname($FASTA_PROT_FILE) . "/" . (split(/\./,basename($FASTA_DNA_FILE)))[0] . "_aln_nj.ph";
	$TREE_GRAPH_FILE = dirname($FASTA_PROT_FILE) . "/" . (split(/\./,basename($FASTA_DNA_FILE)))[0] . "_aln_nj_txt.graph";
	if(! -e $NEWICK_TREE_FILE){ ($NEWICK_TREE_FILE,$TREE_GRAPH_FILE) = run_NEIGHBOR($PROT_DIST_FILE,1) }

	if($user_cluster_boundaries ne '')
	{
		my ($node,@boundaries,$label);

		# 4.1) read newick tree 
		my $newick = new Bio::TreeIO(-file => $NEWICK_TREE_FILE, -format => "newick");
	       	my $tree = $newick->next_tree();
		#printf("# total tree length: %1.1f\n",$tree->total_branch_length)

		# 4.2) find clade boundary nodes
		foreach $node (split(/\,/,$user_cluster_boundaries))
		{
			my @node = $tree->find_node(-id => $node);
			if(@node){ push(@boundaries,@node); }
			else	
			{
				if($node =~ /(\d+)__\d+/)
				{
					$node = $1;
					my @node = $tree->find_node(-id => $node);
					if(@node){ push(@boundaries,@node); }
					else{ die "# $progname : cannot find node >$node< in tree $NEWICK_TREE_FILE\n"; }
				}	
			}
		}
		if(scalar(@boundaries) < 2){ die "# $progname : need at least two valid nodes as clade boundaries\n"; }
		else{ print "# number of valid cluster boundaries = ".scalar(@boundaries)."\n"; } 
	       	
       		# 4.3) find lowest common ancestor for boundaries
	       	my $lca_user_nodes = $tree->get_lca(-nodes => \@boundaries);
	       	
		# 4.4) select all nodes sharing $lca_user_nodes with boundary nodes
		my @treenodes = $tree->get_leaf_nodes();
		foreach $node (@treenodes)
		{
			my ($tmplca,@tmppair);
			
			push(@tmppair,$node);
			push(@tmppair,$boundaries[0]);
			
			$tmplca = $tree->get_lca(  -nodes => \@tmppair);

			if($tmplca == $lca_user_nodes || 
				$tree->is_monophyletic(-nodes => \@tmppair,-outgroup => $lca_user_nodes))
			{
				#print $node->id();
				next if($user_excluded =~ $node->id());
				$label = $node->id();
				$label =~ s/__\d+//;
				push(@{$user_cluster{0}},$label);
			}
			else # add non-members to special -1 cluster 
			{
				push(@{$user_cluster{-1}},$node->id());
			}
		}
		print "# size of user-selected cluster = ".scalar(@{$user_cluster{0}})."\n";
       	}       
}

# 5) cluster sequences according to distance matrix
($dist_matrix,$display_labelled_tree,$display_graph_tree,$scale_string,%CLUSTERS) = cluster_sequences($PROT_DIST_FILE,$maxdist,$filter_redundant,$NEWICK_TREE_FILE,$TREE_GRAPH_FILE,\%FASTAseqs,%user_cluster);

print "# distance matrix = $PROT_DIST_FILE\n\n$mean_distance_string\n$scale_string\n";
print "# NJ tree with cluster labels :\n\n";
print `cat $display_graph_tree`;


if(!$server_output)
{
	print "# labelled tree : $NEWICK_TREE_FILE\n";
	print "# fully labelled tree : $display_labelled_tree\n";
	print "# text formatted tree : $TREE_GRAPH_FILE\n";
	print "# fully labelled text formatted tree : $display_graph_tree\n";
}
else
{
	system("rm -f $display_labelled_tree $display_graph_tree");
}

#../../../bin/simcluster-0.8.14/src/treedraw.bin -i rpoC-Rhizobiales_extended_aln_nj_labelled.ph -l 
# genera rpoC-Rhizobiales_extended_aln_nj_labelled.ph.png     
## esto es muy viejo !!! %CLUSTERS = run_BRANCHCLUST($NEWICK_TREE_FILE,\%TAXA,$filter_paralogs,$clustersize,1); 


# 6) output FASTA cluster file names
foreach $cluster (sort {$a<=>$b} (keys(%CLUSTERS)))
{
	next if($cluster eq '-1'); # orphan sequences with no associated cluster 

	my ($base,$path,$ext) = fileparse($FASTA_PROT_FILE,".faa");
	$clusterfile = "$path/$base\_cluster_$cluster\.fna"; 
	$clusterfile_AA = "$path/$base\_cluster_$cluster\.faa";	

	open(CLUSTFASTA,">$clusterfile") || die "# $progname : cannot write to $clusterfile\n";
	foreach $seq (sort {$a<=>$b} (@{$CLUSTERS{$cluster}}))
	{
		$label = $FASTAseqs{$seq}{'NAME'};
		if(substr($label,0,1) eq '>'){ $label = substr($label,1) }

		if($FASTAseqs{$seq}{'IDENTICALS'} ne '')
		{	
			print CLUSTFASTA ">$seq ==$FASTAseqs{$seq}{'IDENTICALS'}== $label$FASTAseqs{$seq}{'SEQ'}\n";
		}
		else{ print CLUSTFASTA ">$seq $label$FASTAseqs{$seq}{'SEQ'}\n"; }
	}
	close(CLUSTFASTA);
	
	open(CLUSTFASTA,">$clusterfile_AA") || die "# $progname : cannot write to $clusterfile_AA\n";
        foreach $seq (sort {$a<=>$b} (@{$CLUSTERS{$cluster}}))
        {
		$label = $FASTA_AA_seqs{$seq}{'NAME'};
                if(substr($label,0,1) eq '>'){ $label = substr($label,1) }

		if($FASTAseqs{$seq}{'IDENTICALS'} ne '')
		{
                	print CLUSTFASTA ">$seq ==$FASTAseqs{$seq}{'IDENTICALS'}== $label$FASTA_AA_seqs{$seq}{'SEQ'}\n";
		}
		else{ print CLUSTFASTA ">$seq $label$FASTA_AA_seqs{$seq}{'SEQ'}\n"; }
        }
        close(CLUSTFASTA);

	printf("# CLUSTER %d members = %d : %s %s\n",$cluster,scalar(@{$CLUSTERS{$cluster}}),$clusterfile,$clusterfile_AA);
}

if($print_matrix)
{
	print "\n# distance matrix:\n$dist_matrix";
} 



