#!/usr/bin/perl -w
# Bruno Contreras-Moreira, Pablo Vinuesa, ...
# 2005-6 UNAM, Mexico


use strict;
use Getopt::Std;
use FindBin '$Bin';
use lib "$Bin";
use lib "$Bin/bioperl-1.5.2_102/";
use phyTools;
set_phyTools_env();

my $progname = "marfil_prot.pl";
my $REMOVETMPFILES = 1;
my @files2remove;

######################################################

my ($AAfile,$myAAfile,$exclude_file,$do_not_copy_input_file,$seq,$skip_muscle,$skip_ProtTest,,$label_tree,%opts) = ('','','',0,0,0,0);

getopts('hli:e:Dsm:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0)) 
{ 
	print "usage: $progname [options]\n";
	print "-h \t this message\n";
	print "-i \t input file with protein sequences .faa\n";
	print "-e \t exclude file                              (optional)\n";
	print "-s \t skip making alignment, input is aligned   (optional)\n";
	print "-m \t skip ProtTest, run model provided by user (optional,[WAGG|JTTG|...])\n";
	print "-D \t do not copy input file                    (optional, ignored with -e)\n";
	print "-l \t label tree branches with FASTA headers    (optional)\n\n";
	exit; 
}

if(defined($opts{'i'})){ $AAfile = $opts{'i'}; }
else{ die "# $0 : need a valid protein FASTA file\n" }
if(defined($opts{'e'})){ $exclude_file = $opts{'e'}; }
if(defined($opts{'D'} && $exclude_file eq '')){ $do_not_copy_input_file= 1; }
if(defined($opts{'s'})){ $skip_muscle = 1 }
if(defined($opts{'l'})){ $label_tree = 1 }
if(defined($opts{'m'})){ $skip_ProtTest = $opts{'m'} }
else{ $skip_ProtTest = '' }

print "# $progname : -i $AAfile -e $exclude_file -D $do_not_copy_input_file -s $skip_muscle -m $skip_ProtTest -l $label_tree\n";

######################################################

if($do_not_copy_input_file == 0)
{
	# -2) read exclusion file that contains species that should be removed from the analysis
	my @EXCLIST;
	if($exclude_file)
	{
		open(EXCLIST,$exclude_file) || die "#progname : cannot read $exclude_file\n";
		while(<EXCLIST>)
		{
			next if(/^#/);
			push(@EXCLIST,(split)[0]);
		}
		close(EXCLIST);
	}

	# -1) add sequence number to original FASTA sequences
	my %AAfasta = read_FASTA_sequence($AAfile);

	$AAfile = reverse $AAfile;
	$AAfile = (split(/\//,$AAfile))[0];
	$AAfile = reverse $AAfile;

	$myAAfile = "my_" . (split(/\./,$AAfile))[0] . ".faa";

	open(MYAA,">$myAAfile") || die "#progname : cannot write to $myAAfile\n";

	print "\n# FASTA header equivalences (AA):\n";

	foreach $seq (sort {$a<=>$b} (keys(%AAfasta)))
	{
		if($AAfasta{$seq}{'NAME'} =~ /\[(\S+)\]/){ next if(grep(/$1/,@EXCLIST)); }	
		
		$AAfasta{$seq}{'NAME'} =~ s/>//;
		$AAfasta{$seq}{'SEQ'} =~ s/\*//;
		
		print "# $seq => $AAfasta{$seq}{'NAME'}";
		
		print MYAA ">$seq " . $AAfasta{$seq}{'NAME'};
		print MYAA $AAfasta{$seq}{'SEQ'},"\n";
	}

	close(MYAA);
	print "\n";

	push(@files2remove,$myAAfile);
}
else{ $myAAfile = $AAfile }


# 0) run muscle to obtain multiple alignment of protein sequences
my $AA_aln_file;
if($skip_muscle){ $AA_aln_file = $myAAfile }
else
{
	$AA_aln_file = run_muscle($myAAfile,$REMOVETMPFILES); # changed from $myAAfile
	if($AA_aln_file eq $ERROR)
	{
		unlink(@files2remove);
       		die "# $progname : sorry, muscle alignment failed ($AAfile)\n";
	}
}
#push(@files2remove,$myAAfile);
 
# 1) convert DNA alignment to phylip and nexus formats
my $AA_phylip_file = convert_FAS2PHY($AA_aln_file,$REMOVETMPFILES);
if($AA_phylip_file eq $ERROR)
{
	unlink(@files2remove);
	die "# $progname : sorry, FAS2PHY failed ($AAfile)\n";
}

push(@files2remove,$AA_phylip_file);

my ($PROTTEST_bestAICmodel,$PROTTEST_outfile,$PHYMLparams);
unless ($skip_ProtTest)
{       
	# 2) run PROTTEST and get best fit model
	($PROTTEST_bestAICmodel,$PROTTEST_outfile) = run_PROTTEST($AA_phylip_file,$REMOVETMPFILES);
	if($PROTTEST_bestAICmodel eq $ERROR)
	{
		unlink(@files2remove);
		die "# $progname : sorry, run_PROTTEST failed ($AAfile)\n";
	}#print "|$PROTTEST_bestAICmodel|\n";#exit; 
}
else
{
	# 2.1) use model provided by user
	$PROTTEST_bestAICmodel = $skip_ProtTest;
	$PROTTEST_outfile = '-';
}

# 3) select PHYML parameters for this model
$PHYMLparams = getPHYMLparameters_prot($PROTTEST_bestAICmodel);
if($PHYMLparams eq $ERROR)
{
	unlink(@files2remove);
	die "# $progname : sorry, getPHYMLparameters_prot failed with model $PROTTEST_bestAICmodel ($AAfile)\n";
} #print "$PHYMLparams\n";

## 4) run PHYML with this parameters

my ($PHYML_tree,$PHYML_lk,$PHYML_stat,$mean_aLRT,$median_aLRT) = run_PHYML($AA_phylip_file,$PHYMLparams,$REMOVETMPFILES);
if($PHYML_tree eq $ERROR)
{
	unlink(@files2remove);
	die "# $progname : sorry, runPHYML failed\n";
}
elsif($label_tree)
{
        my ($newick,$labelled_newick) = ('','');

        open(TREE,$PHYML_tree) || die "# $0 : cannot read $PHYML_tree\n";
        while(<TREE>){ $newick .= $_; }
        close(TREE);

        if($newick)
        {
                my %ampFASTA = read_FASTA_sequence($AAfile);
                $labelled_newick = add_labels2newick_tree( $newick, \%ampFASTA );
                open(LABELTREE,">$PHYML_tree") || die "#$0 : cannot rewrite $PHYML_tree\n";
                print LABELTREE $labelled_newick;
                close(LABELTREE); #print "#$newick#$labelled_newick#\n";
        }
        else{ print "# $0 : cannot label tree $PHYML_tree\n"; }
}

push(@files2remove,$PHYML_tree,$PHYML_lk,$PHYML_stat);

## 5) create a file with PUZZLE_LM paramerers for protein alignments
my $PUZZLE_param_file = create_PUZZLE_LM_prot_parameters_file($AA_phylip_file,$PHYML_stat,$PROTTEST_bestAICmodel);
if($PUZZLE_param_file eq $ERROR)
{
	unlink(@files2remove);
	die "# $progname : sorry, get_PUZZLE_param_from_PHYML_stat failed ($AAfile)\n";
}

push(@files2remove,$PUZZLE_param_file);

## 6) run PUZZLE TREE for Likelihood Mapping (LM) with ML-param-file returned from get_PUZZLE_param_from_PHYML_stat()
my ($PUZZLE_dist,$PUZZLE_LM_puzzle,$alpha,$aln_stats,$distance_stats) = run_PUZZLE_LM($AA_phylip_file,$PUZZLE_param_file, $REMOVETMPFILES);
if($PUZZLE_LM_puzzle eq $ERROR)
{
	unlink(@files2remove);
	die "# $progname : sorry, runPUZZLE_LM failed ($AAfile)\n";
}

push(@files2remove,$PUZZLE_dist,$PUZZLE_LM_puzzle);

## 7) parse LM_PUZZLE file and print the corresponding MySQL.input table
my $DB_row = parse_LM_PUZZLE_file($PUZZLE_LM_puzzle);
if($DB_row eq $ERROR)
{
	unlink(@files2remove);
	die "# $progname : sorry, parse_LM_PUZZLE_file failed ($AAfile)\n";
} 

print "# output : subs_model = $PROTTEST_bestAICmodel stats = $DB_row\n";
print "# alpha = $alpha alignment stats : $aln_stats\n";
print "# distance stats : $distance_stats\n";
print "# mean aLRT = $mean_aLRT median aLRT = $median_aLRT\n";
print "# outfiles : $AA_phylip_file $PHYML_tree $PUZZLE_param_file $PROTTEST_outfile $PUZZLE_dist $PUZZLE_LM_puzzle $PHYML_stat $AA_aln_file\n";

if($REMOVETMPFILES)
{ 
	`rm -f $PHYML_lk`;
	#if(!$skip_muscle){ unlink($AA_aln_file) }
}


#open ( OUT, ">>Table_puzzle_LM_DNA_alphas.input" ) || die "cannot write Table_puzzle_LM_DNA_alphas.input";
#print OUT "$DB_row\n";
#close(OUT);

# 8) Upload data into MySQL DB
#if(!upload_LM_data2DB(%LM_puzzle_data))
#{
#  die "# $progname : sorry, upload_LM_data2DB failed\n";
#} 





