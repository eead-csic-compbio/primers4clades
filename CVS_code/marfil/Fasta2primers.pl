#!/usr/bin/perl -w
# Bruno Contreras-Moreira, Pablo Vinuesa, ...
# 2005-7 UNAM, Mexico (maya)

# Fasta2primers.pl

# program that takes nucleotide sequence sets and produces pairs of evaluated primers
# for amplicons found in those sequences

use strict;
use FindBin '$Bin';
use lib "$Bin";
use lib "$Bin/bioperl-1.5.2_102/";
use Getopt::Std;
use File::Basename;
use phyTools;

set_phyTools_env();

my $progname = $0; # Fasta2primers.pl

my $CUTAMPLICONSEXE    = "./cut_amplicons.pl";
my $EVALUATEPRIMERSEXE = "./evaluate_primers.pl";
my $RANKPRIMERSEXE     = "./rank_primers.pl";
my $MARFILPROTEXE      = "./marfil_prot.pl -D -l ";
my $MARFILAMPSEXE      = "./marfil_amps.pl -l ";
my $CATEXE             = "/bin/cat";

my $DEFDNAMODEL = 'GTRG';
my $DEFPROTMODEL= 'WAGG';
my $DUMMYCUTOFF = 999;

########################################################################

my (%opts,$marfilEXE,$cwd);
my ($fnafile,$faafile,$exclude_list,$phylo_evaluation,$keep_files,$use_marfil_prot);
my ($only_unique_seqs,$server_output,$bootstrap,$branchLRT,$corrected_primers_file,$user_CUTs_path);
my ($skip_muscle,$limit_evaluations,$limit_runtime,$n_of_evaluations,$amplicon_length_range,$calc_CUT,$default_model,$user_model) = (0,0,0,0,'',0,0,'');
my $codehopT = $ENV{'DEF_CODEHOP_T'};
my $oligo4thermo = 'c';

getopts('hfCSskBLMRlm:c:P:n:p:e:T:o:r:D:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0)) 
{ 
	print   "\nusage: $progname [options]\n\n";
	print   "-h \t this message\n";
	print   "-n \t file with nucleotide sequences in FASTA format              (required unless -c is used)\n";
	print   "-p \t file with protein sequences in FASTA format                 (required unless -c is used)\n";
	print   "-e \t exclude sequence identifiers contained in this file         (optional, ignored if -c is used)\n";
	print   "-u \t use only unique sequences                                   (optional, ignored if -c is used)\n";
	print   "-B \t calculate WAG bootstrap support for amplicon topologies     (optional)\n";
	print   "-L \t calculate WAG+G LRT support for amplicon topologies         (optional, much faster than -B, ignored with -P)\n";
	print   "-c \t use outfile from correct_primers4taxa.pl                    (optional, skips cut_amplicon.pl)\n";
	print   "-P \t phylogenetically evaluate amplicons with quality score\n". 
	        "   \t >= cutoff in the range [best=100,worst=0]                   (optional, example -P 80)\n";
	print   "-M \t use marfil_prot.pl for phylogenetical evaluation            (optional, default marfil_amps.pl, needs -P)\n";
	print   "-f \t fast evaluation using default model (WAGG/GTRG)             (optional, applies to -P/-M)\n"; 
	print   "-m \t use provided model for evaluation                           (optional, example -m WAGG)\n"; 
	print   "-T \t temperature for codehop primers                             (optional, default=$codehopT)\n";
	print   "-o \t oligos to be thermo-evaluated (c=corr,r=relax,d=deg)        (optional, default=$oligo4thermo)\n";
	print   "-r \t range of desired amplicon length (format: -r 400:900)       (optional, by default all amplicons are taken\n";
	print   "-C \t calculate codon usage table from input sequences            (optional)\n";  
	print   "-l \t limit max number of phyloevaluations                        (optional, max=$ENV{'MAXPHYLOEVALS'})\n";
	print   "-R \t limit runtime for phyloevaluations                          (optional, max=$ENV{'MAXRUNTIME'}s)\n";
	print   "-S \t skip multiple alignment if sequences are aligned            (optional)\n";
	print   "-k \t keep all temporary files                                    (optional)\n";
	print   "-s \t format ouput for web server                                 (optional)\n";
	print   "-D \t temporary directory to place CUTs                           (optional, default=$ENV{'BLIMPS_DOCS'})\n";
	exit; 
}

if(defined($opts{'c'}))
{ 
	$corrected_primers_file = $opts{'c'}; 
	if(!-s $corrected_primers_file){ die "# $0 : need a valid outfile from correct_primers4taxa.pl\n"; }
}
else{ $corrected_primers_file = '' }

if(defined($opts{'n'})){ $fnafile = $opts{'n'}; }
else
{ 
	$fnafile = '';
	if(!$corrected_primers_file){ die "# $0 : need a nucleotide file in FASTA format\n"; } 
}
if(defined($opts{'p'})){ $faafile = $opts{'p'}; }
else
{
	$faafile = '';
	if(!$corrected_primers_file){ die "# $0 : need a protein file in FASTA format\n" }
}
if(defined($opts{'e'})){ $exclude_list = $opts{'e'}; }
else{ $exclude_list = '' }
if(defined($opts{'u'})){ $only_unique_seqs = 1; }
else{ $only_unique_seqs = 0; }

if(defined($opts{'r'})){ $amplicon_length_range = $opts{'r'}; }

if(defined($opts{'P'}))
{ 
	$phylo_evaluation = $opts{'P'}; 
	if($phylo_evaluation > 100 || $phylo_evaluation < 0)
	{
		die "# $progname : -P cutoff must be in the range [best=100,worst=0]\n";
	}		  
}
else{ $phylo_evaluation = $DUMMYCUTOFF; } #print "###phylo_evaluation is = $phylo_evaluation####\n";
if(defined($opts{'T'})){ $codehopT = $opts{'T'}; }
if(defined($opts{'S'})){ $skip_muscle = 1; }
if(defined($opts{'k'})){ $keep_files = 1; }
else{ $keep_files = 0 }
if(defined($opts{'s'})){ $server_output = 1; }
else{ $server_output = 0; }
if(defined($opts{'f'})){ $default_model = 1; }
if(defined($opts{'m'})){ $user_model = $opts{'m'}; }
if(defined($opts{'M'}))
{ 
	$use_marfil_prot = 1; 
	$marfilEXE = $MARFILPROTEXE;
	if($default_model){ $marfilEXE .= " -m $DEFPROTMODEL "; }
	elsif($user_model){ $marfilEXE .= " -m $user_model "; }
        $marfilEXE .= ' -i ';
}
else
{ 
	$use_marfil_prot = 0; 
	$marfilEXE = $MARFILAMPSEXE;
	if($default_model){ $marfilEXE .= " -m $DEFDNAMODEL "; }
	elsif($user_model){ $marfilEXE .= " -m $user_model "; }
        $marfilEXE .= ' -n ';
}
if(defined($opts{'B'})){ $bootstrap = 1; }
else{ $bootstrap = 0; }
if(defined($opts{'L'}) && !defined($opts{'P'})){ $branchLRT = 1; }
else{ $branchLRT = 0; }
if(defined($opts{'C'})){ $calc_CUT = 1; }
else{ $calc_CUT = 0; }
if(defined($opts{'o'}))
{
        $oligo4thermo = $opts{'o'};
        if($oligo4thermo !~ m/[c|r|d]/){ $oligo4thermo = 'c' }
}
if(defined($opts{'l'})){ $limit_evaluations = 1 }
if(defined($opts{'R'})){ $limit_runtime = 1 }
if( defined($opts{'D'}) && -d $opts{'D'}){ $user_CUTs_path = $opts{'D'}; }
else{ $user_CUTs_path = '' }

print "# $progname -n $fnafile -p $faafile -e $exclude_list -P $phylo_evaluation -M $use_marfil_prot -k $keep_files -r $amplicon_length_range ";
print "-S $skip_muscle -u $only_unique_seqs -T $codehopT -o $oligo4thermo -s $server_output -B $bootstrap -L $branchLRT -c $corrected_primers_file -C $calc_CUT -f $default_model -m $user_model -l $limit_evaluations -R $limit_runtime -D $user_CUTs_path\n\n";

$cwd = `pwd`; chomp $cwd;

print "# $progname : cwd = $cwd\n\n";

############################################################################

my ($mapsize,@primerfiles,@nt_amp_files,@aa_amp_files,@marfil_files,$listfile,$tabprimersfile,$n_of_sequences,$command,$outfiles,$amp_alignfile) = (0);
my (%table_stats,$table,$args,$nt_amp,$aa_amp,$primer_quality,$amplicon_quality,$codon_table,$marfil_command,$marfil_input,$marfil_no_input,$runtime);
my ($substitution_model,$full_marfil_stats,$phyfile,$PHYMLtree,$puzzleparamsfile,$primersfile,$modelfile,$puzzledistfile,$puzzleLMfile,$PHYMLstatsfile,$tarfile);

# 1) cut possible amplicons ################################################

print "# cutting amplicons...\n";
$args = "-T $codehopT";
if($exclude_list ne ''){ $args .= " -e $exclude_list"; }
if($phylo_evaluation != $DUMMYCUTOFF || $bootstrap || $branchLRT){ $args .= " -A"; }
if($only_unique_seqs == 1){ $args .= " -u"; }
if($server_output == 1){ $args .= " -s"; }
if($amplicon_length_range ne ''){ $args .= " -r $amplicon_length_range " }
if($calc_CUT){ $args .= " -c " }
if($skip_muscle){ $args .= " -S " } 
if($user_CUTs_path){ $args .= " -P $user_CUTs_path " }
if($corrected_primers_file)
{
	$command = "$CATEXE $corrected_primers_file ";
}
else
{
	$command = "$CUTAMPLICONSEXE -n $fnafile -a $faafile $args";
}

open(CUTAMP,"$command |") || die "# $progname : cannot execute $command\n";
while(<CUTAMP>)
{
	#print;
	if(/^# rank_codon_usage_tables : (\S+)\.codon\.use (\d\.\d+)/) 
	{
		$table = $1;
		$table_stats{$table}{'redundancy'} = $2;
		$table_stats{$table}{'n_of_pairs'} = 0;
	}		
	elsif(/^- /)
	{ 
		push(@primerfiles,(split)[1]); 
		$table = (split(/\.primers/,(split(/\__/,(split)[1]))[1]))[0];  
		$table_stats{$table}{'n_of_pairs'}++;
	}
	elsif(/^# reference protein coordinates: \d+ .. (\d+)/){ $mapsize = $1 }
	elsif(/need at least 2 sequences to run this program/)
	{
		print;
		exit(-1);
	}
	elsif(/^# WARNING : failed parsing taxa in FASTA headers/ || /^# number of intron-spanning primers skipped/)
	{
		print;
	}
	elsif(/^# input sequences/){ print }
}
close(CUTAMP); 

print "## table redundancy vs amplicons stats:\n";
foreach $table (sort { $table_stats{$b}{'n_of_pairs'}<=>$table_stats{$a}{'n_of_pairs'} || $table_stats{$a}{'redundancy'}<=>$table_stats{$b}{'redundancy'}}(keys(%table_stats)))
{
	if(!$table_stats{$table}{'redundancy'}){ $table_stats{$table}{'redundancy'} = '?' } # for -c cases
	print "# $table redundancy = $table_stats{$table}{'redundancy'} amplicons = $table_stats{$table}{'n_of_pairs'}\n";
} print "\n";

if(!@primerfiles)
{
	print "# $progname : sorry, could not find pairs of primers for these sequences...\n"; exit;
}
else
{
	if($corrected_primers_file eq ''){ $listfile = basename($fnafile)."_primers.list"; }
	else{ $listfile = basename($corrected_primers_file)."_primers.list"; }
	open(PRIMERSLIST,">$listfile") || die "# $progname : cannot create $listfile\n";
	foreach my $file (@primerfiles)
	{
		print PRIMERSLIST "$file\n";
	}
	close(PRIMERSLIST);
}

# 2) evaluate pairs of primers ##############################################

$tabprimersfile = '';
print "# evaluating primers...\n";
open(EVAL,"$EVALUATEPRIMERSEXE -t -f $listfile -T $oligo4thermo |") || die "# $progname : cannot execute $EVALUATEPRIMERSEXE -t -f $listfile -T $oligo4thermo\n";
while(<EVAL>)
{
	if(/^## summary: crosspot/)
	{ 
		## summary: crosspot=0.40 minsize=910 maxsize=910 tabfile=list.oligo_eval.tab
		$tabprimersfile = (split)[5];
		$tabprimersfile = (split(/=/,$tabprimersfile))[1];
	}
}
close(EVAL);

if(!-s $tabprimersfile)
{
	die "# $progname : sorry, could not evaluate the primer set...\n";
}

print "# primer file.tab = $tabprimersfile\n";

# 3) rank primers ############################################################

if($server_output){ print "\n# min quality for phylogenetic evaluations = $phylo_evaluation\%\n\n"; }

$args = '';
if($server_output){ $args = ' -s ' } 
if($mapsize){ $args .= " -m $mapsize " }

print "# ranking pairs of primers...\n";
open(RANK,"$RANKPRIMERSEXE -i $tabprimersfile -d $cwd $args |") || die "# $progname : cannot execute $RANKPRIMERSEXE -i $tabprimersfile -d $cwd\n";
while(<RANK>)
{
	next if(/$RANKPRIMERSEXE/);
	
	if(/^# source = (\S+)/)
	{
		$primersfile = $1;
		$nt_amp = $primersfile;
		$aa_amp = $primersfile;
		$nt_amp =~ s/\.primers/\.dna_amp/;
		$aa_amp =~ s/\.primers/\.aa_amp/; #print "# $nt_amp | $aa_amp\n";
		push(@nt_amp_files,$nt_amp);
		push(@aa_amp_files,$aa_amp);
	}
	elsif(/^## Amplicon (\d+) \S+ \d+:/)
	{
		$codon_table = (split(/\.primers/,(split(/__/,$primersfile))[1]))[0];
		print "## Amplicon $1 codon_usage_table = $codon_table :\n";
	}
	else
	{
		print; # don't forget to place output
		
		if(/^# primer pair quality = (\-?\d+)/){ $primer_quality = $1; }
		elsif(/^#--/)
		{
			my @times = times(); 
			$runtime = $times[0] + $times[2]; # process + child processes user times 

			# run marfil_amp/marfil_prot as required and print summary
			if($primer_quality >= $phylo_evaluation && 
				(!$limit_evaluations || $n_of_evaluations < $ENV{'MAXPHYLOEVALS'}) &&
				(!$limit_runtime || !$n_of_evaluations || ($n_of_evaluations && $runtime < $ENV{'MAXRUNTIME'}))) 
			{	
				$amplicon_quality = $marfil_command = '';
				
				if($use_marfil_prot)
				{ 
					$marfil_command = "$marfilEXE $aa_amp"; 
					$marfil_input = $aa_amp;
					$marfil_no_input = $nt_amp;
				}
				else
				{ 
					$marfil_command = "$marfilEXE $nt_amp"; 
					$marfil_input = $nt_amp;
					$marfil_no_input = $aa_amp;
				}#die"# $marfil_command\n";
				
				$substitution_model = $full_marfil_stats = $outfiles = '';
				open(MARFIL,"$marfil_command |") || die "# $progname : cannot execute $marfilEXE $nt_amp\n";
				while(<MARFIL>)
				{
					if(/^# output : subs_model = (\S+) stats = (.+)/)
					{
						$substitution_model = $1;
						$full_marfil_stats = $2; #print "$substitution_model $full_marfil_stats\n";
						
						# alignID \t alignname \t n_of_sequences \t nt_sites \t %const_sites \t n_of_site_patterns \t %const_site_patterns \t 
						# percent_gaps \t ChiSQ_failed \t Acc_No_ChiSQ_perc \t %fullyresolved_quartets \t %partiallyresolved_quartets \t
						# %nonresolved_quartets \t normalized_FRQ
						# wanted : 2,3,8,10,11,12,13?
				
						$amplicon_quality .= "# subs.model = $substitution_model\n";
						$amplicon_quality .= "# n_of_seqs = ".(split(/\t/,$full_marfil_stats))[2]." n_alignments_sites = ".(split(/\t/,$full_marfil_stats))[3]."\n";
						$amplicon_quality .= "# n_of_seqs_with_composition_bias = ".(split(/\t/,$full_marfil_stats))[8]."\n";
						$amplicon_quality .= "# \%fully resolved quartets = ".(split(/\t/,$full_marfil_stats))[10]." \%partially resolved quartets = ".
						(split(/\t/,$full_marfil_stats))[11]." \%non resolved quartets = ".(split(/\t/,$full_marfil_stats))[12]."\n";
					}	
					elsif(/^# outfiles : (\S+) (\S+) (\S+) (\S+) (\S+) (\S+) (\S+) (\S+)/)
					{
						#print;
						#guardar para posteriores analisis (al igual que .tab), posiblemente con tar
						($phyfile,$PHYMLtree,$puzzleparamsfile,$modelfile,
						$puzzledistfile,$puzzleLMfile,$PHYMLstatsfile,$amp_alignfile) = 
							($1,$2,$3,$4,$5,$6,$7,$8); 
						$outfiles = "$phyfile $PHYMLtree $puzzleparamsfile $puzzledistfile $puzzleLMfile $PHYMLstatsfile $amp_alignfile";
						if($modelfile ne '-'){ $outfiles .= " $modelfile" }
						#$amplicon_quality .= "# outfiles = $outfiles\n";

						if($server_output){ $amplicon_quality .= "# alignment file = $amp_alignfile\n"; }
					}
					elsif(/^# mean aLRT = \d+\.\d+ median aLRT = \d+\.\d+/)
					{
						$amplicon_quality .= $_;
					}
					elsif(/^# alpha = (\S+)/){ $amplicon_quality .= "# alpha = $1\n"; }
					#elsif(/^: sorry,/){}
				}
				close(MARFIL);
		
				if($full_marfil_stats eq ''){ $amplicon_quality = "# marfil $marfil_input failed\n# end_of_amplicon\n";  }
				else
				{
					$n_of_evaluations++;
					if($server_output){ $amplicon_quality .= "# aLRT likelihood tree file = $PHYMLtree\n"; }
					else
					{
						# tar outfiles
						$tarfile = $primersfile;
						$tarfile =~ s/\.primers/\.tgz/;
						system("tar cfz $tarfile $marfil_input $marfil_no_input $outfiles");
						push(@marfil_files,$phyfile,$PHYMLtree,$puzzleparamsfile,$puzzledistfile,$puzzleLMfile,$PHYMLstatsfile);
						if($modelfile ne '-'){ push(@marfil_files,$modelfile); }
						$amplicon_quality .= "# compressed outfile = $tarfile\n";
					}
				}
				
				print "# phylogenetic amplicon evaluation:\n$amplicon_quality";
				if(!$bootstrap && !$branchLRT){ print "# end_of_amplicon\n"; }		
			}
			elsif($primer_quality < $phylo_evaluation && !$bootstrap && !$branchLRT){ print "# end_of_amplicon\n"; }
			elsif($primer_quality >= $phylo_evaluation && $limit_evaluations && $n_of_evaluations > $ENV{'MAXPHYLOEVALS'})
			{
				print "# phylogenetic evaluation skipped (max $ENV{'MAXPHYLOEVALS'} evaluations)\n# end_of_amplicon\n";
			}
			elsif($primer_quality >= $phylo_evaluation && $limit_runtime && $runtime > $ENV{'MAXRUNTIME'})
			{ 
				print "# phylogenetic evaluation skipped (max runtime $ENV{'MAXRUNTIME'}s)\n# end_of_amplicon\n"; 
			}
			
			$n_of_sequences = 0;	
			if($bootstrap || $branchLRT)
			{			
				open(AMP,$aa_amp) || die "# $0 : cannot read $aa_amp\n";
				while(<AMP>){ if(/^>/){ $n_of_sequences++ } }
				close(AMP);
			}
			
			if($n_of_sequences > 3) # must have 4+ sequences at least to build trees
			{
				if($bootstrap)
				{
					# hacer .phy y quizás en un futuro calcular alpha, de momento solo WAG
					my $ampphyfile = convert_FAS2PHY($aa_amp,!$keep_files); #print "# mira: $ampphyfile\n";
					my %ampFASTA = read_FASTA_sequence($aa_amp);				
					my $bootphyfile = $ampphyfile;
					$bootphyfile =~ s/\.phy/_boot.phy/;
					rename($ampphyfile,$bootphyfile);

					my ($mean_boot_support,$median_boot_support,$bootreefile) = 
						get_mean_boot4WAG_dist_tree($bootphyfile,\%ampFASTA,!$keep_files);
				
					print "# mean bootstrap value = $mean_boot_support median bootstrap value = $median_boot_support\n";
					if($server_output){ print "# bootstrap tree file = $bootreefile\n"; }
					else{ push(@marfil_files,$bootreefile) }
					
					push(@marfil_files,$bootphyfile);
				
					if(!$branchLRT){ print "# end_of_amplicon\n"; }	
				}
			
				if($branchLRT)
				{
					my $ampphyfile = convert_FAS2PHY($aa_amp,!$keep_files); 
					my %ampFASTA = read_FASTA_sequence($aa_amp);			       
					my $aLRTphyfile = $ampphyfile;
					$aLRTphyfile =~ s/\.phy/_aLRT.phy/;
					rename($ampphyfile,$aLRTphyfile);

					my ($mean_aLRT,$median_aLRT,$aLRTtreefile) = 
					       get_mean_aLRT4WAGG_tree($aLRTphyfile,\%ampFASTA,!$keep_files);
				
					print "# mean aLRT = $mean_aLRT median aLRT = $median_aLRT\n";
					if($server_output){ print "# aLRT likelihood tree file = $aLRTtreefile\n"; }
					else{ push(@marfil_files,$aLRTtreefile) }
					
					push(@marfil_files,$aLRTphyfile);
					
					print "# end_of_amplicon\n";
				}
			}
			elsif($bootstrap || $branchLRT){ print "# end_of_amplicon (less than 4 sequences for -L/-B)\n"; }		
		}
	}
}
close(RANK);


# 4) clean tmp files ################################################################
if(!$keep_files)
{
	system("rm -f $listfile");
	unlink(@marfil_files);
	
	if(!$corrected_primers_file)
	{
		unlink(@primerfiles); 
		if(!$server_output){ unlink(@nt_amp_files); }
		unlink(@aa_amp_files);
	}
}

