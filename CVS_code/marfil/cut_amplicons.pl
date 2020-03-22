#!/usr/bin/perl -w
# Bruno Contreras-Moreira, Pablo Vinuesa, ...
# 2005-7 UNAM, Mexico

# finds pairs of primers/amplicons using different codon usage tables until convergence,
# ie, 3 tables in a row must produce no amplicons in order to converge

use strict;
use Digest::MD5 qw(md5_hex);
use FindBin '$Bin';
use lib "$Bin";
use lib "$Bin/bioperl-1.5.2_102/";
use Getopt::Std;
use phyTools;
set_phyTools_env();

my $progname = $0;
my $REMOVETMPFILES = 1;
my $MAXSERVERHEADER = $ENV{'MAXPRIMERHEADER'};

######################################################

my (%opts,$DNAfile,$AAfile,$exclude_file,$user_codon_table,$calc_CUT,$full_CUT_spsum,$CUT_spsum,$user_CUT_path);
my ($print_amplicons,$force_all_tables,$only_unique_seqs,$server_output,$CUT_name,$CUT_table,$skip_muscle);
my $codehopT = $ENV{'DEF_CODEHOP_T'};
my $amplicon_length_range = ''; 

getopts('hsSAun:a:t:e:T:r:FcP:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0)) 
{ 
	print   "usage: $progname [options]\n";
	print   "-h \t this message\n";
	print   "-n \t nucleotide sequences in FASTA format <file.fna>\n";
	print   "-a \t aminoacid sequences in FASTA format <file.faa>\n";
	print   "-A \t print amplicon files                              (optional)\n";
	print   "-e \t exclude list                                      (optional)\n";
	print   "-u \t use only unique sequences                         (optional)\n";
	print   "-F \t force use of all available codon usage tables     (optional, should make it slower)\n";
	print   "-T \t temperature for codehop primers                   (optional, default=$codehopT)\n";
	print   "-r \t range of desired amplicon length (-r 400:900)     (optional, by default all amplicons are taken\n";
	print   "-s \t format output for the web server                  (optional)\n";
	print   "-S \t skip multiple alignment if sequences are aligned  (optional)\n";
	print   "-P \t path to place codon usage tables                  (optional)\n";
	print   "-c \t calculate codon usage tables from input -n        (optional, ignored with -t)\n";
	print   "-t \t user provided codon usage comma-separated file(s) without spaces, e.g. t1,t2,t3...  (optional, *.codon.use)\n";
	exit; 
}

if(defined($opts{'n'})){ $DNAfile = $opts{'n'}; }
else{ die "# $0 : need a valid nucleotide file in FASTA format\n" }
if(defined($opts{'a'})){ $AAfile = $opts{'a'}; }
else{ die "# $0 : need a valid aminoacid file in FASTA format\n" }
if(defined($opts{'e'})){ $exclude_file = $opts{'e'}; }
else{ $exclude_file = ''; }
if(defined($opts{'t'})){ $user_codon_table = $opts{'t'}; }
else{ $user_codon_table = ''; }
if(defined($opts{'S'})){ $skip_muscle = 1; }
else{ $skip_muscle = 0; }
if(defined($opts{'A'})){ $print_amplicons = 1; }
else{ $print_amplicons = 0; }
if(defined($opts{'F'})){ $force_all_tables = 1; }
else{ $force_all_tables = 0; }
if(defined($opts{'u'})){ $only_unique_seqs = 1; }
else{ $only_unique_seqs = 0; }
if(defined($opts{'s'})){ $server_output = 1; }
else{ $server_output = 0; }
if(defined($opts{'T'})){ $codehopT = $opts{'T'}; }
if(defined($opts{'r'})){ $amplicon_length_range = $opts{'r'}; }
if(defined($opts{'c'})){ $calc_CUT = 1; }
else{ $calc_CUT = 0; }
if(defined($opts{'P'})){ $user_CUT_path = $opts{'P'}; }
else{ $user_CUT_path = '';  }

print "# $progname -n $DNAfile -a $AAfile -t $user_codon_table -e $exclude_file -A $print_amplicons -F $force_all_tables -u $only_unique_seqs -T $codehopT -s $server_output -r $amplicon_length_range -c $calc_CUT -S $skip_muscle -P $user_CUT_path\n\n";

######################################################

my ($seq,$ampl,$isSTOPcodon,$lastcodon,$header,$exclterm,@trash,$CUT_path);

if($user_CUT_path){ $CUT_path = $user_CUT_path }
else{ $CUT_path = $ENV{'BLIMPS_DOCS'};  }

# -1) read exclusion file that contains taxa that should be removed from the analysis ###############
my @EXCLIST;
if($exclude_file)
{
	print "# Reading exclude list...\n";

	open(EXCLIST,$exclude_file) || die "# $progname : cannot read $exclude_file\n";
	while(<EXCLIST>)
	{
		next if(/^#|^$/);
		chomp;
		push(@EXCLIST,$_);
		print "# to be excluded: $EXCLIST[$#EXCLIST]\n"; 
	}
	close(EXCLIST);
}


# 0) make full list of taxa to retrieve CUTs before excluding any redundant/to_be_excluded sequences 
if($user_codon_table eq '')
{
	my $uniqueid = substr(md5_hex(time().$$),0,8);
	my %tmpfasta = read_FASTA_sequence($DNAfile);
	my %codon_tables = select_codon_tables( \%tmpfasta );
	
	if(!%codon_tables)
        {
                # 0.1) agregar tabla universal
                print "# WARNING : failed parsing taxa in FASTA headers, using a set of representative codon usage tables\n";
                %codon_tables = split(/,/,$ENV{"GB_UNIVERSAL_CODON_TABLES"});
        }

	if($calc_CUT)
        {
		# 0.0) calculate codon usage frequencies from input sequences
		$CUT_name  = $uniqueid  . ":input_derived_" . $uniqueid;
		$CUT_spsum = $ENV{'GB_SPSUM_DIR'}.'/'.$CUT_name;
		$CUT_table = $CUT_path.'/'.$CUT_name . '.codon.use';
	
		extract_codon_freqs_from_fasta2spsumfile(\%tmpfasta,$CUT_name,$CUT_spsum);

		if(-s $CUT_spsum)
		{
			$codon_tables{$CUT_name} = $CUT_name;
			push(@trash,$CUT_spsum,$CUT_path.'/input_derived_'. $uniqueid . '.codon.use');
		}
		else{ print "# $0 : cannot create input codon table\n" }
        }

	foreach my $table (keys(%codon_tables))
	{
		my $tablename = $table; #print "#$tablename\n";
		$tablename = (split(/\w:/,$tablename))[1];
		$tablename =~ s/\s+|\/|\(|\)/\_/g; # remove chars that cause trouble in the shell
		$tablename =~ s/\.//g;
		
		$tablename .= ".codon.use";
			
		if(-s "$ENV{'GB_SPSUM_DIR'}/$codon_tables{$table}")
		{
			if(!-s $CUT_path.'/'.$tablename)
			{
				open(CODUSE,"$ENV{'BLIMPS_DIR'}/bin/linux/coduse $ENV{'GB_SPSUM_DIR'}/$codon_tables{$table} \"$table\" $CUT_path/$tablename 2>&1|") 
				|| die "# $0 : cannot run $ENV{'BLIMPS_DIR'}/bin/linux/coduse $ENV{'GB_SPSUM_DIR'}/$codon_tables{$table} \"$table\" $CUT_path/$tablename\n";	
				while(<CODUSE>){} #print;}
				close(CODUSE);
			}

			# version antigua
			#sleep(1);system("ls $tablename > /dev/null "); # refresh NFS		
			#if(-s $tablename){ system("mv -f $tablename $CUT_path"); }
			#else{ print "# $0 : cannot create codon usage table $tablename ($table)\n"; }
			
			if(!-s $CUT_path.'/'.$tablename)
			{ 
				print "# $0 : cannot create codon usage table $CUT_path/$tablename ($table)\n"; 
			}
			#comentado para run_Fasta2primers, para evitar que un trabajo borre la tabla que lee otro
			#elsif($user_CUT_path) # reomve CUTs if -P was used in cut_amplicons.pl
			#{
			#	push(@trash,$CUT_path.'/'.$tablename);
			#}

			$user_codon_table .= "$tablename,"; #print "#$table#$tablename\n";
		}
		else{ print "# $progname : cannot locate $ENV{'GB_SPSUM_DIR'}/$codon_tables{$table}, skip it\n"; }
	}
}

# 1) remove identical and sequences included in @EXCLIST #########################################

my %fullDNAfasta = read_FASTA_sequence($DNAfile, 0, 0, $only_unique_seqs); # 0 ,0 = remove gaps, skip bad CDSs
my %fullAAfasta = read_FASTA_sequence($AAfile);
my (%DNAfasta,%AAfasta);

SEQS: foreach $seq (sort {$a<=>$b} (keys(%fullDNAfasta)))
{
	if(@EXCLIST && $fullDNAfasta{$seq}{'NAME'} =~ /\[((\S+\s?){1,4})/)   
	{ 
		$header = $1;  #print ">$header<\n";
		foreach $exclterm (@EXCLIST){ next SEQS if($header =~ /$exclterm/); }
	}	

	# 1.1) this sequence should remain
	$DNAfasta{$seq} = $fullDNAfasta{$seq};
	$AAfasta{$seq} = $fullAAfasta{$seq};

}

if(scalar(keys(%DNAfasta)) != scalar(keys(%AAfasta)))
{
	die "# $progname : input FASTA files contain different numbers of sequences\n";	
}
else
{ 
	if(scalar(keys(%DNAfasta)) < 2)
	{
		die "# $progname : need at least 2 sequences to run this program...\n";
	}
	else { print "# $progname : number of sequences considered = ".scalar(keys(%DNAfasta))."\n"; }
}



# 2) check sequence quality and length ###########################################################
$DNAfile = reverse $DNAfile;
$DNAfile = (split(/\//,$DNAfile))[0];
$DNAfile = reverse $DNAfile;

$AAfile = reverse $AAfile;
$AAfile = (split(/\//,$AAfile))[0];
$AAfile = reverse $AAfile;

my $myDNAfile = "my_" . (split(/\./,$DNAfile))[0] . ".fna";
my $myAAfile  = "my_" . (split(/\./,$AAfile))[0]  . ".faa";

open(MYDNA,">$myDNAfile") || die "# $progname : cannot write to $myDNAfile\n";
open(MYAA,">$myAAfile") || die "# $progname : cannot write to $myAAfile\n";

print "# FASTA header equivalences (DNA):\n";

foreach $seq (sort {$a<=>$b} (keys(%DNAfasta)))
{	
	$DNAfasta{$seq}{'NAME'} =~ s/>//;
	
	print "# $seq => $DNAfasta{$seq}{'NAME'}";

	$DNAfasta{$seq}{'SEQ'} = lc($DNAfasta{$seq}{'SEQ'});	

	$isSTOPcodon = 0;
	$lastcodon = substr($DNAfasta{$seq}{'SEQ'},length($DNAfasta{$seq}{'SEQ'})-3,3);
	if($lastcodon eq 'tga' || $lastcodon eq 'taa' || $lastcodon eq 'tag'){ $isSTOPcodon = 1 } #print "#$lastcodon#$isSTOPcodon#\n";

	# remove STOP codon if found
	if($isSTOPcodon){ $DNAfasta{$seq}{'SEQ'} = substr($DNAfasta{$seq}{'SEQ'},0,-3); }

	if($server_output)
	{ 
		if(length($DNAfasta{$seq}{'NAME'}) > $MAXSERVERHEADER)
		{
			# no cambia nada de momento Nov2008, complica las cosas despues
			print MYDNA ">" . $DNAfasta{$seq}{'NAME'};
		}
		else{ print MYDNA ">" . $DNAfasta{$seq}{'NAME'}; }
	}
	else{ print MYDNA ">$seq " . $DNAfasta{$seq}{'NAME'}; }
	print MYDNA $DNAfasta{$seq}{'SEQ'},"\n"; 
}

print "\n# FASTA header equivalences (AA):\n";

foreach $seq (sort {$a<=>$b} (keys(%AAfasta)))
{
	$AAfasta{$seq}{'NAME'} =~ s/>//;
	$AAfasta{$seq}{'SEQ'} =~ s/\*|#|\+//; ## OPA,OCH,AMB in SEQ.pm: not required if STOP codons are present in nt sequence
	
	print "# $seq => $AAfasta{$seq}{'NAME'}";
	
	if($server_output)
	{ 
		if(length($AAfasta{$seq}{'NAME'}) > $MAXSERVERHEADER)
		{ 
			# no cambia nada de momento Nov2008, complica las cosas despues
			print MYAA ">" . $AAfasta{$seq}{'NAME'};
		}
		else{ print MYAA ">" . $AAfasta{$seq}{'NAME'}; }
	}
	else{ print MYAA ">$seq " . $AAfasta{$seq}{'NAME'}; }
	print MYAA $AAfasta{$seq}{'SEQ'},"\n";
}

close(MYDNA);
close(MYAA);
print "\n";

foreach $seq (sort {$a<=>$b} (keys(%AAfasta)))
{
	if(length($AAfasta{$seq}{'SEQ'})*3 != length($DNAfasta{$seq}{'SEQ'}))
	{
		print "# ERROR: sequence $seq : amino acid and nt lengths do not match (".(length($AAfasta{$seq}{'SEQ'})*3).
		") and (".length($DNAfasta{$seq}{'SEQ'}).") ...\n";
		die "#$AAfasta{$seq}{'SEQ'}\n$DNAfasta{$seq}{'SEQ'}\n"; 
	}
}

print "# reference protein coordinates: 1 .. ".length($AAfasta{(sort {$a<=>$b} (keys(%AAfasta)))[0]}{'SEQ'})."\n\n";

if($skip_muscle)
{
        if(!check_aligned_FASTA_sequences(%AAfasta)){ $skip_muscle = 0; print "# input sequences are not aligned\n"; }
        else{ print "# input sequences seem to be aligned\n"; }
}

# 3) run muscle to obtain multiple alignment of protein sequences (if required)
my $AA_aln_file;
if($skip_muscle){ $AA_aln_file = $myAAfile }
else
{
	$AA_aln_file = run_muscle($myAAfile,$REMOVETMPFILES); 
	if($AA_aln_file eq $ERROR)
	{
	       die "# $progname : sorry, muscle alignment failed\n";
	}
}

# 4) run mablock to obtain AA sequence blocks in BLOCKS format for later use in codehop
my $BLOCKS_file = run_mablock($AA_aln_file,$REMOVETMPFILES); 
if($BLOCKS_file eq $ERROR)
{
	die "# $progname : sorry, run_mablock failed\n";
} 
push(@trash,$BLOCKS_file);


# 5) get a nucleotide alignment based on the one produced by muscle
my $DNA_aln_file; 
if($skip_muscle){ $DNA_aln_file = $myDNAfile }
else
{
	$DNA_aln_file = convert_AAaln2NTaln($myDNAfile,$AA_aln_file);
	if($DNA_aln_file eq $ERROR)
	{
       		die "# $progname : sorry, convert_AAaln2NTaln failed\n";
	}
}

# 6) obtain codehops (COnsensus DEgenerate Hybrid Oligonucleotide Primers) from blocks
my @CODEHOPs_files = run_codehop($BLOCKS_file,\%AAfasta,$user_codon_table,$codehopT,$user_CUT_path); 
if(!@CODEHOPs_files)
{
	die "# $progname : sorry, run_codehop failed\n";
}
#push(@trash,@CODEHOPs_files);

# 7) cut amplicons for each codehop
my ($codehop,$n_of_primers,$n_of_failures,%amp_coordinates) = ('',0,0);
foreach $codehop (@CODEHOPs_files)
{
	my %amplicons;

	print "\n\n> $codehop\n";
	
	# 7.1) parse mablocks & codehop output
	if($amplicon_length_range){ %amplicons = select_amplicons($AA_aln_file,$BLOCKS_file,$codehop,\%amp_coordinates,$amplicon_length_range); }
	else{ %amplicons = select_amplicons($AA_aln_file,$BLOCKS_file,$codehop,\%amp_coordinates); }
	
	#print "amps = ".scalar(keys(%amplicons))."\n";

	# 7.2) cut FASTA alignments in amplicons
	my ($ref_primers_filenames, $ref_amplicon_filenames) = 
		cut_amplicon_aln_files($DNA_aln_file,$AA_aln_file,$codehop,$print_amplicons,%amplicons);

	foreach $ampl (@{$ref_primers_filenames})
	{	
		print "- $ampl\n";
		$n_of_primers++;
	}

	if($print_amplicons)
	{
		foreach $ampl (@{$ref_amplicon_filenames})
		{	
			print "+ $ampl\n";
		}
	}
	
	# 7.3) check convergence
	if(scalar(@{$ref_primers_filenames}) == 0)
	{
		$n_of_failures++;
		last if(!$force_all_tables && $n_of_primers && $n_of_failures > 2);	
	}
}

# 8) remove unwanted outfiles
unlink(@trash);
 

