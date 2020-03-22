#!/usr/bin/perl -w
# Bruno Contreras-Moreira, Pablo Vinuesa, ...
# 2005-7 UNAM, Mexico

# reads ouput from Fasta2primers and an include list of species in order to correct primers
# designed only for some species

use strict;
use FindBin '$Bin';
use lib "$Bin";
use lib "$Bin/bioperl-1.5.2_102/";
use Getopt::Std;
use File::Basename;
use phyTools;
set_phyTools_env();

my $progname = $0;
my $REMOVETMPFILES = 1;

######################################################

my (%opts,$f2pfile,$include_file,$exclude_file,$f2pDIR,$outDIR,$verbose);

getopts('hvp:e:i:d:o:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0)) 
{ 
	print   "usage: $progname [options]\n";
	print   "-h \t this message\n";
	print   "-p \t output file from Fasta2primers.pl\n";
	print   "-d \t directory with output files from Fasta2primers.pl     (optional, overrides -p)\n";
	print   "-o \t full path to output directory                         (optional)\n";
	print   "-v \t be verbose                                            (optional)\n";
	print   "-i \t include list                                          (optional, use either -e or -i)\n";
	print   "-e \t exclude list                                          (optional, use either -e or -i)\n";
	exit; 
}

if(defined($opts{'d'})){ $f2pDIR = $opts{'d'}."/"; }
else{ $f2pDIR = ''; }

if(defined($opts{'p'})){ $f2pfile = $opts{'p'}; }
else
{
	if(!defined($opts{'d'})){ die "# $0 : need a valid output file from Fasta2primers.pl\n" }
	else{ $f2pfile = '' }
}

if(defined($opts{'o'})){ $outDIR = $opts{'o'}."/"; }
else{ $outDIR = '' }
if(defined($opts{'v'})){ $verbose = 1 }
else{ $verbose = 0 }
if(defined($opts{'e'})){ $exclude_file = $opts{'e'}; }
else{ $exclude_file = ''; }
if(defined($opts{'i'})){ $include_file = $opts{'i'}; }
else{ $include_file = ''; }
if($include_file eq '' && $exclude_file eq '')
{
	die "# $0 : need either -e or -i in order to correct primers...\n";
}
if($include_file ne '' && $exclude_file ne '')
{
	die "# $0 : need only -e or -i in order to correct primers...\n";
}

print "# $progname -d $f2pDIR -p $f2pfile -i $include_file -e $exclude_file -o $outDIR -v $verbose\n\n";

######################################################

my ($amp,$term,@trash,@inlist,@exlist,$readOK,$sequenceOK,$excludeOK,@primer_block,@excluded_primer_block,$NC,$ampname,$ampnumber);
my ($old_codehdeg,$old_relaxdeg,$old_degendeg,$new_codehdeg,$new_relaxdeg,$new_degendeg,$cutable);
my ($primersblock,$new_primersfile,$new_aamp_file,$new_dnaamp_file,@f2pfiles,$corrected_file,$tgzfile);

# 0) read inclusion file that contains taxa that should be taken ###############
if($include_file)
{
	print "# Reading $include_file ...\n";

	open(EXCLIST,$include_file) || die "# $progname : cannot read $include_file\n";
	while(<EXCLIST>)
	{
		next if(/^#|^$/);
		$_ =~ s/[\r|\n]//g;
		$_ =~ s/\[//g;
		if(/\[(.*?)\]/){ push(@inlist,$1); }
		elsif(/\[/){ push(@inlist,(split(/\[/,$_))[1]) } 
		else{ push(@inlist,$_); }
		print "# to be included: $inlist[$#inlist]\n";
	}
	close(EXCLIST);
}

# 1) read exclusion file that contains taxa that should be removed from the analysis ###############
if($exclude_file)
{
	print "# Reading $exclude_file ...\n";

	open(EXCLIST,$exclude_file) || die "# $progname : cannot read $exclude_file\n";
	while(<EXCLIST>)
	{
		next if(/^#|^$/);
		$_ =~ s/[\r|\n]//g;
		if($_ =~ m/\[(.*?)\]/){ push(@exlist,$1); }
		elsif(/\[/){ push(@exlist,(split(/\[/,$_))[1]) }
		else{ push(@exlist,$_); }
		print "# to be excluded: $exlist[$#exlist]\n";
	}
	close(EXCLIST);
}

# 2) parse f2p files  ##############################################################################

if($f2pDIR)
{
	opendir(DIR,$f2pDIR) || die "# $0 : cannot list $f2pDIR\n";
	@f2pfiles = readdir(DIR);
	closedir(DIR);
}
else{ push(@f2pfiles,$f2pfile) }

print "\n\n# parsing Fasta2primers output...\n\n";

foreach $f2pfile (@f2pfiles)
{

$f2pfile = $f2pDIR . $f2pfile; 	

if(!$verbose)
{
	$corrected_file = $outDIR . basename($f2pfile) . ".corr"; #$corrected_file = $f2pfile . "___".$exclude_file.$include_file.".corr";
	open(CORR,">$corrected_file") || die "# $0 : cannot create $corrected_file\n";
	print "# $f2pfile => $corrected_file\n";
}

$readOK=0;
open(F2POUT,$f2pfile) || die "# $0 : cannot read $f2pfile\n";
while(<F2POUT>)
{
	if(/^## Amplicon (\S+) codon_usage_table = (\S+)/)
	{ 
		$ampname = $1;
		$cutable = $2;
		if($ampname =~ /aln_/)
		{
			$ampnumber = (split(/aln_/,$ampname))[1];
			$ampname = (split(/aln_/,$ampname))[0] . "aln"; #die "$ampname\n";
		}
		else{ $ampnumber = $ampname; } # for compatibility
		$amp = $_;
		$tgzfile = '';
	}
	elsif(/^# (\S+) redundancy = (\d+\.\d+)/) #176299:Agrobacterium_tumefaciens_str_C58 redundancy = 0.98
	{
		if(!$verbose){ print CORR "# rank_codon_usage_tables : ".$1." ".$2."\n"; }
	}
	elsif(/\w+ 5'\->3' (\w)/)
	{ 
		$readOK = 1; 
		$NC = $1; 
		push(@primer_block,$_);
		next;
	}
	elsif(/^$/ && @primer_block)  # correct previous block
	{ 
		#my @corrected_primers = check_primer_match($ampname,'',$NC,@primer_block);
		my @corrected_primers = check_primer_match_excluded($ampname,'',$NC,\@primer_block,\@excluded_primer_block);
		
		if($NC eq 'N')
		{ 
			if($ampname eq $ampnumber){ $primersblock = "\n## Amplicon $ampname :\n"; }
			else{ $primersblock = "\n## Amplicon $ampname $ampnumber:\n"; } 
		}
		
		foreach my $line (@corrected_primers)
		{
			if($line =~ /(\w+) codeh_corr/)
			{
				$new_codehdeg = calc_oligo_degeneracy($1);
				chomp $line;
				$primersblock .= $line . " ( $old_codehdeg -> $new_codehdeg )\n"; 
			}
			elsif($line =~ /(\w+) relax_corr/)
			{
				$new_relaxdeg = calc_oligo_degeneracy($1);
				chomp $line;
				$primersblock .= $line . " ( $old_relaxdeg -> $new_relaxdeg )\n"; 
			}
			elsif($line =~ /(\w+) degen_corr/)
			{	
				$new_degendeg = calc_oligo_degeneracy($1);
				chomp $line;
				$primersblock .= $line . " ( $old_degendeg -> $new_degendeg )\n";
			}
			else{ $primersblock .= $line; }
		}
		
		$primersblock .= "\n";
		
		$readOK = 0; 
		@primer_block = @excluded_primer_block = (); 
		$old_codehdeg=$old_relaxdeg=$old_degendeg=-1;
	}
	elsif(/^# compressed outfile = (\S+)/) # si hay un tgz abrirlo, buscar los .aa_amp y .dna_amp y ajustarlos
	{
		$tgzfile = $f2pDIR . $1;
	}
	elsif(/^# end_of_amplicon/)
	{
		if($verbose){ print $primersblock; }
		else # big block
		{

		$new_primersfile = $outDIR . basename($ampname) . "_corr" . $ampnumber . "__" . $cutable . ".primers";
                open(NEWPRIMERS,">$new_primersfile") || die "# $0 : cannot create $new_primersfile\n";
                print NEWPRIMERS $primersblock;
                close(NEWPRIMERS);		

		print CORR $primersblock;
                print CORR "- $new_primersfile\n";

                if(-s $tgzfile) # unpack tgz file if required
                {
			my ($aa_ampfile,$dna_ampfile,@trash);
	
			open(TAR,"tar xvfz $tgzfile |");
			while(<TAR>)
			{	
				if(/aa_amp/){ $aa_ampfile = (split)[0];  push(@trash,$aa_ampfile); }
				elsif(/dna_amp/){ $dna_ampfile = (split)[0];  push(@trash,$dna_ampfile); }
				else{ push(@trash,(split)[0]) }
			}
			close(TAR);
			
			if($aa_ampfile)
			{ 
				$new_aamp_file = $outDIR . basename($ampname) . "_corr" . $ampnumber . "__" . $cutable . ".aa_amp";
				open(NEWAA,">$new_aamp_file") || die "# $0 : cannot create $new_aamp_file\n";
				
				my %FASTA = read_FASTA_sequence($aa_ampfile);
				
				foreach my $seq (sort {$a<=>$b} (keys(%FASTA)))
				{
					$sequenceOK = 0;
	
					if($include_file)
					{
						foreach $term (@inlist){ if($FASTA{$seq}{'NAME'} =~ /$term/i){ $sequenceOK=1; last; } }
					}		
					if($exclude_file && !$sequenceOK)
					{
						$excludeOK=0;
						foreach $term (@exlist){ if($FASTA{$seq}{'NAME'}=~ /$term/i){ $excludeOK=1; last } }
						if(!$excludeOK){ $sequenceOK=1; }
					}
					next if(!$sequenceOK);
					
					print NEWAA $FASTA{$seq}{'NAME'}; 
					print NEWAA "$FASTA{$seq}{'SEQ'}\n"; 
				}
				
				close(NEWAA);
			}
			elsif($dna_ampfile)
			{
				$new_dnaamp_file = $outDIR . basename($ampname) . "_corr" . $ampnumber . "__" . $cutable . ".dna_amp";
				
				open(NEWDNA,">$new_dnaamp_file") || die "# $0 : cannot create $new_dnaamp_file\n";
				
				my %FASTA = read_FASTA_sequence($dna_ampfile);
				
				foreach my $seq (sort {$a<=>$b} (keys(%FASTA)))
				{
					$sequenceOK = 0;
	
					if($include_file)
					{
						foreach $term (@inlist){ if($FASTA{$seq}{'NAME'} =~ /$term/i){ $sequenceOK=1; last; } }
					}		
					if($exclude_file && !$sequenceOK)
					{
						$excludeOK=0;
						foreach $term (@exlist){ if($FASTA{$seq}{'NAME'}=~ /$term/i){ $excludeOK=1; last } }
						if(!$excludeOK){ $sequenceOK=1; }
					}
					next if(!$sequenceOK);
					
					print NEWDNA $FASTA{$seq}{'NAME'}; 
					print NEWDNA "$FASTA{$seq}{'SEQ'}\n"; 
				}
				
				close(NEWDNA);
			}
			
			unlink(@trash);
		}
		#else{ print "# $0 : cannot find $tgzfile\n"; } #print "+ $new_aamp_file|$new_dnaamp_file\n";
		
		} # else
	} 
		
	if($readOK)
	{
		$sequenceOK = 0;
	
		if($include_file)
		{
			foreach $term (@inlist){ if($_ =~ /$term/i){ $sequenceOK=1; last; } }
		}
		
		if($exclude_file && !$sequenceOK)
		{
			$excludeOK=0;
			#foreach $term (@exlist){ print "#$_|$term|\n"; if($_ =~ /$term/i){ print "|$term|\n"; $excludeOK=1; last } }
			foreach $term (@exlist){ if($_ =~ /$term/i){ $excludeOK=1; last } }
			if(!$excludeOK){ $sequenceOK=1; }
		}

		# primera linea de codigo desde zaragoza 16112007
		if(/^\w+ \S+/ && !/_corr/)
		{
			if($sequenceOK){ push(@primer_block,$_); }
			else{ push(@excluded_primer_block,$_); }
		}
		elsif(/^(\w+) codeh_corr/)
		{
			$old_codehdeg = calc_oligo_degeneracy($1);
		}
		elsif(/^(\w+) relax_corr/)
		{
			$old_relaxdeg = calc_oligo_degeneracy($1);
		}
		elsif(/^(\w+) degen_corr/)
		{
			$old_degendeg = calc_oligo_degeneracy($1);
		}
	}
}
close(F2POUT);

close(CORR);

}

### Amplicon my_test_cluster_0_aln_1 codon_usage_table = Mesorhizobium_loti_MAFF303099 :
#ATGACCAGTTCCTGatggtnganga 5'->3' N 35 189 (aligned residues)
#....?..?..?..?.....?..?..
#ATGACCAGTTCCTGatggtbgasga codeh_corr test_cluster_0_amp1_N35
#ATGACCAGttyctsatggtbgasga relax_corr test_cluster_0_amp1_N35
#atgaycarttyctsatggtbgasga degen_corr test_cluster_0_amp1_N35


