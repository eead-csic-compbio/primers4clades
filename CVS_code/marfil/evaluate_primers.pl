#!/usr/bin/perl -w
# Bruno Contreras-Moreira, Pablo Vinuesa, ...
# 2006-7 UNAM, Mexico
# evaluate_primers.pl

use strict;
use Getopt::Std;
use FindBin '$Bin';
use lib "$Bin";
use lib "$Bin/bioperl-1.5.2_102/";
use phyTools;
use File::Basename;
set_phyTools_env();

my $progname = $0;
my $REMOVETMPFILES = 1;

#########################################################

my %opts;
my ($INP_primerfile,$INP_tabular,$INP_megablast,$INP_csv,$INP_primers4thermo) = ('',0,0,0,'c');

getopts('htmcT:f:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0)) 
{ 
	print   "usage: $progname [options]\n";
	print   "-h \t this message\n";
	print   "-f \t file containing *.primers (cut_amplicons.pl) filenames\n";
	print   "-t \t .tab output file                                              (optional)\n";
	print   "-c \t .csv output file                                              (optional)\n";
	print   "-T \t pair of primers to be thermo-evaluated (c=corr,r=relax,d=deg) (optional, default=$INP_primers4thermo)\n";
	print   "-m \t run megablast searching for genomic matches                   (optional)\n\n";
	exit; 
}

if(defined($opts{'f'})){ $INP_primerfile = $opts{'f'}; }
else{ die "# $0 : need a valid file containing primed genes produced by parse_AU_result.pl\n" }
if(defined($opts{'t'})){ $INP_tabular = $opts{'t'}; }
if(defined($opts{'c'})){ $INP_csv = $opts{'c'}; }
if(defined($opts{'m'})){ $INP_megablast = 1 }
if(defined($opts{'T'}))
{ 
	$INP_primers4thermo = $opts{'T'};
	if($INP_primers4thermo !~ m/[c|r|d]/){ $INP_primers4thermo = 'c' }
}

print "# $progname -f $INP_primerfile -t $INP_tabular -c $INP_csv -m $INP_megablast -T $INP_primers4thermo\n\n";

###############################################################################

my (@primerfiles,$file,$ampname,$ampnumber,@ampnumbers,$n_of_files,$f);
my ($fdegen,$fTm,$flowTm,$fhpinpot,$fselfpot,$rdegen,$rTm,$rlowTm,$rhpinpot,$rselfpot,$crosspot);
my ($fseqfile,$rseqfile,$sp,$spDB,$ffulldegen,$rfulldegen,$righthit,$lefthit,$frelaxdegen,$rrelaxdegen);
my ($fno_overlaps,$rno_overlaps,$no_overlaps,$replicon,$evalue,$percentid);
my ($freplicon,$fright,$fleft,$fevalue,$fpercentid,$fname,$ampsize);
my ($rreplicon,$rright,$rleft,$revalue,$rpercentid,$rname,$selectedOK);
my ($tabfilename,$csvfilename) = ('','');

# -1) lee archivo de nombres de primer generado por parse_AU_results.pl
#open(FILE,basename($INP_primerfile)) || die "#$0: cannot open file $INP_primerfile $!\n";
open(FILE,$INP_primerfile) || die "#$0: cannot open file $INP_primerfile $!\n";
while(<FILE>)
{
	next if(/^#/);
	push(@primerfiles,(split)[0]);
	if((split)[1]){ push(@ampnumbers,(split)[1]); } #strings such as 1,3, 
}

$n_of_files = scalar(@primerfiles);

# 0) define output filenames and open files
if($INP_tabular)
{
	$tabfilename = $INP_primerfile . ".oligo_eval.tab";

	my $opgbsheader = 'opgbs';
	if(!$INP_megablast){ $opgbsheader = 'void'; }
        open(TAB,">$tabfilename") || die "# $progname : cannot write to $tabfilename\n";
	printf TAB (".primers\tfwd_primer\toligo\trelax_oligo\tdeg\trelaxdeg\tfulldeg\tminTm\tmaxTm\thpinpot\tselfpot\topgbs\trev_primer\toligo\trelax_oligo\tdeg\trelaxdeg\tfulldeg\tminTm\tmaxTm\thpinpot\tselfpot\t$opgbsheader\tpaircrosspot\tminsize\tmaxsize\n");
}
if($INP_csv)
{
	$csvfilename = $INP_primerfile . ".oligo_eval.csv";
	
	my $opgbsheader = 'opgbs';
        if(!$INP_megablast){ $opgbsheader = 'void'; }
	open(CSV,">$csvfilename") || die "# $progname : cannot write to $csvfilename\n";
	printf CSV (".primers,fwd_primer,oligo,relax_oligo,deg,relaxdeg,fulldeg,minTm,maxTm,hpinpot,selfpot,opgbs,rev_primer,oligo,relax_oligo,deg,relaxdeg,fulldeg,minTm,maxTm,hpinpot,selfpot,$opgbsheader,paircrosspot,minsize,maxsize\n");
}


for($f=0;$f<$n_of_files;$f++)
{
        my (%codeh_primers,%degen_primers,%relaxed_primers,%species,%selamps,%genomseqs,@fundegseqs,@rundegseqs);
	
	$selectedOK = 0;
	if(@ampnumbers){ foreach $ampnumber (split(/\,/,$ampnumbers[$f])){ $selamps{$ampnumber}++; $selectedOK++; }}
        
	$file = $primerfiles[$f];
	
        # 1) obtener primers para cada amplicon
        open(PRIMFILE,$file) || die "# $0 : cannot read $file\n";
        while(<PRIMFILE>)
        {
                if(/^## Amplicon \S+ (\d+):/){ $ampnumber = $1;  }
                elsif(/(\w+) codeh_corr (\S+)/ && (!$selectedOK || $selamps{$ampnumber}))
                {
                        push(@{$codeh_primers{$ampnumber}},[$1,$2]);
		}
		elsif(/(\w+) degen_corr (\S+)/ && (!$selectedOK || $selamps{$ampnumber}))
                {
                        push(@{$degen_primers{$ampnumber}},[$1,$2]);
                }
		elsif(/(\w+) relax_corr (\S+)/ && (!$selectedOK || $selamps{$ampnumber}))
                {
                        push(@{$relaxed_primers{$ampnumber}},[$1,$2]); # cool primers!
                }
		# one of two possible formats
		elsif(/^(\w+)\s+>\d+\s+\d?\s?\S+\s+\[(\S+)\]/ && (!$selectedOK || $selamps{$ampnumber}))
		{
			$sp = $2;
                        push(@{$species{$ampnumber}},$sp);
                        push(@{$genomseqs{$ampnumber}{$sp}},$1);
		}
		elsif(/^(\w+)\s+>\d+\s+\d?\s?\[(\S+)\]/ && (!$selectedOK || $selamps{$ampnumber}))
                {
			$sp = $2;
			push(@{$species{$ampnumber}},$sp);
			push(@{$genomseqs{$ampnumber}{$sp}},$1);
		}
	}
	close(PRIMFILE);

	# 2) para cada pareja de primers estimar Tm, hairpins y dimers y dianas en su genoma diana
	foreach $ampnumber (keys(%codeh_primers))
	{
		$fno_overlaps = $rno_overlaps = 0;
		my ($minampsize,$maxampsize,$primerF,$primerR) = (10000000,0,'','');

		if($INP_primers4thermo eq 'r')
		{
			$primerF = $relaxed_primers{$ampnumber}[0][0];
			$primerR = $relaxed_primers{$ampnumber}[1][0];
		}
		elsif($INP_primers4thermo eq 'd')
		{
			$primerF = $degen_primers{$ampnumber}[0][0];
                        $primerR = $degen_primers{$ampnumber}[1][0];
		}
		else
		{
			$primerF = $codeh_primers{$ampnumber}[0][0];
                        $primerR = $codeh_primers{$ampnumber}[1][0];
		}

		open(CHECKPRIMERS,"$ENV{'EXE_CHECKPRIMERS'} $primerF $primerR |") 
		|| die "# $0 : cannot run $ENV{'EXE_CHECKPRIMERS'} $primerF $primerR\n";
		while(<CHECKPRIMERS>)
		{
			if(/^# fdegen/){ $frelaxdegen = (split)[2] }
 			elsif(/^# flowTm/){ $flowTm = (split)[2] }
			elsif(/^# fhighTm/){ $fTm = (split)[2] }
			elsif(/^# fhairpin/){ $fhpinpot = (split)[3] }
			elsif(/^# fselfcompl/){ $fselfpot = (split)[3] }
			elsif(/^# rdegen/){ $rrelaxdegen = (split)[2] }
                        elsif(/^# rlowTm/){ $rlowTm = (split)[2] }			
			elsif(/^# rhighTm/){ $rTm = (split)[2] }
			elsif(/^# rhairpin/){ $rhpinpot = (split)[3] }
			elsif(/^# rselfcompl/){ $rselfpot = (split)[3] }
			elsif(/^# cross-compl/){ $crosspot = (split)[3] }
			#elsif(/^# fundegseq/){ push(@fundegseqs,">".(split)[2],"\n".(split)[2]."\n") }
                        #elsif(/^# rundegseq/){ push(@rundegseqs,">".(split)[2],"\n".(split)[2]."\n") }
		}
		close(CHECKPRIMERS);
	
		$fdegen = calc_oligo_degeneracy($codeh_primers{$ampnumber}[0][0]);
                $rdegen = calc_oligo_degeneracy($codeh_primers{$ampnumber}[1][0]);

		$frelaxdegen = calc_oligo_degeneracy($relaxed_primers{$ampnumber}[0][0]);
		$rrelaxdegen = calc_oligo_degeneracy($relaxed_primers{$ampnumber}[1][0]);

		$ffulldegen = calc_oligo_degeneracy($degen_primers{$ampnumber}[0][0]);
		$rfulldegen = calc_oligo_degeneracy($degen_primers{$ampnumber}[1][0]);

		# 2.1) busca matches en genoma de referencia

		# 2.2.1) escribe archivos FASTA con secuencia de primers
		$fseqfile = $relaxed_primers{$ampnumber}[0][1] . "_f";
		$rseqfile = $relaxed_primers{$ampnumber}[1][1] . "_r";
		open(FSEQ,">$fseqfile") || die "# $0 : cannot write to $fseqfile\n";
		foreach $sp (keys(%{$genomseqs{$ampnumber}})){ print FSEQ ">$sp\n$genomseqs{$ampnumber}{$sp}[0]\n" }
		close(FSEQ);

		open(RSEQ,">$rseqfile") || die "# $0 : cannot write to $rseqfile\n";
                foreach $sp (keys(%{$genomseqs{$ampnumber}})){ print RSEQ ">$sp\n$genomseqs{$ampnumber}{$sp}[1]\n"; }
                close(RSEQ);		
 		
		# 2.2.2) comprueba si hay hits no solapantes, posibles amplicones alternativos en el genoma $spDB
		$revalue = $fevalue = 100;
		foreach $sp (@{$species{$ampnumber}})
		{
			next if(!$INP_megablast);
			$spDB = $ENV{'BLASTNDBPATH'} . "/" . $sp . ".fna";
			if(!-e $spDB.".nsq")
			{ 	
				print "# $progname : cannot find $spDB\n";
			}
			else
			{
				($no_overlaps,$freplicon,$fleft,$fright,$fevalue,$fpercentid) = 
					oligo_megablast($spDB,$fseqfile);
				$fno_overlaps += $no_overlaps;
				
				($no_overlaps,$rreplicon,$rleft,$rright,$revalue,$rpercentid) =                                                       
					oligo_megablast($spDB,$rseqfile);  
				$rno_overlaps += $no_overlaps;
				
				$ampsize = $rright-$fleft;
				if($ampsize < 0)
				{
					$ampsize = $fright-$rleft;
				}
				
				$fname = $codeh_primers{$ampnumber}[0][1];
				$rname = $codeh_primers{$ampnumber}[1][1];	
				
				if($ampsize < $minampsize){ $minampsize = $ampsize }
				if($ampsize > $maxampsize){ $maxampsize = $ampsize }
	
				print "> $sp $freplicon $fname $fleft $fright $fevalue $fpercentid " .
                		"# $rreplicon $rname $rleft $rright $revalue $rpercentid " .
				"# $ampsize\n";
			}
		}	
		
		if($INP_tabular)
		{
			printf TAB ("%s\t%s\t%s\t%s\t%d\t%d\t%d\t%1.1f\t%1.1f\t%1.2f\t%1.2f\t%d\t%s\t%s\t%s\t%d\t%d\t%d\t%1.1f\t%1.1f\t%1.2f\t%1.2f\t%d\t%1.2f\t%d\t%d\n",
			$file,
			# fwd
			$codeh_primers{$ampnumber}[0][1],$codeh_primers{$ampnumber}[0][0],$relaxed_primers{$ampnumber}[0][0],$fdegen,$frelaxdegen,$ffulldegen,$flowTm,$fTm,$fhpinpot,$fselfpot,$fno_overlaps,
			# rev
			$codeh_primers{$ampnumber}[1][1],$codeh_primers{$ampnumber}[1][0],$relaxed_primers{$ampnumber}[1][0],$rdegen,$rrelaxdegen,$rfulldegen,$rlowTm,$rTm,$rhpinpot,$rselfpot,$rno_overlaps,
			# pair
			$crosspot,$minampsize,$maxampsize);
		}
		
		if($INP_csv)
		{
			printf CSV ("%s,%s,%s,%s,%d,%d,%d,%1.1f,%1.1f,%1.2f,%1.2f,%d,%s,%s,%s,%d,%d,%d,%1.1f,%1.1f,%1.2f,%1.2f,%d,%1.2f,%d,%d\n",
			$file,
			# fwd
			$codeh_primers{$ampnumber}[0][1],$codeh_primers{$ampnumber}[0][0],$relaxed_primers{$ampnumber}[0][0],$fdegen,$frelaxdegen,$ffulldegen,$flowTm,$fTm,$fhpinpot,$fselfpot,$fno_overlaps,
			# rev
			$codeh_primers{$ampnumber}[1][1],$codeh_primers{$ampnumber}[1][0],$relaxed_primers{$ampnumber}[1][0],$rdegen,$rrelaxdegen,$rfulldegen,$rlowTm,$rTm,$rhpinpot,$rselfpot,$rno_overlaps,
			# pair
			$crosspot,$minampsize,$maxampsize);  
		}

		print "\n$file\n";
		printf("\n## summary: crosspot=%1.2f minsize=%d maxsize=%d tabfile=%s csvfile =%s\n",$crosspot,$minampsize,$maxampsize,$tabfilename,$csvfilename);
		printf("strand\tprimer\toligo\trelax_oligo\tdeg\trelaxdeg\tfulldeg\tminTm\tmaxTm\thpinpot\tselfpot\topgbs\n");
		printf("fwd\t%s\t%s\t%s\t%d\t%d\t%d\t%1.1f\t%1.1f\t%1.2f\t%1.2f\t%d\n", 
		$codeh_primers{$ampnumber}[0][1],$codeh_primers{$ampnumber}[0][0],$relaxed_primers{$ampnumber}[0][0],$fdegen,$frelaxdegen,$ffulldegen,$flowTm,$fTm,$fhpinpot,$fselfpot,$fno_overlaps);
		printf("rev\t%s\t%s\t%s\t%d\t%d\t%d\t%1.1f\t%1.1f\t%1.2f\t%1.2f\t%d\n", 
		$codeh_primers{$ampnumber}[1][1],$codeh_primers{$ampnumber}[1][0],$relaxed_primers{$ampnumber}[1][0],$rdegen,$rrelaxdegen,$rfulldegen,$rlowTm,$rTm,$rhpinpot,$rselfpot,$rno_overlaps);

		system("rm -f $fseqfile $rseqfile");
	}
}

## print legends #####################
my $bottom_legend = 
"\n\n## abbreviations:\n".
"# .primers    = file containing this primer\n".
"# crosspot    = potential of cross-hybridization [0-1]\n".
"# minsize     = minimum expected amplicon size\n".
"# maxsize     = maximum expected amplicon size\n".
"# tabfile     = tabular output file\n".
"# csvfile     = comma-separated output file\n".
"# primer      = primer name\n".
"# oligo       = primer nucleotide sequence\n".
"# relax_oligo = relaxed primer nucleotide sequence\n".
"# deg         = primer degeneracy in 3' degenerated segment\n".
"# relaxdeg    = relaxed (3' extended segment) degeneracy\n".
"# fulldeg     = full primer degeneracy \n".
"# minTm       = minimum Tm for the pool of codehop primers\n".
"# maxTm       = maximum Tm for the pool of codehop primers\n".
"# hpinpot     = potential of primer hairpin [0-1]\n".
"# selfpot     = potential of primer self-priming [0-1]\n";
if($INP_megablast){ $bottom_legend .= "# opgbs       = other potential genomic binding sites for\n#               the oligo as found with MEGABLAST\n"; }
print $bottom_legend;		

if($INP_tabular)
{ 
	print TAB $bottom_legend;
	close(TAB); 
}

if($INP_csv)
{ 
	print CSV $bottom_legend;
	close(CSV); 
}		
		
