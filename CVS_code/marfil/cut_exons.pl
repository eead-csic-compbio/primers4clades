#! /usr/bin/perl -w

use strict;

# written from Zaragoza in April 2008
# extracts exon segments from genomic sequences that have <<coord..inates,>> in their FASTA headers
# and produces CDS sequences with exon ==boundaries(intronsize),== in header

use FindBin '$Bin';
use lib "$Bin";
use lib "$Bin/bioperl-1.5.2_102/";
use phyTools;
set_phyTools_env();

my $progname = $0;

########################################################################

if(!$ARGV[0] || $ARGV[0] !~ /\.fna$/){ die "# $progname : need an input.fna file\n"; }

my $intronfile = $ARGV[0];
my $exonfile   = $intronfile;
$exonfile =~ s/\.fna/_intronless\.fna/;

print "# $progname -i $intronfile -o $exonfile\n";

###########################################################################

# 0) parse FASTA file, detect and count present exon boundaries
my %FASTAseqs = read_FASTA_sequence( $intronfile, 0, 0, 1 ); # do not make CDSs checks
my %ExBounds  = find_exon_boundaries_FASTA_headers(\%FASTAseqs);

my $n_of_seqs_with_exons = scalar(keys(%ExBounds));
my $n_of_sequences = scalar(keys(%FASTAseqs));

print  "# number of sequences read = $n_of_sequences\n";
print  "# number of sequences with exons boundaries = $n_of_seqs_with_exons\n";

open(EXON,">$exonfile") || die "# $0 : cannot create $exonfile\n";

foreach my $seq (sort {$a<=>$b} (keys(%FASTAseqs)))
{
	my ($CDS,$boundaries,$header) = ('','','');
	
	if($ExBounds{$seq})
	{
		my ($n_of_iebs,$i,$ieb,$ilength,$iebaa) = (scalar(@{$ExBounds{$seq}}));
		for($i=0;$i<$n_of_iebs;$i++) 
		{
			$ieb = $ExBounds{$seq}->[$i];
			if($CDS && $i-1 >= 0)
			{ 
				$iebaa = length($CDS); 
				$ilength = $ExBounds{$seq}->[$i]->[0] - $ExBounds{$seq}->[$i-1]->[1] - 1;
				$boundaries .= "$iebaa($ilength),";  
			}
			$CDS .= substr($FASTAseqs{$seq}{'SEQ'},$ieb->[0]-1,$ieb->[1]-$ieb->[0]+1);	
		}
		if($boundaries){ $boundaries = "==$boundaries==" }
	}
	else
	{
		$CDS = $FASTAseqs{$seq}{'SEQ'};
	}

	$header = $FASTAseqs{$seq}{'NAME'};

	$header =~ s/[\n|\r]//g;

	print EXON "$header $boundaries\n".$CDS."\n";
}

close(EXON);
