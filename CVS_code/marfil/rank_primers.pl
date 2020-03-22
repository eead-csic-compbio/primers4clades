#!/usr/bin/perl -w
# Bruno Contreras-Moreira, Pablo Vinuesa, ...
# 2007 UNAM, Mexico
# rank_primers.pl

use strict;
use Getopt::Std;
use FindBin '$Bin';
use lib "$Bin";
use lib "$Bin/bioperl-1.5.2_102/";
use Bio::Graphics; #requires GD
use Bio::SeqFeature::Generic;
use phyTools;
set_phyTools_env();

# ordena los primers de calidad para un replicon especificado por el usuario

my $progname = $0;

my $MAXQUALITY = 100;
my $MINQUALITY = 0; # the are 11 quality checks, each takes -10
my $QUALITYSTEP = 10;
my $MAPWIDTH = 800;
my $MINPAD = 10;
 
###############################################################################

my %opts;
my ($evaluated_primers_file,$replicon,$minquality,$primersDIR,$server_output,$mapfile,$mapsize) = ('','',$MINQUALITY,'./',0,'',0);

getopts('hsr:i:d:q:m:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0))
{
        print   "usage: $progname [options]\n";
        print   "-h \t this message\n";
        print   "-i \t file.tab containing evaluated primers \n";
	print   "-q \t minimum quality estimation for primer pairs to be printed    (optional, range [best=$MAXQUALITY,worst=$MINQUALITY], default = $MINQUALITY)\n";
        #print   "-r \t print results for this replicon                            (optional)\n";
	print   "-s \t format ouput for web server                                  (optional)\n";
	print	"-m \t generate physical map of amplicons of this aminoacid length  (optional, format=PNG, example: -m 350)\n";
        print   "-d \t directory containing .primers files                          (optional, default ./)\n\n";
        exit;
}

if(defined($opts{'i'})){ $evaluated_primers_file = $opts{'i'}; }
else{ die "# $0 : need a valid file.tab containing primers evaluated by evaluate_primers.pl\n" }
if(defined($opts{'r'})){ $replicon = $opts{'r'}; }
if(defined($opts{'d'})){ $primersDIR = $opts{'d'}; }
if(defined($opts{'q'}) && $opts{'q'} >= $MINQUALITY && $opts{'q'} <= $MAXQUALITY){ $minquality = $opts{'q'}; }
if(defined($opts{'s'})){ $server_output = 1; }
if(defined($opts{'m'}))
{ 
	$mapsize = $opts{'m'}; 
	$mapfile = $evaluated_primers_file . '.png';
}

print "# $progname -i $evaluated_primers_file -q $minquality -d $primersDIR -s $server_output -m $mapsize\n\n";

###############################################################################

# 1) lee archivo de primers evaluados y los ordena por calidad
my ($printOK,$primerblock,$quality_summary,$quality_control);
my ($dum,$fdegen,$ffulldegen,$fhpinpot,$fselfpot,$rdegen,$rfulldegen,$rhpinpot);
my ($rselfpot,$crosspot,$fno_overlaps,$rno_overlaps,$dotprimersfile,$fname,$rname); 
my ($foligo,$frelax_oligo,$frelaxdegen,$roligo,$rrelax_oligo,$rrelaxdegen,$minsize,$maxsize);
my ($fTMmin,$fTMmax,$rTMmin,$rTMmax,%quality_stats,%ranked_pairs,%numbered_pairs,$quality,$ampsize);
my (%coords,$N,$C,$map,$amptrack);
my ($good_primer_pairs,$total_pairs,$pair_number) = (0,0,0);


open(TAB,$evaluated_primers_file) || die "# $0 : cannot read $evaluated_primers_file\n";
while(<TAB>)
{
	next if(/^\.primers/ || /^#/ || /^$/);	
                
	$quality_control = $quality_summary = '';
	$quality = $MAXQUALITY;
	$N = $C = 0;

	# file,fwd_primer,oligo,relax_oligo,deg,relaxdeg,fulldeg,minTm,maxTm,hpinpot,selfpot,opgbs,rev_primer,oligo,relax_oligo,
	# deg,fulldeg,minTm,maxTm,hpinpot,selfpot,opgbs,paircrosspot,minsize,maxsize
	# file.primers,lacasas_univII_cluster_0_amp1_N124,ACCAGGCGGGCACntwytggtayc,ACCAGGCGGGCacntwytggtayc,32,32,4096,64.9,
	# 70.9,0.55,0.67,0,lacasas_univII_cluster_0_amp1_C287,CGCGCGGATCcartagttryc,CGCGcgrayccartagttryc,8,32,288,62.4,69.8,
	# 0.50,0.44,0,0.40,10000000,0
		
	($dotprimersfile,$fname,$foligo,$frelax_oligo,$fdegen,$frelaxdegen,$ffulldegen,$fTMmin,$fTMmax,$fhpinpot,$fselfpot,$fno_overlaps,
	$rname,$roligo,$rrelax_oligo,$rdegen,$rrelaxdegen,$rfulldegen,$rTMmin,$rTMmax,$rhpinpot,$rselfpot,$rno_overlaps,
	$crosspot,$minsize,$maxsize) = split(/\t/); chomp $maxsize;
	
	if($fname =~ /_N(\d+)/){ $N = $1 }
	if($rname =~ /_C(\d+)/){ $C = $1 } #print "$N $C\n";

	if($frelaxdegen >= $ENV{'MAX_RELAXED_DEGEN'})
	{        
		$quality_control .= "# fwd 3' relaxdeg ($frelaxdegen >= $ENV{'MAX_RELAXED_DEGEN'})\n"; 
		$quality_stats{'MAX_RELAXED_DEGEN'}++;
		$quality -= $QUALITYSTEP;
	}
	if($ffulldegen >= $ENV{'MAX_FULLDEGEN'})
	{     
		$quality_control .= "# fwd fulldeg ($ffulldegen  >= $ENV{'MAX_FULLDEGEN'})\n"; 
		$quality_stats{'MAX_FULLDEGEN'}++;
		$quality -= $QUALITYSTEP;
	}
	if($fhpinpot >= $ENV{'MAX_HAIRPINPOT'})
	{      
		$quality_control .= "# fwd hpinpot ($fhpinpot >= $ENV{'MAX_HAIRPINPOT'})\n"; 
		$quality_stats{'MAX_HAIRPINPOT'}++;
		$quality -= $QUALITYSTEP;
	}
	if($fselfpot >= $ENV{'MAX_SELFPOT'})
	{         
		$quality_control .= "# fwd selfpot ($fselfpot >= $ENV{'MAX_SELFPOT'})\n"; 
		$quality_stats{'MAX_SELFPOT'}++;
		$quality -= $QUALITYSTEP;
	}
	if($rrelaxdegen >= $ENV{'MAX_RELAXED_DEGEN'})
	{        
		$quality_control .= "# rev 3' relaxdeg ($rrelaxdegen >= $ENV{'MAX_RELAXED_DEGEN'})\n"; 
		$quality_stats{'MAX_RELAXED_DEGEN'}++;
		$quality -= $QUALITYSTEP;
	}
	if($rfulldegen >= $ENV{'MAX_FULLDEGEN'})
	{     
		$quality_control .= "# rev fulldeg ($rfulldegen >= $ENV{'MAX_FULLDEGEN'})\n"; 
		$quality_stats{'MAX_FULLDEGEN'}++;
		$quality -= $QUALITYSTEP;
	}
	if($rhpinpot >= $ENV{'MAX_HAIRPINPOT'})
	{      
		$quality_control .= "# rev hpinpot ($rhpinpot >= $ENV{'MAX_HAIRPINPOT'})\n"; 
		$quality_stats{'MAX_HAIRPINPOT'}++;
		$quality -= $QUALITYSTEP;
	}
	if($rselfpot >= $ENV{'MAX_SELFPOT'})
	{         
		$quality_control .= "# rev selfpot ($rselfpot >= $ENV{'MAX_SELFPOT'})\n"; 
		$quality_stats{'MAX_SELFPOT'}++;
		$quality -= $QUALITYSTEP;
	}
	if($crosspot >= $ENV{'MAX_F_R_CROSSPOT'})
	{      
		$quality_control .= "# crosspot ($crosspot >= $ENV{'MAX_F_R_CROSSPOT'})\n"; 
		$quality_stats{'MAX_F_R_CROSSPOT'}++;
		$quality -= $QUALITYSTEP;
	}
	if($fno_overlaps >= $ENV{'MAX_NO_OVERLAPS'})
	{ 
		$quality_control .= "# fwd opgbs ($fno_overlaps >= $ENV{'MAX_NO_OVERLAPS'})\n"; 
		$quality_stats{'MAX_NO_OVERLAPS'}++;
		$quality -= $QUALITYSTEP;
	}
	if($rno_overlaps >= $ENV{'MAX_NO_OVERLAPS'})
	{ 
		$quality_control .= "# rev opgbs ($rno_overlaps >= $ENV{'MAX_NO_OVERLAPS'})\n"; 
		$quality_stats{'MAX_NO_OVERLAPS'}++;
		$quality -= $QUALITYSTEP;
	}
	if($quality < $MINQUALITY){ $quality = $MINQUALITY }	

	# 1.1) toma para cada primer seleccionado el alineamiento de hits de primers para el grupo filogenetico estudiado
        $primerblock = "# source = $dotprimersfile \n";
	open(DOTPRIMERS,$dotprimersfile) || die "# $0 : cannot read $dotprimersfile\n";
        while(<DOTPRIMERS>)
	{ 
		$primerblock .= $_; 
		if(/5'\-\>3' N (\d+) (\d+)/)
		{
			$ampsize = ($2 - $1 + 1) * 3;
		}
	}
        close(DOTPRIMERS);
	
	$quality_summary .= "\n# primer pair quality = $quality\%\n";
	
	if($maxsize == 0)
	{
		$quality_summary .= "# expected PCR product length (nt) = $ampsize\n";
	}
	else
	{
		$quality_summary .= "# expected PCR product length (nt) = [$minsize,$maxsize]\n";
	}	
	
	$quality_summary .= "# fwd: minTm = $fTMmin maxTm = $fTMmax\n";
	$quality_summary .= "# rev: minTm = $rTMmin maxTm = $rTMmax\n";
	if($quality_control ne ''){ $quality_summary .= "# quality warnings:\n$quality_control"; }
	$quality_summary .= "#--\n";
	
	$primerblock = $primerblock . "$quality_summary\n\n"; #print "$primerblock\n";	
	
	push(@{$ranked_pairs{$quality}},$primerblock);
	$coords{$primerblock}{'N'} = $N;
	$coords{$primerblock}{'C'} = $C; 

	if($quality_control eq ''){ $good_primer_pairs++ }
	$total_pairs++;

	#COMENTADO HASTA QUE DECIDAMOS QUE HACER CON LOS REPLICONES: SON PARA HACER MAPAS GENOMICOS DE OLIGOS!!!
        #if($replicon eq ''){ print "$primerblock\n" }
        #else
        #{
        # 	if($primerblock =~ $replicon){ print "$primerblock\n" }
        #}		
}
close(TAB);

if($mapfile)
{
	# create new panel 
      	$map = Bio::Graphics::Panel->new( -length => $mapsize,
                                       -key_style => 'between',
                                       -width  => $MAPWIDTH,
                                       -pad_top => $MINPAD,-pad_bottom => $MINPAD,-pad_left => $MINPAD,-pad_right => $MINPAD );

	# add ruler
      	$map->add_track( Bio::SeqFeature::Generic->new( -start=>1,-end=>$mapsize),
                                                         -glyph   => 'arrow',
                                                         -tick    => 1,
                                                         -fgcolor => 'black',
                                                         -double  => 1 );	

	# add amplicon track
	$amptrack = $map->add_track(
                              -glyph     => 'graded_segments',
                              -label     => 1,
                              -bgcolor   => 'green',
                              -min_score => $MINQUALITY,
                              -max_score => $MAXQUALITY);
}

foreach $quality (sort {$b<=>$a} (keys(%ranked_pairs)))
{
        last if($quality < $minquality);
        foreach $primerblock (@{$ranked_pairs{$quality}})
        {
		$pair_number++;

		if($mapfile)
		{
			my $amplicon = Bio::SeqFeature::Generic->new(
                                              -display_name => "$pair_number (quality=$quality\%)",
                                              -score        => $quality,
                                              -start        => $coords{$primerblock}{'N'},
                                              -end          => $coords{$primerblock}{'C'});
  			$amptrack->add_feature($amplicon);
		}

                $primerblock =~ s/## Amplicon/## Amplicon $pair_number/;
		push(@{$numbered_pairs{$quality}},$primerblock);
        }
}

if($mapfile)
{
	open(PNG,">$mapfile") || die "# $0 : cannot write to $mapfile\n";
	print PNG $map->png();
	close(PNG);
	print "# mapfile : $mapfile\n\n";
}

# 2) print ranked pairs of primers ##########################################

foreach $quality (sort {$b<=>$a} (keys(%numbered_pairs)))
{
	last if($quality < $minquality);
	foreach $primerblock (@{$numbered_pairs{$quality}})
	{
		print $primerblock;
	}
}

# 3) print overall stats #####################################################

print "\n\n#################################################################\n"; 

if(!$server_output)
{
	print "## overall quality stats:\n";
	print "# good/bad pairs ratio = $good_primer_pairs/$total_pairs\n";
	print "# individual quality checks stats:\n";

	foreach $quality_control (keys(%quality_stats)){ print "# $quality_control = $quality_stats{$quality_control}\n"; }
}

my $bottom_legend = 
"\n\n## abbreviations:\n".
"# quality     = quality estimation [best = $MAXQUALITY, worst = $MINQUALITY] %\n".
"# crosspot    = potential of cross-hybridization [0-1]\n".
"# relaxdeg    = relaxed (3' extended segment) degeneracy\n".
"# fulldeg     = full primer degeneracy \n".
"# minTm       = minimum Tm for the pool of relaxed primers\n".
"# maxTm       = maximum Tm for the pool of relaxed primers\n".
"# hpinpot     = potential of primer hairpin [0-1]\n".
"# selfpot     = potential of primer self-priming [0-1]\n";
if($fno_overlaps||$rno_overlaps){ $bottom_legend .= "# opgbs       = other potential genomic binding sites for\n#               the oligo as found with MEGABLAST\n"; }
$bottom_legend .= "# aLRT        = approximate likelihood-ratio test [0-1],\n#               overall confidence on trees built using this amplicon\n"; 
print $bottom_legend;	

