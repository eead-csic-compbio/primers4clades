# primers4clades

Primers4clades provides a fully automated pipeline for the design of PCR primers for cross-species amplification of novel sequences from metagenomic DNA or from uncharacterized organisms belonging to user-specified phylogenetic lineages. It implements an extended CODEHOP strategy based on both DNA and protein multiple alignments of coding genes and evaluates thermodynamic properties of the oligonucleotide pairs, as well as the phylogenetic information content of predicted amplicons, computed from the branch support values of maximum likelihood phylogenies. Trees displayed on screen make it easy to target primers to interactively selected clades. 

## Motivation

This repo contains the source code of the Primers4clades server, which has been running since 2009 on two mirrors: http://maya.ccg.unam.mx/primers4clades and http://floresta.eead.csic.es/primers4clades. A tutorial is available at http://maya.ccg.unam.mx/primers4clades/tutorial.html

## Credits

Primers4clades (primers for clades) is developed and maintained by Pablo Vinuesa at CCG/UNAM (Mexico) and Bruno Contreras-Moreira at EEAD/CSIC (Spain), with technical support provided by Romualdo Zayas-Laguna and Victor del Moral. This project is funded by grant DGAPA IN201806-2 from DGAPA/PAPIIT-UNAM, grant P1-60071 from CONACyT-Mexico, ARAID and by grant 200720I038 from CSIC. If you use this resource we ask you to cite the main paper describing this work:

Contreras-Moreira, B., Sachman-Ruiz, B., Figueroa-Palacios, I., and Vinuesa, P. (2009). primers4clades: a web server that uses phylogenetic trees to design lineage-specific PCR primers for metagenomic and diversity studies. Nucl.Acids Res., 37(Web Server issue):W95-W100 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2703966

The server relies on several resources: CODEHOP, MUSCLE, Codon Usage Database, TREE-PUZZLE, PhyML, PHYLIP, Amplicon, NJplot and Bioperl. You should cite the following references to give them credit:

1.  Anisimova, M. and Gascuel, O. (2006) Approximate likelihood-ratio test for branches: A fast, accurate, and powerful alternative. Syst. Biol., 55, 539-352.

2. Edgar, R.C. (2004) MUSCLE: multiple sequence alignment with high accuracy and high throughput. Nucleic Acids Res., 32, 1792-1797.

3. Felsenstein, J. (2004). PHYLIP (Phylogeny Inference Package) v3.6. Distributed by the author. Department of Genetics, University of Washington, Seattle.

4. Guindon, S. and Gascuel, O. (2003) A simple, fast, and accurate algorithm to estimate large phylogenies by maximum likelihood. Syst. Biol., 52, 696-704.

5. Jarman, S.N. (2004) Amplicon: software for designing PCR primers on aligned DNA sequences. Bioinformatics, 20, 1644-1645.

6. Rose, T.M., Henikoff, J.G. and Henikoff, S. (2003) CODEHOP (COnsensus-DEgenerate Hybrid Oligonucleotide Primer) PCR primer design. Nucleic Acids Res., 31, 3763-3766.

7. Rose, T.M., Schultz, E.R., Henikoff, J.G., Pietrokovski, S., McCallum, C.M. and Henikoff, S. (1998) Consensus-degenerate hybrid oligonucleotide primers for amplification of distantly related sequences. Nucleic Acids Res., 26, 1628-1635.

8. Schmidt, H.A., Strimmer, K., Vingron, M. and von Haeseler, A. (2002) TREE-PUZZLE: maximum likelihood phylogenetic analysis using quartets and parallel computing. Bioinformatics, 18, 502-504.

9. Stajich, J.E., Block, D., Boulez, K., Brenner, S.E., Chervitz, S.A., Dagdigian, C., Fuellen, G., Gilbert, J.G., Korf, I., Lapp, H. et al. (2002) The Bioperl toolkit: Perl modules for the life sciences. Genome Res., 12, 1611-1618.

