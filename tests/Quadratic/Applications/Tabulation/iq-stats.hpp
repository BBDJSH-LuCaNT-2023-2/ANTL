#ifndef IQ_STATS_H 
#define IQ_STATS_H 


#define cycprob         .97757481164098618137
#define gamma           0.57721566490153286061
#define fprefix         "latex_data/"


// *** 1. L(1,X) bounds ***

// temp files
#define datafile        "temp/tempdata"
enum {minL0_c, minL1_c, minL5_c, maxL0_c, maxL1_c,
      maxL5_c, minLLI0_c, minLLI1_c, minLLI5_c, 
      maxULI0_c, maxULI1_c, maxULI5_c};

// tex files
#define minL0file	"minL0.tex"
#define minL1file	"minL1.tex"
#define minL5file	"minL5.tex"
#define maxL0file	"maxL0.tex"
#define maxL1file	"maxL1.tex"
#define maxL5file	"maxL5.tex"
#define minLLI0file	"minLLI0.tex"
#define minLLI1file	"minLLI1.tex"
#define minLLI5file	"minLLI5.tex"
#define maxULI0file	"maxULI0.tex"
#define maxULI1file	"maxULI1.tex"
#define maxULI5file	"maxULI5.tex"
#define aveLfile        "aveL.tex"

// graph files
#define localLLIgraph   "localLLIgraph.dat"
#define localULIgraph   "localULIgraph.dat"



// *** 2. Divisibility of h by odd primes ***

// tex files
#define pdivsfile	"pdivs.tex"
#define pdivs0file	"pdivs0.tex"
#define pdivs1file	"pdivs1.tex"
#define pdivsperfile	"pdivs_percent.tex"
#define pdivsper0file	"pdivs_percent0.tex"
#define pdivsper1file	"pdivs_percent1.tex"
#define pdivsratfile    "pdivs_ratio.tex"
#define pdivsrat0file   "pdivs_ratio0.tex"
#define pdivsrat1file   "pdivs_ratio1.tex"
#define pdivssqrfile    "pdivs_square.tex"
#define pdivssqr0file   "pdivs_square0.tex"
#define pdivssqr1file   "pdivs_square1.tex"
#define pdivscubfile    "pdivs_cube.tex"
#define pdivscub0file   "pdivs_cube0.tex"
#define pdivscub1file   "pdivs_cube1.tex"

// graph files
#define pdivs3graph     "pdivs3graph.dat"
#define pdivs5graph     "pdivs5graph.dat"
#define pdivs7graph     "pdivs7graph.dat"


// *** 3a. Noncyclic odd parts of class group ***

// tex files
#define maxhfile        "maxh.tex"
#define maxhoddfile     "maxhodd.tex"
#define ncycfile        "noncyc.tex"
#define ncyc0file       "noncyc0.tex"
#define ncyc1file       "noncyc1.tex"

// graph files
#define ncycgraph       "noncycgraph.dat"

// *** 3b. P-rank probabilities ***

// tex files
#define prankfile       "prank.tex"
#define prank0file      "prank0.tex"
#define prank1file      "prank1.tex"
#define prankperfile    "prank_percent.tex"
#define prankper0file   "prank_percent0.tex"
#define prankper1file   "prank_percent1.tex"
#define prankratfile    "prank_ratio.tex"
#define prankrat0file   "prank_ratio0.tex"
#define prankrat1file   "prank_ratio1.tex"

// graph files
#define prank22graph      "prank22graph.dat"
#define prank32graph      "prank32graph.dat"
#define prank52graph      "prank52graph.dat"
#define prank72graph      "prank72graph.dat"
#define prank23graph      "prank23graph.dat"
#define prank33graph      "prank33graph.dat"
#define prank53graph      "prank53graph.dat"
#define prank73graph      "prank73graph.dat"
#define prank24graph      "prank24graph.dat"

// *** 4. First occurrences of p-sylow groups ***

// tex files
#define first22file        "first22rank.tex"
#define firstp2file        "firstp2rank.tex"
#define first23file        "first23rank.tex"
#define firstp3file        "firstp3rank.tex"
#define first24file        "first24rank.tex"
#define firstp4file        "firstp4rank.tex"
#define first25file        "first25rank.tex"
#define firstp5file        "firstp5rank.tex"
#define first26file        "first26rank.tex"
#define firstp6file        "firstp6rank.tex"
#define twononcycfile      "twononcyc.tex"
#define threenoncycfile    "threenoncyc.tex"
#define fournoncycfile     "fournoncyc.tex"

// *** 5. Number of generators ***

// tex files
#define maxpfile     "maxp.tex"
#define maxratfile   "maxrat.tex"
#define psplitfile   "psplit.tex"
#define kpsplitfile  "kpsplit.tex"

// graph files
#define maxpgraph    "maxpgraph.dat"
#define maxpsuccfile "maxpsucc.dat"


// *** 6. abc conjecture ***

// temp files
#define smallhfile   "smallh.dat"

#endif
