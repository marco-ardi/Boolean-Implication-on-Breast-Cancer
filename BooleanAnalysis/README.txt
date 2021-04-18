Java implementation of StepMiner and BooleanNet.

Author: Giovanni Micale


Requirements:

- Java: https://www.oracle.com/it/java/technologies/javase-downloads.html


References:

- StepMiner:
Debashis Sahoo, David L. Dill, Robert Tibshirani and Sylvia K. Plevritis: "Extracting binary signals from microarray time-course data", 
Nucleic Acids Research, 2007, Vol. 35, No. 11, pp. 3705–3712

- BooleanNet:
Debashis Sahoo, David L. Dill, Andrew J. Gentles, Robert Tibshirani and Sylvia K. Plevritis: "Boolean implication networks derived from large scale, whole genome microarray datasets",
Genome Biology 2008, 9:R157 (doi:10.1186/gb-2008-9-10-r157)


There are two java programs included:

1) StepMiner: run StepMiner algorithm to fit a step function to a series of gene expression data and produce a (quasi-)boolean expression matrix with low (-1), intermediate (0) and high values (1).
2) BooleanNet: find a set of boolean implications between two genes.


BEFORE STARTING: within the folder containing ".java" sources compile java programs

javac -cp . *.java



---  StepMiner ----

Usage: java StepMiner -i <expressionFile> [-d <delta> -o <outputFile>]

REQUIRED PARAMETERS:
-i      Expression matrix file

OPTIONAL PARAMETERS:
-d      Delta threshold around center for low and high values (default=0.5)
-o      Output file where discretized expression matrix will be saved (default=expr.txt)

Delta is a threshold for distinguishing among low, intermediate and high values.
More specifically, if F is a step function defined as:

F(x)=1  iff x>=k
F(x)=-1  iff x<k

Then expression value "e" is set as 
- a low value iff e < k-delta
- an intermediate value iff k-delta <= e <= k+delta
- a high value iff e > k+delta


INPUT:

The input file is an expression matrix with the following format:
- The first row is a list of samples separated by tab ("\t") character
- All next rows contain the name of a gene followed by a list of expression values in the order by which samples have been reported in the first line

Example:

Calu3_totalRNA.S1.4h.A	Calu3_totalRNA.S1.4h.B	Calu3_totalRNA.S1.12h.A	Calu3_totalRNA.S1.12h.B
ENSG00000000003.14	8.9885012255632	9.13263728922046	8.64737959802135	8.7428564738853
ENSG00000000005.6	0	0	1.92445504558834	0.804823140670435	0
ENSG00000000419.12	8.9126319176414	8.85984965190772	8.83307046689681	8.58214897106431


OUTPUTS:

The output file is a (quasi)-boolean expression matrix which has the same structure of the input file and three kinds of expression values:
- "-1" for low values
- "0" for intermediate values
- "1" for high values

Example:

Calu3_totalRNA.S1.4h.A	Calu3_totalRNA.S1.4h.B	Calu3_totalRNA.S1.12h.A	Calu3_totalRNA.S1.12h.B
ENSG00000000003.14	-1	0	0	1
ENSG00000000005.6	-1	-1	1	-1
ENSG00000000419.12	0	0	0	0

By default, the output file is saved in the same folder as "expr.txt"


EXAMPLES:

java StepMiner -i Example/bladder_expr.txt -o Example/bladder_expr_disc.txt



--- BooleanNet ---

Usage: java BooleanNet -i <discretizedMatrixFile> [-g1 <geneA> -g2 <geneB> -rt <relType> -s <statisticThreshold> -p <pvalThreshold> -o <outputFile>]

REQUIRED PARAMETERS:
-i      Discretized expression matrix file

OPTIONAL PARAMETERS:
-g1    First gene (default=all genes)
-g2    Second gene (default=all genes)
-rt     Relationship type (default='any', possible values='any','low-low','low-high','high-low','high-high','equivalent','opposite')
-s      Minimum value of the statistic S of an implication to be considered significant (default=6.0)
-p      Maximum p-value P of an implication to be considered significant (default=0.01)
-o      Output file where implications found will be saved (default=implications.txt)

Extract boolean implications for pairs of genes (e.g. g1 low => g2 high)

The parameter "-rt" determines the type of boolean implications to look for.
Type of implications:
- "low-high": implications of form "g1 low => g2 high"
- "low-low": implications of form "g1 low => g2 low"
- "high-high": implications of form "g1 high => g2 high"
- "high-low": implications of form "g1 high => g2 low"
- "equivalent": implications of form "g1 low => g2 low AND g1 high => g2 high"
- "opposite" : implications of form "g1 low => g2 high AND g1 high => g2 low"

If "-g1" and "-g2" are not specified, then implications for all pairs of genes are searched.
If only "-g1" is specified, then only implications involving gene g1 are searched.
If both "-g1" and "-g2" are specified, then only implications between g1 and g2 are searched.

"-s" and "-p" parameters determine how significant must be the relationship found among genes.
The higher is the value of "-s" and the lower is the value of "-p" the more significant is the implication.


INPUT:

The input file is a discretized expression matrix file, like the one produced by StepMiner java program.


OUTPUT:

The output is a list of boolean implications.

Example:

BSG low => SLC39A13 low
BSG low => RDH10 low
BSG low => ABHD11 low
BSG low => IRF6 low
BSG low => CYP11B1 low

By default, the output file is saved in the same folder as "implications.txt"



EXAMPLES:

#Extract implications for gene A1BG of type high-low and save results to file "bladder_expr_A1BG_hl.txt"
java BooleanNet -i Example/bladder_expr_disc.txt -g1 A1BG -rt high-low -o Example/tcga_A1BG_hl.txt

#Extract implications for gene A1BG of any type having statistic S>=7.0 and p-value P<=0.1 and save results to file "bladder_expr_A1BG_any.txt"
java BooleanNet -i Example/bladder_expr_disc.txt -g1 A1BG -s 7.0 -p 0.1 -o Example/bladder_expr_A1BG_any.txt

#Extract implications for all genes of type equivalent
java BooleanNet -i Example/bladder_expr_disc.txt -rt equivalent

#Extract implications for all genes of all types and save results to "Example/bladder_expr_all_any.txt"
java BooleanNet -i Example/bladder_expr_disc.txt -o Example/bladder_expr_all_any.txt

