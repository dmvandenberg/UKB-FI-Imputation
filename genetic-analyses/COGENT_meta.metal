# === DECLARE SOME PARAMETERS ===
SCHEME STDERR
CUSTOMVARIABLE N
AVERAGEFREQ ON
MINWEIGHT 10000
TRACKPOSITIONS ON

# === DESCRIBE AND PROCESS INPUT 1 ===
CHROMOSOME CHR
POSITION POS
MARKER SNP
FREQLABEL AF1
ALLELE A1 A2
EFFECT BETA
PVALUE P
STDERR  SE
LABEL N AS N
PROCESS ./mega_average_FIS.HRC.completeBeta.fastGWA.gz

# === DESCRIBE AND PROCESS INPUT 2 ===
CHROMOSOME CHR
POSITION POS
MARKER rsID
FREQLABEL EAF
ALLELE EFFECT_ALLELE OTHER_ALLELE
EFFECT BETA
PVALUE PVAL
STDERR SE
LABEL N AS N
PROCESS ./CLEANED.COGENT.withChrPos.txt.gz

# === FINALIZE ===
OUTFILE ./FIS_Meta .txt
ANALYZE
CLEAR
QUIT  