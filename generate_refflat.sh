set -e
if [ $# -ne 1 ];
then
   echo "bash generate_refflat.sh"
   exit -1
fi 

gtf=$1
gtfToGenePred -genePredExt -geneNameAsName2 $gtf annotation.refflat.txt
paste <(cut -f 12 annotation.refflat.txt) <(cut -f 1-10 annotation.refflat.txt) > refFlat.txt


paftools.js gff2bed $gtf >  junction.bed
