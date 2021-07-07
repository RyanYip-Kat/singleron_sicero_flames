set -e
if [ $# -ne 3 ];
then
    echo "usage : bash split_fastq.sh in.fastq num_split split_dir"
    exit -1
fi

fspec=$1
num_files=$2
split_dir=$3

filename=`basename ${fspec}`
extention="${filename##*.}"
echo "file ext : $extention"

if [ ! -d $split_dir ];
then 
mkdir -p $split_dir
fi

# Work out lines per file.
if [ "$extention" == "gz" ];
then
total_lines=$(zcat $fspec | wc -l )
else
total_lines=$(wc -l <${fspec})
fi

((lines_per_file = (total_lines + num_files - 1) / num_files))

# Split the actual file, maintaining lines.
if [ "$extention" == "gz" ];
then
zcat ${fspec} | split --lines=${lines_per_file} - ${split_dir}/out.
else
split --lines=${lines_per_file} ${fspec} ${split_dir}/out.
fi
# Debug information

echo "Total lines     = ${total_lines}"
echo "Lines  per file = ${lines_per_file}"    
wc -l ${split_dir}/out.*
