import argparse
import os
import pysam

def add_bcumi(bam,outfile):
    samfile = pysam.AlignmentFile(bam, "rb")
    outfile = pysam.AlignmentFile(outfile, "wb", template=samfile)
    for read in samfile:
        seqid=read.query_name 
        bcumi,_=seqid.split("#")
        bc,umi=bcumi.split("_")
        
        read.tags += [('U8',umi)]
        read.tags += [('BC',bc)]
        outfile.write(read)

    samfile.close()
    outfile.close()



if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Add BC and UMI By Self With Match_cell_barcode's Fastq")        
    parser.add_argument("--bam",type=str,default=None,help="bam file from minimap2 ...")
    parser.add_argument("--outfile",type=str,default=None,help="out bam file")
    args=parser.parse_args()
    
    outfile=args.outfile
    infile=args.bam
    if outfile is None:
       outfile=os.path.dirname(infile)+"/"+"addBCUMI.bam"
    add_bcumi(infile,outfile)
