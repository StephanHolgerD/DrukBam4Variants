from multiprocessing.managers import BaseManager
from variantfinder.ReadFiles.VcfReader import ReadVcfFile
from variantfinder.WriteFiles.VcfWriter import WriteVcf_cls
from variantfinder.WriteFiles.TsvWriter import WriteTSV_cls
from variantfinder.WriteFiles.BamWriter import WriteBam_cls

from variantfinder.ParseAlignments.VariantCounter import VariantCounter_cls
from time import time
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description='read in a vcf file and count alltAll supoorting alignemnts')
    subparsers = parser.add_subparsers(help='write results to info field of input vcf')

    vcf_parser = subparsers.add_parser("vcf")
    tsv_parser = subparsers.add_parser("tsv")
    

    requiredNamed = vcf_parser.add_argument_group('required arguments')
    requiredNamed.add_argument('-b','--bam', help='Pos. sorted and indexed bam file', required=True)
    requiredNamed.add_argument('-v','--vcf', help='vcf file with variants of interest', required=True)
    requiredNamed.add_argument('-o','--vcfout', help='vcf output', required=True)
    requiredNamed.add_argument('-f','--fasta',help='fasta file of reference',required=True)
    

    requiredNamed = tsv_parser.add_argument_group('required arguments')
    requiredNamed.add_argument('-b','--bam', help='Pos. sorted and indexed bam file', required=True)
    requiredNamed.add_argument('-v','--vcf', help='vcf file with variants of interest', required=True)
    requiredNamed.add_argument('-o','--tsvout', help='vcf output', required=True)
    requiredNamed.add_argument('-f','--fasta',help='fasta file of reference',required=True)


    optArguments = vcf_parser.add_argument_group('optional arguments')
    optArguments.add_argument('--cpu',default=1, help="number of cpu's  to run in paralell",type=int)
    optArguments.add_argument('--bamsout',default=False, help='create bams with variant supporting alignments and non supporting aligments',action='store_true')
    optArguments.add_argument('--bamsdir',default=False, help='folder for bams',type=str)
    
    
    
    
    optArguments = tsv_parser.add_argument_group('optional arguments')
    optArguments.add_argument('--cpu',default=1, help="number of cpu's  to run in paralell",type=int)
    optArguments.add_argument('--bamsout',default=False, help='create bams with variant supporting alignments and non supporting aligments',action='store_true')
    optArguments.add_argument('--bamsdir',default=False, help='folder for bams',type=str)

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args=parser.parse_args()

    
    if args.__contains__('tsvout'):
        fasta = args.fasta
        bam = args.bam
        vcf = args.vcf
        proc = args.cpu
        out = args.tsvout
        writebams = args.bamsout
        bamsdir = args.bamsdir
        x = ReadVcfFile(fasta,bam,vcf)
        VC_object = VariantCounter_cls(x,processes=proc)
        Vcf_Writer_object = WriteTSV_cls(VC_object,x,out)
        if bamsdir and not writebams:
            print('bams dir option is not required without bamsout')
        if writebams:
            WriteBam_cls(VC_object,x,out,bamsdir)
        

    if args.__contains__('vcfout'):
        fasta = args.fasta
        bam = args.bam
        vcf = args.vcf
        out = args.vcfout
        proc = args.cpu
        writebams = args.bamsout
        bamsdir = args.bamsdir
        x = ReadVcfFile(fasta,bam,vcf)
        VC_object = VariantCounter_cls(x,processes=proc)
        Vcf_Writer_object = WriteVcf_cls(VC_object,x,out)
        if bamsdir and not writebams:
            print('bams dir option is not required without bamsout')
        if writebams:
            WriteBam_cls(VC_object,x,out,bamsdir)
        
if __name__ == "__main__":
    main()

