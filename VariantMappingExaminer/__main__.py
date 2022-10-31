from multiprocessing.managers import BaseManager
from VariantMappingExaminer.ReadFiles.VcfReader import ReadVcfFile
from VariantMappingExaminer.WriteFiles.VcfWriter import WriteVcf_cls
from VariantMappingExaminer.ParseAlignments.VariantCounter import VariantCounter_cls
from time import time
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description='read in a vcf file and count alltAll supoorting alignemnts')
    subparsers = parser.add_subparsers(help='write results to info field of input vcf')

    vcf_parser = subparsers.add_parser("vcf")

    requiredNamed = vcf_parser.add_argument_group('required arguments')
    requiredNamed.add_argument('-b','--bam', help='Pos. sorted and indexed bam file', required=True)
    requiredNamed.add_argument('-v','--vcf', help='vcf file with variants of interest', required=True)
    requiredNamed.add_argument('-f','--fasta',help='fasta file of reference',required=True)

    optArguments = vcf_parser.add_argument_group('optional arguments')
    optArguments.add_argument('--cpu',default=1, help="number of cpu's  to run in paralell",type=int)

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args=parser.parse_args()

    if args.__contains__('vcf'):
        t=time()
        fasta = args.fasta
        bam = args.bam
        vcf = args.vcf
        proc = args.cpu
        x = ReadVcfFile(fasta,bam,vcf)
        VC_object = VariantCounter_cls(x,processes=proc)
        Vcf_Writer_object = WriteVcf_cls(VC_object,x)
        print(time()-t)

if __name__ == "__main__":
    main()

