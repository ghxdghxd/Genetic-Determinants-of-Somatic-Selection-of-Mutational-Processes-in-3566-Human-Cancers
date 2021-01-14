#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# GenoImpute.py
# @Author : JT Guo
# @Email  : guojt-4451@163.com
# @Date   : 2018-4-28 09:09:45

# 
# 1. genotypes phasing with shapeit
#     BRCA, 800 samples, chr1, 70868 SNPs, 62 minutes
# 2. imputation with impute2
#     BRCA, 800 samples, chr1, 1-5000000, 40 minutes
# 3. 22 chromsomes, 576(5Mb)，567*40/60*24）= 15.75 days
# 

__version__ = "%(prog)s 3.0"

import os
import tempfile
import subprocess
import sys


def get_args3():
    """Get arguments from commond line."""
    try:
        import argparse as ap
    except ImportError as imerr:
        print("\033[1;31m" + str(imerr) + " \033[0m")
        sys.exit()

    parser = ap.ArgumentParser(fromfile_prefix_chars='@',
                               description=__doc__,
                               formatter_class=ap.RawTextHelpFormatter)
    subparsers = parser.add_subparsers(
        title='Sub Command', description='', help='', dest='subcom')
    shapeit_parser = subparsers.add_parser(
        'phase', help='genotype phasing with shapeit')
    shapeit_parser.add_argument("-g", "--geno",
                                metavar="FILE",
                                dest="geno",
                                required=True,
                                help="")
    shapeit_parser.add_argument("-c", "--chrom",
                                metavar="FILE",
                                dest="chrom",
                                required=True,
                                help="")
    shapeit_parser.add_argument("-o", "--outdir",
                                metavar="DIR",
                                dest="outdir",
                                required=True,
                                help="")
    shapeit_parser.add_argument("-s", "--sample",
                                metavar="FILE",
                                dest="sample",
                                required=True,
                                help="")
    impute_parser = subparsers.add_parser(
        'impute', help='genotype imputation with impute2')
    impute_parser.add_argument("-g",
                               "--geno",
                               metavar="FILE",
                               dest="geno",
                               required=True,
                               help="")
    impute_parser.add_argument("-c",
                               "--chrom",
                               metavar="STR",
                               dest="chrom",
                               required=True,
                               help="")
    impute_parser.add_argument("-s",
                               "--start",
                               metavar="INT",
                               dest="start",
                               required=True,
                               help="")
    impute_parser.add_argument("-e",
                               "--end",
                               metavar="INT",
                               dest="end",
                               required=True,
                               help="")
    impute_parser.add_argument("-o",
                               "--outdir",
                               metavar="DIR",
                               dest="outdir",
                               required=True,
                               help="")
    parser.add_argument('-v', '--version',
                        action='version',
                        version=__version__)
    clustered_group = parser.add_argument_group("Clustered arguments")
    clustered_group.add_argument("--qsub",
                                 action="store_true",
                                 default=False,
                                 dest="qsub",
                                 help="run in cluster [False]")
    clustered_group.add_argument("--nodes",
                                 metavar="STR",
                                 dest="node",
                                 type=str,
                                 help="name of nodes (e.g: n1,n2,...)")
    clustered_group.add_argument("-p",
                                 metavar="INT",
                                 dest="process_number",
                                 type=str,
                                 default='1',
                                 help="number of process [1]")
    if len(sys.argv) == 1:
        parser.print_help()
    else:
        args = parser.parse_args()
        # check_dependencies(["group_name"])
        return args


def check_dependencies(tools):
    """Ensure required tools are present.
    """
    print("Checking required dependencies......")
    if isinstance(tools, list):
        pass
    else:
        tools = [tools]
    for t in tools:
        subp = subprocess.Popen(
            ["which", t], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        if subp.stderr.read():
            # raise OSError("\n\033[1;31m" + "OSError: \033[1;33m" + __file__ + " requires " + tools + "\033[0m\n")
            print("\033[1;31m" + "OSError: \033[1;33m" +
                  __file__ + " requires " + t + "\033[0m")
            sys.exit()
        else:
            print(subp.stdout.read().strip())


def makedir(new_dir, exist_dir=None):
    """Make a directory. If it doesn't exist, handling concurrent race conditions.
    """
    if exist_dir:
        new_dir = os.path.join(exist_dir, new_dir)
    if os.path.exists(new_dir):
        print("The " + new_dir + " is already exist")
    else:
        print("Make " + new_dir)
        os.makedirs(new_dir)


def phase(args):
    name = os.path.basename(args.geno).split(".")[0]
    map = "/share/data0/reference/Genome/ALL_1000G_phase3integrated_v3_impute/" + \
        "genetic_map_" + args.chrom + "_combined_b37.txt"
    cmd = "source /etc/profile.d/set.sh\n"
    cmd = cmd + "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/share/apps/glibc.2.14/lib\n"
    cmd = cmd + "shapeit "
    cmd = cmd + "--input-gen " + args.geno + " " + args.sample + " "
    cmd = cmd + "--input-map " + map + " "
    cmd = cmd + "--output-max " + args.outdir + "/" + \
        name + "_" + args.chrom + "_phased.haps "
    cmd = cmd + "--thread " + args.process_number
    return(cmd)


def impute(args):
    name = os.path.basename(args.geno).split(".")[0]
    name = name + "_" + args.start + "-" + args.end + "_impute.gen"
    cmd = "source /etc/profile.d/set.sh\n"
    cmd = cmd + "impute2 "
    cmd = cmd + "-known_haps_g " + args.geno
    cmd = cmd + " -h /share/data0/reference/Genome/ALL_1000G_phase3integrated_v3_impute/1000GP_Phase3_" + args.chrom + ".hap.gz "
    cmd = cmd + "-l /share/data0/reference/Genome/ALL_1000G_phase3integrated_v3_impute/1000GP_Phase3_" + \
        args.chrom + ".legend.gz "
    cmd = cmd + "-m /share/data0/reference/Genome/ALL_1000G_phase3integrated_v3_impute/genetic_map_" + \
        args.chrom + "_combined_b37.txt "
    cmd = cmd + "-int " + args.start + " " + args.end + " -Ne 20000 "
    cmd = cmd + "-o_gz -o " + args.outdir + "/" + name
    return(cmd)


def qsub(args, cmd):
    ftmp = tempfile.NamedTemporaryFile(mode='w+')
    ftmp.write("#!/bin/bash\n")
    ftmp.write("#PBS -N phase\n")
    ftmp.write("#PBS -l nodes=1:ppn=" +
               args.process_number + ",walltime=10:00:00\n")
    ftmp.write("#PBS -j oe\ncd $PBS_O_WORKDIR\n")
    ftmp.write(cmd)
    ftmp.seek(0)
    # print(ftmp.read())
    os.system("qsub " + ftmp.name)
    ftmp.close()


if __name__ == '__main__':
    args = get_args3()
    makedir(args.outdir)
    if args.subcom == "phase":
        cmd = phase(args)
    elif args.subcom == "impute":
        cmd = impute(args)
    if args.qsub:
        qsub(args, cmd)
    else:
        print(cmd)
