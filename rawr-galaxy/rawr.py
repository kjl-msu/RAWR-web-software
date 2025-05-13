# coding:utf-8
# Available under GPLv3 license
import os, sys, errno
import warnings
import argparse
from rawr_src import sampler, supportEstimator
import shutil 

class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')
    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout
def runShellCmd(self, cmd):        
    ### Alternative subprocess running in commented code below. This is useful for developers because you will see longer subprocesses pr>
    print(cmd)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) # we're not using wait() or communicate() so th>
    sleep(1)
    rc = p.poll()
    while rc != 0:      # only allow programs to exit with valid return code or EOF 
        while True:
            line = p.stdout.readline()
            if not line:   # terminates with EOF
                break
            print(line.decode())
            sys.stdout.flush()   # may not be necessary as we are reading from stdout which prevents it from filling up to buffer size an>
        rc = p.poll()   
    assert rc == 0

def dir_path(string):
    abs_path=os.path.abspath(string)
    if os.path.isdir(abs_path):
       print('Warning: using the same folder for multiple RAWR runs could be problematic because RAxML will not run if existing logs with the same filenames are there. Strongly recommend letting RAWR create folders for you.')
    #else:
    #    raise NotADirectoryError(string)
    return abs_path
def file_path(string):
    if os.path.exists(string):
        return string
    else:
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), string)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='RAndom Walk Resampling command line software.')

    parser.add_argument('alnfile', type=file_path, help='<input> FASTA alignment file.')
    parser.add_argument('treefile', type=file_path,help='<input> Newick phylogeny tree file.')

    parser.add_argument('--algorithm', nargs='?', const=1, help='<optional> resampling algorithm choices [rawr, seres], default rawr.', default="rawr")
    parser.add_argument('--task', nargs='?', const=1, help='<optional> output options: [tree, msa] Either outputs a phylogenetic tree using RAxML or a multiple sequence alignment support estimation, defaults to tree.', default="tree")
    parser.add_argument('--path', nargs='?', const=1, type=dir_path, help="<optional> path to results directory, default output_{algorithm}.", default=os.path.join(os.path.abspath(os.getcwd()),"output_RAWR"))
    parser.add_argument('--n', nargs='?', const=1, type=int, help='<optional> total number of resampled sequences (min 2), defaults to 10.', default=10)
    parser.add_argument('--rate', nargs='?', const=1, type=float, help='<optional> reversal rate between [0, 1], defaults to 0.1.', default=0.1)
    parser.add_argument('--numanchor', nargs='?', const=1, type=int, help='<optional> SERES only. Represents number of anchors, defaults to 20.', default=20)
    parser.add_argument('--lenanchor', nargs='?', const=1, type=int, help='<optional> SERES only. Represents length of anchors, defaults to 5.', default=5)
    parser.add_argument('-seed', type=int, help='<optional> set random seed', default=None)

    parser.add_argument('--mafft', help='<optional> you can supply path to your own mafft if you do not want to use the standalone versions we supply.', default=None)
    parser.add_argument('--raxml', help='<optional> you can supply path to you own RAxML if you do not want to use the standalone versions we supply.', default=None)
    parser.add_argument('--verbose', action='store_true', help='<optional> prints verbose messages')
    args = parser.parse_args()

    if not args.verbose: 
        warnings.simplefilter(action='ignore', category=FutureWarning)
        with HiddenPrints(): 
            if (args.algorithm.lower() == "rawr"):
                # RAWR resample the alignment file into /samples/ directory
                mysampler = sampler.rawrSampler(args.alnfile, args.path,args.rate,args.n) 
                mysampler.sampleSeqs()
            elif (args.algorithm.lower() == "seres"): 
                # SERES resample the alignment file into /samples/ directory
                mysampler = sampler.seresSampler(args.alnfile, args.path,args.rate,args.n,args.numanchor, args.lenanchor)
                mysampler.sampleSeqs()
            else:
                sys.exit('Please enter one of the available algorithm options: rawr or seres.')
    
            if (args.task.lower() == "msa"):
                # use the resampled files to do MSA support estimation
                estimator = supportEstimator.msaSupportEstimator(args.alnfile, args.path, args.n)
                estimator.calculateSupport()
            elif (args.task.lower() == "tree"):
                # use the resampled files to do phylogenetic tree support estimation
                estimator = supportEstimator.treeSupportEstimator(args.alnfile, args.path, args.treefile,args.n) 
                estimator.calculateSupport()
            else:
                sys.exit('Please enter one of the available support estimations: msa or tree.')
   
    else:
        if (args.algorithm.lower() == "rawr"):
            # RAWR resample the alignment file into /samples/ directory
            mysampler = sampler.rawrSampler(args.alnfile, args.path,args.rate,args.n) 
            mysampler.sampleSeqs()
        elif (args.algorithm.lower() == "seres"): 
            # SERES resample the alignment file into /samples/ directory
            mysampler = sampler.seresSampler(args.alnfile, args.path,args.rate,args.n,args.numanchor, args.lenanchor)
            mysampler.sampleSeqs()
        else:
            sys.exit('Please enter one of the available algorithm options: rawr or seres.')
        if (args.task.lower() == "msa"):
            # use the resampled files to do MSA support estimation
            estimator = supportEstimator.msaSupportEstimator(args.alnfile, args.path, args.n)
            estimator.calculateSupport()
        elif (args.task.lower() == "tree"):
            # use the resampled files to do phylogenetic tree support estimation
            estimator = supportEstimator.treeSupportEstimator(args.alnfile, args.path, args.treefile,args.n) 
            estimator.calculateSupport()
        else:
            sys.exit('Please enter one of the available support estimations: msa or tree.')

    # print to standard output so that galaxy can find it
    if (args.task.lower() == "msa"):
        #with open(os.path.join(args.path,"MSA.support.csv"), "r") as result:
        #    sys.stdout.write(result.read())
        with open(os.path.join(args.path,"MSA.support.csv_jalview_annotation.txt"), "r") as result:
            sys.stdout.write(result.read())

    elif (args.task.lower() == "tree"):
        with open(os.path.join(args.path,"tree.support.txt"), "r") as result:
            sys.stdout.write(result.read())

    shutil.rmtree(args.path)
    # remove sampling directories for galaxy
"""
    outputpath=os.path.abspath(os.getcwd())
    if (args.task.lower() == "msa"):
        shutil.copy(os.path.join(args.path,"MSA.support.csv"), outputpath)

    elif (args.task.lower() == "tree"):
        shutil.copy(os.path.join(args.path,"tree.support.txt"), outputpath)

""" 
