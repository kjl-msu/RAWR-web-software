RAWR is a non-parametric resampling method written in Python. We provide two resampling algorithms here, RAndom Walk Resampling (RAWR) and SEquential
RESampling (SERES), for performing support estimation for multiple sequence alignments and phylogenetic trees. RAWR and SERES conceptually work by random walking along biological molecular sequences and retaining the sequence dependence during the resampling process. We present this tutorial for people who prefer to use a GUI application or a website interface. Picture tutorials are available in the supplementary document of the associated publication. The algorithms implemented in RAWR are freely available under the GPLv3 License. 

* [Installation](#installation)
  * [Software Install](#software-install)
  * [Website Install](#website-install)
  * [Local Galaxy Install](#local-galaxy-install)
    * [Compile RAWR](#compile-rawr)
    * [Galaxy Virtualenv](#galaxy-virtualenv)
      * [Run Galaxy in Virtualenv](#running-galaxy-virtualenv)
    * [Galaxy Conda](#galaxy-conda)
      * [Run Galaxy in Conda](#running-galaxy-conda)
    * [Debugging Galaxy](#debugging-galaxy)
  * [MAFFT and RAxML versions](#mafft-and-raxml-versions)
* [RAWR API](#rawr-api)
  * [Examples](#rawr-api-examples)
* [Input](#input)
* [Output](#output)
  * [Sampled sequence output](#sampled-sequence-output)
  * [MSA support estimation output](#msa-support-estimation-output)
  * [Tree support estimation output](#tree-support-estimation-output)
* [Visualize Output](#visualize-output)
  * [MSA support visualizer](#msa-support-visualizer)
  * [Tree support visualizer](#tree-support-visualizer)
  * [Galaxy tree support visualizer](#galaxy-tree-support-visualizer)
* [Related Publications](#related-publications)

Installation
---------------

### Software Install

RAWR has a GUI interface tested to work on Apple macOS Sequoia 15.0, Linux Ubuntu 22.04 LTS, and Windows 11. Please make sure you have set `python` to link to Python3. We have last checked on Python 3.10.15 that the following python dependencies are compatible.

Open a terminal and navigate to the `rawr-web-software` folder. We provide code for Anaconda and Virtualenv to manage your python packages.
- Virtualenv v20.26.6 commands:
```
cd rawr-software
python -m venv .rawr
. .rawr/bin/activate 

pip install "ete3>=3.1.3"
pip install "PyQt5>=5.15.10"
pip install "numpy>=1.26.4"
pip install "scipy>=1.14.1"
pip install "sklearn>=1.5.2"
pip install "Bio>=1.84"
pip install "flask>=3.0.3"
pip install "pathlib>=1.0.1"
pip install "pytz>=2024.1"
pip install "python-dateutil>=2.9.0post0"
pip install "pandas>=2.2.2"
pip install "joblib>=1.4.2"
```
- If you prefer Anaconda package manager, use the command below to install dependencies. We used Anaconda v22.9.0.
```
conda create -n python3.10 python=3.10
conda activate python3.10

conda install -c conda-forge -c etetoolkit ete3 pyqt numpy scipy scikit-learn biopython flask pathlib pytz python-dateutil pandas
```
To run RAWR GUI software, run the following command in `rawr-software` folder:
```
python main.py     
```

### Website Install

We provide the codebase to set up a RAWR web server to run on your local host. Users can directly use the web server for small dataset analysis. The limitation of the web server is 50Mb of the input alignment with no more than 50 taxa.

First, please make sure your local host's operating system is one of the following: Linux or macOS. Please make sure you have set `python` to link to Python3.

If you are using macOS, please change your terminal shell to Bash.

If you're using Windows, please install Windows Subsystem for Linux and reboot prior to setting up the RAWR web service. If you have access to Hyper-V, you can alternatively create a Linux virtual machine.

The web server has been tested on Linux macOS Sequoia 15.0 and Ubuntu 22.04 LTS.

Open a terminal and navigate to the `rawr-web-software` folder. We provide code for Anaconda and Virtualenv to manage your python packages. 
- Virtualenv v20.26.6 commands:
```
cd rawr-web
python -m venv .rawr
. .rawr/bin/activate 

pip install "ete3>=3.1.3"
pip install "PyQt5>=5.15.10"
pip install "numpy>=1.26.4"
pip install "scipy>=1.14.1"
pip install "sklearn>=1.5.2"
pip install "Bio>=1.84"
pip install "flask>=3.0.3"
pip install "pathlib>=1.0.1"
pip install "pytz>=2024.1"
pip install "python-dateutil>=2.9.0post0"
pip install "pandas>=2.2.2"
pip install "joblib>=1.4.2"
```
- Conda v22.9.0 commands:
```
conda create -n python3.10 python=3.10
conda activate python3.10

conda install -c conda-forge -c etetoolkit ete3 pyqt numpy scipy scikit-learn biopython flask flask-mail pathlib pytz python-dateutil pandas
```

Second, set up your server e-mail. 

Its purpose is to notify the user when and where the results can be found, in case they time out or exit from the homepage, they will still get the result safely delivered via email when it's ready. 

We have included instructions below to set up a Gmail account as a server sender. Setting up Flask to use an existing emailing service is the easiest way to use the asynchronous emailing option, but feel free to adapt the code to your needs.

1. Turn on 2FA for your specific gmail account. 

2. Generate an "App password" with Gmail.

3. Edit the file `app.py` in the following app.config.update() section, updating `MAIL_USERNAME` and `MAIL_PASSWORD`. 
```
 app.config.update(
    MAIL_SERVER='smtp.gmail.com',
    MAIL_PORT=465,          
    MAIL_USE_SSL=True,      
    MAIL_USERNAME='username@gmail.com',
    MAIL_PASSWORD='app_password'
)
```
To run RAWR web server, run the following command in `rawr-web` folder and open the website link (http://10.0.2.15:5000/) with any browser
```
python app.py

#You will find this line in the terminal or something similar:
* Running on http://10.0.2.15:5000/
```


### Local Galaxy Install
Galaxy is available for Linux and macOS. This tutorial has been tested on Linux macOS Sequoia 15.0 and Ubuntu 22.04 LTS. 

Please make sure you have set `python` to link to Python3.

If you are using macOS, please change your terminal shell to Bash.

If you're using Windows, please install Windows Subsystem for Linux and reboot prior to setting up Galaxy. If you have access to Hyper-V, you can alternatively create a Linux virtual machine.

Navigate to `rawr-galaxy` folder and download galaxy following instructions from official website: https://galaxyproject.org/admin/get-galaxy/ 

Folder structure after Galaxy is set up properly:
```
rawr-galaxy/
├── galaxy/              # created when you download galaxy
├── conda_stuff/
├── dist/
├── rawr_src/
│   ├── mafft/
│   └── raxmlHPC/
```
Three steps to set up Galaxy:
- (1) compile RAWR 
- (2) install Planemo and set up Galaxy
- (3) use Planemo to serve Galaxy locally (after initial setup, only use this step)

## Compile RAWR
This tutorial used Python v3.10.15, but feel free to use a newer version.

1. navigate to `rawr-galaxy` folder, create rawr environment, and compile RAWR binary executable

- Virtualenv v20.26.6 commands:
```
cd rawr-galaxy
python -m venv .rawr
. .rawr/bin/activate 

pip install Bio>=1.84
pip install "python-dateutil>=2.9.0"
pip install "ete3>=3.1.3"
pip install "numpy>=2.1.2"
pip install "scipy>=1.14.1"
pip install "sklearn>=1.5.2"
pip install "pandas>=2.2.3"
pip install "pyinstaller>=6.10.0"

pyinstaller --hidden-import="sklearn.utils._cython_blas" --hidden-import="sklearn.neighbors.typedefs" rawr.py

deactivate
```
- Conda v22.9.0 commands:
```
cd rawr-galaxy

conda create -n rawr python=3.10
conda activate rawr

conda config --show channels
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --remove channels defaults

conda install -c conda-forge -c etetoolkit ete3 numpy scipy scikit-learn biopython pandas pyinstaller

pyinstaller --hidden-import="sklearn.utils._cython_blas" --hidden-import="sklearn.neighbors.typedefs" rawr.py

conda deactivate
```
- if you would prefer, we provide precompiled RAWR tarballs at https://gitlab.msu.edu/liulab/rawr-web-software/-/tree/master/rawr-galaxy

2. modify the `rawr-galaxy/rawr.xml` file on line 3 with full path to your rawr executable. To find the full path in, navigate to `dist/rawr/` and type `pwd` into your terminal.
```
1 <tool id="rawr" name="RAWR random sampler" version="1.0" python_template_version="3.5">
2    <command detect_errors="exit_code"><![CDATA[
3     /full/path/to/dist/rawr/rawr
```

3. Copy over MAFFT and RAxML folders from `rawr_src`
```
cp -r rawr_src/mafft* dist/rawr/_internal/
cp -r rawr_src/raxml* dist/rawr/_internal/
```

4. double check that rawr is working
```
/full/path/to/dist/rawr/rawr
#returns:
usage: rawr [-h] [--algorithm [ALGORITHM]] [--task [TASK]] [--path [PATH]] [--n [N]] [--rate [RATE]] [--numanchor [NUMANCHOR]]
            [--lenanchor [LENANCHOR]] [-seed SEED] [--mafft MAFFT] [--raxml RAXML] [--verbose]
            alnfile treefile
rawr: error: the following arguments are required: alnfile, treefile

# small example
/full/path/to/dist/rawr/rawr example/dataset-10taxa/alignment.fasta example/dataset-10taxa/infer.tree --task tree --algorithm rawr --n 2 --rate 0.1
```

## Galaxy Virtualenv
This Virtualenv tutorial is for installing Planemo and Galaxy. We used Python v3.10.15 in Virtualenv v20.26.6 with Galaxy v24.1.2. 

1. Please download Galaxy from official website: https://galaxyproject.org/admin/get-galaxy/. This creates the `galaxy` folder.

2. navigate to `rawr-galaxy/galaxy` and create virtualenv
```
cd galaxy
python -m venv .venv
. .venv/bin/activate
```

3. install Galaxy dependencies. The `common_startup.sh` script will set up the Galaxy server for the first time. 
```
PYTHONPATH= sh scripts/common_startup.sh --no-create-venv
```

- Package versions used in tutorial, installed by the `common_startup.sh` script:
```
virtualenv                20.26.6      
python                    3.10.15                   
pip                       22.0.2
planemo                   0.75.26
galaxy-tool-util          24.1.2
galaxy-util               24.1.2
cwltool                   3.1.20241007082533
yarn                      1.22.22           
```

4. locate the modified `rawr.xml` file. Make sure it's in `rawr-galaxy` folder (or wherever you call planemo from, make sure rawr.xml is in the same folder)
```
cd ..
ls       # rawr-galaxy folder
you should see the rawr.xml file you edited earlier
dist  example  galaxy  rawr.py  rawr.spec  rawr_src  rawr.xml 
```

5. install planemo and serve Galaxy
```
pip install planemo
planemo serve --galaxy_root galaxy 
```

If you run into trouble, look at the `debugging Galaxy` section.

#### Running Galaxy Virtualenv
In the future, you do not need to set up again. Simply do the following:
```
cd rawr-galaxy
. galaxy/.venv/bin/activate
planemo serve --galaxy_root galaxy

#You should see this line:
Serving on http://127.0.0.1:9090
# Click to open the link with a browser
```
## Galaxy Conda
This Conda tutorial teaches you how to install Planemo, set up Galaxy, and create RAWR conda package. This was tested with Python v3.10.15 in Anaconda v24.5.0 with Galaxy v24.1.2. 

Note: you can set up Galaxy Conda in one of two ways: 
- (I) manually, as described in the instructions below. This is the recommended method.
- (II) allow Planemo to handle Galaxy set up, it will pull from the latest github repository. If you wish to let Planemo do it, simply install planemo, complete step 1, and skip to running galaxy conda.

Manual set up (recommended)

1. Please download Galaxy from official website: https://galaxyproject.org/admin/get-galaxy/. This creates the `galaxy` folder.

2. create local conda package for RAWR
```
conda activate rawr
cp -r conda_stuff/* dist 
cd dist
```
- Note: folder `conda_stuff` contains a different rawr.xml file intended for Conda.
- Modify `dist/rawr.xml` on line 6 with full path to your rawr executable. To find the full path, navigate to `dist/rawr/` and type `pwd` into your terminal.

```
1 <tool id="rawr" name="RAWR random sampler" version="1.0" python_template_version="3.5">
2    <requirements>
3        <requirement type="package" version="1.0">rawr</requirement>
4    </requirements>
5    <command detect_errors="exit_code"><![CDATA[
6     /full/path/to/dist/rawr/rawr
```
- build conda package
```
conda install conda-build
conda build conda-recipe
conda deactivate
```
3. ensure your environment has the correct version of python you want to use with Galaxy (code has been tested on Python 3.8 as well, and you can try newer Python versions).
```
python --version
Python 3.10.15
```

4. navigate to the "galaxy" folder and start up galaxy with "sh run.sh".
```
cd galaxy
sh run.sh
```

5. The Galaxy server will then be set up for the first time, followed immediately by starting the Galaxy server. If you run into issues, help is available in the `debugging Galaxy` section.

If you see the following services starting, Galaxy has finished set up and is trying to start the server. It is safe to stop the process here.
```
celery                           STARTING  
celery-beat                      STARTING  
gunicorn                         STARTING  
```
Or if you see this line, Galaxy server is already up and running, indicating both server set up and launch has completed successfully. It is safe to stop the process here.
```
Serving on http://127.0.0.1:8080
```

You can close the Galaxy server process.
```
# Break the server process with (macOS) command+. or (Linux) ctrl+c
```

- During setting up Galaxy for the first time, a `_galaxy_` conda environment will be created automatically for you. 

- Package versions used in tutorial for `_galaxy_` (this conda environment is created for you when setting up Galaxy for the first time).
```
python                    3.10.15
galaxy-tool-util          24.1.2
galaxy-util               24.1.2
planemo                   0.75.26
cwltool                   3.1.20240708091337
stdlib-list               0.10.0
virtualenv                20.26.6
yarn                      1.22.22
```

6. install Planemo for Galaxy
```
conda deactivate
conda activate _galaxy_   # created for us by Galaxy

conda config --show channels
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --remove channels defaults

conda install -c bioconda -c conda-forge cwltool stdlib-list planemo
```

7. use Planemo to serve Galaxy
```
cd ../dist        # where rawr.xml is located
planemo serve --conda_use_local --galaxy_root ../galaxy

#You should see this line:
Serving on http://127.0.0.1:9090
# Click to open the link with a browser
```

If the installation does not go smoothly, some tips are available in the `debugging Galaxy` section.

#### Running Galaxy Conda
In the future, you do not need to set up again. Simply do the following:
```
cd rawr-galaxy/dist      
conda activate _galaxy_
planemo serve --conda_use_local --galaxy_root ../galaxy

#You should see this line:
Serving on http://127.0.0.1:9090
# Click to open the link with a browser
```

## Debugging Galaxy

### Could not solve environment:
```
Could not solve for environment specs
The following packages are incompatible
```
During Galaxy server set up with `sh run.sh`, you may run into problem where the environment cannot be solved. 
Check that your environment has the same package versions as tutorial with:
```
# go through the list of packages, check all packages are the same version or higher
conda list | grep galaxy-tool-util    
```

### Socket error:
```
Error: Cannot open an HTTP server: socket.error reported AF_UNIX path too long
```
This error occurs when setting up Galaxy server for the first time. It means your filepath is too long for sockets, which are limited to 104 characters. Just move the folder containing galaxy to a directory with shorter path. For example:
```
cp -r rawr-galaxy ~/
```
Then navigate to inside the galaxy folder and try setting up the server.
```
cd ~/rawr-galaxy/galaxy
sh run.sh
```

### Error 137:
```
Browserslist: caniuse-lite is outdated. Please run:
  npx update-browserslist-db@latest
Killed
error Command failed with exit code 137.
```
Code 137 is usually an indication of out of memory. Building the Galaxy client is very memory intensive. Try separating the server set up from the client set up:

1) navigate inside Galaxy folder and build the server with the following command:
```
GALAXY_SKIP_CLIENT_BUILD=1 sh run.sh
```
Once you see the following text, you know that the Galaxy server is being initiated. You can stop the process and proceed to build the client.
```
celery                           STARTING  
celery-beat                      STARTING  
gunicorn                         STARTING  
```

2) ensure you have yarn v1 (at the time of this writing, Galaxy depends on yarn v1)
```
npm install -g yarn@1
```
3) then build Galaxy client:
```
make client
```
4) (optional) if you're still running into code 137, try editing the Makefile to allow Galaxy to use more memory.
Increase the number after "max-old-space-size=". For example:
```
max-old-space-size=10240 # MB. This lets Galaxy use up to 10 GB of RAM, which is definitely more than enough.
```

### Yarn error:
```
Usage: yarn [options]
yarn: error: no such option: --network-timeout
```
During Galaxy set up, you may randomly run into yarn issue. Double check your yarn version is correct. Yarn v1 is what Galaxy requires at the time of this writing. 
```
npm install -g yarn@1
yarn --version
sh run.sh                        # retry setting up Galaxy 
```
 
### Galaxy server fails to launch on missing tool id
```
Exception: Missing tool 'id' for tool at 'XmlToolSource[/path/to/rawr-galaxy/galaxy/packages/app/galaxy/tools/bundled/expression_tools/expression_macros.xml]'
```
If you see errors about missing id, it's likely that Planemo was launched from an unexpected folder. Planemo must be launched in the folder containing the correct version of `rawr.xml`. 
- Conda Galaxy: navigate to `rawr-galaxy/dist`, ensure `rawr.xml` is present, and call `
planemo serve --galaxy_root ../galaxy`
- Virtualenv Galaxy: navigate to `rawr-galaxy`, ensure `rawr.xml` is present, and call `
planemo serve --galaxy_root ../galaxy`

### Virtualenv Galaxy -- Unclickable RAWR:
This is for using Galaxy with Virtualenv. If everything looks fine but you can't use RAWR (and you see error messages about package dependencies or `h11._util.LocalProtocolError: Too much data for declared Content-Length`), you may be using the wrong version of `rawr.xml`.
- Make sure the file `rawr-galaxy/rawr.xml` does **not** have the following `<requirements>` code on lines 2 - 4:
```
2    <requirements>
3        <requirement type="package" version="1.0">rawr</requirement>
4    </requirements>
```
- If you see these three lines in `rawr-galaxy/rawr.xml`, delete them. 


### Conda Galaxy -- Unclickable RAWR:

This is for using Galaxy with Conda. If everything looks like it's working but you can't click on RAWR inside of Galaxy (plus error messages with loading tools in the terminal), there may be an issue with automatically creating a conda environment for RAWR. Make sure of the following:
- `rawr-galaxy/dist/rawr.xml` file has the correct path to RAWR
- `rawr-galaxy/dist` is the folder you're inside when calling Planemo
- `_galaxy_` has most updated version of `conda-build`
-  `/path/to/anaconda3/conda-bld/rawr-1.0-0.tar.bz2` exists -- this is the result of building compiled RAWR into a conda package.

### MAFFT and RAxML versions
RAWR source comes with standalone versions of MAFFT v7.487, RAxML v8.2.12 for macOS and Linux, and RAxML v8.2.10 for Windows. 

If you want updated versions of these software, you can download or compile the standalone versions from their respective websites and update their respective folders inside `src` folder (please make sure to double check that the binary filenames remain the same).
  - MAFFT: https://mafft.cbrc.jp/alignment/software/.
  - RAxML: https://github.com/stamatak/standard-RAxML.

## RAWR API
RAWR comes with application programming interface (API) which you can adapt into your own pipeline. 

Examples of how to implement RAWR API functions are in `main.py`, which you can run:
```
python main.py
```
The current implementation of RAWR API is sequential, but you can parallelize it by modifying `src/supportEstimator.py`.

### RAWR API Examples
- How to use RAWR algorithm to resample an alignment file:
```
rawr = sampler.rawrSampler(alnFile, outputDir)  
rawr.sampleSeqs()
```

- How to use SERES algorithm to resample an alignment file:
```
seres = sampler.seresSampler(alnFile, OutputDir)  
seres.sampleSeqs()
```

- How to adjust parameters for RAWR resampling:
```
rawr = sampler.rawrSampler(alnFile, outputDir, reverseRate, samplenum)  
rawr.sampleSeqs()
```

- How to adjust parameters for SERES resampling:
```
seres = sampler.seresSampler(alnFile, OutputDir, reverseRate,samplenum, anchorLen, anchorNum)  
seres.sampleSeqs()
``` 

- How to use RAWR algorithm with MSA support estimation:
```
rawr = sampler.rawrSampler(alnFile, outputDir)  
rawr.sampleSeqs()
rawr_msa = supportEstimator.msaSupportEstimator(alnFile, outputDir) 
rawr_msa.calculateSupport()
```

- How to use SERES algorithm with phylogeny support estimation:
```
seres = sampler.seresSampler(alnFile, OutputDir)  
seres.sampleSeqs()
seres_tree = supportEstimator.treeSupportEstimator(alnFile, outputDir, treeFile, samplenum) 
seres_tree.calculateSupport()
``` 

- How to use RAWR algorithm with MSA support estimation and custom MAFFT software:
```
rawr = sampler.rawrSampler(alnFile, outputDir)  
rawr.sampleSeqs()
rawr_msa = supportEstimator.msaSupportEstimator(alnFile, outputDir, samplenum, "/path/to/mafft") 
rawr_msa.calculateSupport()
``` 

- How to use SERES algorithm with phylogeny support estimation and custom MAFFT & RAxML software:
```
seres = sampler.seresSampler(alnFile, OutputDir)  
seres.sampleSeqs()
seres_tree = supportEstimator.treeSupportEstimator(alnFile, outputDir, treeFile, samplenum, treeFile,samplenum,"/path/to/mafft","/path/to/raxml") 
seres_tree.calculateSupport()
```

- How to output JalView format MSA support file to import into and visualize in JalView
```
from src.MSA_support_csv_2_jalview_sequence_annotation import support_csv_2_jalview

rawr = sampler.rawrSampler(alnFile, outputDir)  
rawr.sampleSeqs()
rawr_msa = supportEstimator.msaSupportEstimator(alnFile, outputDir) 
rawrMSAoutput = rawr_msa.calculateSupport()
support_csv_2_jalview(rawrMSAoutput,"pink")
```

---------------
Input
---------------
If you are running RAWR for MSA estimation, please prepare a multiple alignment file in FASTA format. 
If you are running RAWR for phylogeny alignment estimation, please prepare a multiple alignment file and a phylogeny in Newick tree format.

In this tutorial, we will not go over preparing the inputs. Instead, we provide a small dataset to try out RAWR. The example dataset has been prepared from sequence files and aligned with MAFFT. Then, we used the alignment file as input into RAxML and reconstructed a phylogeny under the General Time Reversal (GTR) model.

1) Find the example datasets in your RAWR directory under the examples folder. Here, you should see dataset-10taxa. 

2)	Select either RAWR or SERES as resampling algorithm.

3)	Select either multiple sequence alignment (MSA) estimation or phylogenetic tree estimation.

4)	Optionally, you can edit the default parameter values.
```
By default:
Sample Number: 10 (minimum 2. Represents total number of sampled sequences.)
Reverse Rate: 0.1 (0 to 1 not including 0 or 1. Represents the reverse rate of the random walk.)
Anchor Number: 20 (SERES only. Represents the number of anchors.)
Anchor Length: 5 (SERES only. Represents the length of the anchors.)
```

5)	Select a file input for the FASTA alignment file. In our example datasets, choose alignment.fasta

6)	Select a file input for the Newick phylogeny tree file. In out example dataset, choose infer.tree.

7)	Click OK or Submit and let the program run to completion. 

Please note that if you are using the GUI software version of RAWR v1.0, it is best to provide an empty results folder to avoid any problems with RAxML refusing to run due to previous runs' files present or improper software termination.

---------------
Output
---------------
### Sampled sequence output
The Sequence sampler will create a tar compressed directory called `samples`. The sampled sequences and indices can be found here.

```
*.seq.fasta: the unaligned sampled sequence in FASTA format.

*.index: this file contains a series of column indices of the input alignment. The column indices represent sampled columns of input alignment in order.
```

### MSA support estimation output

The MSA support estimator will calculate the support value for input alignment and save the support value to `MSA.support.csv`. 

`MSA.support.csv` contains four columns. The first column locates the position of a residue pair. The second column and third column are the row index of two residues in this residue pair. Only residue pairs with no gaps are contained in this result. The last column is the support value, which ranges from 0 to 1.
`MSA.support.csv_jalview_annotation.txt` contains a JalView annotation file for visualization. (Note: the local Galaxy client is implemented such that running MSA support estimation with RAWR v1.0 will only output JalView annotation file.)

### Tree support estimation output

The tree support estimator will generate annotated tree in Newick format, which is saved to `tree.support.txt`. It will also generate a figure called `tree.support.png` to visualize the tree structure and support values. 

---------------
Visualize Output
---------------

### MSA support visualizer
We used JalView 2.11.4.0 in the following example.
1. Click “File” > “Input Alignment” > “From File” to enter your original “alignment.fasta” input file. 
This will open an alignment viewer inside JalView. 

2. Inside the alignment viewer, click “File” > “Load Features/ Annotations” and select the file “MSA.support.csv_jalview_annotation.txt”. 
 
### Tree support visualizer
We provide a visualization of the Newick tree with confidence support. You can also choose other phylogeny visualizers such as Interactive Tree Of Life (iTOL). 

### Galaxy tree support visualizer
Galaxy can access a suit of visualization and other tools, including Newick Display for visualizing tree support. Please look at the RAWR 1.0 documentation for help and instructions regarding how to download the Newick Display tool from Galaxy Tool Shed.

---------------
Related Publications
---------------
Wang, W., Hejasebazzi, A., Zheng, J., & Liu, K. J. (2021). Build a better bootstrap and the RAWR shall beat a random path to your door: phylogenetic support estimation revisited. Bioinformatics, 37(Supplement_1), i111-i119.

Wang, W., Smith, J., Hejase, H. A., & Liu, K. J. (2020). Non-parametric and semi-parametric support estimation using SEquential RESampling random walks on biomolecular sequences. Algorithms for Molecular Biology, 15(1), 1-15.

Wang, W., Wuyun, Q., & Liu, K. J. (2019, November). An application of random walk resampling to phylogenetic HMM inference and learning. In 2019 IEEE International Conference on Bioinformatics and Biomedicine (BIBM) (pp. 44-51). IEEE.
