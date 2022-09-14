# VirBot: a protein-based RNA virus detector for metagenomic data
VirBot is an RNA virus detection tool, which allows accurate and sensitive identificaton for RNA viral contigs. Currently we are still optimizating the codes. But this version is ready to test/use and we welcome any suggestions.

## Dependency:
* Prodigal
* HMMER3
* DIAMOND
* python 3.x

### Quick install

1. We highly recommend using `conda` to install all the dependencies.
```bash
# create the environment and install the dependencies
conda create -n virbot -c bioconda prodigal hmmer diamond python=3
# activate the environment
conda activate virbot
```

2. We currently do not support download VirBot from conda. Please download VirBot by "git lfs clone" instread of "git clone".
```
git lfs clone https://github.com/GreyGuoweiChen/RNA_virus_detector
```

3. You may want to add permissions to all database file.
```
chmod -R 777 RNA_virus_detector/ref
```

4. You may want to add VirBot to your environment variables path.
For examples (if your .bashrc file and RNA_virus_detector are under the user directory):
    
    For macOS:
    ```
    vi ~/.bashrc
    export PATH=$PATH:~/RNA_virus_detector
    source ~/.bashrc
    ```
    
    For Linux:
    ```
    vi ~/.bashrc
    export PATH="$PATH:~/RNA_virus_detector"
    source ~/.bashrc
    ```


## Usage:

```
# If the environment variables path is correctly set.
VirBot.py [--input INPUT_CONTIG] [--output PREDICTION] [--temp_dir TEMPORAY FOLDER]

# If you run VirBot.py in the software folder.
python VirBot.py [--input INPUT_CONTIG] [--output PREDICTION] [--temp_dir TEMPORAY FOLDER]
```

### Options 

```
--input: The input contig file in fasta format.
--output: The output fasta file containing all the predicted RNA virus contigs.
--temp_dir (Optional): The temporary directory used to hold temporary files (default: VB_result).
--sen: Use the sensitive mode for VirBot.
```

### Example:
  
```
VirBot.py --input test/test_input.fa --output test_output.fna

VirBot.py --input test/test_input.fa --output test_output.fna --sen
```
