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

2. Please download VirBot by "git clone"
```
git clone https://github.com/GreyGuoweiChen/RNA_virus_detector
cd RNA_virus_detector
```

3. Please download the database file from [OneDrive](https://portland-my.sharepoint.com/:f:/g/personal/gwchen3-c_my_cityu_edu_hk/EufG0D1CYLREg_7K1UgMvpwBg6bbBIJSM0vdV5udvw1k_w?e=nOJo3G) and uncompress them in the same directory with VirBot.py .
```
unzip ref.zip
```

4. You may want to add permissions to all database file.
```
chmod -R 777 ref
```

5. (Optional) If you want to use VirBot as an executable file, please add VirBot to your environment variables path.
For examples (if your .bashrc file is under the user directory):
    
    For macOS and Linux:
    ```
    chmod +x VirBot.py
    echo 'PATH=$PATH:YOUR_DIC/RNA_virus_detector' >> ~/.bashrc
    ```
    Please replace YOUR_DIC by the path of RNA_virus_detector.
    
    Please remember to activate the PATH by:
    ```
    source ~/.bashrc
    ```
    and check its avaiability by:
    ```
    echo $PATH
    ```
    
## Usage:

```
# If you run VirBot.py as python script.
python VirBot.py [--input INPUT_CONTIG] [--output OUTPUT_DIRECTORY] [--sen]

# If the environment variables path is correctly set and you run VirBot.py as an executable file.
VirBot.py [--input INPUT_CONTIG] [--output OUTPUT_DIRECTORY] [--sen]
```

### Options 

```
--input: The input contig file in fasta format.
--output: The output directory (default: VB_result).
--sen (Optional): Use the sensitive mode for VirBot.
```

### Example:
  
```
python3 VirBot.py --input test/test_input.fa

VirBot.py --input test/test_input.fa --output VB_result --sen
```
