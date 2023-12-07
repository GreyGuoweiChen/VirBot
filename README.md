# VirBot: an RNA viral contig detector for metagenomic data
Without relying on cultivation, metagenomic sequencing greatly accelerated the novel RNA virus detection. However, it is not trivial to accurately identify RNA viral contigs from a mixture of species. The low content of RNA viruses in metagenomic data requires a highly specific detector, while new RNA viruses can exhibit high genetic diversity, posing a challenge for alignment-based tools.
    
Here, we develop VirBot, an easy-to-use yet effective RNA virus detection tool from metagenomic data. It takes assembled contigs as input and detect ones from RNA viruses. For each detected contig, we assign it the likely taxon label according to the encoded protein that has the highest score against our pHMMs.
    
VirBot is deisgned based on a comprehensive RNA viral protein family database. Besides the hallmar gene in RNA virus, RdRps, we also leverage other essential proteins, including capsid proteins, envelop proteins, viral auxiliary proteins, etc. By using the adaptative bit score cutoff, VirBot shows its higher specificity in metagenomic dataset and sensitivity in novel RNA virus dataset. VirBot supports identifying contigs as short as 500bp. The construction of RNA viral pHMMs and key components of VirBot are shown as Fig (a,b). For each RNA viral pHMM, we use voting strategy among the clustered viral proteins to determine the taxon, which is used for the taxonomic assignment of the detected contigs. 

We validated VirBot in various scenarios and benchmarked it with other 7 RNA virus detection tools. VirBot achieves higher recall than other tools. And it also demonstrates high specificity in metagenomic data that only contains a small number of RNA viruses. Here, we briefly show the result of VirBot on metagenomics data (Fig (c)) and pure RNA virus dataset (Fig (d)).

| ![Image](images/github.png) |
|:--:|
| (a) Construction of the RNA viral pHMMs database. (b) Sketch of the key components of VirBot. (c) Detection performance on simulated data: ERR1992810 (left) and ERR2185279 (right). (d) Recall on RNA viral datasets: RNA phages dataset (left) and marine water RNA virome dataset (right). 
where the metagenomics dataset contains sequences from 82 eukaryotes, 365 prokaryotes, and DNA/RNA viruses; and for the RNA virus datasets, one comprises 8,849 RNA phages that were barely detected before, while eukaryotic viruses dominate known RNA viruses; another dataset is an RNA virome sample sequenced from marine water containing 114,139 RNA viral seqeunces.| 

## Dependency:
* Prodigal 2.6.3
* HMMER3 3.3.2
* DIAMOND 2.0.15
* pandas 1.5.2
* python 3.x

### Quick install

We highly recommend using `conda` to install all the dependencies.
To install, please download VirBot by "git clone"
```
git clone https://github.com/GreyGuoweiChen/RNA_virus_detector
cd RNA_virus_detector

# create the environment and install the dependencies using conda or mamba
mamba env create -f environment.yml

# activate the environment
conda activate virbot
```

Next you need to get the reference file `ref.zip`. If git-lfs is installed globally, this will have been downloaded already into `virbot/data/ref.zip`. If not you will now need to retrieve them with:
```
git lfs install
git lfs fetch
git lfs pull
```
or alternatively, it is available from [OneDrive](https://portland-my.sharepoint.com/:f:/g/personal/gwchen3-c_my_cityu_edu_hk/EufG0D1CYLREg_7K1UgMvpwBg6bbBIJSM0vdV5udvw1k_w?e=nOJo3G) and can be placed into the virbot/data directory.

Finally, unzip the reference files and install
```
# unzip reference files
cd virbot/data; gunzip ref.zip; cd ../..

# install
pip install .
```

## Usage:
The VirBot is friendly to use. It requires input as fasta format and will return the fasta format output file with the identified RNA virus sequences. 

```
# If the environment variables path is correctly set and you run VirBot.py as an executable file.
virbot [--input INPUT_CONTIG] [--output OUTPUT_DIRECTORY] [--sen] [--threads]

# If you run VirBot.py as python script.
python VirBot.py [--input INPUT_CONTIG] [--output OUTPUT_DIRECTORY] [--sen] [--threads]
```

### Options 

```
--input: The input contig file in fasta format.
--output: The output directory (default: VB_result).
--sen (Optional): Use the sensitive mode for VirBot.
--threads (Optional): The threads number run for HMMER and DIAMOND (default: 8)
```

### Example:
  
```
virbot --input test/test_input.fa

virbot --input test/test_input.fa --output VB_result --sen --threads 8
```

### Benchmark datasets:
Please check the link in [OneDrive](https://portland-my.sharepoint.com/:f:/g/personal/gwchen3-c_my_cityu_edu_hk/EufG0D1CYLREg_7K1UgMvpwBg6bbBIJSM0vdV5udvw1k_w?e=nOJo3G).
