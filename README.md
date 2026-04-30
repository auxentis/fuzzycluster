# fuzzycluster
Web-server for FuzzyClusTeR available at https://utils.researchpark.ru/bio/fuzzycluster

DNA repeats constitute a large fraction of eukaryotic genomes and play important roles in genome stability and evolution. While tandem repeats such as microsatellites have been extensively studied, the genomic organization and potential functions of dispersed or loosely organized repeat patterns remain poorly understood. Here we present FuzzyClusTeR, a web server for the identification, visualization and enrichment analysis of DNA repeat clusters in genomic sequences. Using parameterized metrics, FuzzyClusTeR detects both classical tandem repeats and regions where related motifs occur in proximity without forming perfect tandem arrays, which we term diffuse (or fuzzy) repeat clusters. The server supports analysis of user-defined sequences as well as genome-scale datasets, including the T2T-CHM13 and GRCh38 human genome assemblies, and provides interactive visualization and statistical tools for assessing the genomic distribution of repetitive motifs and corresponding clusters. As a demonstration, we analyzed telomeric-like repeats in the T2T-CHM13v2.0 genome and identified families of diffuse clusters enriched in these motifs. Comparison with simulated sequences suggests that these clusters represent non-random genomic patterns with potential evolutionary and functional significance. FuzzyClusTeR enables systematic exploration of repeat clustering across genomic regions or entire genomes. 


All backend with added sqlite3 jobs:

```
for management of fasta files uploads and chroomosome (worker.py)
for management of multifasta files uploads and entire genomes (heavy_worker.py)
for management of simple sequence requests (light_worker.py)
a common entry.dp is added for all job_ids
```

```
fuzzycluster-web/
|- app.py                  # Main Flask application file
|- config.py                # Config file for SQLite DB and jobs
|- worker.py                # SQLite Worker for fasta processing
|- heavy_worker.py          # SQLite Heavy_worker for multifasta processing 
|- light_worker.py          # SQLite Light_worker for simple requests processing
|- fuzzy_functions.py          # All essential functions
|- index.html             # HTML template (rendered with Jinja2)
|- about.html              # HTML home template (rendered with Jinja2)
|- waiting.html             # HTML template for fasta queues (rendered with Jinja2)
|- heavy_waiting.html       # HTML template for multifasta queues (rendered with Jinja2)
|- light_waiting.html       # HTML template for light queues (rendered with Jinja2)
|- static/ 
    |- shared/ # shared between containers folder for sessions        
    |- css/ # style.css is here
    |- genomes/ # folder for genomes

```
