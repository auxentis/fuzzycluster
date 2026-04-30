# fuzzycluster
Web-server for FuzzyCluster 

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
