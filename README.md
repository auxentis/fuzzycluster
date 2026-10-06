# FuzzyClusTeR

Web server for **FuzzyClusTeR**: https://utils.researchpark.ru/bio/fuzzycluster

DNA repeats constitute a large fraction of eukaryotic genomes and play important roles in genome stability and evolution. While tandem repeats such as microsatellites have been extensively studied, the genomic organization and potential functions of dispersed or loosely organized repeat patterns remain poorly understood.

**FuzzyClusTeR** is a web server for the identification, visualization, and enrichment analysis of DNA repeat clusters in genomic sequences. Using parameterized metrics, FuzzyClusTeR detects both classical tandem repeats and regions where related motifs occur in proximity without forming perfect tandem arrays, which we term **diffuse (or fuzzy) repeat clusters**.

The server supports analysis of user-defined sequences as well as genome-scale datasets, including the T2T-CHM13 and GRCh38 human genome assemblies, and provides interactive visualization and statistical tools for assessing the genomic distribution of repetitive motifs and corresponding clusters.

As a demonstration, we analyzed telomeric-like repeats in the T2T-CHM13v2.0 genome and identified families of diffuse clusters enriched in these motifs. Comparison with simulated sequences suggests that these clusters represent non-random genomic patterns with potential evolutionary and functional significance.

FuzzyClusTeR enables systematic exploration of repeat clustering across genomic regions or entire genomes.

## Repository structure

```text
fuzzycluster-web/
├── app.py                    # Main Flask application
├── config.py                 # Configuration, SQLite databases, and jobs
├── worker.py                 # Processing of FASTA sequence jobs
├── heavy_worker.py           # Processing of multi-FASTA and genome-scale jobs
├── light_worker.py           # Processing of simple sequence requests
├── fuzzy_functions.py        # Core FuzzyClusTeR functions
├── wsgi.py                   # WSGI entry point
├── gunicorn.conf.py          # Gunicorn configuration
├── index.html                # Main HTML template
├── about.html                # About page
├── waiting.html              # FASTA job queue template
├── heavy_waiting.html        # Multi-FASTA job queue template
├── light_waiting.html        # Simple request queue template
├── static/
│   ├── shared/               # Shared files and session data
│   ├── css/                  # Stylesheets
│   └── genomes/              # Genome datasets
└── requirements.txt          # Python dependencies
```

## SQLite databases

FuzzyClusTeR uses SQLite databases to manage submitted jobs:

* `worker.py` — FASTA sequence and chromosome processing
* `heavy_worker.py` — multi-FASTA and genome-scale processing
* `light_worker.py` — simple sequence requests
* `entry.db` — common job-entry database

The required databases are initialized automatically when the Gunicorn master process starts.

## Installation

### 1. Create a Python virtual environment

```bash
python3 -m venv fuzzy
source fuzzy/bin/activate
```

### 2. Install dependencies

```bash
pip install -r requirements.txt
pip install gunicorn
```

## Running FuzzyClusTeR locally

### 1. WSGI entry point

The `wsgi.py` file should contain:

```python
from app import app
```

### 2. Start the application with Gunicorn

```bash
gunicorn -c gunicorn.conf.py -w 4 -b 127.0.0.1:8000 wsgi:app
```

Here, `8000` is the local port used by Gunicorn and can be changed if necessary.

The application can then be accessed in a web browser at:

```text
http://127.0.0.1:8000/
```

### 3. Optional: test the Flask application directly

For development and debugging, the application can also be started using Flask:

```bash
python app.py
```

It will normally be available at:

```text
http://127.0.0.1:5000/
```

## Configuration

The Flask application uses the following configuration:

```python
ROUTE_PREFIX = os.environ.get("ROUTE_PREFIX", "")
app = Flask(
    __name__,
    static_folder="static",
    static_url_path=f"{ROUTE_PREFIX}/static",
    template_folder="."
)

app.config["ROUTE_PREFIX"] = ROUTE_PREFIX
app.config["SECRET_KEY"] = SECRET_KEY
app.config["UPLOAD_FOLDER"] = "static/upload"
app.config["CSV_FOLDER"] = "static/csvfiles/"
```

Client-side API requests use the configured API prefix, for example:

```javascript
const response = await fetch(
    `${window.API_PREFIX}/job_status/${jobId}`
);
const data = await response.json();
```

## Public instance

A publicly accessible instance of FuzzyClusTeR is available at:

https://utils.researchpark.ru/bio/fuzzycluster

*FuzzyClusTeR: a web server for analysis of tandem and diffuse DNA repeat clusters with application to telomeric-like repeats*
https://www.biorxiv.org/content/10.64898/2026.03.19.712643v1
