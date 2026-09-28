# fuzzycluster
Web-server for FuzzyClusTeR available at https://utils.researchpark.ru/bio/fuzzycluster

DNA repeats constitute a large fraction of eukaryotic genomes and play important roles in genome stability and evolution. While tandem repeats such as microsatellites have been extensively studied, the genomic organization and potential functions of dispersed or loosely organized repeat patterns remain poorly understood. Here we present FuzzyClusTeR, a web server for the identification, visualization and enrichment analysis of DNA repeat clusters in genomic sequences. Using parameterized metrics, FuzzyClusTeR detects both classical tandem repeats and regions where related motifs occur in proximity without forming perfect tandem arrays, which we term diffuse (or fuzzy) repeat clusters. The server supports analysis of user-defined sequences as well as genome-scale datasets, including the T2T-CHM13 and GRCh38 human genome assemblies, and provides interactive visualization and statistical tools for assessing the genomic distribution of repetitive motifs and corresponding clusters. As a demonstration, we analyzed telomeric-like repeats in the T2T-CHM13v2.0 genome and identified families of diffuse clusters enriched in these motifs. Comparison with simulated sequences suggests that these clusters represent non-random genomic patterns with potential evolutionary and functional significance. FuzzyClusTeR enables systematic exploration of repeat clustering across genomic regions or entire genomes. 


## FuzzyClusTeR structure
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
    |- genomes//2025_Genomes_for_fuzzycluster/ # folder for genomes

```

```
sqlite jobs
for management of fasta files uploads and chroomosome (worker.py)
for management of multifasta files uploads and entire genomes (heavy_worker.py)
for management of simple sequence requests (light_worker.py)
a common entry.dp is added for all job_ids
```

## FuzzyClusTeR installation

### 1. Set Up a Python Virtual Environment
Create project `fuzzy`:

```bash
python3 -m venv fuzzy
source fuzzy/bin/activate
```
### 2. Install dependencies and servers (in the activated environment)

```bash
pip install -r requirements.txt  # requirements.txt contains all the libraries required for the application work
pip install gunicorn  # Gunicorn is WSGI HTTP server for running Python web apps (like Flask and Django). It sits between the application code and the internet, managing requests efficiently.
sudo apt install nginx  # Nginx is a high-performance web server and reverse proxy.
```
### 3. Test the Flask App (this is just on Flask's development server)

```bash
python app.py
```
The output should look like this:

```
* Running on http://127.0.0.1:5000/
* Restarting with stat
* Debugger is active!
* Debugger PIN: 123-456-789
```
In web-browser it can be seen at [http://127.0.0.1:5000/](http://127.0.0.1:5000/)

### 4. Set Up Gunicorn (Production server)
In `wsgi.py` file create one line:
```python
from app import app
```
Production-ready WSGI deployment:
```bash
gunicorn -w 4 -b 127.0.0.1:<desired_port> wsgi:app  # Here we need to connect to a desired port, usually 8000
```
### 5. Create a systemd Service File
Create `/etc/systemd/system/fuzzycluster.service` (fuzzycluster is the name of the service)

File body:
```
[Unit]
Description=Gunicorn instance to serve fuzzycluster
After=network.target

[Service]
User=username
Group=www-data
WorkingDirectory=/home/username/fuzzycluster-web
Environment="PATH=/home/username/fuzzycluster-web/fuzzy/bin"
ExecStart=/home/username/fuzzycluster-web/fuzzy/bin/gunicorn -w 4 -b 127.0.0.1:<desired_port> wsgi:app

[Install]
WantedBy=multi-user.target
```

Reload and start:
```bash
sudo systemctl daemon-reexec
sudo systemctl daemon-reload
sudo systemctl start fuzzycluster
sudo systemctl enable fuzzycluster
sudo systemctl status fuzzycluster

sudo curl http://127.0.0.1:8000
```
Outputs html page like

### 6. Set Up Nginx
Create a config file: `/etc/nginx/sites-available/fuzzycluster`
Domain: # domain address here

File body:
```nginx
server {
    listen 80;
    fuzzycluster domain.ru; # Enter the domain name here

    location / {
        proxy_pass http://127.0.0.1:<desired_port>; # Gunicorn is listening here, usually port 8000
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
    }
}
```
Enable the config by creating a link:
```bash
sudo ln -s /etc/nginx/sites-available/fuzzycluster /etc/nginx/sites-enabled
sudo nginx -t
sudo systemctl restart nginx
sudo systemctl status nginx
sudo systemctl stop nginx # To stop serving
```

Set Up Firewall:
```bash
sudo ufw allow 'Nginx Full'
```
