from flask import Flask, request, render_template,session,make_response, flash
from flask import get_flashed_messages
from flask import send_file, abort
from flask_wtf import FlaskForm
from wtforms import FileField, StringField, SubmitField, BooleanField, DecimalField, DateTimeField, IntegerField, RadioField, SelectField, TextAreaField
from wtforms.validators import InputRequired, Length, Optional, NumberRange, Regexp
from werkzeug.utils import secure_filename
from random import choice
from shutil import rmtree
import datetime
import os
import re
import string
import csv
import threading
import time
import pandas as pd
from glob import glob
from io import BytesIO

from zipfile import ZipFile
import sqlite3
from flask import jsonify
import json
from flask import url_for
import uuid

from config import setup_logging, SECRET_KEY

setup_logging()

import logging
logger = logging.getLogger("Flask")


# Create the app object
ROUTE_PREFIX = os.environ.get("ROUTE_PREFIX", "")
app = Flask(
    __name__,
    static_folder="static",
    static_url_path=f"{ROUTE_PREFIX}/static",
    template_folder="."
)

app.config["ROUTE_PREFIX"] = ROUTE_PREFIX

app.config["SECRET_KEY"] = SECRET_KEY
app.config['UPLOAD_FOLDER'] = 'static/upload'
app.config["CSV_FOLDER"] = 'static/csvfiles/'

#We decided to limit sizes of upload files
app.config['MAX_CONTENT_LENGTH'] = 250 * 1000000 # Max upload in bytes
ALLOWED_EXTENSIONS = set(['fasta', 'txt'])

# importing queue functions
from config import LIGHT_DB_PATH, DB_PATH, HEAVY_DB_PATH, ENTRY_DB_PATH, enable_wal_mode, with_retry, get_light_connection, get_connection, get_heavy_connection
# importing fuzzy function for calculations
from fuzzy_functions import search_functionGC
from fuzzy_functions import loop_function
from fuzzy_functions import cluster_functionGC
from fuzzy_functions import image_functionGC
from fuzzy_functions import cluster_loopsGC
from fuzzy_functions import cluster_plot
from fuzzy_functions import highlighted_clusterGC
from fuzzy_functions import fasta_reader
from fuzzy_functions import highlight_functionGC
from fuzzy_functions import repeats_dataframeGC
from fuzzy_functions import save_repeats_csvGC
from fuzzy_functions import save_clusters_csvGC
from fuzzy_functions import csvGC_repeats_open
from fuzzy_functions import csvGC_clusters_open
from fuzzy_functions import save_image
from fuzzy_functions import image_loops_multifasta
from fuzzy_functions import image_hist_loops_multifasta
from fuzzy_functions import image_clusters_multifasta
from Bio.Seq import Seq
from fuzzy_functions import pattern_element_length



class MyForm(FlaskForm):
    file = FileField('Upload File')
    seqname = StringField('Name', validators=[Optional(), Length(min=0, max=50, message="Sequence too long!")])
    dna_pattern = re.compile(r'''
        ^                # Start of string
        [
            ATGC         # Standard bases
            NRYBDKMHVSW  # Ambiguity codes
            \s           # Whitespace (newlines, spaces, tabs)
        ]+               # Match one or more
        $                # End of string
    ''', re.IGNORECASE | re.VERBOSE)

    seq = TextAreaField('Sequence', validators=[
    Optional(), 
    Length(min=30, max=10000, message="Sequence too long!"), Regexp(dna_pattern, message="Sequence contains non-IUPAC symbols")],
    filters=[
        # This cleans the data BEFORE validation
        lambda x: re.sub(r'\s+', '', x).upper().replace('U', 'T') if x else x
    ])
    clusterscore = DecimalField('Score', validators=[Optional(strip_whitespace=True), NumberRange(min=0, max=100000)])
    input_looplength = IntegerField('Loop', validators=[Optional(strip_whitespace=True),  NumberRange(min=0, max=100000)])
    ssr = DecimalField('SSR', validators=[Optional(strip_whitespace=True), NumberRange(min=0.0, max=1.0)])
    submit = SubmitField('Submit')
    printind = BooleanField('Print indexes of repeats')
    printseq = BooleanField('Print sequence with repeats')
    printclust = BooleanField('Print sequence with clusters of repeats')
    csvrepeats = BooleanField('Save found repeats in a file')
    csvclusters = BooleanField('Save found clusters in a file')   
    pattern = StringField('Search custom pattern', validators=[
        Optional(),
        Length(max=200),
        pattern_element_length(3, 30, 3),
        Regexp(r'^[ACGTacgt|]+$', message="Only A, C, G, T and '|' are allowed, avoid empty spaces")
        ])
    
    tandem_repeats = BooleanField(' Tandem repeats only')
    repeats_number = IntegerField('Repeats', validators=[Optional(strip_whitespace=True), NumberRange(min=2)])
    choice = RadioField('Search patterns', choices=[('TTAGGG', 'TTAGGG'), ('FuzzyTel', 'FuzzyTel'), ('Custom pattern', 'Custom pattern')])
    choicegenome = RadioField('Genomes', choices=[('CHM13-T2T', 'CHM13-T2T'), ('GRCh38', 'GRCh38')], default='default', validate_choice=False)
    chromosome = IntegerField('Select chromosome', validators=[Optional(strip_whitespace=True), NumberRange(min=1, max=24, message="Please, choose the chromosome")])
    region = StringField('Enter region')


@app.route(f"{app.config['ROUTE_PREFIX']}/health", methods=['GET'])
def health_check():
    return jsonify({'status': 'healthy'}), 200

if __name__ == '__main__':
    app.run(debug=True)

################################################################
# 251206 entry.db for all job_id
def init_entry_db():
    conn = sqlite3.connect(ENTRY_DB_PATH)
    c = conn.cursor()

    c.execute("CREATE TABLE IF NOT EXISTS jobs ("
              "job_id TEXT PRIMARY KEY, "
              "output_dir TEXT NOT NULL, "
              "created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP)"
    )

    conn.commit()
    conn.close()

def register_job(job_id, session_dir):
    print('register_job', job_id)
    conn = sqlite3.connect(ENTRY_DB_PATH)
    c = conn.cursor()

    c.execute(
        "INSERT INTO jobs (job_id, output_dir) VALUES (?, ?)",
        (job_id, session_dir)
    )

    conn.commit()
    conn.close()

###############################################################
# 251105 SQLite3 light queues for complete_prediction_simple
@with_retry(retries=5, delay=0.2)
def insert_light_job(params_dict, job_id):
    """Insert a new job with a UUID and return its ID."""
    with get_light_connection(LIGHT_DB_PATH) as conn:
        c = conn.cursor()
        c.execute(
            "INSERT INTO jobs (id, session_id, status, params, created_at) VALUES (?, ?, ?, ?, datetime('now'))",
            (job_id, session.get("id"), "queued", json.dumps(params_dict))         
        )
        conn.commit()
    return job_id


# chenge reset to False later
def init_light_db(reset=True):
    conn = sqlite3.connect(LIGHT_DB_PATH)
    c = conn.cursor()
    if reset:
        c.execute("DROP TABLE IF EXISTS jobs")
    c.execute("""
        CREATE TABLE IF NOT EXISTS jobs (
            id TEXT PRIMARY KEY,
            session_id TEXT,
            status TEXT NOT NULL,
            params TEXT,
            result TEXT,
            error TEXT, 
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            started_at TIMESTAMP,
            completed_at TIMESTAMP
        )
    """)
    conn.commit()
    conn.close()
    logger.info(f"Light jobs database initialized ({LIGHT_DB_PATH})")


@app.route(f"{app.config['ROUTE_PREFIX']}/light_job_status/<job_id>")
def light_job_status(job_id):
    with get_light_connection(LIGHT_DB_PATH) as conn:
        c = conn.cursor()
        c.execute("SELECT status, result, error FROM jobs WHERE id=?", (job_id,))
        row = c.fetchone()

        if not row:
            return jsonify({"error": "job not found"}), 404

        status, result, error = row

        if status in ("queued", "running"):
            return jsonify({"status": status})

        if status == "failed":
            return jsonify({"status": "failed", "error": error})

        # completed
        return jsonify({"status": "completed", "result": result})

@app.route(f"{app.config['ROUTE_PREFIX']}/light_job_result/<job_id>")
def light_job_result(job_id):
    with get_light_connection(LIGHT_DB_PATH) as conn:
        c = conn.cursor()
        c.execute("SELECT status, result FROM jobs WHERE id=?", (job_id,))
        row = c.fetchone()

        if not row:
            return "Job not found", 404

        status, result_json = row
        if status != "completed":
            return f"Job status: {status}", 400
        result_data = json.loads(result_json)
        # print('results_json', result_data)
        form = MyForm()

        return render_template(
            'index.html', 
            form=form, 
            job_id = job_id,
            utc_dt=datetime.datetime.utcnow().replace(microsecond=0).isoformat(' '),
            prediction_text0=result_data["prediction_text0"],
            prediction_text1=result_data["prediction_text1"],
            prediction_text2=result_data["prediction_text2"],
            prediction_text3=result_data["prediction_text3"],
            prediction_text4=result_data["prediction_text4"],
            prediction_text5=result_data["prediction_text5"],
            prediction_text6=result_data["prediction_text6"],
            prediction_text7=result_data["prediction_text7"],
            tablesG=result_data["tablesG"],
            tablesC=result_data["tablesC"],
            titlesG=result_data["titlesG"],
            titlesC=result_data["titlesC"],
            prediction_text8=result_data["prediction_text8"],
            prediction_text9=result_data["prediction_text9"],
            prediction_text10=result_data["prediction_text10"],
            prediction_text11=result_data["prediction_text11"],
            download_files=result_data["download_files"],
            plot_image=result_data["plot_image"],
            plot1_url=url_for('static', filename=f"shared/{result_data['plot1_url']}"),
            plot2_url=url_for('static', filename=f"shared/{result_data['plot2_url']}"),
            plot3_url=url_for('static', filename=f"shared/{result_data['plot3_url']}"),
            prediction_text12=result_data["prediction_text12"],
            session_dir=result_data["session_dir"],
            cluster_data_threshold=result_data["cluster_data_threshold"],
            download_figures=result_data["download_figures"]
        )

#################################################################
# 251105 SQLite3 queues

@with_retry(retries=5, delay=0.2)
def insert_job(params_dict, job_id):
    """Insert a new job with a UUID and return its ID"""
    with get_connection(DB_PATH) as conn:
        c = conn.cursor()
        c.execute(
            "INSERT INTO jobs (id, session_id, status, params, created_at) VALUES (?, ?, ?, ?, datetime('now'))",
            (job_id, session.get("id"), "queued", json.dumps(params_dict))
        )
        conn.commit()
    return job_id

def init_db(reset=True):
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()
    if reset:
        c.execute("DROP TABLE IF EXISTS jobs")
    c.execute("""
        CREATE TABLE IF NOT EXISTS jobs (
            id TEXT PRIMARY KEY,
            session_id TEXT,
            status TEXT NOT NULL,
            params TEXT,
            result TEXT,
            error TEXT, 
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            started_at TIMESTAMP,
            completed_at TIMESTAMP
        )
    """)
    conn.commit()
    conn.close()
    print(f"Standart jobs database initialized ✅ ({DB_PATH})")


@app.route(f"{app.config['ROUTE_PREFIX']}/job_status/<job_id>")
def job_status(job_id):
    with get_connection(DB_PATH) as conn:
        c = conn.cursor()
        c.execute("SELECT status, result, error FROM jobs WHERE id=?", (job_id,))
        row = c.fetchone()

        if not row:
            return jsonify({"error": "job not found"}), 404

        status, result, error = row

        if status in ("queued", "running"):
            return jsonify({"status": status})

        if status == "failed":
            return jsonify({"status": "failed", "error": error})

        # completed
        return jsonify({"status": "completed", "result": result})

@app.route(f"{app.config['ROUTE_PREFIX']}/job_result/<job_id>")
def job_result(job_id):
    with get_connection(DB_PATH) as conn:
        c = conn.cursor()
        c.execute("SELECT status, result FROM jobs WHERE id=?", (job_id,))
        row = c.fetchone()

        if not row:
            return "Job not found", 404

        status, result_json = row
        if status != "completed":
            return f"Job status: {status}", 400

        result_data = json.loads(result_json)
        form = MyForm()

        return render_template(
            'index.html', 
            form=form, 
            job_id=job_id,
            utc_dt=datetime.datetime.utcnow().replace(microsecond=0).isoformat(' '),
            prediction_text0=result_data["prediction_text0"],
            prediction_text1=result_data["prediction_text1"],
            prediction_text2=result_data["prediction_text2"],
            prediction_text3=result_data["prediction_text3"],
            prediction_text4=result_data["prediction_text4"],
            prediction_text5=result_data["prediction_text5"],
            prediction_text6=result_data["prediction_text6"],
            prediction_text7=result_data["prediction_text7"],
            tablesG=result_data["tablesG"],
            tablesC=result_data["tablesC"],
            titlesG=result_data["titlesG"],
            titlesC=result_data["titlesC"],
            prediction_text8=result_data["prediction_text8"],
            prediction_text9=result_data["prediction_text9"],
            prediction_text10=result_data["prediction_text10"],
            prediction_text11=result_data["prediction_text11"],
            download_files=result_data["download_files"],
            plot_image=result_data["plot_image"],
            plot1_url=url_for('static', filename=f"shared/{result_data['plot1_url']}"),
            plot2_url=url_for('static', filename=f"shared/{result_data['plot2_url']}"),
            plot3_url=url_for('static', filename=f"shared/{result_data['plot3_url']}"),
            prediction_text12=result_data["prediction_text12"],
            session_dir=result_data["session_dir"],
            cluster_data_threshold=result_data["cluster_data_threshold"],
            download_figures=result_data["download_figures"]
        )
# end sqlite 
################################################################
# heavy worker for multifasta
@with_retry(retries=5, delay=0.2)
def insert_heavy_job(params_dict, job_id):
    """Insert a new job with a UUID and return its ID"""
    with get_heavy_connection(HEAVY_DB_PATH) as conn:
        c = conn.cursor()
        c.execute(
            "INSERT INTO jobs (id, session_id, status, params, created_at) VALUES (?, ?, ?, ?, datetime('now'))",
            (job_id, session.get("id"), "queued", json.dumps(params_dict))
        )
        conn.commit()
    return job_id


def init_heavy_db(reset=True):
    conn = sqlite3.connect(HEAVY_DB_PATH)
    c = conn.cursor()
    if reset:
        c.execute("DROP TABLE IF EXISTS jobs")
    c.execute("""
        CREATE TABLE IF NOT EXISTS jobs (
            id TEXT PRIMARY KEY,
            session_id TEXT,
            status TEXT NOT NULL,
            params TEXT,
            result TEXT,
            error TEXT, 
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            started_at TIMESTAMP,
            completed_at TIMESTAMP
        )
    """)
    conn.commit()
    conn.close()
    logger.info(f"Heavy jobs database initialized ({HEAVY_DB_PATH})")



@app.route(f"{app.config['ROUTE_PREFIX']}/heavy_job_status/<job_id>")
def heavy_job_status(job_id):
    with get_heavy_connection(HEAVY_DB_PATH) as conn:
        c = conn.cursor()
        c.execute("SELECT status, result, error FROM jobs WHERE id=?", (job_id,))
        row = c.fetchone()

    if not row:
        return jsonify({"error": "heavy_job not found"}), 404

    status, result, error = row

    if status in ("queued", "running"):
        return jsonify({"status": status})

    if status == "failed":
        return jsonify({"status": "failed", "error": error})

    # completed
    return jsonify({"status": "completed", "result": result})

@app.route(f"{app.config['ROUTE_PREFIX']}/heavy_job_result/<job_id>")
def heavy_job_result(job_id):
    with get_heavy_connection(HEAVY_DB_PATH) as conn:
        c = conn.cursor()
        c.execute("SELECT status, result FROM jobs WHERE id=?", (job_id,))
        row = c.fetchone()

        if not row:
            return "Heavy job not found", 404

        status, result_json = row
        if status != "completed":
            return f"Job status: {status}", 400

        result_data = json.loads(result_json)
        form = MyForm()

        return render_template(
            'index.html', 
            form=form, 
            job_id=job_id,
            utc_dt=datetime.datetime.utcnow().replace(microsecond=0).isoformat(' '),
            download_files=result_data["download_files"],
            download_figures=result_data["download_figures"]
        )


################################################################
def random_id(n=32):
    """
    Returns a random ID of a given length
    made of numbers and letters
    """
    return ''.join(choice(string.ascii_letters + string.digits) for _ in range(n))

def vacuum_entry_db():
    conn = sqlite3.connect(ENTRY_DB_PATH)
    c = conn.cursor()

    c.execute("SELECT job_id, output_dir FROM jobs")
    rows = c.fetchall()

    removed = 0
    for job_id, output_dir in rows:
        if not os.path.isdir(output_dir):
            c.execute("DELETE FROM jobs WHERE job_id = ?", (job_id,))
            removed += 1

    conn.commit()
    conn.close()

    if removed:
        logger.info(f"[VACUUM] Removed aged jobs {output_dir} from entry")

def vacuum_thread():
    """
    Cleans up all aged user session directories
    """
    while True:
        for dirname in os.listdir(temp_root):
            if dirname.startswith(temp_signature):
                session_dir=f'{temp_root}/{dirname}'
                start_time_path=f'{session_dir}/{start_time_filename}'
                if os.path.isfile(start_time_path):
                    start_time=datetime.datetime.fromisoformat(open(start_time_path).read())
                    session_age_seconds=(datetime.datetime.utcnow()-start_time).seconds
                    if session_age_seconds>session_duration_minutes*60:
                        rmtree(session_dir)
                        logger.info(f"[VACUUM] Removed aged session {session_dir}")

        # clean entry.db after filesystem cleanup
        vacuum_entry_db()
        time.sleep(vacuum_period_seconds)


@app.route(f"{app.config['ROUTE_PREFIX']}/remove")
def delete_files_in_directory(directory_path_images, directory_path_uploads, directory_path_csvfiles):
   try:
     files_images = os.listdir(directory_path_images)
     files_uploads = os.listdir(directory_path_uploads)
     files_csvfiles = os.listdir(directory_path_csvfiles)
     for file in files_images:
       file_path_images = os.path.join(directory_path_images, file)
       if os.path.isfile(file_path_images):
         os.remove(file_path_images)
     for file in files_uploads:
        file_path_uploads = os.path.join(directory_path_uploads, file)
        if os.path.isfile(file_path_uploads):
         os.remove(file_path_uploads)
     for file in files_csvfiles:
        file_path_csvfiles = os.path.join(directory_path_csvfiles, file)
        if os.path.isfile(file_path_csvfiles):
         os.remove(file_path_csvfiles)
         
   except OSError:
     print("Error occurred while deleting files.")


@app.route(f"{app.config['ROUTE_PREFIX']}/download/<uuid:job_id>")
def download(job_id):
    job_id = str(job_id)

    conn = sqlite3.connect(ENTRY_DB_PATH)
    conn.row_factory = sqlite3.Row
    c = conn.cursor()

    row = c.execute(
        "SELECT output_dir FROM jobs WHERE job_id = ?",
        (job_id,)
    ).fetchone()

    conn.close()
    # print('entrydb_row', row)

    if not row:
        return ('Not found', 404)

    target = row["output_dir"]

    stream = BytesIO()
    if os.path.isdir(target):
        with ZipFile(stream, 'w') as zf:
            for filename in glob(os.path.join(target, '*.csv')):
                zf.write(filename, arcname=os.path.join(os.path.basename(target),
                                                    os.path.basename(filename)))
        stream.seek(0)
        download_name='archive.zip'

    else:
        logger.info(f'EE download(): Session dir not found: {target}')
        stream=BytesIO()
        download_name='archive-empty.zip'

    #print('DD: download(): bytes ',stream.getbuffer().nbytes)
    logger.info(f'download job_id {job_id}')

    return send_file(
        stream,
        as_attachment=True,
        download_name=download_name)

    

@app.route(f"{app.config['ROUTE_PREFIX']}/download_figs/<uuid:job_id>")
def download_figs(job_id):
    job_id = str(job_id)

    conn = sqlite3.connect(ENTRY_DB_PATH)
    conn.row_factory = sqlite3.Row
    c = conn.cursor()

    row = c.execute(
        "SELECT output_dir FROM jobs WHERE job_id = ?",
        (job_id,)
    ).fetchone()

    conn.close()

    if not row:
        return ('Not found', 404)

    target = row["output_dir"]

    stream = BytesIO()
    if os.path.isdir(target):
        with ZipFile(stream, 'w') as zf:
            for ext in ('*.pdf', '*.png'):  # include both png and pdf
                for filename in glob(os.path.join(target, ext)):
                    zf.write(filename, arcname=os.path.join(os.path.basename(target),
                                                    os.path.basename(filename)))

        stream.seek(0)
        download_name = 'archive-figures.zip'

    else:
        logger.info(f'EE download_figs(): Session dir not found: {target}')
        stream=BytesIO()
        download_name='archive-empty.zip'

    # print('DD: download(): bytes ',stream.getbuffer().nbytes)
    logger.info(f'download_figs job_id {job_id}')
    return send_file(
        stream,
        as_attachment=True,
        download_name=download_name
        )


# Root of temporary directories (for user sessions)
# Make sure it's writeable for the Flask process
temp_root='static/shared'
temp_signature='fuzzy'
session_duration_minutes=120
vacuum_period_seconds=60
start_time_filename='start_time.txt'

if not os.path.isdir(temp_root):
    os.mkdir(temp_root)

# Start the vacuum cleaner to watch and
# clean up user session data on the server
vacuum=threading.Thread(target=vacuum_thread)
vacuum.start()

@app.before_request
def make_session_permanent():
    # Set the sessions maximal duration of 60 minutes
    # Within this time, the user can log in and out
    # of the session and keep the data preserved.
    session.permanent = True
    app.permanent_session_lifetime = datetime.timedelta(minutes=session_duration_minutes)


@app.route(f"{app.config['ROUTE_PREFIX']}/about", methods=['GET', 'POST'])
def about():
    return render_template("about.html")


def complete_prediction_fasta(filename, patternG_obj, patternC_obj, clusterscore, 
                              input_looplength, print_indexes, print_sequence, print_clusters, 
                              csv_file_repeats, csv_file_clusters, csv_file_repeatsG, csv_file_repeatsC, 
                              csv_file_clustersG, csv_file_clustersC, input_ssr, input_repeats, tandem_select, selected_option, fasta_id, session_id, session_dir):    
    # print('Welcome to complete_prediction_fasta')
    # session_id=session['id']
    # session_dir = f'{temp_root}/{temp_signature}-{session_id}'
    session_id=session_id
    session_dir=session_dir

    # We need to initialize the values, in case if
    # the user does not select any outputs,
    # and we don't want to trigger a "Variable undefined' error:
        
    print_indexesG = 'Printing is turned off or the sequence is too large'
    print_indexesC = 'Printing is turned off or the sequence is too large'
    print_sequenceG = 'Printing is turned off or the sequence is too large'
    print_sequenceC = 'Printing is turned off or the sequence is too large'
    print_clustersG = 'Printing is turned off or the sequence is too large'
    print_clustersC = 'Printing is turned off or the sequence is too large'
    titlesG='No titles'
    titlesC='No titles'
    if selected_option == 'Custom pattern':
        orientationF = 'F'
        orientationR = 'R'
    else: 
        orientationF = 'G'
        orientationR = 'C'
    plot_image='No plot'
    plot=None
    cluster_fig=None
    loop_medianG=0.0
    loop_medianC=0.0
    clusterC=pd.DataFrame()
    clusterG=pd.DataFrame()
    # print('collected initial variables')
    for entry in fasta_reader(filename):
        # print('fasta entry', entry)
        seqname = entry.id
        seq = entry.seq
        description = entry.description
        
        if fasta_id is None or seqname == fasta_id:
            # slice_start = 200465808
            # slice_end = 200469048
            # seq = seq[slice_start:slice_end]
            seq = str(seq).replace(' ', '').upper().replace('U', 'T')
            lenseq = len(seq)
            resultG = search_functionGC(seq, patternG_obj)
            resultC = search_functionGC(seq, patternC_obj)

            patternG_span_open='<span style="background-color: #79B6FF">'
            patternC_span_open='<span style="background-color: #8cfb9b">'
            span_close='</span>'
            clusterG_span = '<span style="background-color: #79B6FF">\g<0></span>'
            clusterC_span = '<span style="background-color: #8cfb9b">\g<0></span>'
            cluster_span_open = '<span style="color:red;font-weight:bold">'
            indG = resultG[0]
            indC = resultC[0]
            ind_startG = resultG[1]
            ind_endG = resultG[2]
            ind_startC = resultC[1]
            ind_endC = resultC[2]
            spans_medianG = resultG[3]
            spans_medianC = resultC[3]
            highlightedG = highlight_functionGC(seq, patternG_span_open, span_close, indG)
            highlightedC = highlight_functionGC(seq, patternC_span_open, span_close, indC)

            
            #######################################################################################
            loop_plotG = loop_function(seq, ind_startG, ind_endG)
            loop_plotC = loop_function(seq, ind_startC, ind_endC)
            loop_indG = loop_plotG[0]
            loop_listG = loop_plotG[1]
            loop_list_arrG = loop_plotG[2]
            loop_indC = loop_plotC[0]
            loop_listC = loop_plotC[1]
            loop_list_arrC = loop_plotC[2]
            loop_medianG = loop_plotG[3]
            loop_medianC = loop_plotC[3]

            clusterG_result = cluster_functionGC(seq, seqname, clusterscore, ind_startG, ind_endG, loop_list_arrG, input_looplength, loop_medianG, spans_medianG, input_ssr, input_repeats, tandem_select)
            clusterC_result = cluster_functionGC(seq, seqname, clusterscore, ind_startC, ind_endC, loop_list_arrC, input_looplength, loop_medianC, spans_medianC, input_ssr, input_repeats, tandem_select)

            clusterG = clusterG_result[0]
            clusterC = clusterC_result[0]
            loop_estimateG = clusterG_result[1]
            loop_estimateC = clusterC_result[1]

            plot = image_functionGC(seq, loop_indG, loop_listG, loop_list_arrG, loop_indC, loop_listC, loop_list_arrC, description, loop_medianG, loop_medianC,session_dir, selected_option)
            plot_image = 'plot'
            image1 = plot[0]
            image2 = plot[1]

            cluster_arrayG = clusterG.iloc[:, [1,2,3]]
            cluster_arrayC = clusterC.iloc[:, [1,2,3]]

            clusterG_loop = cluster_loopsGC(seq, cluster_arrayG)
            clusterC_loop = cluster_loopsGC(seq, cluster_arrayC)

            clusterG_list = clusterG_loop[0]
            clusterG_ind = clusterG_loop[1]
            clusterC_list = clusterC_loop[0]
            clusterC_ind = clusterC_loop[1]
            cluster_fig = cluster_plot(seq, clusterG_list, clusterG_ind, clusterC_list, clusterC_ind, description,session_dir, selected_option)
            clusterG_indexes = clusterG_loop[2]
            clusterC_indexes = clusterC_loop[2]
            highlighted_clusterG = highlighted_clusterGC(seq, clusterG_indexes, clusterG_list, clusterG_span, cluster_span_open, span_close, patternG_obj)
            highlighted_clusterC = highlighted_clusterGC(seq, clusterC_indexes, clusterC_list, clusterC_span, cluster_span_open, span_close, patternC_obj)


            titlesG = (f'Clusters {orientationF}-orientation, cluster score = {clusterscore}, loop length = {loop_estimateG}, SSR = {input_ssr}, Number of repeats = {input_repeats}')
            titlesC = (f'Clusters {orientationR}-orientation, cluster score = {clusterscore}, loop length = {loop_estimateC}, SSR = {input_ssr}, Number of repeats = {input_repeats}')

            
            if print_indexes and lenseq <= 10000:
                print_indexesG = indG
                print_indexesC = indC
            else:
                print_indexesG = 'Printing is turned off or the sequence is too large'
                print_indexesC = 'Printing is turned off or the sequence is too large'

                
            if print_sequence and lenseq <= 10000:
                print_sequenceG = highlightedG
                print_sequenceC = highlightedC
            else:
                print_sequenceG = 'Printing is turned off or the sequence is too large'
                print_sequenceC = 'Printing is turned off or the sequence is too large'

            if print_clusters and lenseq <= 10000:
                print_clustersG = highlighted_clusterG
                print_clustersC = highlighted_clusterC
            else:
                print_clustersG = 'Printing is turned off or the sequence is too large'
                print_clustersC = 'Printing is turned off or the sequence is too large'


            # print('clustersG', len(clusterG))
            
            # 250523 limit rendering tables for data with many clusters
            # Labelling of G/F and C/R in csv files is F/R did't input condition on telomeric
            cluster_threshold = len(clusterG) >=200 or len(clusterC) >=200
            # print('cluster_threshold', cluster_threshold)
        
            repeatsG = repeats_dataframeGC(seqname, resultG[1], resultG[2])
            repeatsC = repeats_dataframeGC(seqname, resultC[1], resultC[2])

            if csv_file_repeats is True:
                saved_repeats = save_repeats_csvGC(repeatsG, repeatsC, csv_file_repeatsG, csv_file_repeatsC)
            else:
                pass
    
            #if should_use_clusters_csv is True:
            if csv_file_clusters is True:
                saved_clusters = save_clusters_csvGC(clusterG, clusterC, csv_file_clustersG, csv_file_clustersC)
                        
            else:
                pass


            if csv_file_repeats is True or csv_file_clusters is True or cluster_threshold is True:
            #if csv_file_repeats is True or should_use_clusters_csv is True:
                # download()
                download_files = 1
            else:
                download_files = 0
        
            return  (print_indexesG, print_indexesC, print_sequenceG, print_sequenceC, clusterG, clusterC, 
                titlesG, titlesC, print_clustersG, print_clustersC, download_files, plot_image, image1, image2, cluster_fig, loop_medianG, loop_medianC, cluster_threshold)
        else:
            pass
###################################################################################################################################
def complete_prediction_multifasta(filename, patternG_obj, patternC_obj, clusterscore, 
                              input_looplength, csv_file_repeats, csv_file_clusters, 
                              csv_file_repeatsG, csv_file_repeatsC, csv_file_clustersG, csv_file_clustersC, input_ssr, input_repeats, tandem_select, selected_option, session_id, session_dir):    
    loops_map = []
    loops_hist = []
    clusters_map = []

    # session_id=session['id']
    # session_dir = f'{temp_root}/{temp_signature}-{session_id}'
    session_id=session_id
    session_dir=session_dir

    pdf_loops_map = f'{session_dir}/loops_map.pdf'
    pdf_loops_hist = f'{session_dir}/loops_hist.pdf'
    pdf_clusters = f'{session_dir}/clusters_map_cluster_{clusterscore}_loop_{input_looplength}_ssr_{input_ssr}_nrepeats_{input_repeats}.pdf'
    

    for entry in fasta_reader(filename):
        seqname = entry.id
        seq = entry.seq
        description = entry.description
        # print(description)
        seq = str(seq).replace(' ', '').upper().replace('U', 'T')
        resultG = search_functionGC(seq, patternG_obj)
        resultC = search_functionGC(seq, patternC_obj)

        indG = resultG[0]
        indC = resultC[0]
        ind_startG = resultG[1]
        ind_endG = resultG[2]
        ind_startC = resultC[1]
        ind_endC = resultC[2]
        spans_medianG = resultG[3]
        spans_medianC = resultC[3]
        
        
        
        #######################################################################################
        loop_plotG = loop_function(seq, ind_startG, ind_endG)
        loop_plotC = loop_function(seq, ind_startC, ind_endC)
        loop_indG = loop_plotG[0]
        loop_listG = loop_plotG[1]
        loop_list_arrG = loop_plotG[2]
        loop_indC = loop_plotC[0]
        loop_listC = loop_plotC[1]
        loop_list_arrC = loop_plotC[2]
        loop_medianG = loop_plotG[3]
        loop_medianC = loop_plotC[3]

        clusterG_result = cluster_functionGC(seq, seqname, clusterscore, ind_startG, ind_endG, loop_list_arrG, input_looplength, loop_medianG, spans_medianG, input_ssr, input_repeats, tandem_select)
        clusterC_result = cluster_functionGC(seq, seqname, clusterscore, ind_startC, ind_endC, loop_list_arrC, input_looplength, loop_medianC, spans_medianC, input_ssr, input_repeats, tandem_select)

        clusterG = clusterG_result[0]
        clusterC = clusterC_result[0]
        loop_estimateG = clusterG_result[1]
        loop_estimateC = clusterC_result[1]


        cluster_arrayG = clusterG.iloc[:, [1,2,3]]
        cluster_arrayC = clusterC.iloc[:, [1,2,3]]

        clusterG_loop = cluster_loopsGC(seq, cluster_arrayG)
        clusterC_loop = cluster_loopsGC(seq, cluster_arrayC)

        clusterG_list = clusterG_loop[0]
        clusterG_ind = clusterG_loop[1]
        clusterC_list = clusterC_loop[0]
        clusterC_ind = clusterC_loop[1]
        repeatsG = repeats_dataframeGC(seqname, resultG[1], resultG[2])
        repeatsC = repeats_dataframeGC(seqname, resultC[1], resultC[2])
        saved_repeats = save_repeats_csvGC(repeatsG, repeatsC, csv_file_repeatsG, csv_file_repeatsC)
        saved_clusters = save_clusters_csvGC(clusterG, clusterC, csv_file_clustersG, csv_file_clustersC)

        figure1 = image_loops_multifasta(seq, loop_indG, loop_listG, loop_indC, loop_listC, description, loops_map, selected_option)
        figure2 = image_hist_loops_multifasta(loop_list_arrG, loop_list_arrC, description, loops_hist, loop_medianG, loop_medianC, selected_option)
        figure3 = image_clusters_multifasta(seq, clusterG_list, clusterG_ind, clusterC_list, clusterC_ind, description, clusters_map, selected_option)
        # print('len figure1', len(figure1))

    # download()
    download_files = 1 

    save_image(pdf_loops_map, figure1)
    save_image(pdf_loops_hist, figure2)
    save_image(pdf_clusters, figure3)


    # download_figs()
    download_figures = 1

    return  download_files, download_figures, session_dir
####################################################################################################################################

@app.route(f"{app.config['ROUTE_PREFIX']}/complete_prediction_simple")
def complete_prediction_simple(seqname, seq, patternG_obj, patternC_obj, clusterscore, 
                               input_looplength, print_indexes, print_sequence, print_clusters, 
                               csv_file_repeats, csv_file_clusters, csv_file_repeatsG, csv_file_repeatsC, 
                               csv_file_clustersG, csv_file_clustersC, input_ssr, input_repeats, tandem_select, selected_option, session_id, session_dir):    
    # session_id=session['id']
    # session_dir = f'{temp_root}/{temp_signature}-{session_id}'
    session_id=session_id
    session_dir=session_dir

    
    seqname = seqname
    seq = seq
    description = seqname
            
    print_indexesG = 'Printing is turned off or the sequence is too large'
    print_indexesC = 'Printing is turned off or the sequence is too large'
    print_sequenceG = 'Printing is turned off or the sequence is too large'
    print_sequenceC = 'Printing is turned off or the sequence is too large'
    print_clustersG = 'Printing is turned off or the sequence is too large'
    print_clustersC = 'Printing is turned off or the sequence is too large'
    titlesG='No titles'
    titlesC='No titles'
    if selected_option == 'Custom pattern':
        orientationF = 'F'
        orientationR = 'R'
    else: 
        orientationF = 'G'
        orientationR = 'C'
    plot_image='No plot'
    plot=None
    cluster_fig=None
    loop_medianG=0.0
    loop_medianC=0.0
    clusterG=pd.DataFrame()
    clusterC=pd.DataFrame()

    seq = str(seq).replace(' ', '').upper().replace('U', 'T')
    lenseq = len(seq)
    resultG = search_functionGC(seq, patternG_obj)
    resultC = search_functionGC(seq, patternC_obj)

    patternG_span_open='<span style="background-color: #79B6FF">'
    patternC_span_open='<span style="background-color: #8cfb9b">'
    span_close='</span>'
    clusterG_span = '<span style="background-color: #79B6FF">\g<0></span>'
    clusterC_span = '<span style="background-color: #8cfb9b">\g<0></span>'
    cluster_span_open = '<span style="color:red;font-weight:bold">'
    indG = resultG[0]
    indC = resultC[0]
    ind_startG = resultG[1]
    ind_endG = resultG[2]
    ind_startC = resultC[1]
    ind_endC = resultC[2]
    spans_medianG = resultG[3]
    spans_medianC = resultC[3]
    highlightedG = highlight_functionGC(seq, patternG_span_open, span_close, indG)
    highlightedC = highlight_functionGC(seq, patternC_span_open, span_close, indC)

    #######################################################################################
    loop_plotG = loop_function(seq, ind_startG, ind_endG)
    loop_plotC = loop_function(seq, ind_startC, ind_endC)
    loop_indG = loop_plotG[0]
    loop_listG = loop_plotG[1]
    loop_list_arrG = loop_plotG[2]
    loop_indC = loop_plotC[0]
    loop_listC = loop_plotC[1]
    loop_list_arrC = loop_plotC[2]
    loop_medianG = loop_plotG[3]
    loop_medianC = loop_plotC[3]


    clusterG_result = cluster_functionGC(seq, seqname, clusterscore, ind_startG, ind_endG, loop_list_arrG, input_looplength, loop_medianG, spans_medianG, input_ssr, input_repeats, tandem_select)
    clusterC_result = cluster_functionGC(seq, seqname, clusterscore, ind_startC, ind_endC, loop_list_arrC, input_looplength, loop_medianC, spans_medianC, input_ssr, input_repeats, tandem_select)

    clusterG = clusterG_result[0]
    clusterC = clusterC_result[0]
    loop_estimateG = clusterG_result[1]
    loop_estimateC = clusterC_result[1]

    #######################################################################################

    plot = image_functionGC(seq, loop_indG, loop_listG, loop_list_arrG, loop_indC, loop_listC, loop_list_arrC, description, loop_medianG, loop_medianC,session_dir, selected_option)
    plot_image = 'plot'
    image1 = plot[0]
    image2 = plot[1]

    cluster_arrayG = clusterG.iloc[:, [1,2,3]]
    cluster_arrayC = clusterC.iloc[:, [1,2,3]]

    clusterG_loop = cluster_loopsGC(seq, cluster_arrayG)
    clusterC_loop = cluster_loopsGC(seq, cluster_arrayC)

    clusterG_list = clusterG_loop[0]
    clusterG_ind = clusterG_loop[1]
    clusterC_list = clusterC_loop[0]
    clusterC_ind = clusterC_loop[1]
    cluster_fig = cluster_plot(seq, clusterG_list, clusterG_ind, clusterC_list, clusterC_ind, description,session_dir, selected_option)
    clusterG_indexes = clusterG_loop[2]
    clusterC_indexes = clusterC_loop[2]
    highlighted_clusterG = highlighted_clusterGC(seq, clusterG_indexes, clusterG_list, clusterG_span, cluster_span_open, span_close, patternG_obj)
    highlighted_clusterC = highlighted_clusterGC(seq, clusterC_indexes, clusterC_list, clusterC_span, cluster_span_open, span_close, patternC_obj)
    
    titlesG = (f'Clusters {orientationF}-orientation, cluster score = {clusterscore}, loop length = {loop_estimateG}, SSR = {input_ssr}, Number of repeats = {input_repeats}')
    titlesC = (f'Clusters {orientationR}-orientation, cluster score = {clusterscore}, loop length = {loop_estimateC}, SSR = {input_ssr}, Number of repeats = {input_repeats}')
    cluster_threshold = len(clusterG) >=200 or len(clusterC) >=200
    
    if print_indexes is True and lenseq <= 10000:
        print_indexesG = indG
        print_indexesC = indC

    else:
        print_indexesG = 'Printing is turned off or the sequence is too large'
        print_indexesC = 'Printing is turned off or the sequence is too large'
        
    if print_sequence is True and lenseq <= 10000:
        print_sequenceG = highlightedG
        print_sequenceC = highlightedC

    else:
        print_sequenceG = 'Printing is turned off or the sequence is too large'
        print_sequenceC = 'Printing is turned off or the sequence is too large'

    if print_clusters is True and lenseq <= 10000:
        print_clustersG = highlighted_clusterG
        print_clustersC = highlighted_clusterC

    else:
        print_clustersG = 'Printing is turned off or the sequence is too large'
        print_clustersC = 'Printing is turned off or the sequence is too large'

    if csv_file_repeats is True:
        repeatsG = repeats_dataframeGC(seqname, resultG[1], resultG[2])
        repeatsC = repeats_dataframeGC(seqname, resultC[1], resultC[2])
        saved_repeats = save_repeats_csvGC(repeatsG, repeatsC, csv_file_repeatsG, csv_file_repeatsC)
    else:
        pass

    if csv_file_clusters is True:
        saved_clusters = save_clusters_csvGC(clusterG, clusterC, csv_file_clustersG, csv_file_clustersC)
                    
    else:
        pass

    if csv_file_repeats is True or csv_file_clusters is True:
        # download()
        download_files = 1
    else:
        download_files = 0
    

    return  (print_indexesG, print_indexesC, print_sequenceG, print_sequenceC, clusterG, clusterC, 
        titlesG, titlesC, print_clustersG, print_clustersC, download_files, plot_image, image1, image2, cluster_fig, loop_medianG, loop_medianC, cluster_threshold)


@app.route(f"{app.config['ROUTE_PREFIX']}", methods=['GET', 'POST'])
@app.route(f"{app.config['ROUTE_PREFIX']}/", methods=['GET', 'POST'])
@app.route(f"{app.config['ROUTE_PREFIX']}/predictGC", methods=['GET', 'POST'])
def predictGC():
    
    form = MyForm()
    
    if form.validate_on_submit():
        logger.info("Form validated")

        session['id'] = random_id()
        job_id = str(uuid.uuid4()) 

        session_id = session['id']
        session_dir = f"{temp_root}/{temp_signature}-{session_id}"
        os.mkdir(session_dir)
        session['start']=datetime.datetime.utcnow()
        open(f'{session_dir}/{start_time_filename}','w').write(session['start'].strftime('%Y-%m-%dT%H:%M:%S'))
        register_job(job_id, session_dir) 

        is_seq_filled = bool(form.seq.data)
        is_choicegenome_filled = bool(form.choicegenome.data and form.choicegenome.data != 'default') 
        is_file_filled = bool(form.file.data)
        filled_fields = sum([is_seq_filled, is_choicegenome_filled, is_file_filled])
        #filled_fields = sum(bool(field) for field in [form.seq.data, form.choicegenome.data, form.file.data])
        # print('filled_fields', filled_fields)

        #if filled_fields != 1:
        if filled_fields > 1:
            flash("Please, choose only one option: enter the sequence into the text field, select a genome, or upload a file", "warning")
            return render_template('index.html', utc_dt=datetime.datetime.utcnow().replace(microsecond=0).isoformat(' '), form=form)
        # Check that at least one of the fields has data
        if filled_fields < 1:
            flash("Please, enter the sequence into the text field, select a genome, or upload a file", "warning")
            return render_template('index.html', utc_dt=datetime.datetime.utcnow().replace(microsecond=0).isoformat(' '), form=form)

    
        else:      
            file = form.file.data if is_file_filled else None
            seq = form.seq.data if is_seq_filled else None
            selected_genome = form.choicegenome.data if is_choicegenome_filled else None
            seqname = form.seqname.data
            chromosome = form.chromosome.data
            region = form.region.data
            
            csv_file_repeats = form.csvrepeats.data
            csv_file_clusters = form.csvclusters.data
            print_sequence = form.printseq.data
            print_clusters = form.printclust.data
            print_indexes = form.printind.data
            
            user_clusterscore = form.clusterscore.data
            user_looplength = form.input_looplength.data
            tandem_repeats = form.tandem_repeats.data
            user_ssr = form.ssr.data
            repeats_number = form.repeats_number.data


            if not user_clusterscore:
                clusterscore = 'undefined'
            else:
                clusterscore = float(user_clusterscore)

            if not user_ssr:
                input_ssr = 'undefined'
            else:
                input_ssr = float(user_ssr)

            if not repeats_number:
                input_repeats = 'any'
            else:
                input_repeats = repeats_number
            
            if tandem_repeats is True:
                input_looplength = 0
                input_repeats = 2
                tandem_select = True
            
            else:
                if not user_looplength:
                    input_looplength = 'estimate'
                    tandem_select = False
                else:
                    # input_looplength = int(input_looplength.strip())
                    input_looplength = user_looplength
                    tandem_select = False
            

            selected_option = form.choice.data

            if selected_option == 'TTAGGG':
                patternG = "T{2}A{1}G{3}"
                patternC = "C{3}T{1}A{2}"
                pattern_printG = 'TTAGGG'
                pattern_printC = 'CCCTAA'
                
            elif selected_option == 'FuzzyTel':
                patternG = "T{1,2}A{0,1}G{3,5}|T{1,2}\D{1}A{1}G{3,5}|T{2}A{1}G{2}"
                patternC = "C{3,5}T{0,1}A{1,2}|C{3,5}T{1}\D{1}A{1,2}|C{2}T{1}A{2}"
                pattern_printG = 'FuzzyTel G-orientation'
                pattern_printC = 'FuzzyTel C-orientation' 

            elif selected_option == 'Custom pattern':
                pattern = form.pattern.data
                if not pattern:
                    enter_the_pattern = 'Please, enter the pattern'
                    return render_template('index.html', utc_dt=datetime.datetime.utcnow().replace(microsecond=0).isoformat(' '), enter_the_pattern = enter_the_pattern, form=form)
                else:
                    patternG = str(form.pattern.data).strip().upper()
                    patternC = str(patternG[::-1].translate(str.maketrans('ATCG','TAGC')))
                    pattern_printG = patternG
                    pattern_printC = patternC

            ###############################################################################

            
            # print('selected_option', selected_option)
            patternG_obj = re.compile(patternG)
            patternC_obj = re.compile(patternC)

            ######################################################################################
            chrom_names_T2T = {'chr1': 'NC_060925.1', 
                              'chr2': 'NC_060926.1', 
                              'chr3': 'NC_060927.1', 
                              'chr4': 'NC_060928.1',
                              'chr5': 'NC_060929.1', 
                              'chr6': 'NC_060930.1', 
                              'chr7': 'NC_060931.1', 
                              'chr8': 'NC_060932.1',
                              'chr9': 'NC_060933.1', 
                              'chr10': 'NC_060934.1', 
                              'chr11': 'NC_060935.1', 
                              'chr12': 'NC_060936.1',
                              'chr13': 'NC_060937.1', 
                              'chr14': 'NC_060938.1', 
                              'chr15': 'NC_060939.1', 
                              'chr16': 'NC_060940.1',
                              'chr17': 'NC_060941.1', 
                              'chr18': 'NC_060942.1', 
                              'chr19': 'NC_060943.1', 
                              'chr20': 'NC_060944.1',
                              'chr21': 'NC_060945.1', 
                              'chr22': 'NC_060946.1', 
                              'chr23': 'NC_060947.1', 
                              'chr24': 'NC_060948.1'}
            
            chrom_names_h38 = {'chr1': 'NC_000001.11', 
                              'chr2': 'NC_000002.12', 
                              'chr3': 'NC_000003.12', 
                              'chr4': 'NC_000004.12',
                              'chr5': 'NC_000005.10', 
                              'chr6': 'NC_000006.12', 
                              'chr7': 'NC_000007.14', 
                              'chr8': 'NC_000008.11',
                              'chr9': 'NC_000009.12', 
                              'chr10': 'NC_000010.11', 
                              'chr11': 'NC_000011.10', 
                              'chr12': 'NC_000012.12',
                              'chr13': 'NC_000013.11', 
                              'chr14': 'NC_000014.9', 
                              'chr15': 'NC_000015.10', 
                              'chr16': 'NC_000016.10',
                              'chr17': 'NC_000017.11', 
                              'chr18': 'NC_000018.10', 
                              'chr19': 'NC_000019.10', 
                              'chr20': 'NC_000020.11',
                              'chr21': 'NC_000021.9', 
                              'chr22': 'NC_000022.11', 
                              'chr23': 'NC_000023.11', 
                              'chr24': 'NC_000024.10'}
            #######################################################################################
            ######################################################################################
            #The main procedure
            
            if not seq:
            
                if form.file.data:
                    filename = form.file.data
                    file.save(os.path.join(f'{session_dir}', secure_filename(f'{file.filename}.fasta')))

                    filenames = glob(f'{session_dir}/*.fasta')
                    filenames_ctime_sorted=sorted( [ (f,os.path.getctime(f)) for f in filenames], key=lambda x: x[1], reverse=True)
                    # Select the latest file
                    latest_filename_ctime=filenames_ctime_sorted[0]
                    # Take only the filename and discard the ctime
                    filename=latest_filename_ctime[0]
                    # print('Selected FASTA:',filename)

                    # Verify that the file is in FASTA format
                    with open(filename) as f:
                        lines = f.read().strip().splitlines()
                        # print(lines)

                    if len(lines) < 2 or not lines[0].startswith(">"):
                        flash("Uploaded file doesn’t seem to be in FASTA format (missing '>' header line)", "warning")
                        return render_template('index.html', utc_dt=datetime.datetime.utcnow().replace(microsecond=0).isoformat(' '), form=form)
                    
                    elif len(lines[1]) < 30:
                        flash("Uploaded nucleotide sequence is too short, at least 30 is required", "warning")
                        return render_template('index.html', utc_dt=datetime.datetime.utcnow().replace(microsecond=0).isoformat(' '), form=form)

                    else:
                        num = len([1 for line in lines if line.startswith(">")])
                   
                    if num > 17:
                        if seqname == 'codetel':
                            num = num
                        else:
                            # flash("Please, upload only one sequence: your file contains several", "warning")
                            flash("Analysis of up to 17 sequence entries is supported at present, your file contains more", "warning")
                            return render_template('index.html', utc_dt=datetime.datetime.utcnow().replace(microsecond=0).isoformat(' '), form=form)
                    else:
                        # num = 1
                        num = num


                else:
                    is_choicegenome_filled = bool(form.choicegenome.data not in [None, '', 'default'])
                    # print("Choice_genome", is_choicegenome_filled)
                    if is_choicegenome_filled:
                    #selected_genome = form.choicegenome.data
                        if selected_genome == 'CHM13-T2T':
                            #filename = 'static/genomes/human_chr1_fragment.fasta'
                            filename = 'static/genomes/2025_Genomes_for_fuzzycluster/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna'
            
                        elif selected_genome == 'GRCh38':
                            #filename = 'static/genomes/influenza_short_h4_Telinsert.fasta'
                            filename = 'static/genomes/2025_Genomes_for_fuzzycluster/GCF_000001405.40_GRCh38.p14_genomic.fna'  
                        
                        if seqname == 'codetel':
                            num = len([1 for line in open(filename) if line.startswith(">")])
                        else:
                            num = 1
            
        
                if num > 1:
                    csv_file_repeats = True
                    csv_file_clusters = True
                    repeats_files = csvGC_repeats_open(pattern_printG, pattern_printC,session_dir, selected_option)
                    csv_file_repeatsG = repeats_files[0]
                    csv_file_repeatsC = repeats_files[1]
                    cluster_files = csvGC_clusters_open(pattern_printG, pattern_printC, clusterscore, input_looplength, input_ssr, input_repeats,session_dir, selected_option)
                    csv_file_clustersG = cluster_files[0]
                    csv_file_clustersC = cluster_files[1]

                    params_dict = {
                            "filename": filename,
                            "patternG_obj": patternG_obj.pattern,
                            "patternC_obj": patternC_obj.pattern,
                            "clusterscore": clusterscore,
                            "input_looplength": input_looplength,
                            "csv_file_repeats": csv_file_repeats,
                            "csv_file_clusters": csv_file_clusters,
                            "csv_file_repeatsG": csv_file_repeatsG,
                            "csv_file_repeatsC": csv_file_repeatsC,
                            "csv_file_clustersG": csv_file_clustersG,
                            "csv_file_clustersC": csv_file_clustersC,
                            "input_ssr": input_ssr,
                            "input_repeats": input_repeats,
                            "tandem_select": tandem_select,
                            "selected_option": selected_option,
                            "session_id": session_id, 
                            "session_dir": session_dir,
                            "pattern_printG": pattern_printG,
                            "pattern_printC": pattern_printC
                        }
                    # print('filename', filename)
                    job_id = insert_heavy_job(params_dict, job_id)
                    print(f"Queued heavy job {job_id}")
                    logger.info(f"Queued heavy job {job_id}")
                    return render_template("heavy_waiting.html", route_prefix=app.config['ROUTE_PREFIX'], job_id=job_id)
            


                else:
                    if is_choicegenome_filled:
                        if not chromosome:
                            flash("Please, choose the chromosome", "warning")
                            # print(get_flashed_messages(with_categories=True))
                            return render_template('index.html', utc_dt=datetime.datetime.utcnow().replace(microsecond=0).isoformat(' '), form=form)
                        elif chromosome and selected_genome == 'CHM13-T2T':
                            # print('chromosome', chromosome)
                            num_chromosome = str(chromosome)
                            frontend_chrom = 'chr'+ num_chromosome
                            fasta_id = chrom_names_T2T.get(frontend_chrom)
                            # print('fasta_id', fasta_id)
                        elif chromosome and selected_genome == 'GRCh38':
                            # print('chromosome', chromosome)
                            num_chromosome = str(chromosome)
                            frontend_chrom = 'chr'+ num_chromosome
                            fasta_id = chrom_names_h38.get(frontend_chrom)
                            # print('fasta_id', fasta_id)
                    else:
                        fasta_id = None

                    # I force here generation of all files
                    csv_file_repeats = True
                    csv_file_clusters = True

                    if csv_file_repeats is True:
                        repeats_files = csvGC_repeats_open(pattern_printG, pattern_printC,session_dir, selected_option)
                        csv_file_repeatsG = repeats_files[0]
                        csv_file_repeatsC = repeats_files[1]
                        
                    else:
                        csv_file_repeatsG = False
                        csv_file_repeatsC = False 

                    if csv_file_clusters is True:
                        cluster_files = csvGC_clusters_open(pattern_printG, pattern_printC, clusterscore, input_looplength, input_ssr, input_repeats,session_dir, selected_option)
                        csv_file_clustersG = cluster_files[0]
                        csv_file_clustersC = cluster_files[1]

                    else:
                        csv_file_clustersG = False
                        csv_file_clustersC = False 

                    params_dict = {
                            "filename": filename,
                            "patternG_obj": patternG_obj.pattern,
                            "patternC_obj": patternC_obj.pattern,
                            "clusterscore": clusterscore,
                            "input_looplength": input_looplength,
                            "print_indexes": print_indexes,
                            "print_sequence": print_sequence,
                            "print_clusters": print_clusters,
                            "csv_file_repeats": csv_file_repeats,
                            "csv_file_clusters": csv_file_clusters,
                            "csv_file_repeatsG": csv_file_repeatsG,
                            "csv_file_repeatsC": csv_file_repeatsC,
                            "csv_file_clustersG": csv_file_clustersG,
                            "csv_file_clustersC": csv_file_clustersC,
                            "input_ssr": input_ssr,
                            "input_repeats": input_repeats,
                            "tandem_select": tandem_select,
                            "selected_option": selected_option,
                            "fasta_id": fasta_id,
                            "session_id": session_id, 
                            "session_dir": session_dir, 
                            "pattern_printG": pattern_printG,
                            "pattern_printC": pattern_printC
                        }


                    job_id = insert_job(params_dict, job_id)
                    print(f"Queued job {job_id}")
                    logger.info(f"Queued job {job_id}")
                    return render_template("waiting.html", route_prefix=app.config['ROUTE_PREFIX'], job_id=job_id)
                    

            else:
                if len(seq) < 30:
                    enter_the_sequence = 'Please, enter the sequence of at least 30 nt, or choose a genome or upload the file'
                    return render_template('index.html', utc_dt=datetime.datetime.utcnow().replace(microsecond=0).isoformat(' '), enter_the_sequence = enter_the_sequence, form=form)

                else:
                    prediction_allow = False
                    num_gt = seq.strip().count('>')
                    if num_gt > 1:
                        flash("Please, enter only one sequence, letters A, T, G, C", "warning")
                        return render_template('index.html', utc_dt=datetime.datetime.utcnow().replace(microsecond=0).isoformat(' '), form=form)
                    else:
                        prediction_allow = True

                    if prediction_allow == True:

                        if selected_option == 'Custom pattern':
                            orientationF = 'F'
                            orientationR = 'R'
                        else: 
                            orientationF = 'G'
                            orientationR = 'C'

                        csv_file_clustersG = False
                        csv_file_clustersC = False
                        csv_file_repeatsG = False
                        csv_file_repeatsC = False
                        if csv_file_repeats is True:
                            repeats_files = csvGC_repeats_open(pattern_printG, pattern_printC,session_dir, selected_option)
                            csv_file_repeatsG = repeats_files[0]
                            csv_file_repeatsC = repeats_files[1]

                        #Corrected here 240717
                        if csv_file_clusters is True:
                            cluster_files = csvGC_clusters_open(pattern_printG, pattern_printC, clusterscore, input_looplength, input_ssr, input_repeats,session_dir, selected_option)
                            csv_file_clustersG = cluster_files[0]
                            csv_file_clustersC = cluster_files[1]
                        

                        params_dict = {
                            "seqname": seqname,
                            "seq": seq,
                            "patternG_obj": patternG_obj.pattern,
                            "patternC_obj": patternC_obj.pattern,
                            "clusterscore": clusterscore,
                            "input_looplength": input_looplength,
                            "print_indexes": print_indexes,
                            "print_sequence": print_sequence,
                            "print_clusters": print_clusters,
                            "csv_file_repeats": csv_file_repeats,
                            "csv_file_clusters": csv_file_clusters,
                            "csv_file_repeatsG": csv_file_repeatsG,
                            "csv_file_repeatsC": csv_file_repeatsC,
                            "csv_file_clustersG": csv_file_clustersG,
                            "csv_file_clustersC": csv_file_clustersC,
                            "input_ssr": input_ssr,
                            "input_repeats": input_repeats,
                            "tandem_select": tandem_select,
                            "selected_option": selected_option,
                            "session_id": session_id, 
                            "session_dir": session_dir, 
                            "pattern_printG": pattern_printG,
                            "pattern_printC": pattern_printC
                        }

                        job_id = insert_light_job(params_dict, job_id)
                        print(f"Queued light job {job_id}")
                        logger.info(f"Queued light job {job_id}")
                        return render_template("light_waiting.html", route_prefix=app.config['ROUTE_PREFIX'], job_id=job_id)
    else:
        session_id = session.get('id')
        if request.method == 'POST':
            logger.info('Form not validated')
            is_seq_filled = bool(form.seq.data)
            is_choicegenome_filled = bool(form.choicegenome.data and form.choicegenome.data != 'default') 
            is_file_filled = bool(form.file.data)
            is_choice_filled = bool(form.choice.data)
            filled_fields = sum([is_seq_filled, is_choicegenome_filled, is_file_filled])
            filled_fields_and_pattern = sum([filled_fields, is_choice_filled])
            if "Sequence too long!" in form.seq.errors:
                flash("Please enter a sequence between 30 and 10,000 characters, longer sequences can be submitted via Upload File", "warning")
                return render_template('index.html', utc_dt=datetime.datetime.utcnow().replace(microsecond=0).isoformat(' '), form=form)
            if filled_fields_and_pattern == 0:
                flash("Please fill out the text field, select a genome, or upload a file and don't forget to choose a pattern to analyze", "warning")
                return render_template('index.html', utc_dt=datetime.datetime.utcnow().replace(microsecond=0).isoformat(' '), form=form)
            if "Sequence contains non-IUPAC symbols" in form.seq.errors:
                flash("Please enter nucleotide sequence, letters A, C, G, T. Fasta sequences can be submitted via Upload File", "warning")
                return render_template('index.html', utc_dt=datetime.datetime.utcnow().replace(microsecond=0).isoformat(' '), form=form)
            if filled_fields == 0:
                flash("Please fill out the text field, select a genome, or upload a file.")
                return render_template('index.html', utc_dt=datetime.datetime.utcnow().replace(microsecond=0).isoformat(' '), form=form)
            if is_choice_filled !=True:
                flash("Please choose a pattern to analyze", "warning")
                return render_template('index.html', utc_dt=datetime.datetime.utcnow().replace(microsecond=0).isoformat(' '), form=form)
        print(form.errors)
        return render_template('index.html', utc_dt=datetime.datetime.utcnow().replace(microsecond=0).isoformat(' '), form=form)



if __name__ == "__main__":
    init_entry_db()
    init_light_db()
    init_db()
    init_heavy_db()

    enable_wal_mode(DB_PATH)
    enable_wal_mode(HEAVY_DB_PATH)
    enable_wal_mode(LIGHT_DB_PATH)


    app.run(debug=False, threaded=True)
