# worker.py
import time
import sqlite3
import json
import time, os, traceback
import re
import shutil
from datetime import datetime, timedelta
from config import DB_PATH, enable_wal_mode, with_retry, get_connection
from app import complete_prediction_fasta
from config import setup_logging
import logging
from config import SECRET_KEY


setup_logging()

logger = logging.getLogger("Worker")

enable_wal_mode(DB_PATH)  # WAL mode for this DB

@with_retry()      
def cleanup_old_jobs_db(days: int = 1): #JOB_EXPIRATION_DAYS = 1
    """Remove jobs older than `days` from the jobs database."""
    cutoff = datetime.utcnow() - timedelta(days=days)
    with sqlite3.connect(DB_PATH) as conn:
        c = conn.cursor()
        c.execute("DELETE FROM jobs WHERE created_at < ?", (cutoff,))
        conn.commit()
    logger.info(f"[Cleaner] Removed jobs older than {days} days")


def recover_stale_jobs(timeout_minutes=30):
    cutoff = datetime.utcnow() - timedelta(minutes=timeout_minutes)
    with get_connection(DB_PATH) as conn:
        c = conn.cursor()
        c.execute("""
            UPDATE jobs
            SET status='queued', started_at=NULL
            WHERE status='running' AND started_at < ?
        """, (cutoff,))
        conn.commit()

@with_retry()
def mark_running(job_id):
    with get_connection(DB_PATH) as conn:
        c = conn.cursor()
        c.execute(
            "UPDATE jobs SET status='running', started_at=? WHERE id=?",
            (datetime.utcnow(), job_id),
        )
        conn.commit()


@with_retry()
def mark_completed(job_id, result_data):
    with get_connection(DB_PATH) as conn:
        c = conn.cursor()
        c.execute(
            """
            UPDATE jobs
            SET status='completed',
                result=?,
                completed_at=?
            WHERE id=?
            """,
            (json.dumps(result_data), datetime.utcnow(), job_id),
        )
        conn.commit()


@with_retry()
def mark_failed(job_id, error_message):
    with get_connection(DB_PATH) as conn:
        c = conn.cursor()
        c.execute(
            """
            UPDATE jobs
            SET status='failed',
                error=?,
                completed_at=?
            WHERE id=?
            """,
            (error_message, datetime.utcnow(), job_id),
        )
        conn.commit()


def get_next_queued_job():
    """Atomically claim the next queued job (SQLite-safe)."""
    try:
        with sqlite3.connect(DB_PATH, isolation_level=None) as conn:
            conn.execute("BEGIN IMMEDIATE")
            c = conn.cursor()

            c.execute("""
                SELECT id, params
                FROM jobs
                WHERE status='queued'
                ORDER BY created_at
                LIMIT 1
            """)
            row = c.fetchone()

            if not row:
                conn.commit()
                return None

            job_id, params_json = row

            c.execute("""
                UPDATE jobs
                SET status='running',
                    started_at=?
                WHERE id=?
            """, (datetime.utcnow(), job_id))

            conn.commit()

        return {
            "id": job_id,
            "params": json.loads(params_json)
        }

    except Exception:
        logger.exception("Failed to fetch next job")
        return None

def process_job(job, worker_name):
    job_id = job["id"]
    params = job["params"]
    patternG_obj = re.compile(params["patternG_obj"])
    patternC_obj = re.compile(params["patternC_obj"])
    pattern_printG = params["pattern_printG"]
    pattern_printC = params["pattern_printC"]
    logger.info(f"Started complete_prediction_fasta {job_id}")

    try:
        prediction = complete_prediction_fasta(
            params["filename"],
            patternG_obj,
            patternC_obj, 
            params["clusterscore"],
            params["input_looplength"],
            params["print_indexes"],
            params["print_sequence"],
            params["print_clusters"],
            params["csv_file_repeats"],
            params["csv_file_clusters"],
            params["csv_file_repeatsG"],
            params["csv_file_repeatsC"],
            params["csv_file_clustersG"],
            params["csv_file_clustersC"],
            params["input_ssr"],
            params["input_repeats"],
            params["tandem_select"],
            params["selected_option"],
            params["fasta_id"],
            params["session_id"],
            params["session_dir"]
        )
        logger.info(f"Finished complete_prediction_fasta {job_id}")
        print_indexesG = prediction[0]
        print_indexesC = prediction[1]
        print_sequenceG = prediction[2]
        print_sequenceC = prediction[3]
        clusterG = prediction[4]
        clusterC = prediction[5]
        titlesG = prediction[6]
        titlesC = prediction[7]
        print_clustersG = prediction[8]
        print_clustersC = prediction[9]
        download_files = prediction[10]
        plot_image = prediction[11]
        image1 = prediction[12]
        image2 = prediction[13]
        cluster_fig = prediction[14]
        cluster_data_threshold = prediction[17]
        if params["input_repeats"] == 'any':
            input_repeats_num = 2
        else:
            input_repeats_num = params["input_repeats"]
        conditions_summary = f'Entered parameters: <br> Pattern: {pattern_printG}/{pattern_printC} <br> Cluster score: {params["clusterscore"]} <br> Loop size: {params["input_looplength"]} <br> SSR: {params["input_ssr"]} <br> Number of repeats per cluster: {input_repeats_num}'
        if params["selected_option"] == 'Custom pattern':
            orientationF = 'F'
            orientationR = 'R'
        else: 
            orientationF = 'G'
            orientationR = 'C'
        
        download_figures = 1
        session_dir = params["session_dir"]
        logger.info(f"Generated results after complete_prediction_fasta {job_id}")


        result_data = {
            "prediction_text0": f'Indexes of repeats {orientationF}-orientation (0-based)',
            "prediction_text1": print_indexesG,
            "prediction_text2": f'Indexes of repeats {orientationR}-orientation (0-based)', 
            "prediction_text3": print_indexesC, 
            "prediction_text4": f'Found repeats {orientationF}-orientation', 
            "prediction_text5": print_sequenceG, 
            "prediction_text6": f'Found repeats {orientationR}-orientation', 
            "prediction_text7": print_sequenceC, 
            "tablesG": clusterG.to_html(classes=["table table-bordered table-striped table-sm"]), 
            "titlesG": titlesG, 
            "tablesC": clusterC.to_html(classes=["table table-bordered table-striped table-sm"]), 
            "titlesC": titlesC,
            "prediction_text8": titlesG, 
            "prediction_text9": print_clustersG, 
            "prediction_text10": titlesC, 
            "prediction_text11": print_clustersC, 
            "download_files": download_files, 
            "plot_image": plot_image, 
            "plot1_url": image1,
            "plot2_url": image2,
            "plot3_url": cluster_fig,
            "prediction_text12": conditions_summary,
            "session_dir": session_dir, 
            "cluster_data_threshold": cluster_data_threshold, 
            "download_figures": download_figures
        }

        mark_completed(job_id, result_data)

        logger.info(f"Completed job {job_id}")

    except Exception as e:
        logger.exception(f"Job {job_id} failed")
        mark_failed(job_id, str(e))

def worker_loop(worker_name="Worker"):
    logger.info("Started")

    last_cleanup = datetime.utcnow()
    empty_queue_count = 0

    while True:
        try:
            now = datetime.utcnow()

            if now - last_cleanup > timedelta(days=1):
                logger.info("Worker cleaning old jobs")
                cleanup_old_jobs_db(days=1)
                last_cleanup = now

            job = get_next_queued_job()

            if job:
                empty_queue_count = 0
                process_job(job, worker_name)
            else:
                empty_queue_count += 1
                sleep_time = min(empty_queue_count * 2, 30)
                logger.debug(f"Empty queue, I sleep sweetly {sleep_time} secs")
                time.sleep(sleep_time)

        except KeyboardInterrupt:
            logger.info("Stopped by user")
            break
        except Exception:
            logger.exception("Crashed")
            time.sleep(5)


if __name__ == "__main__":
    print("Worker started...")
    recover_stale_jobs()
    worker_loop()
