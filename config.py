# config.py
import os
import sqlite3
import time
import functools
from pathlib import Path
import logging
from logging.handlers import RotatingFileHandler
import sys

BASE_DIR = os.path.dirname(os.path.abspath(__file__))  # project folder
DB_PATH = os.path.join(BASE_DIR, "db", "jobs.db")
HEAVY_DB_PATH = os.path.join(BASE_DIR, "db", "heavyjobs.db")
LIGHT_DB_PATH = os.path.join(BASE_DIR, "db", "lightjobs.db")
ENTRY_DB_PATH = os.path.join(BASE_DIR, "db", "entry.db")

LOG_DIR = Path(BASE_DIR) / "logs"
LOG_DIR.mkdir(exist_ok=True)
LOG_FILE = LOG_DIR / "app.log"
LOG_FORMAT = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"

SECRET_KEY_FILE = Path(BASE_DIR) / "secrets" / "secret_key.txt"

if os.path.exists(SECRET_KEY_FILE):
    with open(SECRET_KEY_FILE) as f:
        SECRET_KEY = f.read().strip()
else:
    SECRET_KEY = os.environ.get("SECRET_KEY", "dev-secret-change-me")

if not SECRET_KEY:
    raise RuntimeError("SECRET_KEY not set")

# Logging with rotation
def setup_logging(level=logging.INFO):
    root = logging.getLogger()
    if root.handlers:
        return  # already configured
    root.setLevel(level)

    formatter = logging.Formatter(LOG_FORMAT)
    root.handlers.clear()

    # File only and rotating
    file_handler = RotatingFileHandler(
        LOG_FILE,
        maxBytes=10 * 1024 * 1024,  # 10 MB
        backupCount=5
    )
    file_handler.setFormatter(formatter)
    root.addHandler(file_handler)
    # We may switch it on later to silence 
    # logging.getLogger("werkzeug").setLevel(logging.WARNING)
    # logging.getLogger("gunicorn.access").setLevel(logging.WARNING)
    # logging.getLogger("gunicorn.error").setLevel(logging.INFO)



# --- 1. Enable WAL mode (run once) ---
def enable_wal_mode(DB_PATH):
    with sqlite3.connect(DB_PATH) as conn:
        conn.execute("PRAGMA journal_mode=WAL;")
        conn.execute("PRAGMA synchronous=NORMAL;")  # optional, slightly faster
        conn.commit()
    print("SQLite WAL mode enabled ✅")

# --- 2. Retry decorator for locked writes ---
def with_retry(retries=5, delay=0.1):
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            for attempt in range(retries):
                try:
                    return func(*args, **kwargs)
                except sqlite3.OperationalError as e:
                    if "locked" in str(e).lower():
                        time.sleep(delay)
                    else:
                        raise
            raise RuntimeError(f"Database write failed after {retries} retries")
        return wrapper
    return decorator


def get_connection(db_path=DB_PATH, timeout=5.0):
    """Return a connection with a short timeout for locking situations."""
    conn = sqlite3.connect(db_path, timeout=timeout)
    conn.execute("PRAGMA journal_mode=WAL;")
    return conn

def get_heavy_connection(db_path=HEAVY_DB_PATH, timeout=5.0):
    """Return a connection with a short timeout for locking situations."""
    conn = sqlite3.connect(db_path, timeout=timeout)
    conn.execute("PRAGMA journal_mode=WAL;")
    return conn

def get_light_connection(db_path=LIGHT_DB_PATH, timeout=5.0):
    """Return a connection with a short timeout for locking situations."""
    conn = sqlite3.connect(db_path, timeout=timeout)
    conn.execute("PRAGMA journal_mode=WAL;")
    return conn