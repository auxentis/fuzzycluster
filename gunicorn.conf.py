def on_starting(server):
    from app import init_entry_db, init_light_db, init_db, init_heavy_db
    from config import DB_PATH, HEAVY_DB_PATH, LIGHT_DB_PATH, enable_wal_mode

    init_entry_db()
    init_light_db()
    init_db()
    init_heavy_db()

    enable_wal_mode(DB_PATH)
    enable_wal_mode(HEAVY_DB_PATH)
    enable_wal_mode(LIGHT_DB_PATH)

    print("Databases initialized (master process)")