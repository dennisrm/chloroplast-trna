#
"""
For creating, adding, accessing sqlite database of genome minimizer maps
"""


import sqlite3

class ConnectDB:
    def __init__(self, dbfilename):
        self.mapdb = MinimapDB(dbfilename)

    def __enter__(self):
        return self.mapdb

    def __exit__(self, exception_type, exception_value, traceback):
        from traceback import print_tb
        print_tb(traceback)
        print("exit", exception_type, exception_value, traceback)
        self.mapdb.close()
        return True


class MinimapDB:

    def __init__(self, dbfilename):
        self.conn = sqlite3.connect(dbfilename)
        self.c = self.conn.cursor()
        self.tables = []

    def create(self,tabname,colname = "minimer"):
        self.c.execute(f"CREATE TABLE {tabname} ({colname} text)")


    def addmap2(self,tabname,colname,minimap,addrows=True):
        self.c.execute(f"ALTER TABLE {tabname} ADD {colname} text NULL")
        if addrows:
            self.c.executemany(f"INSERT OR IGNORE INTO {tabname}(minimer) VALUES (?)", [(m,) for m in minimap.keys()])
        self.c.execute(f"CREATE INDEX minimers ON {tabname}(minimer)")
        for minimer,positions in minimap.items():
            positions = " ".join([str(pos) for pos in positions])
            self.c.execute(f"UPDATE {tabname} SET {colname} = '{positions}' WHERE minimer = '{minimer}'")
        print("UPDATED")
        self.c.execute("DROP INDEX minimers")
        self.conn.commit()



    def select(self, columns="minimer, positions", table="test", fetchall=False, printselect=False):
        self.c.execute("SELECT " + columns + " FROM " + table)
        if fetchall:
            fetch = self.c.fetchall()
        else:
            fetch = self.c.fetchone()
        if printselect:
            print("select:", fetch)
        return fetch

    def selectJoin(self, table1, table2):
        c = self.conn.cursor()
        c.execute("SELECT {0}.minimer, {0}.positions,{1}.positions FROM {0},{1} where {0}.minimer = {1}.minimer".format(
            table2, table2))
        print(c.fetchall())

    def close(self):
        self.conn.commit()
        self.conn.close()



