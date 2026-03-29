# Import required modules
import csv
import sqlite3

# Assume database already exists and is created with DB Browser for SQLite or another browser you prefer
connection = sqlite3.connect('location to you sqlite3 db')

cursor = connection.cursor()

file = open('location to your sample data from report_summary_CryptoSSU_v3.0 csv files')

contents = csv.reader(file)
# Insert data into table
insert_records = "INSERT INTO sample_name (sample, single_end, fastq_1,fastq_2) VALUES (?,?,?,?)"

cursor.executemany(insert_records, contents)

select_all = "SELECT * FROM sample_name"
rows = cursor.execute(select_all).fetchall()

# Output to the console screen
for r in rows:
    print(r)

# Committing the changes
connection.commit()

# closing the database connection
connection.close()
