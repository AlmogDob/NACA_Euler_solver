clc; clear;

db_file_str = 'NACA_str.db';
db_table_str = sqlread(sqlite(db_file_str),'NACA_data');

db_file_int = 'NACA_int.db';
db_table_int = sqlread(sqlite(db_file_int),'NACA_data');


