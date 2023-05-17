import tabula
import pandas as pd

# Replace 'file_path.pdf' with the path to your PDF file
file_path = '/Users/zivseker/Desktop/Projects/bio-project/Analysis of correlated mutations.pdf'

headers = ['Residue i', 'Residue j', 'df','x^2','P']

# Read the table from the PDF file
tables = tabula.read_pdf(file_path, pages='all', multiple_tables=True, pandas_options={'header': None})
print("after read pdf")

non_empty_tables = [table for table in tables if not table.empty]
data_frame_tables = [pd.DataFrame(table) for table in non_empty_tables]
first_table, other_tables = data_frame_tables[0], data_frame_tables[1:]
first_table = first_table.dropna(axis=1, how='all')
first_table = first_table.drop(0)
first_table.columns = headers
final_tables = []
for table in other_tables:
    try:
        new_table = table.dropna(axis=1, how='all')
        new_table = new_table.drop(0)
        new_table.columns = headers
        final_tables.append(new_table)
    except:

        print(table)
final_tables = [first_table, *final_tables]
print(len(final_tables))
final_df = pd.concat(final_tables, axis=0, ignore_index=True)
final_df.to_csv('my-csvvvvv.csv')
