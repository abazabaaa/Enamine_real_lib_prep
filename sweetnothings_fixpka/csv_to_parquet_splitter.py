# You can install pyarrow with conda install -c conda-forge pyarrow
# Thats it!! 
from pyarrow import Table
from pyarrow.parquet import ParquetWriter
import pyarrow as pa
import pandas as pd
from pyarrow import csv
 
 
class InputStreamReader:
    def __init__(self, file_stream):
        self.file_stream = file_stream
        self._stream = None
 
    def batches(self):
        i = tries = 0
        while True:
            try:
                batch = self.__next_batch()
                i += 1
                yield i, batch
            except StopIteration:
                break
# This is the real workhorse. This function allows us to read the next
# recordbatch from the CSVStreamReader instance that we create.
    def __next_batch(self):
        return self.stream.read_next_batch()
 
# These params include some of the things that we will need to know as
# we stream the CSV one record batch at a time.
    @property
    def stream(self):
        if not self._stream:
 
# This is very important.. it will dictate how many lines are read
# within each recordbatch. The default settings currently allow for
# 6.3M smiles strings to be read into memory per chunk. My workstation
# has around 64GB of RAM, but for this exercise I want to keep things
# below 16GB. This took some trial and error, but it works well.
            read_options = pa.csv.ReadOptions(block_size=chunksize)
 
# You will need to set this for document. The enamine REAL is tab- 
# separated.
            parse_options = pa.csv.ParseOptions(delimiter='\t')
 
# This is part of the magic of Arrow. We can choose which columns we want
# to read. Now, for some reason.. the folks at Enamine decided to use
# True/False for some columns. This is regrettable. It breaks Arrow 
# (atleast in my feeble understanding of the process) , so we
# need to exclude all columns but the SMILES and the id.
 
            convert_options = pa.csv.ConvertOptions(include_columns=include_columns)
 
# This is where the CSVStreamReader instance will begin.
            self._stream = pa.csv.open_csv(
                self.file_stream, read_options=read_options,
                parse_options=parse_options,
                convert_options=convert_options
                 
            )
        return self._stream

# Here is where we set the needed columns. If you decompress the file
# from enamine and open it with a text editor, you will see the first 
# few lines. There are many other ways to get the header, and I will leave
# that up to you to find. 
 
include_columns = ['smiles', 'idnumber']
 
# This gives us our chunksize, and the InputStreamReader instance
# will use this to decide how many rows to read.
 
chunksize = 1048576*1000
 
# Additional Arrow magic.. It will automatically decide on the compression
# found within the document and decompress as it reads.
 
file_stream = '/data/Enamine_REAL_HAC_6_20_CXSMILES.cxsmiles.bz2'
 
# This is our input stream reader instance.
 
input_stream_reader = InputStreamReader(file_stream)

# Since the stream is essentially a list of batches, we can create a for 
# loop that allows us read each batch for i batches in 
# input_stream_reader.batches(). 
 
for i, batch in input_stream_reader.batches():
 
    # We will convert our first batch into a pandas dataframe
    # This makes dealing with schema easier.
    df = batch.to_pandas()
 
    # Now we convert our df into an Arrow table. This will allow
    # for us to write as a feather, or parquet.
    table = pa.Table.from_pandas(df)
 
    # Since we will be writing this to parquet, lets  pull the schema
    # and the smiles. We can use the smiles to just tell us
    # how many we are writing. More clever readers will see that 
    # we could easily do many other things with RDKit here. : D
 
    schema = table.schema
    smiles = list(df['smiles'])
    print(f'Writing a total of {len(smiles)} to disk.')
     
    # Lets write to disk. Each parquet will take up 130MB and you 
    # should end with about 28 of them.
    # I needed about 3GB RAM for each chunk. If readers try this and 
    # get kernel crashes, then lower the chunksize.
 
    ParquetWriter(f'/pathtoparquet/Enamine_REAL_{i}.parquet', \
        schema).write_table(table)
 