# You can install pyarrow with conda install -c conda-forge pyarrow
# Thats it!! 
from pyarrow import Table
from pyarrow.parquet import ParquetWriter
import pyarrow as pa
import pandas as pd
from pyarrow import csv
import sys

file_stream = sys.argv[1]
output_dir = sys.argv[2]
chunksize = int(sys.argv[3])


 
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
 
# chunksize = 1048576*1000
 
# Additional Arrow magic.. It will automatically decide on the compression
# found within the document and decompress as it reads.
 

 
# This is our input stream reader instance.
 
input_stream_reader = InputStreamReader(file_stream)

# Since the stream is essentially a list of batches, we can create a for 
# loop that allows us read each batch for i batches in 
# input_stream_reader.batches(). 
 
for i, batch in input_stream_reader.batches():
 
    # We will convert our batch into a pandas dataframe

    df = batch.to_pandas()
 
    smiles = list(df['smiles'])
    print(f'Writing a total of {len(smiles)} to disk.')
     
    # Lets write to disk.
    
    df.to_csv(f'{output_dir}/Enamine_REAL_{i}.smi', sep='\t', index=False, escapechar='n')

 
