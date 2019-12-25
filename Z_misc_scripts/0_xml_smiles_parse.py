#!/usr/bin/env python3

import sys,os,re
import pandas as pd
import xml.etree.ElementTree as ET
from CommonUtility import file_handle

# do this if running in jupyter
# pd.set_option('display.max_columns', None)

# convert XML to dataframe (assumes only one layer of nesting)
def xml2df(xml_data):
    root = ET.XML(xml_data) # element tree
    all_records = []
    for i, child in enumerate(root):
        record = {}
        for subchild in child:
            record[subchild.tag] = subchild.text
        all_records.append(record)
    df = pd.DataFrame(all_records)

    # how to make datetimes from unix epoch ints
    df['CreatedTimestamp'] = pd.to_datetime(df['CreatedDate'], unit='s')
    df['ModifiedTimestamp'] = pd.to_datetime(df['ModifiedDate'], unit='s')

    return df

# load XML to dataframe (gotta be small)
xml_data = file_handle(sys.argv[0]).read()
df = xml2df(xml_data)

print(xml_doc)