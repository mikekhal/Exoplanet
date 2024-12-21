import os
import json
import requests
from urllib.parse import quote as urlencode
from astropy.table import Table
import numpy as np
import pprint
#This investigation was done on solar data with temperatures in between 5000-6000 K, metallicities of -0.5 and 0.5, and solar masses of 0.9999-1.0001 on the data 
#This range resulted in 21344 objects, from which 200 objects were used as part of the TESS mission and had a light curve that could be analysed.
pp = pprint.PrettyPrinter(indent=4)

def mast_query(request):
    """Perform a MAST query."""
    request_url = 'https://mast.stsci.edu/api/v0/invoke'
    headers = {
        "Content-type": "application/x-www-form-urlencoded",
        "Accept": "text/plain",
        "User-agent": "python-requests/3.9"
    }
    req_string = urlencode(json.dumps(request))
    resp = requests.post(request_url, data="request=" + req_string, headers=headers)
    return resp.headers, resp.content.decode('utf-8')

def download_files(science_products, base_dir):
    """Download files for the given science products."""
    download_url = 'https://mast.stsci.edu/api/v0.1/Download/file/'
    for row in science_products:
        out_path = os.path.join(base_dir, row['obs_collection'], row['obs_id'])
        os.makedirs(out_path, exist_ok=True)
        out_path = os.path.join(out_path, os.path.basename(row['productFilename']))
        payload = {"uri": row['dataURI']}
        resp = requests.get(download_url, params=payload)
        with open(out_path, 'wb') as FLE:
            FLE.write(resp.content)
        if not os.path.isfile(out_path):
            print("ERROR: " + out_path + " failed to download.")
        else:
            print("COMPLETE:", out_path)
import pandas as pd
TIC_objs = pd.read_csv('CTL_v8_01_Advanced_Search_1.csv', sep=',', header=0, usecols=[0], skiprows=5)
TIC_objs=(TIC_objs)
objects_of_interest = []
objects_of_interest.extend(['TIC ' + str(tic) for tic in TIC_objs.iloc[:, 0]])



# Base directory for all files
base_output_dir = "mastFiles"

for obj in objects_of_interest:
    print(f"Processing object: {obj}")
    
    # Resolve object coordinates
    resolver_request = {'service': 'Mast.Name.Lookup', 'params': {'input': obj, 'format': 'json'}}
    headers, resolved_object_string = mast_query(resolver_request)
    resolved_object = json.loads(resolved_object_string)
    
    if not resolved_object['resolvedCoordinate']:
        print(f"ERROR: Object {obj} could not be resolved.")
        continue
    
    obj_ra = resolved_object['resolvedCoordinate'][0]['ra']
    obj_dec = resolved_object['resolvedCoordinate'][0]['decl']
    
    # Query observations
    mashup_request = {
        "service": "Mast.Caom.Filtered",
        "format": "json",
        "params": {
            "columns": "*",
            "filters": [{"paramName": "target_name", "values": [obj.split()[-1]]}] #https://www.geeksforgeeks.org/how-to-substring-a-string-in-python/
        }
    }
    
    headers, out_string = mast_query(mashup_request)
    count = json.loads(out_string)
    
    if not count.get('data'):
        print(f"ERROR: No observations found for {obj}.")
        continue
    
    mast_data_table1 = Table()
    for col, atype in [(x['name'], x['type']) for x in count['fields']]:
        if atype == "string":
            atype = "str"
        if atype == "boolean":
            atype = "bool"
        mast_data_table1[col] = np.array([x.get(col, None) for x in count['data']], dtype=atype)
    
    interesting_observation = mast_data_table1[0]
    obsid = interesting_observation['obsid']
    
    # Query products
    product_request = {'service': 'Mast.Caom.Products', 'params': {'obsid': obsid}, 'format': 'json', 'pagesize': 100, 'page': 1}
    headers, obs_products_string = mast_query(product_request)
    obs_products = json.loads(obs_products_string)
    
    sci_prod_arr = [x for x in obs_products['data'] if x.get("productType", None) == 'SCIENCE']
    science_products = Table()
    for col, atype in [(x['name'], x['type']) for x in obs_products['fields']]:
        if atype == "string":
            atype = "str"
        if atype == "boolean":
            atype = "bool"
        if atype == "int":
            atype = "float"  # array may contain nan values, and they do not exist in numpy integer arrays
        science_products[col] = np.array([x.get(col, None) for x in sci_prod_arr], dtype=atype)
    
    print(f"Number of science products for {obj}: {len(science_products)}")
    obj_dir = os.path.join(base_output_dir, obj.replace(" ", "_"))
    download_files(science_products, obj_dir)
