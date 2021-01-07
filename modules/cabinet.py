# Interface to cabinet
# Raymond Tsang, Oct 15, 2015
#
# Modifed to run as standalone script
# Raymond Tsang, Mar 6, 2017

import requests
from elasticsearch import Elasticsearch 

# Rootfiles will be saved to this directory
localdir = '.'

# ================================
es = Elasticsearch(['http://xenon:pi=3.14159@nexo.ph.ua.edu:80/elastic'])
dburl = 'http://xenon:136@nexo.ph.ua.edu/material_database/'

# key = 'X-AAA'
def get_doc(key):
  results = es.search(q='key:'+key)
  for r in results['hits']['hits']:
    doc = r['_source']
    if 'doctype' in doc.keys() and doc['key'] == key:
      return doc

# k = 'R-AAA.B.C.D'
def get_measurement(k,iso):
  if k[0] == 'P': #retro
    doc = get_doc(k)
    ms = doc['measurements']
    for m in ms:
      if iso == m['isotope'].replace('-','').lower():
        return { 
          'specific_activity': m['value'], 
          'error': m['error'], 
          'error_type': 'Symmetric error (68% C.L.)' if m['type'] == 'obs' else 'Upper limit (90% C.L.)', 
          'altlimit': '', 
          'unit': m['unit'],
          'asym_lower_error': '',
          'systerr': ''
        }
  else: # regular
    did, sid, cid, mid = k.split('.')
    doc = get_doc(did)

    meas = doc['samples'][int(sid)-1]['counting'][int(cid)-1]['analysis_result'][int(mid)-1]
    usual_isotopes = ['u238','th232','k40','co60','cs137']
    if iso.replace('-','').lower() in usual_isotopes:
      return meas['results'][iso.replace('-','').lower()]
    else:
      for oth in meas['results']['otherisotopes']:
        if iso.replace('-','').lower() == oth['isotope'].replace('-','').lower():
          return oth
      return {}
  return {}

# Get montecarlo doc with key="key".
# key = MC-XXX.Y
def get_mc(key,localfilename):
  isokey = key.split('.')[-1]
  dockey = '.'.join(key.split('.')[:-1])
  doc = get_doc(dockey)
  did = doc['_id']
  filename = doc['rootfiles'][int(isokey)-1]['filename'] 
  
  return download_file(dburl+did+'/'+filename, localdir+'/'+localfilename)

def get_mc_meta(dockey,iso):
  doc = get_doc(dockey)
  for f in doc['rootfiles']:
    if f['isotope'] == iso:
      return f
  return {}

# url: url of file to be downloaded
# path: path of file on local disk 
# return: path 
def download_file(url,path):
  r = requests.get(url, stream=True)
  with open(path, 'wb') as f:
    for chunk in r.iter_content(chunk_size=1024): 
      if chunk: # filter out keep-alive new chunks
          f.write(chunk)
  return path
 
# ================================
if __name__ == '__main__':

  # The default detector document 
  detector_doc = get_doc('D-005')

  # Loop over all parts
  for comp in detector_doc['components']:
    print comp['name'], comp['material'],comp['montecarloid'], comp['radioassayid']
   
    # Get MC document (One MC doc has the rootfiles for all isotopes associated with this part)
    mcdoc = get_doc(comp['montecarloid'])

    # Download all rootfiles in this MC doc
    for i in range(len(mcdoc['rootfiles'])):
      key = comp['montecarloid']+'.'+str(i+1)
      iso = mcdoc['rootfiles'][i]['isotope']
      filename = 'nEXO_Histos_%s_%s.root' % (comp['name'].replace(' ',''), iso)
      path = get_mc(key,filename)
      print 'MC file (%s) for %s is saved to %s' % (key, iso, path)
      
