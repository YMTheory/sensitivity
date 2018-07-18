# Interface to cabinet
# Raymond Tsang, Oct 15, 2015
from __future__ import print_function
import sys
import requests
from elasticsearch import Elasticsearch 
from settings import *

#es = Elasticsearch(['http://xenon:pi=3.14159@particle1.ph.ua.edu:3306'])
#dburl = 'http://xenon:136@particle1.ph.ua.edu:5984/material_database/'

#es = Elasticsearch(['http://xenon:pi=3.14159@nexo.ph.ua.edu/elastic/'])
#dburl = 'http://xenon:136@nexo.ph.ua.edu/material_database/'

es = Elasticsearch(['http://xenon:pi=3.14159@127.0.0.1:9200'])
dburl = 'http://xenon:136@127.0.0.1:5984/material_database/'

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

# id = 'X-AAA'
def get_doc(key):
  results = es.search(q='key:'+key)
  #print 'HITS:',results['hits']['hits']
  for r in results['hits']['hits']:
    doc = r['_source']
    #if key=='R-002': print 'KEY is',key, doc['key']
    if 'doctype' in doc.keys() and doc['key'] == key:
      #print "Found:",key,doc['key']
      #if key == 'MC-012': print doc['rootfiles']
      return doc

# k = 'R-AAA.B.C.D'
def get_measurement(k,iso):
  eprint('get_measurement',k,iso)
  if k[0] == 'P': #retro
    doc = get_doc(k)
    #print 'retro'
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

    #print 'DEBUG', did, sid, cid, mid
    meas = doc['samples'][int(sid)-1]['counting'][int(cid)-1]['analysis_result'][int(mid)-1]
    usual_isotopes = ['u238','th232','k40','co60','cs137']
    #print 'iso',iso
    if iso in usual_isotopes:
      return meas['results'][iso]
    else:
      #print 'meas:',meas['results'].keys()
      for oth in meas['results']['otherisotopes']:
        if iso == oth['isotope'].replace('-','').lower():
          return oth
      return {}
  return {}

# Get montecarlo doc with key="key".
def get_mc(key):
  isokey = key.split('.')[-1]
  dockey = '.'.join(key.split('.')[:-1])
  #print dockey,isokey
  doc = get_doc(dockey)
  did = doc['_id']
  #rev = doc['_rev'].split('-')[0]
  filename = doc['rootfiles'][int(isokey)-1]['filename'] 
  
  #return download_file(dburl+did+'/'+filename,localdir+'/rootfiles/'+did+'_'+rev+'_'+filename)
  return download_file(dburl+did+'/'+filename,localdir+'/rootfiles/'+did+'_'+filename)

def get_mc_meta(dockey,iso):
  #isokey = key.split('.')[-1]
  #dockey = '.'.join(key.split('.')[:-1])
  doc = get_doc(dockey)#['rootfiles'][int(isokey)-1]
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

#doc = get_doc('R-001')
#print doc['title']
def main():
  #print 'Import this with: "import cabinet".'
  #meas = get_measurement('R-045.1.2.1','k40')
  #print meas
  print(get_mc('MC-012.2'))
 
if __name__ == '__main__':
  main()
