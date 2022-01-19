import os
import uproot as up
import numpy as np
import pandas as pd


datadir = '/p/vast1/nexo/data/merged-v11-mcid-labels/'

fileslist = os.listdir(datadir)

parts_dict = dict()

mc_id_map = pd.read_csv('mc_id_map.txt',delimiter=' ',header=None)
print(mc_id_map)

for filename in fileslist:

  component_name = filename.split('.')[0].split('_')[-1]
  isotope = filename.split('.')[0].split('_')[-2]  
  if component_name not in mc_id_map[1].values:
     continue

  rootfile = up.open(datadir+filename)

  mask = mc_id_map[1].values == component_name
  #print(mask)
  print('{}\t{}'.format(filename,mc_id_map[0].loc[mask].values[0]))

  if 'MC-' not in filename:
     new_filename = 'MC-{:03}_{}'.format(mc_id_map[0].loc[mask].values[0],filename)
     print('Renaming {}'.format(filename))
     os.rename(datadir+filename,datadir+new_filename)
     print('------> renamed to {}'.format(new_filename))

  try:
      num_primaries = rootfile['NPrimaries'].values[0]
  except KeyError as e:
      print(e)
      continue

  if not component_name in parts_dict.keys():
     parts_dict[component_name] = dict()
#  else:
#     if parts_dict[component_name] != num_primaries:
#        print('ERROR: file {} has the wrong number of primaries for the {} component'.format(filename,component_name))
#        print('Currently {} has {}, but the new file has {}'.format(component_name,parts_dict[component_name],num_primaries))
#        if parts_dict[component_name] == 0:
#           parts_dict[component_name] = num_primaries      
  parts_dict[component_name][isotope] = num_primaries



for component, this_dict in parts_dict.items():
  print('{}:'.format(component))
  for isotope, num_primaries in this_dict.items():
    print('\t{}\t{}'.format(isotope,num_primaries))

