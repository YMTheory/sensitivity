
import os

dir = '../results'
name = 'disc_hamamatsu_0nu_eff_tpc_rdm_%0.1f_years_%0.1f_counts'
years = [5.0] #[0.5,1.0,2.5,5.0,10.,12.5,15.0,20.0] #[0.5,1.0,2.5,5.0,10.0]
counts = range(11)

for year in years:
    for count in counts:
        cmd = 'hadd -f %s/full/all%s.root %s/done/%s_*.root' % (dir,name%(year,count),dir,name%(year,count))
        print cmd
        os.system(cmd)
