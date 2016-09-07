
import os

dir = '../results'
name = 'fits_hamamatsu_0nu_eff_tpc_rdm_%s_%0.2f_%0.1f_years_0.0_counts'
years = [5.0] #[0.5,1.0,2.5,5.0,10.,12.5,15.0,20.0] #[0.5,1.0,2.5,5.0,10.0]
bkgs = [0.25,0.5,1.0,2.0,4.0] #[1.0,2.0,2.5,3.0,4.0]
group = 'InternalTh232' #'VesselTh232' #'InternalTh232' #'InternalTh232' #'Far' #'LXeXe137' #'VesselU238' #'LXeRn222' #'InternalU238'

for year in years:
    for bkg in bkgs: 
        outname = name%(group,bkg,year)
        #outname = outname.replace('_rdm_','_ul_')
        inname = name%(group,bkg,year)
        cmd = 'hadd -f %s/full/all%s.root %s/done/%s_*.root' % (dir,outname,dir,inname)
        print cmd
        os.system(cmd)
