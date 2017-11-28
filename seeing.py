#!/usr/bin/python
import pandas as pd
import matplotlib.pyplot as plt
import sys
import matplotlib

"""
The script plots the seeing conditions from ESO night seeing log.
Author: E. Borisova

"""
C12 = 0.74
C23 = 0.82
C34 = 0.92

y1=0.5
y2=1.6

if len(sys.argv) < 2:
    print '*****'
    sys.exit('Usage: %s logfile_name' % sys.argv[0])
    print '*****'

print '-------'
print 'Starting the seeing script'
print '-------'
head_names = ['Date','Time','Dimm','FWHM_Autoguiding_Tel','FWHM_Autoguiding_Zenith','FWHM_Im_Analysis_Normal','FWHM_Im_Analysis_LinearObs','FWHM_Im_Analysis_LinearFit','Airmass','Humidity','Wind_Direction','Top_Wind_Speed','ASM_Wind_Speed']

plt.ioff()


#seeing_list = pd.read_csv(argv[1],header=(0), sep="\t",parse_dates=[[0,1]])
seeing_list = pd.read_csv(sys.argv[1],header=None, sep="\t",parse_dates=[[0,1]], names=head_names)
seeing_list.plot('Date_Time' , 'FWHM_Autoguiding_Tel', grid=1, figsize=(12,5),  lw=1.5)

plt.ylabel('Seeing', fontsize=16)
plt.title('GTO-09   %s' % seeing_list.Date_Time[0].date())
plt.axhline(C12, ls='--', c='green', lw = 2)
plt.axhline(C23, ls='--', c='orange', lw = 2)
plt.axhline(C34, ls='--', c='darkred', lw = 2)
plt.ylim(y1,y2)
#plt.tight_layout()
plt.text(seeing_list.Date_Time[1], 0.60, 'C1', weight='bold', fontsize=14)
plt.text(seeing_list.Date_Time[1], 0.75, 'C2', weight='bold', fontsize=14)
plt.text(seeing_list.Date_Time[1], 0.85, 'C3', weight='bold', fontsize=14)
plt.text(seeing_list.Date_Time[1], 1.25, 'C4', weight='bold', fontsize=14)

plotname='seeing_%s.png' % seeing_list.Date_Time[0].date()
plt.savefig(plotname, dpi=200)
print '----->  Figure was saved: %s' % plotname
