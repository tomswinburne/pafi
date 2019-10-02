""" 
Sample plotting of PAFI data
(c) TD Swinburne 2019
swinburne@cinam.univ-mrs.fr
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumtrapz
from scipy import interpolate

#from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
#from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
try:
  from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition)
  has_inset=True
except ImportError:
  has_inset=False




""" 
KEY PARAMETERS- 
looks for file of form dumps/free_energy_profile_TTTK_EEE, 
where TTT=temperature in K and EEE=suffix
"""

epoch = 6 # EEE
Tlist = [100,200,300,400,500] # TTT
file_prefix = "dumps/free_energy_profile" 
plot_title = "SIA dumbell in EAM bcc Fe"
output_file = "PAFI_SIA_Fe.pdf"

fig = plt.figure(figsize=(8,6))
ax = plt.subplot('111')

xf = np.linspace(0.,1.0,101)

if has_inset:
  ax2 = plt.axes([0,0,1,1])
  ip = InsetPosition(ax, [0.3,0.1,0.4,0.4])
  ax2.set_axes_locator(ip)

inset_data =[]

for i,T in enumerate(Tlist):
  fn = file_prefix+"_%dK_%d" % (T,epoch)
  data = np.loadtxt(fn)
  cF = cumtrapz(-data[:,2],x=data[:,0])
  F = data[:,1]
  eb = cumtrapz(data[:,3],x=data[:,0])
  ax.fill_between(data[:-1,0],y1=F[:-1]-eb,y2=F[:-1]+eb,color='C'+str(i),alpha=0.5)
  ax.plot(data[:-1,0],data[:-1,1],'C'+str(i)+'-',label="%dK" % T,lw=2)
  inset_data.append([T,F.max(),eb[F[:-1].argmax()]])


inset_data = np.r_[inset_data]

if has_inset:
  m = 0
  for i in [1,2]:
    m += 0.5*(inset_data[:,1][i]-inset_data[:,1][0])/inset_data[:,0][i]
  ax2.plot(inset_data[:,0],inset_data[:,1][0]+m*inset_data[:,0],'k--',lw=2,label=r"$\Delta$F(0)+cT")
  ax2.plot(inset_data[:,0],inset_data[:,1],lw=2,label=r"$\Delta$F(T)")
  ax2.fill_between(inset_data[:,0],y1=inset_data[:,1]-inset_data[:,2],y2=inset_data[:,1]+inset_data[:,2],color='C0',alpha=0.5)


ax.set_title(plot_title)
ax.legend()
ax.set_ylabel("Free Energy Profile [eV]")
ax.set_xlabel("Reaction Coordinate")
ax2.set_title("Free Energy Barrier [eV]")
ax2.set_xlabel("Temperature [K]")
ax2.legend()
#plt.tight_layout()
#plt.show()
plt.savefig(output_file)
