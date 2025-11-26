import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as scistats

#npts = 50
npts = 100

nmems = 40

modmean =0
modstd = 1

obsmean = 0
obsstd = 1

repeats = 1

obspciles = []
all_obs = []

for n in np.arange(npts):

    percentiles = []
    
    for x in np.arange(repeats):
    
        #Generate 'model' distribution at this point:
        mods = np.random.normal(modmean,modstd,nmems)

        obs = np.random.normal(obsmean,obsstd)
        all_obs.append(obs)
    
    
        thispcile = scistats.percentileofscore(mods,obs)

        #plt.hist(mods,bins=15)
        #ax = plt.gca()
        #ax.axvline(obs)
        #plt.title('pcile='+str(thispcile))
        #plt.show()

        percentiles.append(thispcile)
        
    mean = np.mean(percentiles)
    
    obspciles.append(mean)


obspciles = np.array(obspciles)
all_obs = np.array(all_obs)

print(obspciles)

#plt.hist(obspciles,bins=10)
##plt.hist(all_obs,bins=10)
#plt.show()

testpts = np.linspace(0,100,101)

#At each testpt, find out how many obs were at that percentile or above:
fracpts = []
for t in testpts:
    n_gt_t = len(obspciles[obspciles>=t])

    fracpts.append(n_gt_t/npts)

fracpts = np.array(fracpts)

plt.plot(testpts,fracpts*100.0)
#plot perfect line:
xpts = np.linspace(0,100,10)
ypts = 100 - xpts
plt.plot(xpts,ypts,color='black')
plt.show()
