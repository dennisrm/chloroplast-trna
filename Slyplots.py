
#########
from cpgenome.soltools import *
from cpgenome.plotmethods import *

#########



trnaList = gettRNAdb()

ts = [t for t in trnaList if t.species == "lycopersicum"]

cpts = [t for t in ts if t.chromosome == "chloroplast"]

tDic = {t.anticodon:[] for t in cpts}

chr8ts = [t for t in ts if t.chromosome == "chromosome8"]

for t in chr8ts:
    if t.anticodon in tDic:
        tDic[t.anticodon].append(t)



for a in tDic.values():
    begins = []
    ends = []
    for t in a:
        begins.append(t.begin)
        ends.append(t.end)


cpCan = faReader("./genomes/Solanum_cp_genbank.fasta",True)


chr8Can = faReader("./genomes/Solanum_lycopersicum_chr8.fasta",True)[0]


cp = [cp for cp in cpCan if cp.species == "lycopersicum"][0]
xs,ys = cluster(cp,chr8Can)




plt.figure(figsize=(4,4))
panel = plt.axes([0,0,1,1],frameon=False)

drawCircle(panel,1)
#panel.add_patch(mplpatches.Rectangle([-0.5,-0.5],1,1,facecolor="gray"))

panel.plot((-.75,.75,.75,-.75),(1/20,1/20,-1/20,-1/20),color="black",linewidth=1)

xmax = len(cp.sequence)
ymax = len(chr8Can.sequence)
randis = np.random.randint(0,len(xs),100) 
for i in range(len(xs)):
#for i in randis:
    lin = ys[i]/ymax
    rad = xs[i]/xmax

    theta = 2*np.pi*rad
    x1 = np.cos(theta)
    y1 = np.sin(theta)
    level = abs(0.5-rad)
    c = (level,level,1-level)
    #c = (abs(0.5-rad**2)%1,abs(0.5-rad**2)%1,abs(0.5+rad**2)%1)

    if theta < np.pi and theta >= 0:
        x2 = lin-1/2
        y2 = 1/50
        x3 = lin-1/2
        y3 = -1/50
        
        x2 *= 1.5
        x3 *= 1.5
        
        #m = (x1**2 - x2**2 - (1/20)**2 + y1**2)/(2*(x1-x2))
        #R = np.sqrt(y1**2+(x1-m)**2)
        #theta1 = np.arccos((x1-m)/R)
        #theta2 = np.arccos((x2-m)/R)
        #drawCircle(panel,R,midpoint =(m,0), thetaRange=(min(theta1,theta2),max(theta1,theta2)),color=c,alpha=0.1,linewidth=0.2)
                   
    elif theta >= np.pi and theta < 2*np.pi:
        x2 = lin-1/2
        y2 = -1/50
        x3 = lin-1/2
        y3 = 1/50
        
        x2 *= 1.5
        x3 *= 1.5
        
        #m = (x1**2 - x2**2 - (1/20)**2 + y1**2)/(2*(x1-x2))
        #R = np.sqrt(y1**2+(x1-m)**2)
        #theta1 = np.arccos((x1-m)/R)
        #theta2 = np.arccos((x2-m)/R)
        #drawCircle(panel,R,midpoint =(m,0), thetaRange=(min(-theta1,-theta2),max(-theta1,-theta2)),color=c,alpha=0.1,linewidth=0.2)
        

    panel.plot((x1,x2,x3),(y1,y2,y3),color=c,alpha=0.2,linewidth=0.2)
    
circleSwarm(panel,xs,down=False)
panel.set_xlim(-2,2)
panel.set_ylim(-2,2)
    

plt.savefig("teststraight2.png",dpi=600)


# In[25]:


print(xmax,ymax)
plt.figure(figsize=(5,5))
panel = plt.axes([0,0,1,1],frameon=False)
#ys = list(np.random.normal(0,0.4,1500))+list(np.random.normal(6,1,100))
#ys = []
#while len(ys) < 5000:
#    ys.extend(np.random.normal(np.random.uniform(0,10),np.random.uniform(0,2),np.random.randint(250)))
circleSwarm(panel,xs)
panel.set_xlim(-1.5,1.5)
panel.set_ylim(-1.5,1.5)

plt.savefig("testcpswarm.png",dpi=600)


# In[10]:



#with open("output-xs-ys-sly-chr8-cp.txt","w") as outfile:
#    for i in range(len(xs)):
#        line = str(xs[i]) + "\t" + str(ys[i]) + "\n"
#        outfile.write(line)


# In[4]:


# Read
xs=[]
ys=[]
with open("output-xs-ys-sly-chr8-cp.txt","r") as infile:
    for line in infile:
        x,y = [int(a) for a in line.split()]
        xs.append(x)
        ys.append(y)


# In[27]:


plt.figure(figsize=(10,5))
plt.scatter(ys,xs,s=5,color="red",alpha=0.9)

tlines = [t.begin for t in chr8ts if t.anticodon in tDic]
print(len(tlines))
plt.vlines(tlines,0,max(xs),color="blue",alpha=0.2)

#plt.savefig("sly-cp-chr8-seeds-tlines.svg")


# In[28]:


plt.figure(figsize=(10,5))
plt.axes(xlim=(.1e7,.2e7))
plt.scatter(ys,xs,s=5)
plt.vlines(tlines,0,15e4,color="blue",alpha=0.2)

plt.show()


# In[295]:


plt.figure(figsize=(5,5))
plt.axes(xlim=(1.143e7,1.15e7))
plt.scatter(ys,xs,s=10,color="black")
plt.vlines(tlines,0,1.6e5,color="blue",alpha=0.2)


plt.show()


# In[296]:


plt.figure(figsize=(10,5))
plt.axes(xlim=(1.125e7,1.15e7))
plt.scatter(ys,xs,s=5,color="black")
plt.vlines(tlines,0,1.6e5,color="blue",alpha=0.7)


plt.show()


# In[28]:


plt.figure(figsize=(10,5))
#plt.axes(ylim=(-18,150))
plt.hist(xs,bins=100,color="purple")
plt.hlines(0,0,max(xs))
plt.vlines([t.begin for t in cpts if t.accession == cp.accession],-15,0,color="green",alpha=0.9)

#plt.savefig("xxsly-cp-hist-seeds-tlines.svg")


# In[32]:


print(len(xs))


# In[13]:


plt.figure(figsize=(10,5))
plt.axes(xlim=(3e4,4e4))
plt.hist(xs,bins=2000,color="gray")
plt.vlines([t.begin for t in cpts if t.accession == cp.accession],-2,0,color="green",alpha=1)
plt.hlines(0,0,max(xs))


# In[14]:


plt.figure(figsize=(10,5))
plt.axes(xlim=(60000,74000))
plt.hist(xs,bins=2000,color="gray")
plt.vlines([t.begin for t in cpts if t.accession == cp.accession],-2,0,color="green",alpha=1)
plt.hlines(0,0,max(xs))


# In[202]:


plt.figure(figsize=(10,5))
plt.axes(xlim=(30000,50000))
plt.hist(xs,bins=1000,color="green")
plt.vlines([t.begin for t in cpts if t.accession == "MH559323.1"],0,30,alpha=0.8)


# In[7]:


cpCan = [cp for cp in cpCan if cp.species == "lycopersicum"]


# In[8]:


minimerset = set(cpCan[0].miniMap)
for cp in cpCan:
    minimerset = minimerset.intersection(cp.miniMap)
    
rs = []
thetas = []
for i,cp in enumerate(cpCan):
    for min in minimerset:
        for a in cp.miniMap[min]:
            rs.append(i)
            thetas.append(2*np.pi*a/len(cp.sequence))
    

print(len(thetas))


# In[12]:


ts = []
r = 0
crs = []
cols = []
for cp in cpCan:
    cp.trnas = []
    for t in cpts:
        if t.accession == cp.accession:
            cp.trnas.append(t)

for cp1 in [cp for cp in cpCan if True]:
    length = len(cp1.sequence)
    
    ts.extend([2*np.pi*t.begin/length for t in cp1.trnas])
    crs.extend([r]*len(cp1.trnas))
    
    ts.extend([2*np.pi*t.end/length for t in cp1.trnas])
    crs.extend([r]*len(cp1.trnas))
    
    r += 1
    
    if cp1.genus == "Capsicum":
        col = "green"
    elif cp1.genus == "Solanum":
        col = "red"
    else:
        col = "black"
    cols.extend([col]*len(cp1.trnas))
    cols.extend([col]*len(cp1.trnas))
    
    
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, projection='polar')
for i,t in enumerate(ts):
    plt.polar((t,t),(crs[i]-0.3,crs[i]+.3),c=cols[i])
print("mid")
plt.scatter(thetas,rs,alpha=0.008,c="blue",s=50)


ax.set_rmax(r-.5)
ax.set_rorigin(-2*r)
ax.set_theta_zero_location('N')
ax.set_theta_direction(-1)
ax.set_yticklabels([])
ax.set_xticklabels([0])
ax.set_xticks([0])

