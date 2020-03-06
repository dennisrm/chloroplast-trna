
import numpy as np
from minimapdb import *
from soltools import *
from plotmethods import *
import matplotlib.patheffects as path_effects
import matplotlib.patches as mplpatches



with ConnectDB("cpgenome.db") as mapdb:
    c = mapdb.conn.cursor()
    #c.execute("CREATE INDEX testindex ON sly(minimer)")
    columns = {}
    tablist = "can cba sly spe".split()
    for tabname in tablist:
        c.execute(f"PRAGMA table_info({tabname})")
        fetched = c.fetchall()
        columns[tabname] = [x[1] for x in fetched]
        



with ConnectDB("cpgenome.db") as mapdb:
    c = mapdb.conn.cursor()
    #c.execute("CREATE INDEX testindex ON sly(minimer)")
    columns = {}
    tablist = "can cba sly spe".split()
    for tabname in tablist:
        c.execute(f"PRAGMA table_info({tabname})")
        fetched = c.fetchall()
        columns[tabname] = [x[1] for x in fetched]


# In[5]:


def swarmplot(panel,y_values,xcoor,panel_width,panel_height,xmin,xmax,ymin,ymax,plot_width,minimum_distance,shift,size,color):
        
    def tooClose(point1,point2):
        x1,y1,c = point1
        x2,y2,c = point2
        if (x2-x1)**2+(y2-y1)**2 < minimum_distance**2:
            return True
        return False
    
    panel.plot(4,80)
    xscale = (xmax-xmin)/panel_width
    yscale = (ymax-ymin)/panel_height
    points = []
    for i,y in enumerate(y_values): 
        points.append([xcoor/xscale ,y/yscale,color[i]])
    points.sort()
    
    
    overlaps = {}
    for i,p in enumerate(points):
        overlaps[i] = []
        j = i+1
        while j<len(points) and abs(p[1]-points[j][1]) < minimum_distance:
            overlaps[i].append(j)
            j += 1
        j = i-1
        while j >= 0 and abs(p[1]-points[j][1]) < minimum_distance:
            overlaps[i].append(j)
            j -= 1
    
    
    indices = np.arange(len(points))
    np.random.shuffle(indices)

    plotted = set()
    
    for i in indices:
        
        leftRight = shift * (1 - 2*(i%2))
        point1 = points[i]

        js = [j for j in overlaps[i] if j in plotted]
        move = True
        while move:
            move = False
            for j in js:
                point2 = points[j]
                while tooClose(point1,point2):
                    point1[0] += leftRight
                    move = True
            
            
        
        plotted.add(i)
    
    xs = [x*xscale for x,y,c in points]
    ys = [y*yscale for x,y,c in points]
    colors = [c for x,y,c in points]

    #xlow = xcoor - plot_width/2
    #xhigh = xcoor + plot_width/2
    #for i,x in enumerate(xs):
    #    if x < xlow:
    #        xs[i] = xlow
    #    elif x > xhigh:
    #        xs[i] = xhigh

    panel.scatter(xs,ys,s=size,c=colors,linewidth = 0)


# In[6]:


trnaList = gettRNAdb()


# In[7]:


cptrnas = []
chtrnas = []
for trna in trnaList:
    if trna.species in "lycopersicum annuum pennellii baccatum".split():
        if trna.chromosome == "chloroplast":
            cptrnas.append(trna)
        elif "chromosome" in trna.chromosome:
            chtrnas.append(trna)


# In[8]:


def circleSwarm(panel, ys, r0=1, minimum_distance=1 / 100, shift=1 / 360, size=0.8, colors="black", ymax=None,
                up=True, down=True, jitter = False, clockOrigin=False):
    def tooClose(point1, point2):
        r1, t1, txx = point1
        r2, t2, txx = point2
        x1 = r1 * np.cos(t1)
        y1 = r1 * np.sin(t1)
        x2 = r2 * np.cos(t2)
        y2 = r2 * np.sin(t2)
        if (x2 - x1) ** 2 + (y2 - y1) ** 2 < minimum_distance ** 2:
            return True
        return False

    if ymax is None:
        ymax = max(ys)

    toTheta = 2 * np.pi / ymax
    points = []
    for i, y in enumerate(ys):
        if clockOrigin:
            theta = np.pi/2 - y * toTheta
        else:
            theta = y * toTheta
        points.append([r0 + 2 * shift, theta, None])

    points.sort()

    overlaps = {}
    for i, p in enumerate(points):
        overlaps[i] = []
        j = i + 1

        while j < len(points) and abs(p[1] - points[j][1]) / (r0) < minimum_distance:
            overlaps[i].append(j)
            j += 1
        j = i - 1
        while j >= 0 and abs(p[1] - points[j][1]) / (r0) < minimum_distance:
            overlaps[i].append(j)
            j -= 1

    indices = np.arange(len(points))
    np.random.shuffle(indices)

    plotted = set()
    clusters = []
    for i in indices:

        point1 = points[i]

        js = [j for j in overlaps[i] if j in plotted]
        move = True
        while move:
            move = False
            for j in js:
                point2 = points[j]
                while tooClose(point1, point2):
                    if point1[2] is None and point1[0] > 1.35:
                        #clusters.append([i])
                        #center = sum([points[i][1] for x in overlaps[i]])/len(overlaps[i])
                        #point1[2] = center + (point1[1]-center)*1000
                        point1[2] = np.random.normal(point1[1],shift*12)
                        #point1[2] = np.random.uniform(point1[1]-shift*30 / (point1[0]), point1[1]+shift*30 / (point1[0]))
                        point1[1] = point1[2]
                    point1[0] += shift
                    move = True
        plotted.add(i)
        
    """
    def checkConnection(clust1, clust2):
        for a in clust1:
            for b in clust2:
                if b in overlaps[a]:
                    return True
        return False

    while clusters:
        clust1 = clusters.pop()

        for clust2 in clusters:
            if checkConnection(clust1, clust2):
                clust2.extend(clust1)
                break
        else:
            components.append(clust1)
            """
    



    # drawCircle(panel,r0)
    xs = [r * np.cos(theta) for r, theta,txx in points]
    ys = [r * np.sin(theta) for r, theta,txx in points]
    # colors = [c for x,y,c in points]
    panel.scatter(xs, ys, s=size, c=colors, linewidth=0)


# In[9]:


cpmax = {"can":156781,"cba":157145, "sly":155461 ,"spe":155254}

canmax = {1:3.01019445  , 2:163962470 , 3:261510930 , 4:215701946 , 5:217274494 , 6:219521584 , 7:222112641 , 8:153299543  , 9:238794889 , 10:205736368 , 11:220335243 , 12:229934170 }
cbamax = {1:210184557  , 2:169143200  , 3:297848814 , 4:219314553  , 5:227523855  , 6:253157891  , 7:254917563 , 8: 210184557 , 9:196862403 , 10:229738584 , 11:263656806 , 12:236504799 }
slymax = {1:90311507 , 2:49921527 , 3:64845585 , 4:64069666 , 5:65026213  , 6:46045610 , 7:65272350  , 8:63035830 , 9:67665719 , 10:64837920 , 11:53391715 , 12:65488263 }
spemax = {1:109333515 , 2:59803892 , 3:75414019 , 4:77197300 , 5:77991103  , 6:60730942 , 7:79292169 , 8:70546659 , 9:84057508  , 10:82529941, 11:66223686 , 12:83305730 }

nucmax = {"can":canmax,"cba":cbamax,"sly":slymax,"spe":spemax}


# In[10]:


superymax = 307401564
with ConnectDB("cpgenome.db") as mapdb:
    c = mapdb.conn.cursor()
    for tabname in "sly spe can cba".split()[:]:
        cpcols = [col for col in columns[tabname] if f"_cp_" in col]
        for cpcol in cpcols[0:1]:
            for number in range(4,5):
                
                chromoNum = [col for col in columns[tabname] if f"chr{number}" in col][0] 
                fetched = mapdb.select(f"{cpcol},{chromoNum}",tabname,fetchall=True)
                
                trnas = [trna for trna in cptrnas if trna.accession.replace(".","") in cpcol]
                nuctrnas = [trna for trna in chtrnas if trna.chromosome == f"chromosome{number}" and tabname[1:] in trna.species]
                
                aacodlist = [f"{trna.aminoacid}-{trna.anticodon}" for trna in sorted(trnas,key = lambda x: x.begin)]
                aacodset = set(aacodlist)
                ntlabs = [(trna.begin,f"{trna.aminoacid}-{trna.anticodon}") for trna in nuctrnas if f"{trna.aminoacid}-{trna.anticodon}" in aacodset]
                nonclabs = [(trna.begin,f"{trna.aminoacid}-{trna.anticodon}") for trna in nuctrnas if f"{trna.aminoacid}-{trna.anticodon}" not in aacodset]
                xs = []
                ys = []
                for cp,nuc in fetched:
                    if cp and nuc:
                        for c in cp.split():
                            for n in nuc.split():
                                if n and c:
                                    xs.append(int(c))
                                    ys.append(int(n))

                plt.figure(figsize=(8,4))
                panel = plt.axes([0,0,0.5,1],frameon=False)
                panel2 = plt.axes([0.5,0,1/8,1],frameon=False)
                panel3 = plt.axes([0.65,0.1,0.35,.8],frameon=False)
       
                drawCircle(panel,1)

                panel.plot((1/50,1/50,-1/50,-1/50,1/50),(-.75,.75,.75,-.75,-0.75),color="black",linewidth=0.75)
                print("lenx:",len(xs))
                print("//75: ",len(xs)//75)
                sampleIndexes = np.random.choice(list(range(len(xs))),len(xs)//75)
                

                xmax =  max(max(xs),cpmax[tabname])
                ymax = max(max(ys),nucmax[tabname][number])
                
                superymax = max(ymax,superymax)
                print("superymax:",superymax)
            
                drawCircle(panel,1.3,color="gray",alpha=0.5)
                
                lasttheta = None
                trnaLocations = []
                colordict = {}
                for trna in sorted(trnas,key=lambda x:x.begin,reverse=False):
                    trnaLocations.append(trna.begin)
                    theta = np.pi*5/2 - 2*np.pi*trna.begin/xmax
                    theta2 = np.pi*5/2 - 2*np.pi*trna.end/xmax
                    
                    
                    if trna.begin < trna.end:
                        direction = "f"
                        c = (0.2,0.2,0.5)
                        radcolor = c
                        colordict[f"{trna.aminoacid}-{trna.anticodon}"] = radcolor
         
                    else:
                        c = (0.8,0.4,0.2)
                        direction = "r"
                        radcolor = c

                    colordict[f"{trna.aminoacid}-{trna.anticodon}"] = radcolor
                    if lasttheta and lasttheta-theta < 0.04:
                        texttheta = lasttheta - 0.04
                    else:
                        texttheta = theta
                    drawRadial(panel,theta,1.26,1.4,color=radcolor,linewidth=.25)
                    drawRadial(panel,theta2,1.28,1.32,color=(1/3,1/3,1/3),linewidth=.25)
    
                    textx = np.cos(texttheta)*1.8
                    texty = np.sin(texttheta)*1.8
                    
                    panel.plot((1.4*np.cos(theta),1.4*np.cos(texttheta), 1.48*np.cos(texttheta)),(1.4*np.sin(theta),1.4*np.sin(texttheta), 1.48*np.sin(texttheta)),
                               color=radcolor,linewidth=0.25)
                  
                    rotation = theta*180/np.pi
                    label =  f"{trna.begin} : {trna.aminoacid}-{trna.anticodon}"
                    text = panel.text(textx,texty,label,rotation = rotation,va='center',ha='center',color=c,fontsize=4)
                    text.set_path_effects([path_effects.Stroke(linewidth=2, foreground='white',alpha=0.5),
                       path_effects.Normal()])

                    lasttheta = texttheta


                sampleSet = set(sampleIndexes)
                cs = [(ys[i]/ymax,ys[i]/ymax,1-ys[i]/ymax) for i in sampleIndexes]
                cs = [(0.25,0.25,0.25) for i in sampleIndexes]
                print(tabname,number)
                ####
            
                color=[]
                ####
                for i in sampleIndexes:
                    lin = ys[i]/ymax
                    rad = xs[i]/xmax

                    theta = np.pi*(5/2) - 2*np.pi*rad
                    x1 = np.cos(theta)
                    y1 = np.sin(theta)
                                        
                    level = min(abs(0.5-rad),0.5)
                    if tabname == "can":
                        c = (0,0, 2*level)
                    elif tabname == "cba":
                        c = (level,0,2*level)
                    elif tabname == "sly":
                        c = (2*level,0, 0)
                    elif tabname == "spe":
                        c = ( 2*level, .75*level,0)

                    color.append(c)
            
                    if x1 < -1/50:
                        x2 = -1/50
                        x3,x4 = x2,x2
                        y2 = 1.5*(lin-1/2)
                        y3,y4,y5 = y2,y2,y2
                    elif x1 > 1/50:
                        x2 = 1/50
                        x3,x4 = x2,x2
                        y2 = 1.5*(lin-1/2)
                        y3,y4,y5 = y2,y2,y2
                    elif x1 > 0:
                        x2 = x1
                        x3 = 1/50
                        x4 = x3
                        if y1 > 0:
                            y2 = .75
                        else:
                            y2 = -.75
                        y3 = y2
                        y4 = 1.5*(lin-1/2)
                        y5 = y4
                    else:
                        x2 = x1
                        x3 = 1/50
                        x4 = x3
                        if y1 > 0:
                            y2 = .75
                        else:
                            y2 = -.75
                        y3 = y2
                        y4 = 1.5*(lin-1/2)
                        y5 = y4
                    
                    x5 = -x4

                        ###
                    panel.plot((x1,x2,x3,x4,x5),(y1,y2,y3,y4,y5),color=c,alpha=0.2,linewidth=0.2)
                circleSwarm(panel,[xs[i] for i in sampleIndexes],down=False,jitter=True,colors=cs,clockOrigin=True)
                
                
                panel_width = 1
                panel_height = 4
                swxmin,swxmax = -1,1
                swymin,swymax = -(superymax-ymax)/2,ymax+(superymax-ymax)/2
                plot_width = 1
                minimum_distance = .012
                shift = 0.0006
                size = .5
                
                swarmplot(panel2,[ys[i] for i in sampleIndexes],0,panel_width,panel_height,swxmin,swxmax,swymin,swymax,plot_width,minimum_distance,shift,size,color)
                scale = ymax/superymax
                #for trna in nuctrnas:
                #    panel2.plot((.1,.2),(trna.begin*scale,trna.begin*scale),alpha=0.1,color="blue")
                
                gap = (superymax-ymax)/2
           
                for nuc in nuctrnas:
                    
                    #panel3.plot((0,.2,.75,1),( gap*ymax/superymax + ymax*nuc.begin/superymax, gap*ymax/superymax + ymax*nuc.begin/superymax,nuc.begin,nuc.begin),color="black",alpha=0.1)
                    #panel3.plot((0, 0.1, 0.75, 1),( 1.1*gap*ymax/superymax + 1.1*ymax*nuc.begin/superymax,1.1*gap*ymax/superymax + 1.1*ymax*nuc.begin/superymax,nuc.begin,nuc.begin),color="black",alpha=0.1,linewidth=.25)
                    panel3.plot(( 0.75, 1),(nuc.begin,nuc.begin),color="black",alpha=0.3,linewidth=.25)
                    
                print("ntlabs",len(ntlabs))
                print("nonclabs",len(nonclabs))
                

                histx = 1.25
                hist,bins = np.histogram([x[0] for x in ntlabs],bins = np.linspace(0,ymax,200))
                binsize = bins[1]-bins[0]
                
                
                for i in range(len(hist)):
                    maxhist = max(hist)
                    if hist[i] > 0:
                        level = 0.25 + 0.75*hist[i]/maxhist
                    else: 
                        level = 0
                    rectangle = mplpatches.Rectangle([histx,bins[i]],.15,binsize,facecolor=(1-level,1-level/4,1-level),
                                 linewidth=0)
                    panel3.add_patch(rectangle)
                    

                
                panel3.text(.75,ymax+binsize,f"All tRNAs in\nchromosome{number}",rotation=90,fontsize=4,color="black",ha='left',va='bottom')

                panel3.text(histx,ymax+binsize,"Chloroplast anticodon\nn : {}".format(sum(hist)),rotation=90,fontsize=4,color=(0,.7,0),ha='left',va='bottom')
                #panel3.text(histx+.05,ymax*1.15,f"n : {sum(hist)}",rotation=90,fontsize=4,color=(0,.7,0),ha='left',va='center')
                #panel3.text(histx,-ymax*0.03-(maxhist+1)*binsize,f"0 to\n{maxhist}",fontsize=4,color="black",ha='left',va='top')

                    
                rectangle = mplpatches.Rectangle([histx-.05,0],.5, ymax,fill=False,
                                 linewidth=.5,color="black") 
                panel3.add_patch(rectangle)
                
                histx += 0.25
                
                hist,bins = np.histogram([x[0] for x in nonclabs],bins = np.linspace(0,ymax,200))
                binsize = bins[1]-bins[0]
                panel3.text(histx,ymax+binsize,"Non-chloroplast\nn : {}".format(sum(hist)),rotation=90,fontsize=4,color="black",ha='left',va='bottom')
                #panel3.text(histx+.05,ymax*1.15,f"n : {sum(hist)}",rotation=90,fontsize=4,color="black",ha='left',va='center')

                for i in range(len(hist)):
                    maxhistnon = max(hist)
                    if hist[i] > 0:
                        level = 0.25 + 0.75*hist[i]/maxhistnon
                    else:
                        level = 0
                    rectangle = mplpatches.Rectangle([histx,bins[i]],.15,binsize,facecolor=(1-level,1-level,1-level),
                                 linewidth=0)
                    panel3.add_patch(rectangle)
                    """xcursor += .1
                # 
                xcursor = xcsave
                for i in range(1,maxhistnon+1):
                    level = i/maxhistnon
                    rectangle = mplpatches.Rectangle([xcursor,-ymax*0.03-2*binsize],.1,binsize,facecolor=(1-level,1-level,1-level),linewidth=0)
                    panel3.add_patch(rectangle)
                    
                    if i == maxhistnon:
                        panel3.text(xcursor,-ymax*0.03-7*binsize,i,fontsize=4,ha="left",va="top")
                    xcursor += .11"""
                    
                #panel3.text(histx,-ymax*0.03-(maxhistnon+1)*binsize,f"0 to\n{maxhistnon}",fontsize=4,color="black",ha='left',va='top')

                    
                
                    
                rectangle = mplpatches.Rectangle([histx-.3,0],.5, ymax,fill=False,
                                 linewidth=.5,color="black") 
                panel3.add_patch(rectangle)
                
                histx += 0.5
                
                rectleft = histx-.05
                
                xcursor = rectleft+.05
                xcsave = xcursor
                for i in range(1,maxhist+1):
                    level = 0.25 + 0.75*i/maxhist
                    rectangle = mplpatches.Rectangle([xcursor,-ymax*0.03-0*binsize],.1,binsize,facecolor=(1-level,1-level/4,1-level),linewidth=0)
                    panel3.add_patch(rectangle)
                    
                    if i == maxhist:
                        panel3.text(xcursor,-ymax*0.03-7*binsize,i,fontsize=4,ha="left",va="top")
                        
                    xcursor += .11
                    
                xcursor = xcsave
                for i in range(1,maxhistnon+1):
                    level = i/maxhistnon
                    rectangle = mplpatches.Rectangle([xcursor,-ymax*0.03-2*binsize],.1,binsize,facecolor=(1-level,1-level,1-level),linewidth=0)
                    panel3.add_patch(rectangle)
                    
                    if i == maxhistnon:
                        panel3.text(xcursor,-ymax*0.03-7*binsize,i,fontsize=4,ha="left",va="top")
                    xcursor += .11
                                    
                
                aahist = {}
                maxhistaa = 0
                for aacod in aacodlist:
                    hist,bins = np.histogram([x[0] for x in ntlabs if x[1] == aacod],bins = np.linspace(0,ymax,200))
                    aahist[aacod] = hist
                    maxhistaa = max(maxhistaa,max(hist))
                    binsize = bins[1]-bins[0]
                
                for aacod in aacodlist:
                    hist = aahist[aacod]
                    panel3.text(histx,ymax+binsize,aacod,rotation=90,fontsize=4,color=colordict[aacod],ha='left',va='bottom')
                    for i in range(len(hist)):                        
                        if hist[i] > 0:
                            level = 0.25 + 0.75*hist[i]/maxhistaa
                            facecolor = [1-(1-c)*level for c in colordict[aacod]]
                        else:
                            facecolor = (1,1,1)
                        rectangle = mplpatches.Rectangle([histx,bins[i]],.1,binsize,facecolor=facecolor,linewidth=0)
                        panel3.add_patch(rectangle)
                    
                    panel3.text(histx,ymax*1.13,f"n : {sum(hist)}",rotation=90,fontsize=4,color=colordict[aacod],ha='left',va='center')
                    
                    
                    histx += 0.11
                
                histscalex = xcsave
                for i in range(1,maxhistaa+1):
                    level = i/maxhistaa
                    level = .25+.75*level
                    facecolor = [1-(1-a)*level for a in (0.2,0.2,0.5)]
                    rectangle = mplpatches.Rectangle([histscalex,-ymax*0.03-4*binsize],.1,binsize,facecolor=facecolor,linewidth=0)
                    panel3.add_patch(rectangle)
                    
                    facecolor = [1-(1-a)*level for a in (0.8,0.4,0.2)] 
                    rectangle = mplpatches.Rectangle([histscalex,-ymax*0.03-6*binsize],.1,binsize,facecolor=facecolor,linewidth=0)
                    panel3.add_patch(rectangle)
                    
                    if i in {1,maxhistaa}:
                        panel3.text(histscalex,-ymax*0.03-7*binsize,i,fontsize=4,ha="left",va="top")
                    
                    histscalex += .11
                    

                panel3.plot(( 0.75, 1),(-ymax*0.03,-ymax*0.03),color="black",alpha=0.3,linewidth=.25)
                plt.text(0.9,-ymax*0.065,"1\n tRNA",ha='center',va='center',fontsize=4)

                rectangle = mplpatches.Rectangle([.75,-ymax*0.115],histx-.75+.05, ymax*0.105,fill=False,
                                 linewidth=.5,color="black")                    
                panel3.add_patch(rectangle)
                
                
                panel3.text(rectleft,-ymax*.11,f"tRNAs per {binsize//100/10} kb [~{int(10*binsize/xmax)/10}x chloroplast genome]",va="bottom",ha="left",fontsize=4)
                

                    
                    
                    

                rectangle = mplpatches.Rectangle([rectleft,0],histx-rectleft+.05, ymax,fill=False,
                                 linewidth=.5,color="black")
                panel3.add_patch(rectangle)
   
                gap = (superymax-ymax)/2
    
                nuchigh = gap*ymax/superymax + ymax*ymax/superymax 
                nuclow = (gap*ymax/superymax)
                
                panel3.plot((.75,1,1,.75,.75),(0,0,ymax,ymax,0),color="black",linewidth = 0.5)
                
                
                #panel3.plot((0,.2,.2,0,0),(nuclow,nuclow,nuchigh,nuchigh,nuclow),color="black",linewidth = 0.5)
                    
                panel2.text (0,ymax,"\n".join([f"chromosome{number}",f"{len(sampleIndexes)} of ",f"{len(xs)}","alignment seeds",f"{(ymax//10**5)/10} Mb"]),va='bottom',ha='center',fontsize=5)
                    
                    
                left,right = .8,2.12
                bottom,top = 1.1,1.7
                
                panel.text (left+.05,bottom+.5,"Chloroplast tRNA genes",fontsize=4)
                panel.plot((left,right),(bottom+.45,bottom+.45),linewidth=0.25,color="black")
                panel.text (left+.05,bottom+.25,"Fordward strand tRNA\n  \u2014 Poisition : Amino acid-Anticdon",color=(0.2,0.2,0.5),fontsize=4)
                panel.text (left+.05,bottom+.05,"Reverse strand tRNA\n  \u2014 Poisition : Amino acid-Anticdon",color=(0.8,0.4,0.2),fontsize=4)
                

                panel.plot((left,right,right,left,left),(bottom,bottom,top,top,bottom),linewidth=0.5,color="black")
                
    

                
                
                #panel2.plot((.94,.94),(0,ymax),color="black",linewidth=.75)
                panel2.plot((.94,.94),(0,ymax),color="black",linewidth=.75,linestyle=":")
                panel3.plot((.1,.1),(0,ymax),color="black",linewidth=.75,linestyle=":")
                
                for y in (0,ymax/4,ymax/2,ymax*3/4,ymax):
                    panel2.plot((.9,1),(y,y),color="black",linewidth=.5)
                    panel3.plot((0.05,.15),(y,y),color="black",linewidth=.5)


                    

                
                panel2.plot((-.94,-.94),(0,ymax),color="black",linewidth=.75)
                
                #panel3.plot((.65,.65),(0,ymax),color="black",linewidth=.75)
                
                for y in range(0,ymax,2*10**7):
                    #panel2.plot((.95,1),(y,y),color="black",linewidth=.5)
                    panel2.plot((-.95,-1),(y,y),color="black",linewidth=.5)

                    panel3.plot((.7,.75),(y,y),color="black",linewidth=.5)
                    
                    panel2.text(-1,y,f"{y//10**6} Mb",va='center',ha='right',fontsize=4)
                    panel3.text(.69,y,f"{y//10**6} Mb",va='center',ha='right',fontsize=4)
                    
                

                
                panel.set_xlim(-2.5,2.5)
                panel.set_ylim(-2.5,2.5)
                
                #panel2.set_xlim(-2,2)
                panel2.set_ylim(-(superymax-ymax)/2,ymax+(superymax-ymax)/2)
                panel2.set_xlim(swxmin,swxmax)
                
                panel3.set_ylim(-ymax*0.12,ymax*1.12)
                panel3.set_xlim(0,5.4)
                
                for p in [panel,panel2,panel3]:
                    p.tick_params(axis='both',which='both',bottom=False,labelbottom=False,left=False,labelleft=False,right=False,labelright=False,top=False,labeltop=False)
                    p.ticklabel_format(style='plain')
                plt.savefig(f"./cpgenome_images/2pc_chromoplot_{tabname}{number}_{cpcol}.png",dpi=600)
                plt.close()
            


# In[11]:


print(binsize)
print(xmax)
print(binsize/xmax)


# In[12]:



plt.figure(figsize = (8,2))

gsps = ["can" ,"cba","sly","spe"]
cscales =  [(0,0, 2),(1,0,2),(2,0,0),( 2,.75,0)]
dimensions = [[0,0,.24,1],[.25,0,.24,1],[.5,0,0.24,1],[0.75,0,0.24,1]]

for i in range(4):
    gsp = gsps[i]
    cscale = cscales[i]
    dims = dimensions[i]
    keypanel = plt.axes(dims,frameon=False)
    
    keypanel.set_xlim(-2,2)
    keypanel.set_ylim(-2,2)
    keypanel.tick_params(axis='both',which='both',bottom=False,labelbottom=False,left=False,labelleft=False,right=False,labelright=False,top=False,labeltop=False)


    
    

        
    drawCircle(keypanel,r=1)
    
    for i in np.linspace(0,1,1000):
        
        level = abs(0.5-i)
        c = [k*level for k in cscale]
        theta = np.pi/2 - 2*np.pi*i
        x = np.cos(theta)
        y = np.sin(theta)
        keypanel.plot((x*(.75),x),(y*.75,y),color=c,linewidth=2,alpha=0.2)
        
    drawCircle(keypanel,r=1.01,linewidth=1)
    drawCircle(keypanel,r=0.74,linewidth=1)
    
    for x in range(0,151,25):
        x = x*1000
        xp = x/cpmax[gsp]
        ticktheta = np.pi/2 - 2*np.pi*xp
        drawRadial(keypanel,ticktheta,.9,1.15,linewidth=1,color=(.25,.25,.25))
        
    for x in 50*1000, 100*1000:
        xp = x/cpmax[gsp]
        ticktheta = np.pi/2 - 2*np.pi*xp
        textx = np.cos(ticktheta) * 1.4
        texty = np.sin(ticktheta) * 1.4
        keypanel.text(textx,texty,f"{x//1000} KB",ha='center',va='center',color=(.25,.25,.25),fontsize=7)
    
    
    drawRadial(keypanel,np.pi/2,0.9,1.5,linewidth=1)

    keypanel.text(.025,1.25,f"0 KB",ha='left',va='center',fontsize = 7)
    keypanel.text(-.02,1.4,f"{cpmax[gsp]//1000} KB",ha='right',va='center',fontsize=7)

    
    gspdict = {"can":"Capsicum annuum chloroplast","cba":"Capsicum baccatum chloroplast","sly":"Solanum lycopersicum chloroplast","spe":"Solanum pennellii chloroplast"}
    mainstring = "\n".join(gspdict[gsp].split())
    keypanel.text(0,0,mainstring,ha='center',va='center',fontsize=7)
    
plt.savefig("./figure_images/circle_colorscales.png",dpi=600)

