import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np


def drawCircle(panel, r, midpoint=(0, 0), thetaRange=(0, 2 * np.pi), linewidth=1, linestyle=None, color="black",
               alpha=1):
    x, y = midpoint
    xs, ys = [], []
    for theta in np.arange(thetaRange[0], thetaRange[1], 0.01):
        xs.append(r * np.cos(theta) + x)
        ys.append(r * np.sin(theta) + y)
    panel.plot(xs, ys, linestyle=linestyle, linewidth=linewidth, color=color, alpha=alpha)

def drawRadial(panel, theta, rmin, rmax, color="black", linewidth=1):
    xs = [r * np.cos(theta) for r in (rmin, rmax)]
    ys = [r * np.sin(theta) for r in (rmin, rmax)]
    panel.plot(xs, ys, color=color, linewidth=linewidth)

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
    uppoints = []
    downpoints = []
    for i, y in enumerate(ys):
        if clockOrigin:
            theta = np.pi/2 - y * toTheta
        else:
            theta = y * toTheta
        uppoints.append([r0 + 2 * shift, theta, theta])
        downpoints.append([r0 - 2 * shift, theta,theta])

    pointslist = []
    if up:
        pointslist.append([uppoints, 1])
    if down:
        pointslist.append([downpoints, -1])
    for points, updown in pointslist:

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

        for i in indices:

            point1 = points[i]

            js = [j for j in overlaps[i] if j in plotted]
            move = True
            while move:
                move = False
                for j in js:
                    point2 = points[j]
                    while tooClose(point1, point2):
                        if jitter and point1[0] > 1.5:
                            leftRight = np.random.choice([-1,1])
                            #jitterDistance = np.random.uniform(point1[1]-shift*(point1[0]**3)/2, point1[1]+shift*(point1[0]**3)/2)
                            #jitterDistance = np.random.normal(point1[1], shift)
                            jitterDistance = np.random.uniform(point1[2]-shift*40 / (point1[0]), point1[2]+shift*40 / (point1[0]))
                            #jitterDistance = np.random.normal(point1[2],  shift * 10)



                            point1[1] = jitterDistance

                        point1[0] += updown * shift
                        move = True

            if updown == -1:
                js = [j for j in plotted if points[j][0] < point1[0] + minimum_distance]
                move = True
                while move:
                    move = False
                    for j in js:
                        point2 = points[j]
                        while tooClose(point1, point2):
                            point1[0] += updown * shift
                            move = True

            plotted.add(i)

    points = []
    if up:
        points.extend(uppoints)
    if down:
        points.extend(downpoints)
    # drawCircle(panel,r0)
    xs = [r * np.cos(theta) for r, theta,txx in points]
    ys = [r * np.sin(theta) for r, theta,txx in points]
    # colors = [c for x,y,c in points]
    panel.scatter(xs, ys, s=size, c=colors, linewidth=0)

