#This is the basis for an N-Body Code created by 
#Mykhailo Lobodin and Sagar Ramchandani

import numpy as np
import matplotlib.pyplot as plt

GravitationalConstant = 6.67408e-11 

class body:                             # universal body container
    
    def __init__ (self, mass, position, velocity):  # body creation 
        self.m = mass
        self.x = position
        self.v = velocity
        self.a = np.array([0,0,0])
        
    def distance (self, otherbody):     # distance to another body calculator
        r = self.x - otherbody.x
        r = np.sqrt(np.dot(r, r))
        return r
          
    def gravitate (self, otherbody):    # two body gravitation acceleration calculator
        r= self.x - otherbody.x
        r = r*((np.dot(r,r))**(-1.5))
        self.a += -GravitationalConstant*otherbody.m*r
        otherbody.a += GravitationalConstant*self.m*r
        
    def resetacc(self):                 # acceleration reset procedure
        self.a = [0,0,0]
 
    def velocity2(self):                 # V^2
        return np.dot(self.v,self.v)

def magnitude(vector):
    return np.sqrt(np.dot(vector,vector))

def orbitalcalc(ecc, perihelium, msun): # orbit parameters
    return np.sqrt(GravitationalConstant*msun*(1+ecc)/(perihelium))
 
def initialiser(total_time, steps_per_year, mode, integration_scheme, par1 ='+'): 
    global N,SS,dt,endtime,eccen    #total num of bodies and the solar system
    if mode=='kozai':
        N=3
        vel=orbitalcalc(0.5, 1.496e10, 1.989e30)
        SS=np.array([body(1.989e30,np.array([0.,0.,0.]),np.array([0.,0.,0.])),
                body(1.989e29,np.array([1.496e11,0.,0.]),np.array([0.,2.979e4,0.])),
                body(5.972e24,np.array([-1.496e10,0,0]),np.array([0,-vel*0.5,vel*(np.sqrt(3)/2)]))])
    elif mode=='cluster1':
        N=30
        #Define parsec in metres
        parsec=3.0857e16
        solarmass=1.989e30
        radius=0.25*parsec
        SS=[]
        for i in range(N):
            #Random position in spherical coordinates
            pos=np.random.random_sample((3,))*np.array([radius,np.pi,2*np.pi])
            posx=pos[0]*np.sin(pos[1])*np.cos(pos[2])
            posy=pos[0]*np.sin(pos[1])*np.sin(pos[2])
            posz=pos[0]*np.cos(pos[1])
            pos=np.array([posx,posy,posz])
            #Velocity dispersion in m/s
            dispersion=0.2e3
            velocity=np.array([np.random.normal(0,dispersion),np.random.normal(0,dispersion),np.random.normal(0,dispersion)])
            SS.append(body(solarmass,pos,velocity))
    elif mode=='cluster2':
        N=30
        #Define parsec in metres
        parsec=3.0857e16
        solarmass=5*1.989e30
        radius=0.5*parsec
        SS=[]
        for i in range(N):
            #Random position in spherical coordinates
            pos=np.random.random_sample((3,))*np.array([radius,np.pi,2*np.pi])
            posx=pos[0]*np.sin(pos[1])*np.cos(pos[2])
            posy=pos[0]*np.sin(pos[1])*np.sin(pos[2])
            posz=pos[0]*np.cos(pos[1])
            pos=np.array([posx,posy,posz])
            #Velocity dispersion in m/s
            dispersion=2e3
            velocity=np.array([np.random.normal(0,dispersion),np.random.normal(0,dispersion),np.random.normal(0,dispersion)])
            SS.append(body(solarmass,pos,velocity))    

    else: 
        print('not implemented yet')   
    tscale = 31558150
    dt = tscale / steps_per_year
    endtime = total_time*tscale
    if integration_scheme == 'kdk':
       accelerator()
# initialisation of star system for given 
# mode, integraton method and disretisation
                                        
def accelerator():                      # calculating resulting accelerations
    global N,SS
    for i in range(N):
        SS[i].resetacc()
    for i in range (1,N):
        for j in range(i):
            SS[i].gravitate(SS[j])
             
def euler():                            # euler method
    global N,SS,dt
    accelerator()
    for i in range(N):
        SS[i].x += SS[i].v * dt
        SS[i].v += SS[i].a * dt
            
def kdkleap():                          # kdk method integrator
    global N,SS,dt 
    for i in range(N):
        SS[i].v += SS[i].a * dt * 0.5
        SS[i].x += SS[i].v * dt
    accelerator()
    for i in range(N):
        SS[i].v += SS[i].a * dt * 0.5
        
def dkdleap():                          # dkd method integrator
    global N,SS,dt 
    for i in range(N):
        SS[i].x += SS[i].v * dt * 0.5
    accelerator()
    for i in range(N):
        SS[i].v += SS[i].a * dt
        SS[i].x += SS[i].v * dt * 0.5
   
def workloop(total_time, steps_per_year, mode, integration_scheme, plot=False):
    global N,SS,dt,endtime,current_time
    if integration_scheme == 'Euler':   # assigning choosen integrator 
        integrator = euler
    elif integration_scheme == 'kdk':
        integrator = kdkleap
    elif integration_scheme == 'dkd':
        integrator = dkdleap
    current_time = 0
    
    if mode == 'cluster1' or mode=='cluster2':
       initialiser(total_time, steps_per_year, mode, integration_scheme)
       xPositions=[]
       yPositions=[]
       zPositions=[]
       for i in range(N):
           pos=SS[i].x
           xPositions.append([pos[0]])
           yPositions.append([pos[1]])
           zPositions.append([pos[2]])
       timearr=[]
       halfMRadius=[]
       while current_time < endtime:       # main loop 
           integrator()
           halfMRadius.append(halfMassRadius())
           timearr.append(current_time)
           for i in range(N):
               pos=SS[i].x
               xPositions[i].append(pos[0])
               yPositions[i].append(pos[1])
               zPositions[i].append(pos[2])
           current_time += dt
       if plot==True:
           ranger=1.5*halfMassRadius()
           for i in range(N):
               plt.plot(xPositions[i],yPositions[i])
               plt.xlim(-ranger,ranger)
               plt.ylim(-ranger,ranger)
           plt.show()
           for i in range(N):
               plt.plot(yPositions[i],zPositions[i])
               plt.xlim(-ranger,ranger)
               plt.ylim(-ranger,ranger)
           plt.show()
           for i in range(N):
               plt.plot(xPositions[i],zPositions[i])
               plt.xlim(-ranger,ranger)
               plt.ylim(-ranger,ranger)
           plt.show()
           plt.plot(timearr,halfMRadius)
           plt.show()
    elif mode == 'kozai':
        initialiser(total_time, steps_per_year, mode, integration_scheme)
        xPositions=[]
        yPositions=[]
        zPositions=[]
        for i in range(N):
           pos=SS[i].x
           xPositions.append([pos[0]])
           yPositions.append([pos[1]])
           zPositions.append([pos[2]])
        distanceEP=[]
        distancePS=[]
        velE=[]
        timearr=[]
        ecc=[]
        inclin=[]
        appo=np.array([.0,.0])
        peri=np.array([magnitude(SS[2].x),.0])
        tappo=.0
        tperri=.0
        tecc=[]
        k=1
        while current_time < endtime:       # main loop 
           integrator()
           inclin.append(inclination())
           timearr.append(current_time)
           distanceEP.append(SS[0].distance(SS[2]))
           distancePS.append(SS[0].distance(SS[1]))
           velE.append(magnitude(SS[2].a))
           for i in range(N):
               pos=SS[i].x
               xPositions[i].append(pos[0])
               yPositions[i].append(pos[1])
               zPositions[i].append(pos[2])
           if k==0:
               peri[1]=SS[0].distance(SS[2])
               if peri[0]>peri[1]:
                   peri[0]=peri[1]
               elif peri[0]<peri[1]:
                   tperri=current_time
                   k=1
           elif k==1:
               appo[1]=SS[0].distance(SS[2])
               if appo[1]>appo[0]:
                   appo[0]=appo[1]
               elif appo[1]<appo[0]:
                   k=0
                   tappo=current_time
                   tecc.append((tappo+tperri)/2)
                   a=(peri[0]+appo[0])/2
                   e=1-peri[0]/a
                   ecc.append(e)   
           current_time += dt
        if plot==True:
           for i in range(N):
               plt.plot(xPositions[i],yPositions[i])
           plt.show()
           for i in range(N):
               plt.plot(yPositions[i],zPositions[i])
           plt.show()
           for i in range(N):
               plt.plot(xPositions[i],zPositions[i])
           plt.show()
           plt.plot(distanceEP)
           plt.plot(timearr,distancePS)
           plt.plot(timearr,velE)
           plt.show()
           plt.plot(timearr,inclin)
           plt.show()
           plt.plot(tecc,ecc)
           plt.show()
           fig, ax1 = plt.subplots()
           ax1.plot(tecc,ecc)
           ax2 = ax1.twinx()
           ax2.plot(timearr,inclin)
           fig.tight_layout()
           plt.show()
          
def inclination():
    global N,SS
    d1=np.cross((SS[1].x-SS[0].x),(SS[1].v-SS[0].v))
    d1=d1/magnitude(d1)
    d2=np.cross((SS[2].x-SS[0].x),(SS[2].v-SS[0].v))
    d2=d2/magnitude(d2)
    ang=np.arcsin(magnitude(np.cross(d1,d2)))*180/3.14159265
    return ang
    
def halfMassRadius():
    global N, SS
    COM=np.zeros((3,))
    for i in range(N):
        pos=SS[i].x
        COM+=pos
    COM=COM/N
    distances=[]
    for i in range(N):
        pos=SS[i].x
        relPos=COM-pos
        distances.append(magnitude(relPos))
    sortedDistances=np.sort(distances)
    if N%2==0:
        return (sortedDistances[int(N/2)] + sortedDistances[int(N/2)-1])/2
    else:
        return sortedDistances[int(N/2)+1]

#workloop(1e4,0.01,'cluster1','dkd',True)
workloop(100,1000,'kozai','dkd',True)
#initialiser(1, 1000, 'kozai','dkd')
#print(inclination())
