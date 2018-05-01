from obstacle import * 

class Circle(Obstacle):

    def __init__(self, Rc, dens, R):
        Obstacle.__init__(self,Rc,dens)
        if R <= 0:
            raise Exception('Circle radius has to be positive')
        self.R = R

    def make(self, offset = 0):
        N = int(np.round(2*np.pi*self.R*self.dens))
        phi = np.linspace(0,2*np.pi*(1-1/N),N)
        x = self.Rc[0] + self.R*np.cos(phi)
        y = self.Rc[1] + self.R*np.sin(phi)
        z = np.zeros(N)
        for i in range(N):
            p = Particle(i + offset)
            p.radius = 1.0
            p.type = self.obstacle_type 
            p.r = np.array([x[i],y[i],z[i]])
            p.v = np.array([0,0,0])
            p.n = np.array([0,0,0])
            p.in_tissue = 0
            p.boundary = 0
            self.pos.append(p)
            

    