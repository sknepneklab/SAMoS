from datetime import *
import numpy as np 
import argparse

from particle import *

def projection_matrix(axis):
  n = axis/np.sqrt(np.dot(axis,axis))
  return (np.identity(3) - np.outer(n,n))


def rotation_matrix(axis,theta):
  n = axis/np.sqrt(np.dot(axis,axis))
  a = np.cos(theta/2)
  b,c,d = -n*np.sin(theta/2)
  return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                  [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                  [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])

class Torus:
  
  def __init__(self, R, r, m, n, l):
    self.R = R
    self.r = r
    self.m = m
    self.n = n 
    self.l = l    
    self.N = int(2*self.l*(self.m**2 + self.m*self.n + self.n**2)/np.gcd(self.m + 2*self.n, 2*self.m + self.n))
    self.particles = [Particle(i) for i in range(self.N)]
    self.__generate_pos()
    self.__generate_vel(1.0)
    self.__generate_director()

  def __generate_pos(self):
    e1 = np.array([1, 0])
    e2 = np.array([0, 1])
    a1 = np.array([1, 0])
    a2 = np.array([0.5, 0.5*np.sqrt(3)])
    c = self.m * a1 + self.n * a2 
    t = -self.l * ((self.m + 2*self.n) * a1 - (2*self.m + self.n) * a2)/np.gcd(self.m + 2*self.n, 2*self.m + self.n)
    A = np.linalg.inv(np.array([[np.dot(c,e1), np.dot(t,e1)],[np.dot(c,e2), np.dot(t,e2)]]))
    i = 0
    for p in range(-200,201):
      for q in range(-200,201):
        r = A @ (p * a1 + q * a2)
        if (0 <= r[0] < 1) and (0 <= r[1] < 1) and (i < self.N):
          theta, phi = 2*np.pi * r
          x = (self.R + self.r*np.cos(theta))*np.cos(phi)
          y = (self.R + self.r*np.cos(theta))*np.sin(phi)
          z = self.r*np.sin(theta)
          if not self.__is_close([x,y,z], self.particles[:i]):
            self.particles[i].r = [x,y,z]
            i += 1
    
  def __generate_vel(self,vav=1.0):
    for p in self.particles:
      x, y, z = p.r
      fact1 = np.sqrt(x*x + y*y)
      fact2 = (self.R - fact1)/fact1
      axis = np.array([-2*x*fact2,-2*y*fact2,2*z])
      v = np.array([np.random.uniform(-0.5,0.5), np.random.uniform(-0.5,0.5), np.random.uniform(-0.5,0.5)])
      v = np.dot(projection_matrix(axis),v)
      vlen = np.sqrt(sum(v**2))
      v *= vav/vlen
      theta = np.random.uniform(0,2*np.pi)
      v = np.dot(rotation_matrix(axis,theta),v)
      p.v = v
 
  def __generate_director(self):
    for p in self.particles:
      x, y, z = p.r
      fact1 = np.sqrt(x*x + y*y)
      fact2 = (self.R - fact1)/fact1
      axis = np.array([-2*x*fact2,-2*y*fact2,2*z])
      n = np.array([np.random.uniform(-0.5,0.5), np.random.uniform(-0.5,0.5), np.random.uniform(-0.5,0.5)])
      n = np.dot(projection_matrix(axis),n)
      nlen = np.sqrt(sum(n**2))
      n *= 1.0/nlen
      theta = np.random.uniform(0,2*np.pi)
      n = np.dot(rotation_matrix(axis,theta),n)
      p.n = n

  def __is_close(self, r, particles):
    for p in particles:
      if np.all(np.isclose(np.array(r)-np.array(p.r), 0)):
        print('papaya!')
        return True 
    return False
      
 
  def write(self,outfile):
    gentime = datetime.now()
    with open(outfile,'w') as out:
      print(f'# Total of {self.N} particles', file=out)
      print(f'# Generated on : {gentime}', file=out)
      print(f'# id  type  x   y   z   vx   vy   vz   nx   ny   nz', file=out)
      for p in self.particles:
        x, y, z = p.r
        vx, vy, vz = p.v
        nx, ny, nz = p.n
        print(f'{p.idx} {p.tp} {x} {y} {z} {vx} {vy} {vz}  {nx} {ny} {nz}', file=out)
   

parser = argparse.ArgumentParser()
parser.add_argument("--R", dest='R', type=float, default=10.0, help="large radius of torus")
parser.add_argument("--r", dest='r', type=float, default=4.0, help="small radius of torus")
parser.add_argument("--m", dest='m',  type=int, default=6, help="m parameter")
parser.add_argument("--n", dest='n',  type=int, default=3, help="n parameter")
parser.add_argument("--l", dest='l',  type=int, default=6, help="l parameter")
parser.add_argument("--o", dest='out', type=str, default='torus.dat', help="output file")
args = parser.parse_args()

print()
print("\tSoft Active Matter on Surfaces (SAMoS)")
print("\tBuilding of a random configuration on torus")
print()
print("\tRastko Sknepnek")
print("\tUniversity of Dundee")
print("\t(c) 2023")
print("\t----------------------------------------------")
print()
print("\tLarge radius : ", args.R)
print("\tSmall radius : ", args.r)
print("\tm : ", args.m)
print("\tn : ", args.n)
print("\tl : ", args.l)

N = int(2*args.l*(args.m**2 + args.m*args.n + args.n**2)/np.gcd(args.m + 2*args.n, 2*args.m + args.n))
print("\tExpected number of particles : ", N)
start = datetime.now()

T = Torus(args.R, args.r, args.m, args.n, args.l)
T.write(args.out)

end = datetime.now()

total = end - start

print("\tGenerated number of particles : ", T.N)
print("\tOutput file : ", args.out)
print()
print(f"  *** Completed in {total.total_seconds()} seconds *** ")
print()


