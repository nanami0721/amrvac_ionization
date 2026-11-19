#!/usr/bin/env python3
"""Independent finite-difference curl and geometric checks for the reference."""
import math,argparse
p=argparse.ArgumentParser();p.add_argument('--radius',type=float,default=.35);p.add_argument('--radial-min',type=float,default=1.);p.add_argument('--radial-max',type=float,default=2.);args=p.parse_args()
R=args.radius;eps=.001;a0=eps*R/((8/math.sqrt(7))*(6/7)**3)
n=[math.cos(math.radians(146.25)),math.sin(math.radians(146.25)),0]
def potential(q):
 s=sum(x*x for x in q)/R**2
 return [a0*max(0,1-s)**4*x for x in n]
def field(q):
 s=sum(x*x for x in q)/R**2
 cross=[q[1]*n[2]-q[2]*n[1],q[2]*n[0]-q[0]*n[2],q[0]*n[1]-q[1]*n[0]]
 return [-8*a0/R**2*max(0,1-s)**3*x for x in cross]
def inside(p,patch):
 if patch:p=[-p[0],p[2],p[1]]
 theta=math.acos(p[2]/math.sqrt(sum(x*x for x in p)));phi=math.atan2(p[1],p[0])
 return .109375*2*math.pi<=theta<=.390625*2*math.pi and abs(phi)<=.390625*2*math.pi
h=1e-6
for q in [[.05,.08,.11],[-.1,.04,-.15],[.4,0,0]]:
 deriv=[]
 for k in range(3):
  plus=q.copy();minus=q.copy();plus[k]+=h;minus[k]-=h
  deriv.append([(x-y)/(2*h) for x,y in zip(potential(plus),potential(minus))])
 curl=[deriv[1][2]-deriv[2][1],deriv[2][0]-deriv[0][2],deriv[0][1]-deriv[1][0]]
 assert max(abs(x-y) for x,y in zip(curl,field(q)))<1e-10*eps
assert abs(math.sqrt(sum(x*x for x in field([0,0,R/math.sqrt(7)])))/eps-1)<1e-14
for direction in [1,-1]:
 v=[-direction*n[1],direction*n[0],0]
 def center(t):return [1.5*n[k]+(t-.3)*v[k] for k in range(3)]
 t=.3+1.5*math.tan(math.radians(-5.625))/direction
 assert inside(center(t-1e-4),0)!=inside(center(t+1e-4),0)
 for step in range(601):
  c=center(step/1000);rad=math.sqrt(sum(x*x for x in c))
  assert rad-R>args.radial_min and rad+R<args.radial_max and inside(c,1)
 print('direction',direction,'Yin boundary crossing',t)
print('Analytic curl, peak normalization, crossing and radial support checks passed.')
