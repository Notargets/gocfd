import dolfin
import numpy
#import pylab
from matplotlib import pyplot

# define an exact stream function
psi_exact_str = 'x[1]<=pi ? epsilon*cos(x[0])-(1.0/(cosh((x[1]-0.5*pi)/delta)*cosh((x[1]-0.5*pi)/delta)))/delta : epsilon*cos(x[0]) + (1.0/(cosh((1.5*pi-x[1])/delta)*cosh((1.5*pi-x[1])/delta)))/delta'   
epsilon = 0.05
delta = numpy.pi/15.0

#psi_exact = dolfin.Expression(psi_exact_str,epsilon=epsilon,delta=delta)
psi_exact = dolfin.Expression(psi_exact_str,degree=int(10),epsilon=epsilon,delta=delta)

# define the number of elements and polynomial degree of basis functions
n_elements = [5,10,20,40,80]

pp = numpy.arange(1,5)

# define the mesh
xBounds = numpy.array([0.,2.0*numpy.pi])
yBounds = numpy.array([0.,2.0*numpy.pi])

div_max = numpy.zeros([len(pp),len(n_elements)])

for kp,p in enumerate(pp):
    for kn,n in enumerate(n_elements):
        print ('p = ' + str(p) + '   k = ' + str(n))
        mesh = dolfin.RectangleMesh(dolfin.Point(xBounds[0],yBounds[0],0.0),dolfin.Point(xBounds[1],yBounds[1],0.0),n,n,'crossed')

        # define the function spaces
        U = dolfin.FunctionSpace(mesh,'RT',int(p))      # velocity
        PSI = dolfin.FunctionSpace(mesh,'CG',int(p))    # stream function
        P = dolfin.FunctionSpace(mesh,'DG',int(p-1))    # divergence of velocity


        # compute the finite element approximation of the analytical stream function
        psi_h = dolfin.project(psi_exact,PSI,solver_type='lu')

        # initialize the velocity and its divergence
        u_h = dolfin.Function(U)
        p_h = dolfin.Function(P)

        # compute the velocity 
        u_h.assign(dolfin.project(dolfin.curl(psi_h),U,solver_type='lu'))

        # compute the divergence 
        p_h.assign(dolfin.project(dolfin.div(u_h),P,solver_type='lu'))


        # get the maximum value of the divergence
        #div_max[kp,kn] = numpy.abs(p_h.vector().array()).max()
        div_max[kp,kn] = numpy.abs(p_h.vector().get_local()).max()



pyplot.semilogy(pp,div_max)
pyplot.legend(n_elements,loc=4)
pyplot.xlabel('p')
pyplot.show()
