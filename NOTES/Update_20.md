### Update: (Nov 20 2020):
Graphics :-D
Very happy to see the first [avs renderings of the 2D density!](images/render-mesh-isentropic-vortex-initial-zoom-7.PNG) This is a 7-th order mesh, zoomed in on the center of the initial solution for the isentropic vortex.

For comparison - [here is a first order version of the same view](images/render-mesh-isentropic-vortex-initial-zoom.PNG). You can see the triangulation detail - note that there is a nested triangle configuration, a larger layer of triangles has within it a set of internal triangles. Those are the nodes of the finite elements, and the triangulation inside the larger triangles is a Delaunay triangulation, intended to be used as a linear approximation to the polynomial field. 

