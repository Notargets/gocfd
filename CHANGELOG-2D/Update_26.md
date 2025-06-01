## Update: (Dec 31, 2020):

New Year's Day progress!

I partitioned the 2D solver by elements for enhanced parallelism. The solver now computes the RHS and time stepping fully
parallel, but must still synchronize between each time sub-step (within the Runge-Kutta solver) to exchange data at edges,
which is also done in parallel. So there are now two discrete stages/types of parallelism, one for the full domain of
elements, and one for the edge exchanges and flux computation.

As a result of the partitioning, the parallelism is far better and we can get much faster solution times. The level of
parallelism is well suited to machines with less than 100 processors. For more parallelism, we will need to use a mesh
partitioning algorithm that selects the element partitions so as to minimize the surface area shared between partitions, such
as the commonly used tool "metis". When the elements are partitioned in that way, we can isolate the time spent in
synchronization to exactly the minimum necessary. By comparison, now we're doing all edge flux computation during the
synchronization period.


[Back to Index](../CHANGELOG-2D.md)
