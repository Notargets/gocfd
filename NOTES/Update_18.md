## Update: Nov 2 2020):
Divergence is now tested correct for transformed triangles, including the use of the ||n|| scale factor to carry
((Flux) dot (face normal)) correctly into the RT element degree of freedom for edges.

I'm working now on the actual data structures that will efficiently compute the flux values, etc., with an eye on the
memory / CPU/GPU performance tradeoffs. Contiguous space matrix multiplications are supported well by GPU, so I'm focusing
on making most everything a contiguous space matrix multiply.


[Back to Index](../NOTES_Index.md)
