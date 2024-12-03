## Update: [9/18/21]

Working on an enhanced flux transfer scheme that promises to protect against odd-even decoupling ("wiggles") while
minimizing artificial dissipation that can destroy turbulence fields, etc.

The scheme I've located is by [Xue-Song Li, et al](../research/filters_and_flux_limiters/roe-er-li.pdf) from an AIAA paper
in 2020, who described the "Roe-ER" flux, which combines features of the Harten entropy fix and the rotated Roe flux schemes.
It looks very efficient and in the papers it seems to do a very good job at the wiggle issue while delivering excellent
accuracy.

My remaining concern/question is about whether the use of a flux limiter for transfer of flux at the element faces is
sufficient to damp oscillations, or whether I'll have to couple elements more deeply, thus removing a principal advantage of
Nodal Galerkin methods - that of "compact support". We'll soon see!



[Back to Index](../NOTES_Index.md)
