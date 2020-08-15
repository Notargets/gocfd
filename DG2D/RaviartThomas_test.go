package DG2D

import (
	"fmt"
	"testing"
)

func TestRTElement(t *testing.T) {
	N := 6
	fmt.Printf("N = %d\n", N)
	fmt.Printf("i\tj\n")
	var count int
	for i := 0; i <= N; i++ {
		for j := 0; j <= (N - i); j++ {
			fmt.Printf("%d\t%d\n", i, j)
			//DG1D.JacobiP()
			count++
		}
	}
	fmt.Printf("Expected Count = %d\n", (N+1)*(N+2)/2)
	fmt.Printf("Count = %d\n", count)
}
