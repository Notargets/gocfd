package gonudg

import (
	"math"
	"testing"
)

func TestRSTtoABC(t *testing.T) {
	tests := []struct {
		name string
		r, s, t float64
		wantA, wantB, wantC float64
		tol float64
	}{
		{
			name: "origin",
			r: 0, s: 0, t: 0,
			wantA: -1.0/3.0, wantB: -1.0/3.0, wantC: 0,
			tol: 1e-14,
		},
		{
			name: "vertex1",
			r: -1, s: -1, t: -1,
			wantA: -1, wantB: -1, wantC: -1,
			tol: 1e-14,
		},
		{
			name: "vertex2",
			r: 1, s: -1, t: -1,
			wantA: 1, wantB: -1, wantC: -1,
			tol: 1e-14,
		},
		{
			name: "vertex3",
			r: -1, s: 1, t: -1,
			wantA: -1, wantB: 1, wantC: -1,
			tol: 1e-14,
		},
		{
			name: "vertex4",
			r: -1, s: -1, t: 1,
			wantA: -1, wantB: -1, wantC: 1,
			tol: 1e-14,
		},
		{
			name: "degenerate_case_s_plus_t_zero",
			r: 0.5, s: -0.5, t: 0.5,
			wantA: -1, wantB: -1, wantC: 0.5,
			tol: 1e-14,
		},
		{
			name: "degenerate_case_t_one",
			r: 0, s: 0, t: 1,
			wantA: -1, wantB: -1, wantC: 1,
			tol: 1e-14,
		},
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			// Test single point version
			a, b, c := RSTtoABCSingle(tc.r, tc.s, tc.t)
			if math.Abs(a-tc.wantA) > tc.tol {
				t.Errorf("RSTtoABCSingle: a = %v, want %v", a, tc.wantA)
			}
			if math.Abs(b-tc.wantB) > tc.tol {
				t.Errorf("RSTtoABCSingle: b = %v, want %v", b, tc.wantB)
			}
			if math.Abs(c-tc.wantC) > tc.tol {
				t.Errorf("RSTtoABCSingle: c = %v, want %v", c, tc.wantC)
			}

			// Test array version
			rArr := []float64{tc.r}
			sArr := []float64{tc.s}
			tArr := []float64{tc.t}
			aArr, bArr, cArr := RSTtoABC(rArr, sArr, tArr)
			if math.Abs(aArr[0]-tc.wantA) > tc.tol {
				t.Errorf("RSTtoABC: a[0] = %v, want %v", aArr[0], tc.wantA)
			}
			if math.Abs(bArr[0]-tc.wantB) > tc.tol {
				t.Errorf("RSTtoABC: b[0] = %v, want %v", bArr[0], tc.wantB)
			}
			if math.Abs(cArr[0]-tc.wantC) > tc.tol {
				t.Errorf("RSTtoABC: c[0] = %v, want %v", cArr[0], tc.wantC)
			}
		})
	}
}

