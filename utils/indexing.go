package utils

type Index []int

func (I Index) ApplyFunc(f func(val int) int) (r Index) {
	r = make([]int, len(I))
	for i, val := range I {
		r[i] = f(val)
	}
	return
}
