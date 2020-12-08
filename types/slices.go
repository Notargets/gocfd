package types

func GrowSlice(mysliceI interface{}, newCap int) (biggerSliceI interface{}) {
	switch myslice := mysliceI.(type) {
	case []EdgeKey:
		l := len(myslice)
		var biggerSlice []EdgeKey
		if l >= newCap {
			return myslice
		} else {
			biggerSlice = make([]EdgeKey, newCap)
		}
		for i, val := range myslice {
			biggerSlice[i] = val
		}
		return biggerSlice
	case []EdgeInt:
		l := len(myslice)
		var biggerSlice []EdgeInt
		if l >= newCap {
			return myslice
		} else {
			biggerSlice = make([]EdgeInt, newCap)
		}
		for i, val := range myslice {
			biggerSlice[i] = val
		}
		return biggerSlice
	case []int:
		l := len(myslice)
		var biggerSlice []int
		if l >= newCap {
			return myslice
		} else {
			biggerSlice = make([]int, newCap)
		}
		for i, val := range myslice {
			biggerSlice[i] = val
		}
		return biggerSlice
	case []float64:
		l := len(myslice)
		var biggerSlice []float64
		if l >= newCap {
			return myslice
		} else {
			biggerSlice = make([]float64, newCap)
		}
		for i, val := range myslice {
			biggerSlice[i] = val
		}
		return biggerSlice
	case []string:
		l := len(myslice)
		var biggerSlice []string
		if l >= newCap {
			return myslice
		} else {
			biggerSlice = make([]string, newCap)
		}
		for i, val := range myslice {
			biggerSlice[i] = val
		}
		return biggerSlice
	}
	return
}
