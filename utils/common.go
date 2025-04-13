package utils

const (
	NODETOL = 1.e-12
)

type EvalOp uint8

const (
	Equal EvalOp = iota
	Less
	Greater
	LessOrEqual
	GreaterOrEqual
)

func Split1D(iMax, numThreads, threadNum int) (istart, iend int) {
	var (
		Npart = iMax / numThreads
	)
	istart = threadNum * Npart
	iend = istart + Npart
	if threadNum == numThreads-1 {
		iend = iMax
	}
	return
}

func Split1DMaxChunk(iMax, numThreads int) (maxChunk int) {
	// Finds the max chunk size within split
	var (
		Npart = iMax / numThreads
	)
	maxChunk = Npart + iMax%numThreads
	return
}

type DynBuffer[T any] struct {
	buffer []T
}

// NewDynBuffer creates a structure with an initial reserved capacity.
func NewDynBuffer[T any](cap int) *DynBuffer[T] {
	return &DynBuffer[T]{
		buffer: make([]T, 0, cap),
	}
}

// Add appends a new shocked cell index.
func (s *DynBuffer[T]) Add(datum T) {
	s.buffer = append(s.buffer, datum)
}

// Reset efficiently clears the list without reallocating memory.
func (s *DynBuffer[T]) Reset() {
	s.buffer = s.buffer[:0]
}

// Cells returns the underlying slice.
func (s *DynBuffer[T]) Cells() []T {
	return s.buffer
}
