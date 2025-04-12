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

type DynIntBuffer struct {
	buffer []int
}

// NewDynBuffer creates a structure with an initial reserved capacity.
func NewDynIntBuffer(initialCapacity int) *DynIntBuffer {
	return &DynIntBuffer{
		buffer: make([]int, 0, initialCapacity),
	}
}

// Add appends a new shocked cell index.
func (s *DynIntBuffer) Add(datum int) {
	s.buffer = append(s.buffer, datum)
}

// Reset efficiently clears the list without reallocating memory.
func (s *DynIntBuffer) Reset() {
	s.buffer = s.buffer[:0]
}

// Cells returns the underlying slice.
func (s *DynIntBuffer) Cells() []int {
	return s.buffer
}
