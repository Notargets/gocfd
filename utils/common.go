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
