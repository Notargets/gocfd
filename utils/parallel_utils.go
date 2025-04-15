package utils

import "fmt"

type MailBox[T any] struct {
	NP           int
	MessageChans []chan *DynBuffer[T]    // One for each thread
	PostMsgQs    []map[int]*DynBuffer[T] // One for each thread,
	// key is target thread
	ReceiveMsgQs []*DynBuffer[T] // One for each thread
	MailFlag     []bool          // MyThread receiver has messages in outbox
}

func NewMailBox[T any](NP int) *MailBox[T] {
	mb := &MailBox[T]{
		NP:           NP,
		MessageChans: make([]chan *DynBuffer[T], NP),
		PostMsgQs:    make([]map[int]*DynBuffer[T], NP),
		ReceiveMsgQs: make([]*DynBuffer[T], NP),
		MailFlag:     make([]bool, NP),
	}
	for n := 0; n < NP; n++ {
		mb.MessageChans[n] = make(chan *DynBuffer[T], NP) // Worst case is all-to-all
		mb.PostMsgQs[n] = make(map[int]*DynBuffer[T])
		mb.ReceiveMsgQs[n] = NewDynBuffer[T](0)
	}
	return mb
}

func (mb *MailBox[T]) PostMessage(myThread, targetThread int, msg T) {
	var (
		exists bool
		tgt    *DynBuffer[T]
	)
	if tgt, exists = mb.PostMsgQs[myThread][targetThread]; !exists {
		mb.PostMsgQs[myThread][targetThread] = NewDynBuffer[T](0)
	}
	tgt = mb.PostMsgQs[myThread][targetThread]
	tgt.Add(msg)
	if !mb.MailFlag[myThread] {
		mb.MailFlag[myThread] = true
	}
}

func (mb *MailBox[T]) PostMessageToAll(myThread int, msg T) {
	for k := 0; k < mb.NP; k++ {
		if k != myThread {
			mb.PostMessage(myThread, k, msg)
		}
	}
}
func (mb *MailBox[T]) DeliverMyMessages(myThread int) {
	if mb.MailFlag[myThread] {
		// fmt.Printf("Here in mailbox after MailFlag\n")
		for targetThread, msgBuffer := range mb.PostMsgQs[myThread] {
			if targetThread < 0 || targetThread > mb.NP-1 {
				panic(fmt.Sprintf("Target thread %d out of bounds", targetThread))
			}
			// fmt.Printf("Target thread, #msgs: %d, %d\n", targetThread,
			// 	len(msgBuffer.Cells()))
			// fmt.Printf("Message[%d]: %v\n", i, msg)
			mb.MessageChans[targetThread] <- msgBuffer
		}
		mb.MailFlag[myThread] = false
	}
}

func (mb *MailBox[T]) ReceiveMyMessages(myThread int) {
	for {
		select {
		case msgBuffer := <-mb.MessageChans[myThread]:
			// fmt.Println("Length of msgBuffer = ", len(msgBuffer.Cells()))
			for _, msg := range msgBuffer.Cells() {
				mb.ReceiveMsgQs[myThread].Add(msg)
			}
			msgBuffer.Reset() // Reset the originating buffer
		default:
			return
		}
	}
}
func (mb *MailBox[T]) ClearMyMessages(myThread int) {
	mb.ReceiveMsgQs[myThread].Reset()
}

type PartitionMap struct {
	MaxIndex       int // MaxIndex is partitioned into ParallelDegree partitions
	ParallelDegree int
	Partitions     [][2]int // Beginning and end index of partitions
}

func NewPartitionMap(ParallelDegree, maxIndex int) (pm *PartitionMap) {
	pm = &PartitionMap{
		MaxIndex:       maxIndex,
		ParallelDegree: ParallelDegree,
		Partitions:     make([][2]int, ParallelDegree),
	}
	for n := 0; n < ParallelDegree; n++ {
		pm.Partitions[n] = pm.Split1D(n)
	}
	return
}

func (pm *PartitionMap) GetBucket(kDim int) (bucketNum, min, max int) {
	_, bucketNum, min, max = pm.getBucketWithTryCount(kDim)
	return
}

func (pm *PartitionMap) getBucketWithTryCount(kDim int) (tryCount, bucketNum, min, max int) {
	// Initial guess
	bucketNum = int(float64(pm.ParallelDegree*kDim) / float64(pm.MaxIndex))
	for !(pm.Partitions[bucketNum][0] <= kDim && pm.Partitions[bucketNum][1] > kDim) {
		if pm.Partitions[bucketNum][0] > kDim {
			bucketNum--
		} else {
			bucketNum++
		}
		if bucketNum == -1 || bucketNum == pm.ParallelDegree {
			return 0, -1, 0, 0
		}
		tryCount++
	}
	min, max = pm.Partitions[bucketNum][0], pm.Partitions[bucketNum][1]
	/*
		if tryCount != 0 {
			fmt.Printf("bn, kDim, maxIndex, ParallelDegree, tryCount = %d, %d, %d, %d, %d\n",
				bucketNum, kDim, pm.MaxIndex, pm.ParallelDegree, tryCount)
		}
	*/
	return
}

func (pm *PartitionMap) GetBucketRange(bucketNum int) (kMin, kMax int) {
	kMin, kMax = pm.Partitions[bucketNum][0], pm.Partitions[bucketNum][1]
	return
}

func (pm *PartitionMap) GetLocalK(baseK int) (k, Kmax, bn int) {
	var (
		kmin, kmax int
	)
	bn, kmin, kmax = pm.GetBucket(baseK)
	Kmax = kmax - kmin
	k = baseK - kmin
	return
}

func (pm *PartitionMap) GetGlobalK(kLocal, bn int) (kGlobal int) {
	if bn == -1 {
		kGlobal = kLocal
		return
	}
	var (
		kMin = pm.Partitions[bn][0]
	)
	kGlobal = kMin + kLocal
	return
}

func (pm *PartitionMap) GetBucketDimension(bn int) (kMax int) {
	if bn == -1 {
		kMax = pm.MaxIndex
		return
	}
	var (
		k1, k2 = pm.GetBucketRange(bn)
	)
	kMax = k2 - k1
	return
}

func (pm *PartitionMap) Split1D(threadNum int) (bucket [2]int) {
	// This routine splits one dimension into c.ParallelDegree pieces, with a maximum imbalance of one item
	var (
		Npart            = pm.MaxIndex / (pm.ParallelDegree)
		startAdd, endAdd int
		remainder        int
	)
	remainder = pm.MaxIndex % pm.ParallelDegree
	if remainder != 0 { // spread the remainder over the first chunks evenly
		if threadNum+1 > remainder {
			startAdd = remainder
			endAdd = 0
		} else {
			startAdd = threadNum
			endAdd = 1
		}
	}
	bucket[0] = threadNum*Npart + startAdd
	bucket[1] = bucket[0] + Npart + endAdd
	return
}

type NeighborMsg struct {
	// From the perspective of the receiver
	KNeighborGlobal, KMyGlobal int
}

type NeighborNotifier struct {
	PartitionMap  *PartitionMap
	NeighborsEtoE [][3]int
	mb            *MailBox[*NeighborMsg]
}

func NewNeighborNotifier(PartitionMap *PartitionMap,
	NeighborsEtoE [][3]int) (nen *NeighborNotifier) {
	nen = &NeighborNotifier{
		mb:            NewMailBox[*NeighborMsg](PartitionMap.ParallelDegree),
		PartitionMap:  PartitionMap,
		NeighborsEtoE: NeighborsEtoE,
	}
	return
}

func (nen *NeighborNotifier) PostNotification(myThread, myKLocal int) {
	// The pattern here is:
	// for range messages {Post}; Deliver; blockWait; Receive
	var (
		myK       = nen.PartitionMap.GetGlobalK(myKLocal, myThread)
		Neighbors = nen.NeighborsEtoE[myK]
	)
	// The message to each neighbor is simply to connect to this node via it's
	// local EtoE list of element->face mappings.
	for _, nbrK := range Neighbors { // For each neighboring tri
		if nbrK == -1 {
			continue
		}
		targetThread, tgtKMin, tgtKMax := nen.PartitionMap.GetBucket(nbrK)
		// fmt.Printf("Posting %d to thread %d, Kmin %d, Kmax %d\n", nbrK,
		// 	targetThread, tgtKMin, tgtKMax)
		if nbrK < tgtKMin || nbrK >= tgtKMax {
			panic("Neighbor K out of range")
		}
		nen.mb.PostMessage(myThread, targetThread, &NeighborMsg{
			myK,
			nbrK,
		})
	}
}

func findEdge(NeighborList [3]int, k int) (edge int) {
	for i := 0; i < 3; i++ {
		if NeighborList[i] == k {
			return i
		}
	}
	panic("Edge not found")
	return
}

func (nen *NeighborNotifier) DeliverNotifications(myThread int) {
	// Must be called in myThread before receivers can receive notifications
	nen.mb.DeliverMyMessages(myThread)
}

func (nen *NeighborNotifier) ReadNotifications(myThread int) (faces [][3]int) {
	// Must be called after waiting for the DeliverNotifications is done in a
	// sync
	nen.mb.ReceiveMyMessages(myThread)
	for _, msg := range nen.mb.ReceiveMsgQs[myThread].Cells() {
		// Find which face the remote element connects to via EtoE
		// Neighbors := c.DFR.Tris.EtoE[msg.KMyGlobal]
		Neighbors := nen.NeighborsEtoE[msg.KMyGlobal]
		nEdge := findEdge(Neighbors, msg.KNeighborGlobal)
		faces = append(faces, [3]int{msg.KNeighborGlobal, msg.KMyGlobal, nEdge})
	}
	nen.mb.ClearMyMessages(myThread)
	return
}
