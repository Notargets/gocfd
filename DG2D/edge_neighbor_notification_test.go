package DG2D

import (
	"sync"
	"testing"

	"github.com/notargets/gocfd/utils"
	"github.com/stretchr/testify/assert"
)

func TestNeighborNotify(t *testing.T) {
	var (
		N        = 0
		NThreads = 5
	)
	if !testing.Verbose() {
		return
	}
	meshFile := "test_data/test_10tris_centered.neu"
	dfr := NewDFR2D(N, false, meshFile)
	Partitions := utils.NewPartitionMap(NThreads, dfr.K)
	t.Log(dfr.Tris.EtoE)

	nen := utils.NewNeighborNotifier(Partitions, dfr.Tris.EtoE)
	wg := sync.WaitGroup{}
	PostDeliver := func(myThread int) {
		defer wg.Done()
		// KMax := c.RK.Kmax[myThread]
		KMax := Partitions.GetBucketDimension(myThread)
		for k := 0; k < KMax; k++ {
			nen.PostNotification(myThread, k)
		}
		nen.DeliverNotifications(myThread)
	}
	var myThread int
	for myThread = 0; myThread < NThreads; myThread++ {
		// myThread = 0
		wg.Add(1)
		go PostDeliver(myThread)
	}
	wg.Wait()
	Receive := func(myThread int) {
		defer wg.Done()
		t.Logf("Messages for thread: %d\n", myThread)
		for _, conn := range nen.ReadNotifications(myThread) {
			neighborGlobalK := conn[0]
			myGlobalK := conn[1]
			myFace := conn[2]
			t.Logf("RemoteElement[%d] ConnectsToMy[%d] OnFace[%d]\n",
				neighborGlobalK, myGlobalK, myFace)
			assert.Equal(t, dfr.Tris.EtoE[myGlobalK][myFace], neighborGlobalK)
		}
	}
	for myThread = 0; myThread < NThreads; myThread++ {
		wg.Add(1)
		go Receive(myThread)
	}
	wg.Wait()
}
