package mpc

import (
	"sync"
	"testing"
)

func TestCommunicationStatsAggregateDeltaAndReset(t *testing.T) {
	nets := ParallelNetworks{
		&Network{
			pid:           2,
			SentBytes:     map[int]uint64{0: 10, 1: 20},
			ReceivedBytes: map[int]uint64{0: 7, 1: 11},
			commSent:      map[int]int{0: 1, 1: 2},
			commReceived:  map[int]int{0: 3, 1: 4},
			loggingActive: true,
		},
		&Network{
			pid:           2,
			SentBytes:     map[int]uint64{0: 30},
			ReceivedBytes: map[int]uint64{1: 13},
			commSent:      map[int]int{0: 5},
			commReceived:  map[int]int{1: 6},
			loggingActive: true,
		},
	}

	got := nets.GetCommunicationStats()
	want := CommunicationStats{SentBytes: 60, ReceivedBytes: 31, SentMessages: 8, ReceivedMessages: 13}
	if got != want {
		t.Fatalf("aggregate: got %+v want %+v", got, want)
	}

	start := CommunicationStats{SentBytes: 50, ReceivedBytes: 30, SentMessages: 7, ReceivedMessages: 10}
	wantDelta := CommunicationStats{SentBytes: 10, ReceivedBytes: 1, SentMessages: 1, ReceivedMessages: 3}
	if delta := got.Sub(start); delta != wantDelta {
		t.Fatalf("delta: got %+v want %+v", delta, wantDelta)
	}

	nets.ResetNetworkLog()
	if zero := nets.GetCommunicationStats(); zero != (CommunicationStats{}) {
		t.Fatalf("reset: got %+v want all-zero counters", zero)
	}
}

func TestCommunicationStatsConcurrentAccess(t *testing.T) {
	net := &Network{
		SentBytes:     make(map[int]uint64),
		ReceivedBytes: make(map[int]uint64),
		commSent:      make(map[int]int),
		commReceived:  make(map[int]int),
		loggingActive: true,
	}
	nets := ParallelNetworks{net}

	const workers = 8
	const updatesPerWorker = 1_000
	start := make(chan struct{})
	var wg sync.WaitGroup
	for worker := 0; worker < workers; worker++ {
		wg.Add(1)
		go func(peer int) {
			defer wg.Done()
			<-start
			for update := 0; update < updatesPerWorker; update++ {
				net.UpdateSenderLog(peer, 3)
				net.UpdateReceiverLog(peer, 5)
				_ = nets.GetCommunicationStats()
			}
		}(worker)
	}
	close(start)
	wg.Wait()

	want := CommunicationStats{
		SentBytes:        workers * updatesPerWorker * 3,
		ReceivedBytes:    workers * updatesPerWorker * 5,
		SentMessages:     workers * updatesPerWorker,
		ReceivedMessages: workers * updatesPerWorker,
	}
	if got := nets.GetCommunicationStats(); got != want {
		t.Fatalf("concurrent aggregate: got %+v want %+v", got, want)
	}
}

func TestCommunicationStatsConcurrentResetAndToggle(t *testing.T) {
	net := &Network{
		SentBytes:     make(map[int]uint64),
		ReceivedBytes: make(map[int]uint64),
		commSent:      make(map[int]int),
		commReceived:  make(map[int]int),
		loggingActive: true,
	}
	nets := ParallelNetworks{net}

	const iterations = 1_000
	start := make(chan struct{})
	var wg sync.WaitGroup
	for _, run := range []func(){
		func() {
			for i := 0; i < iterations; i++ {
				net.UpdateSenderLog(1, 3)
				net.UpdateReceiverLog(1, 5)
			}
		},
		func() {
			for i := 0; i < iterations; i++ {
				nets.ResetNetworkLog()
				_ = nets.GetCommunicationStats()
			}
		},
		func() {
			for i := 0; i < iterations; i++ {
				net.DisableLogging()
				net.EnableLogging()
			}
		},
	} {
		wg.Add(1)
		go func(run func()) {
			defer wg.Done()
			<-start
			run()
		}(run)
	}
	close(start)
	wg.Wait()

	// Leave a deterministic final state after the deliberately unordered phase.
	net.EnableLogging()
	nets.ResetNetworkLog()
	net.UpdateSenderLog(1, 3)
	net.UpdateReceiverLog(1, 5)
	want := CommunicationStats{SentBytes: 3, ReceivedBytes: 5, SentMessages: 1, ReceivedMessages: 1}
	if got := nets.GetCommunicationStats(); got != want {
		t.Fatalf("post-reset aggregate: got %+v want %+v", got, want)
	}
}
