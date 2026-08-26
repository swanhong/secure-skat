package protocol

import "time"

const NoTimingPhenotype = -1

type TimingEvent struct {
	Phase          string
	ParentPhase    string
	PhenotypeIndex int
	Elapsed        time.Duration
}

type TimingObserver func(TimingEvent)

func selectTimingObserver(observers []TimingObserver) TimingObserver {
	if len(observers) == 0 {
		return nil
	}
	if len(observers) > 1 {
		panic("at most one timing observer is supported")
	}
	return observers[0]
}

func startTiming(
	observer TimingObserver,
	phase string,
	parentPhase string,
	phenotypeIndex int,
) func() {
	if observer == nil {
		return func() {}
	}

	started := time.Now()
	return func() {
		observer(TimingEvent{
			Phase:          phase,
			ParentPhase:    parentPhase,
			PhenotypeIndex: phenotypeIndex,
			Elapsed:        time.Since(started),
		})
	}
}
