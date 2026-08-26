package main

import (
	"encoding/csv"
	"os"
	"path/filepath"
	"strconv"
	"time"

	"github.com/hhcho/sfgwas/rewrite/protocol"
)

const noTimingIndex = -1

var timingFields = []string{
	"component",
	"scope",
	"chromosome",
	"party",
	"lane",
	"phase",
	"parent_phase",
	"event_index",
	"batch_index",
	"batch_width",
	"batch_gene_count",
	"gene_index",
	"gene_id",
	"phenotype_index",
	"phenotype_name",
	"trace_mode",
	"measurement_kind",
	"elapsed_seconds",
	"sample_count_a",
	"sample_count_b",
	"public_variant_count",
	"private_variant_count",
	"ckks",
	"data_bits",
	"frac_bits",
	"probes",
	"status",
}

type timingRow struct {
	Component           string
	Scope               string
	Chromosome          string
	Party               int
	Lane                int
	Phase               string
	ParentPhase         string
	EventIndex          int
	BatchIndex          int
	BatchWidth          int
	BatchGeneCount      int
	GeneIndex           int
	GeneID              string
	PhenotypeIndex      int
	PhenotypeName       string
	TraceMode           string
	MeasurementKind     string
	Elapsed             time.Duration
	SampleCountA        int
	SampleCountB        int
	PublicVariantCount  int
	PrivateVariantCount int
	CKKS                string
	DataBits            int
	FractionalBits      int
	Probes              int
	Status              string
}

type timingRecorder struct {
	options runOptions
	input   *preprocessedInput
	now     func() time.Time
	rows    []timingRow
}

type timingSpan struct {
	recorder *timingRecorder
	row      timingRow
	started  time.Time
}

type batchTimingSpan struct {
	recorder             *timingRecorder
	input                preprocessedInput
	batchIndex           int
	batch                protocol.GeneBatch
	privateVariantCounts []int
	phase                string
	parentPhase          string
	phenotypeIndex       int
	started              time.Time
}

func newTimingRecorder(options runOptions) *timingRecorder {
	return &timingRecorder{
		options: options,
		now:     time.Now,
		rows:    []timingRow{},
	}
}

func (recorder *timingRecorder) setInput(input preprocessedInput) {
	recorder.input = &input
}

func (recorder *timingRecorder) start(row timingRow) *timingSpan {
	return &timingSpan{
		recorder: recorder,
		row:      row,
		started:  recorder.now(),
	}
}

func (recorder *timingRecorder) startPhase(
	scope string,
	phase string,
	parentPhase string,
) *timingSpan {
	row := newTimingRow(scope, phase, parentPhase)
	return recorder.start(row)
}

func (recorder *timingRecorder) startBatchPhase(
	input preprocessedInput,
	batchIndex int,
	batch protocol.GeneBatch,
	privateVariantCounts []int,
	phase string,
	parentPhase string,
	phenotypeIndex int,
) *batchTimingSpan {
	return &batchTimingSpan{
		recorder:             recorder,
		input:                input,
		batchIndex:           batchIndex,
		batch:                batch,
		privateVariantCounts: privateVariantCounts,
		phase:                phase,
		parentPhase:          parentPhase,
		phenotypeIndex:       phenotypeIndex,
		started:              recorder.now(),
	}
}

func (span *timingSpan) finish(status string) {
	span.row.Elapsed = span.recorder.now().Sub(span.started)
	span.row.Status = status
	span.recorder.add(span.row)
}

func (span *batchTimingSpan) finish(status string) {
	span.recorder.addBatchEvent(
		span.input,
		span.batchIndex,
		span.batch,
		span.privateVariantCounts,
		protocol.TimingEvent{
			Phase:          span.phase,
			ParentPhase:    span.parentPhase,
			PhenotypeIndex: span.phenotypeIndex,
			Elapsed:        span.recorder.now().Sub(span.started),
		},
		status,
	)
}

func (recorder *timingRecorder) batchObserver(
	input preprocessedInput,
	batchIndex int,
	batch protocol.GeneBatch,
	privateVariantCounts []int,
) protocol.TimingObserver {
	return func(event protocol.TimingEvent) {
		recorder.addBatchEvent(
			input,
			batchIndex,
			batch,
			privateVariantCounts,
			event,
			"success",
		)
	}
}

func (recorder *timingRecorder) add(row timingRow) {
	if row.Component == "" {
		row.Component = "secure"
	}
	row.Party = recorder.options.Party
	row.Lane = recorder.options.Lane
	if row.MeasurementKind == "" {
		row.MeasurementKind = "actual"
	}
	row.CKKS = recorder.options.CKKS
	row.DataBits = recorder.options.DataBits
	row.FractionalBits = recorder.options.FractionalBits
	row.Probes = recorder.options.Probes
	if recorder.input != nil {
		if row.Chromosome == "" && len(recorder.input.Genes) > 0 {
			row.Chromosome = recorder.input.Genes[0].Chromosome
		}
		if row.SampleCountA == 0 {
			row.SampleCountA = recorder.input.Metadata.SampleCountA
		}
		if row.SampleCountB == 0 {
			row.SampleCountB = recorder.input.Metadata.SampleCountB
		}
	}
	row.EventIndex = len(recorder.rows)
	recorder.rows = append(recorder.rows, row)
}

func (recorder *timingRecorder) addBatchEvent(
	input preprocessedInput,
	batchIndex int,
	batch protocol.GeneBatch,
	privateVariantCounts []int,
	event protocol.TimingEvent,
	status string,
) {
	publicVariantCount := 0
	privateVariantCount := noTimingIndex
	if len(privateVariantCounts) == len(batch.GeneIndices) {
		privateVariantCount = 0
	}
	for position, geneIndex := range batch.GeneIndices {
		publicVariantCount += input.DataParams.Genes[geneIndex].VariantCount
		if privateVariantCount >= 0 {
			privateVariantCount += privateVariantCounts[position]
		}
	}

	traceMode := "exact"
	if batch.W >= recorder.options.Probes {
		traceMode = "hutchinson"
	}
	base := timingRow{
		Scope:               "batch",
		Phase:               event.Phase,
		ParentPhase:         event.ParentPhase,
		BatchIndex:          batchIndex,
		BatchWidth:          batch.W,
		BatchGeneCount:      len(batch.GeneIndices),
		GeneIndex:           noTimingIndex,
		PhenotypeIndex:      event.PhenotypeIndex,
		TraceMode:           traceMode,
		Elapsed:             event.Elapsed,
		PublicVariantCount:  publicVariantCount,
		PrivateVariantCount: privateVariantCount,
		Status:              status,
	}
	if event.PhenotypeIndex >= 0 {
		base.Scope = "batch_phenotype"
		base.PhenotypeName = input.Metadata.PhenotypeColumns[event.PhenotypeIndex]
	}
	recorder.add(base)

	geneCount := len(batch.GeneIndices)
	if geneCount == 0 {
		return
	}
	for position, geneIndex := range batch.GeneIndices {
		gene := input.Genes[geneIndex]
		geneRow := base
		geneRow.Scope = "gene"
		if event.PhenotypeIndex >= 0 {
			geneRow.Scope = "gene_phenotype"
		}
		geneRow.BatchGeneCount = 1
		geneRow.GeneIndex = geneIndex
		geneRow.GeneID = gene.ID
		geneRow.MeasurementKind = "amortized"
		geneRow.Elapsed = event.Elapsed / time.Duration(geneCount)
		geneRow.PublicVariantCount = input.DataParams.Genes[geneIndex].VariantCount
		geneRow.PrivateVariantCount = noTimingIndex
		if len(privateVariantCounts) == len(batch.GeneIndices) {
			geneRow.PrivateVariantCount = privateVariantCounts[position]
		}
		recorder.add(geneRow)
	}
}

func newTimingRow(scope string, phase string, parentPhase string) timingRow {
	return timingRow{
		Scope:               scope,
		Party:               noTimingIndex,
		Lane:                noTimingIndex,
		Phase:               phase,
		ParentPhase:         parentPhase,
		BatchIndex:          noTimingIndex,
		BatchWidth:          noTimingIndex,
		BatchGeneCount:      noTimingIndex,
		GeneIndex:           noTimingIndex,
		PhenotypeIndex:      noTimingIndex,
		PublicVariantCount:  noTimingIndex,
		PrivateVariantCount: noTimingIndex,
	}
}

func (recorder *timingRecorder) write(path string) error {
	if path == "" {
		return nil
	}
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		return err
	}

	temporary := path + ".tmp"
	file, err := os.Create(temporary)
	if err != nil {
		return err
	}
	writer := csv.NewWriter(file)
	if err := writer.Write(timingFields); err != nil {
		file.Close()
		return err
	}
	for _, row := range recorder.rows {
		if err := writer.Write(row.csvRecord()); err != nil {
			file.Close()
			return err
		}
	}
	writer.Flush()
	if err := writer.Error(); err != nil {
		file.Close()
		return err
	}
	if err := file.Close(); err != nil {
		return err
	}
	return os.Rename(temporary, path)
}

func (row timingRow) csvRecord() []string {
	return []string{
		row.Component,
		row.Scope,
		row.Chromosome,
		optionalTimingInteger(row.Party),
		optionalTimingInteger(row.Lane),
		row.Phase,
		row.ParentPhase,
		strconv.Itoa(row.EventIndex),
		optionalTimingInteger(row.BatchIndex),
		optionalTimingInteger(row.BatchWidth),
		optionalTimingInteger(row.BatchGeneCount),
		optionalTimingInteger(row.GeneIndex),
		row.GeneID,
		optionalTimingInteger(row.PhenotypeIndex),
		row.PhenotypeName,
		row.TraceMode,
		row.MeasurementKind,
		strconv.FormatFloat(row.Elapsed.Seconds(), 'f', 9, 64),
		optionalTimingInteger(row.SampleCountA),
		optionalTimingInteger(row.SampleCountB),
		optionalTimingInteger(row.PublicVariantCount),
		optionalTimingInteger(row.PrivateVariantCount),
		row.CKKS,
		optionalTimingInteger(row.DataBits),
		optionalTimingInteger(row.FractionalBits),
		optionalTimingInteger(row.Probes),
		row.Status,
	}
}

func optionalTimingInteger(value int) string {
	if value < 0 {
		return ""
	}
	return strconv.Itoa(value)
}
