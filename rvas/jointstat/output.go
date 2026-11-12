package jointstat

import (
	"encoding/csv"
	"encoding/json"
	"os"
	"path/filepath"
	"strconv"
)

func atomicWrite(path string, write func(*os.File) error) error {
	if err := os.MkdirAll(filepath.Dir(path), 0700); err != nil {
		return err
	}
	f, err := os.CreateTemp(filepath.Dir(path), ".joint-stat-")
	if err != nil {
		return err
	}
	defer os.Remove(f.Name())
	if err = write(f); err != nil {
		f.Close()
		return err
	}
	if err = f.Close(); err != nil {
		return err
	}
	return os.Rename(f.Name(), path)
}
func writeJSON(path string, value interface{}) error {
	return atomicWrite(path, func(f *os.File) error { enc := json.NewEncoder(f); enc.SetIndent("", "  "); return enc.Encode(value) })
}
func writeTSV(path string, header []string, rows [][]string) error {
	return atomicWrite(path, func(f *os.File) error {
		w := csv.NewWriter(f)
		w.Comma = '\t'
		w.Write(header)
		w.WriteAll(rows)
		return w.Error()
	})
}
func writeResults(path string, results []row) error {
	var rows [][]string
	for _, r := range results {
		rows = append(rows, []string{r.Ancestry, strconv.FormatFloat(r.MAF, 'g', -1, 64), r.Annotation, strconv.Itoa(r.Chromosome), "analysis_genes", r.Panel, r.Series, r.Lower, r.Upper, strconv.Itoa(r.Count), "released"})
	}
	return writeTSV(path, []string{"ancestry", "maf_threshold", "annotation", "chromosome", "gene_set", "panel", "series", "lower", "upper", "count", "status"}, rows)
}
func writeMetrics(path string, metrics []measurement) error {
	var rows [][]string
	for _, m := range metrics {
		rows = append(rows, []string{m.Stage, m.Ancestry, strconv.Itoa(m.Chromosome), strconv.FormatFloat(m.MAF, 'g', -1, 64), strconv.FormatFloat(m.Seconds, 'f', 6, 64), strconv.FormatUint(m.Sent, 10), strconv.FormatUint(m.Received, 10), strconv.FormatUint(m.SentMessages, 10), strconv.FormatUint(m.ReceivedMessages, 10)})
	}
	return writeTSV(path, []string{"stage", "ancestry", "chromosome", "maf_threshold", "seconds", "sent_bytes", "received_bytes", "sent_messages", "received_messages"}, rows)
}
