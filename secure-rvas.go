package main

import (
	"fmt"
	"os"

	"github.com/hhcho/sfgwas/rvas/jointstat"
	"github.com/hhcho/sfgwas/rvas/workflow"
)

const (
	version    = "1.0.2"
	updateDate = "11-11-2026"
)

func main() {
	fmt.Printf("run secure-skat version: %s (update date: %s)\n", version, updateDate)
	var err error
	if len(os.Args) > 1 && os.Args[1] == "joint-stat" {
		err = jointstat.Main(os.Args[2:])
	} else {
		err = workflow.Main(os.Args[1:])
	}
	if err != nil {
		fmt.Fprintf(os.Stderr, "secure-rvas: %v\n", err)
		os.Exit(1)
	}
}
