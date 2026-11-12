package main

import (
	"fmt"
	"os"

	"github.com/hhcho/sfgwas/rvas/jointstat"
)

func main() {
	if err := jointstat.Main(os.Args[1:]); err != nil {
		fmt.Fprintln(os.Stderr, "joint-stat:", err)
		os.Exit(1)
	}
}
