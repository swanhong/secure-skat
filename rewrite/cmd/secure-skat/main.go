package main

import (
	"flag"
	"fmt"
	"os"

	securecrypto "github.com/hhcho/sfgwas/crypto"
)

func main() {
	if len(os.Args) < 2 || os.Args[1] != "run" {
		fmt.Fprintln(os.Stderr, "usage: secure-skat run [options]")
		os.Exit(2)
	}

	flags := flag.NewFlagSet("run", flag.ContinueOnError)
	options := runOptions{}
	flags.IntVar(&options.Party, "party", -1, "party ID: 0, 1, or 2")
	flags.IntVar(&options.Lane, "lane", 0, "Workbench chromosome lane index")
	flags.StringVar(&options.Input, "input", "", "preprocessed input directory")
	flags.StringVar(&options.Output, "output", "", "party A secure result CSV")
	flags.StringVar(
		&options.TimingOutput,
		"timing-output",
		"",
		"optional benchmark timing CSV for this party",
	)
	flags.IntVar(&options.PortBase, "port-base", 18000, "first localhost party port")
	flags.StringVar(&options.SharedKeys, "shared-keys", "", "shared PRG key directory")
	flags.StringVar(
		&options.CKKS,
		"ckks",
		securecrypto.CKKSParamsPN14QP436S45,
		"CKKS parameter name",
	)
	flags.IntVar(&options.DataBits, "data-bits", 0, "MPC data bits (required)")
	flags.IntVar(&options.FractionalBits, "frac-bits", 30, "MPC fractional bits")
	flags.IntVar(&options.Probes, "probes", 50, "Hutchinson probe count")
	flags.Int64Var(&options.Seed, "seed", 42, "trace probe seed")
	if err := flags.Parse(os.Args[2:]); err != nil {
		os.Exit(2)
	}

	if err := runSecure(options); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
}
