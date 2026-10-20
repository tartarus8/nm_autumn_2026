package labio

import (
	"flag"
	"fmt"
	"os"

	"nm/internal/linalg"
)

func Main(run func(bool) error) {
	extra := flag.Bool("extraInfo", false, "print intermediate steps and checks")
	flag.Parse()
	if err := run(*extra); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
}

type Parameters struct {
	Epsilon       float64
	MaxIterations int
}

func (p Parameters) Check() error {
	if p.Epsilon <= 0 || p.MaxIterations <= 0 {
		return fmt.Errorf("epsilon and max_iterations must be positive")
	}
	return nil
}

func PrintHistory(name string, history []linalg.Step, reference []float64) {
	fmt.Print(name + ": iteration | estimate | residual")
	if reference != nil {
		fmt.Print(" | reference error")
	}
	fmt.Println(" | approximation")
	for _, s := range history {
		fmt.Printf("%5d  %12.5g  %12.5g", s.Iteration, s.Estimate, s.Residual)
		if reference != nil {
			fmt.Printf("  %12.5g", linalg.Distance(s.X, reference))
		}
		fmt.Printf("  %.10g\n", s.X)
	}
}
