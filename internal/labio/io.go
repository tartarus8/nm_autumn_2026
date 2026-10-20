package labio

import (
	"encoding/json"
	"flag"
	"fmt"
	"os"
	"path/filepath"
	"runtime"

	"nm/internal/linalg"
)

type Options struct {
	Input     string
	ExtraInfo bool
}

func Main(run func(Options) error) {
	_, source, _, _ := runtime.Caller(1)
	dir := filepath.Dir(source)
	input := flag.String("input", filepath.Join(dir, "tests.json"), "input JSON file")
	extra := flag.Bool("extraInfo", false, "print intermediate steps and checks")
	flag.Parse()
	if err := run(Options{Input: *input, ExtraInfo: *extra}); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
}

func ReadTasks[T any](path string) ([]T, error) {
	data, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	var input struct {
		Tasks []T `json:"tasks"`
	}
	if err = json.Unmarshal(data, &input); err != nil {
		return nil, err
	}
	return input.Tasks, nil
}

type Parameters struct {
	Epsilon       float64 `json:"epsilon"`
	MaxIterations int     `json:"max_iterations"`
}

func (p Parameters) Check() error {
	if p.Epsilon <= 0 || p.MaxIterations <= 0 {
		return fmt.Errorf("epsilon and max_iterations must be positive")
	}
	return nil
}

// PrintHistory keeps convergence analysis in the terminal.
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
