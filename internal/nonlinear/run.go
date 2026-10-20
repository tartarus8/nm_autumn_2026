package nonlinear

import (
	"fmt"
	"nm/internal/labio"
	"nm/internal/linalg"
)

type Config struct {
	labio.Parameters
	ID                 string     `json:"id"`
	Description        string     `json:"description"`
	Functions          []string   `json:"functions"`
	Derivatives        [][]string `json:"jacobian"`
	Mapping            []string   `json:"phi"`
	MappingDerivatives [][]string `json:"phi_jacobian"`
	Initial            []float64  `json:"initial"`
	Lower              []float64  `json:"lower"`
	Upper              []float64  `json:"upper"`
	Q                  float64    `json:"contraction_bound"`
	Reference          []float64  `json:"reference"`
	InitialNote        string     `json:"initial_note"`
}

func (c Config) Build() (Problem, error) {
	p := Problem{Initial: c.Initial, Lower: c.Lower, Upper: c.Upper, Q: c.Q}
	var err error
	n := len(c.Initial)
	if p.F, err = CompileVector(c.Functions, n); err != nil {
		return p, err
	}
	if p.Jacobian, err = CompileJacobian(c.Derivatives, n); err != nil {
		return p, err
	}
	if p.Phi, err = CompileVector(c.Mapping, n); err != nil {
		return p, err
	}
	if p.PhiJacobian, err = CompileJacobian(c.MappingDerivatives, n); err != nil {
		return p, err
	}
	return p, nil
}

func Run(o labio.Options, scalar bool) error {
	tasks, err := labio.ReadTasks[Config](o.Input)
	if err != nil {
		return err
	}
	for _, c := range tasks {
		if err = c.Check(); err != nil {
			return err
		}
		if scalar && len(c.Initial) != 1 {
			return fmt.Errorf("scalar command expects one unknown")
		}
		p, err := c.Build()
		if err != nil {
			return fmt.Errorf("%s: %w", c.ID, err)
		}
		fmt.Printf("Task %s: %s\ninitial = %v; epsilon = %g\nGraphical choice: %s\n", c.ID, c.Description, c.Initial, c.Epsilon, c.InitialNote)
		for _, newton := range []bool{false, true} {
			name, metric := "fixed_point", "error_bound"
			if newton {
				name, metric = "newton", "step_norm"
			}
			history, solveErr := Solve(p, c.Epsilon, c.MaxIterations, newton)
			if solveErr != nil {
				return fmt.Errorf("%s/%s: %w", c.ID, name, solveErr)
			}
			last := history[len(history)-1]
			fmt.Printf("%s: x = %.12g; iterations = %d; %s = %.3e; residual = %.3e\n", name, last.X, last.Iteration, metric, last.Estimate, last.Residual)
			if o.ExtraInfo {
				labio.PrintHistory(name+" (estimate column: "+metric+")", history, c.Reference)
				if newton {
					for i := 1; i < len(history); i++ {
						previous := history[i-1].X
						linalg.PrintMatrix(fmt.Sprintf("Jacobian at iteration %d", i-1), p.Jacobian(previous))
						correction := make([]float64, len(previous))
						for j := range previous {
							correction[j] = previous[j] - history[i].X[j]
						}
						fmt.Printf("F(x) = %.10g; correction delta = %.10g; next x = x - delta\n", p.F(previous), correction)
					}
				} else {
					fmt.Printf("Supplied bound ||Phi'(x)||inf <= %g on the configured region\n", p.Q)
				}
			}
			if c.Reference != nil {
				fmt.Printf("reference error = %.3e\n", linalg.Distance(last.X, c.Reference))
			}
		}
		fmt.Println()
	}
	return nil
}
