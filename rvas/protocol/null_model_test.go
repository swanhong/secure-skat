package protocol

import (
	"crypto/rand"
	"fmt"
	"math"
	"net"
	"os"
	"path/filepath"
	"runtime/debug"
	"sync"
	"testing"
	"time"

	mpc_core "github.com/hhcho/mpc-core"
	"github.com/hhcho/sfgwas/mpc"
	"gonum.org/v1/gonum/mat"
)

// Exercise the real three-party fixed-point operations, including a design
// whose unnormalized Cholesky pivots exceed the 70/40 normalizer's range.
func TestNormalizedNullAndFinalize(t *testing.T) {
	const rows, siteARows, covariates, phenotypes = 32, 11, 5, 2
	x, y := normalizationFixture(rows, siteARows)
	wantBeta, wantOmega, wantRSS := normalizationOracle(t, x, y)
	scales := normalizationScales(x)
	var gram mat.Dense
	gram.Mul(x.T(), x)
	const replicas = 32
	for _, column := range []int{1, 3} { // PC and age squared
		pivot := normalizationPivot(t, x, column)
		if pivot*replicas <= math.Ldexp(1, 30) || pivot/rows >= math.Ldexp(1, 30) {
			t.Fatalf("fixture column %d does not cross the raw-only range boundary: pivot=%g", column, pivot)
		}
	}

	type result struct {
		beta, omega [][]float64
		rss, b, z   []float64
	}
	results := runNormalizationParties(t, func(obj *mpc.MPC) []result {
		var output []result
		for _, copies := range []int{1, replicas} {
			params := DataParams{N: rows * copies, NA: siteARows * copies,
				NB: (rows - siteARows) * copies, C: covariates, PhenotypeCount: phenotypes,
				Genes: []Gene{{GeneID: "synthetic", VariantCount: 4}}}
			var localX, localY *mat.Dense
			if obj.GetPid() != auxiliaryPartyID {
				start, end := 0, siteARows
				if obj.GetPid() == cohortBPartyID {
					start, end = siteARows, rows
				}
				localX = replicateNormalizationRows(x, start, end, copies)
				localY = replicateNormalizationRows(y, start, end, copies)
			}
			observe := func(string) func() { return func() {} }
			beta, omega, rss := SetupNull(obj, params, localX, localY, observe)
			r := result{
				beta:  obj.RevealSymMat(beta).ToFloat(40),
				omega: obj.RevealSymMat(omega).ToFloat(40),
				rss:   obj.RevealSymVec(rss).ToFloat(40),
			}

			// Known normalized moments of a rank-four equal-eigenvalue kernel.
			// Use independent RSS/N inputs here so an N-factor error in SetupNull
			// cannot cancel an opposite error in Finalize.
			share := func(values ...float64) mpc_core.RVec {
				out := mpc_core.InitRVec(obj.GetRType().Zero(), len(values))
				if obj.GetPid() == obj.GetHubPid() {
					out = mpc_core.FloatToRVec(obj.GetRType(), values, 40)
				}
				return out
			}
			b, z := Finalize([]*mpc.MPC{obj}, params,
				mpc_core.RMat{share(8, 2)}, mpc_core.RMat{share(2, 1)},
				mpc_core.RMat{share(4, 1)}, mpc_core.RMat{share(1, 0.5)},
				share(2, 0.5), share(1.5), share(0.25), share(0.25), share(0.0625), observe)
			r.b, r.z = obj.RevealSymVec(b).ToFloat(40), obj.RevealSymVec(z).ToFloat(40)
			output = append(output, r)
		}
		return output
	})

	for pid := cohortAPartyID; pid <= cohortBPartyID; pid++ {
		for index, copies := range []int{1, replicas} {
			t.Run(fmt.Sprintf("party%d/replicas%d", pid, copies), func(t *testing.T) {
				got := results[pid][index]
				for row := 0; row < covariates; row++ {
					for col := 0; col < phenotypes; col++ {
						checkNormalizedClose(t, fmt.Sprintf("scaled beta[%d,%d]", row, col), got.beta[row][col]*scales[row], wantBeta.At(row, col)*scales[row], 2e-4)
					}
					for col := 0; col < covariates; col++ {
						checkNormalizedClose(t, fmt.Sprintf("omega[%d,%d]", row, col), got.omega[row][col], wantOmega.At(row, col), 2e-4)
						// Check the inverse in balanced column units, so small
						// PC/age-squared inverse entries cannot pass merely
						// because their absolute magnitudes are small.
						product := 0.0
						for k := 0; k < covariates; k++ {
							product += gram.At(row, k) / rows * got.omega[k][col] * scales[col] / scales[row]
						}
						identity := 0.0
						if row == col {
							identity = 1
						}
						checkNormalizedClose(t, fmt.Sprintf("balanced Gram*omega[%d,%d]", row, col), product, identity, 2e-3)
					}
				}
				for col := 0; col < phenotypes; col++ {
					checkNormalizedClose(t, fmt.Sprintf("RSS/N[%d]", col), got.rss[col], wantRSS[col], 2e-4)
					rssMean := []float64{2, 0.5}[col]
					alpha := (1 - float64(covariates)/float64(rows*copies)) / (2 * rssMean)
					linear := []float64{3, 1.5}[col]
					wantB := math.Sqrt(alpha * linear * linear / 1.5)
					q := alpha * []float64{12, 3}[col] / 4
					eta := 1.0 / 18
					wantZ := (math.Cbrt(q) - 1 + eta) / math.Sqrt(eta)
					checkNormalizedClose(t, fmt.Sprintf("burden[%d]", col), got.b[col], wantB, 2e-5)
					checkNormalizedClose(t, fmt.Sprintf("WH[%d]", col), got.z[col], wantZ, 2e-5)
				}
			})
		}
	}
}

func normalizationFixture(rows, siteARows int) (*mat.Dense, *mat.Dense) {
	x, y := mat.NewDense(rows, 5, nil), mat.NewDense(rows, 2, nil)
	for i := 0; i < rows; i++ {
		pc := 10000 * float64((i*7)%19-9) / 9
		age := float64(20 + 10*((i*13)%17))
		site := 0.0
		if i >= siteARows {
			site = 1
		}
		x.SetRow(i, []float64{1, pc, age, age * age, site})
		noise := float64((i*11)%13 - 6)
		y.SetRow(i, []float64{100 + 0.001*pc + 0.3*age - 0.0005*age*age + 2*site + noise,
			4 + 0.00002*pc - 0.002*age + 0.00001*age*age - 0.2*site + 0.07*noise})
	}
	return x, y
}

// QR on column-scaled input gives an independent least-squares oracle, and
// SVD singular values give the inverse Gram without forming normal equations.
func normalizationOracle(t *testing.T, x, y *mat.Dense) (*mat.Dense, *mat.Dense, []float64) {
	t.Helper()
	n, c := x.Dims()
	_, p := y.Dims()
	z := mat.DenseCopyOf(x)
	scales := normalizationScales(x)
	for j := 0; j < c; j++ {
		for i := 0; i < n; i++ {
			z.Set(i, j, x.At(i, j)/scales[j])
		}
	}
	var qr mat.QR
	qr.Factorize(z)
	var scaledBeta mat.Dense
	if err := qr.SolveTo(&scaledBeta, false, y); err != nil {
		t.Fatal(err)
	}
	beta := mat.NewDense(c, p, nil)
	for i := 0; i < c; i++ {
		for j := 0; j < p; j++ {
			beta.Set(i, j, scaledBeta.At(i, j)/scales[i])
		}
	}
	var fitted mat.Dense
	fitted.Mul(z, &scaledBeta)
	rss := make([]float64, p)
	for i := 0; i < n; i++ {
		for j := 0; j < p; j++ {
			r := y.At(i, j) - fitted.At(i, j)
			rss[j] += r * r / float64(n)
		}
	}
	var svd mat.SVD
	if !svd.Factorize(z, mat.SVDThin) {
		t.Fatal("oracle SVD failed")
	}
	var v mat.Dense
	svd.VTo(&v)
	singular := svd.Values(nil)
	omega := mat.NewDense(c, c, nil)
	for i := 0; i < c; i++ {
		for j := 0; j < c; j++ {
			value := 0.0
			for k, s := range singular {
				value += v.At(i, k) * v.At(j, k) / (s * s)
			}
			omega.Set(i, j, value*float64(n)/(scales[i]*scales[j]))
		}
	}
	return beta, omega, rss
}

func normalizationScales(x *mat.Dense) []float64 {
	n, c := x.Dims()
	scales := make([]float64, c)
	for j := 0; j < c; j++ {
		for i := 0; i < n; i++ {
			scales[j] = math.Max(scales[j], math.Abs(x.At(i, j)))
		}
	}
	return scales
}

func normalizationPivot(t *testing.T, x *mat.Dense, column int) float64 {
	t.Helper()
	n, _ := x.Dims()
	var qr mat.QR
	qr.Factorize(x.Slice(0, n, 0, column))
	var beta, fitted mat.Dense
	if err := qr.SolveTo(&beta, false, x.Slice(0, n, column, column+1)); err != nil {
		t.Fatal(err)
	}
	fitted.Mul(x.Slice(0, n, 0, column), &beta)
	var ss float64
	for i := 0; i < n; i++ {
		r := x.At(i, column) - fitted.At(i, 0)
		ss += r * r
	}
	return ss
}

func replicateNormalizationRows(source *mat.Dense, start, end, copies int) *mat.Dense {
	_, columns := source.Dims()
	result := mat.NewDense((end-start)*copies, columns, nil)
	for i := 0; i < (end-start)*copies; i++ {
		result.SetRow(i, source.RawRowView(start+i%(end-start)))
	}
	return result
}

func checkNormalizedClose(t *testing.T, name string, got, want, tolerance float64) {
	t.Helper()
	if math.IsNaN(got) || math.IsInf(got, 0) || math.Abs(got-want) > tolerance*math.Max(1, math.Abs(want)) {
		t.Errorf("%s: got %.12g want %.12g (tolerance %g)", name, got, want, tolerance)
	}
}

// This harness is test-only and has no dependency on local diagnostic helpers.
func runNormalizationParties[T any](t *testing.T, body func(*mpc.MPC) T) []T {
	t.Helper()
	keys := t.TempDir()
	for _, name := range []string{"global", "0_1", "0_2", "1_2"} {
		key := make([]byte, 32)
		if _, err := rand.Read(key); err != nil {
			t.Fatal(err)
		}
		if err := os.WriteFile(filepath.Join(keys, "shared_key_"+name+".bin"), key, 0o600); err != nil {
			t.Fatal(err)
		}
	}
	ports := make([]string, 3)
	listeners := make([]net.Listener, len(ports))
	for i := range ports {
		listener, err := net.Listen("tcp", "127.0.0.1:0")
		if err != nil {
			for _, opened := range listeners {
				if opened != nil {
					opened.Close()
				}
			}
			t.Fatalf("three-party integration test needs localhost sockets: %v", err)
		}
		listeners[i] = listener
		ports[i] = fmt.Sprint(listener.Addr().(*net.TCPAddr).Port)
	}
	for _, listener := range listeners {
		listener.Close()
	}
	servers := map[string]mpc.Server{
		"party0": {IpAddr: "127.0.0.1", Ports: map[string]string{"party1": ports[0], "party2": ports[1]}},
		"party1": {IpAddr: "127.0.0.1", Ports: map[string]string{"party2": ports[2]}},
		"party2": {IpAddr: "127.0.0.1", Ports: map[string]string{}},
	}
	results := make([]T, 3)
	errors := make([]any, 3)
	var wait sync.WaitGroup
	for pid := 0; pid < 3; pid++ {
		wait.Add(1)
		go func(pid int) {
			defer wait.Done()
			defer func() {
				if err := recover(); err != nil {
					errors[pid] = fmt.Sprintf("%v\n%s", err, debug.Stack())
				}
			}()
			networks := mpc.InitCommunication("127.0.0.1", servers, pid, 3, 1, keys)
			defer networks[0].CloseAll()
			obj := mpc.InitParallelMPCEnv(networks, mpc_core.LElem256Zero, 70, 40)[0]
			obj.SetHubPid(1)
			obj.SetDivSqrtMaxLen(1024)
			results[pid] = body(obj)
		}(pid)
		time.Sleep(100 * time.Millisecond)
	}
	wait.Wait()
	for pid, err := range errors {
		if err != nil {
			t.Fatalf("party %d panicked: %v", pid, err)
		}
	}
	return results
}
